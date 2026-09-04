#!/usr/bin/env python3
"""Merge per-shard single-cell sparse matrices into the library-wide matrices.

Replaces the single-task build over the merged 49 GB tracking file.  Reads only
the shard artifacts written by tracking_to_shard_sparse.py, and holds no matrix
in memory: RAM is the shards' file handles, the barcode remaps, and one output
row.

Why no summing is ever needed
-----------------------------
A shard is a (cluster, contig) partition.  feature -> contig is a function and
cell -> cluster is a function, so every (feature, cell) pair belongs to exactly
one shard.  Merged nnz is therefore additive, which also means the MatrixMarket
header can be written before any data and no counting pass is required.

  basic mode          features disjoint across shards, cells shared
                      -> each output row comes from exactly one shard
  cluster-guided      features shared across clusters, cells disjoint
                      -> an output row is assembled from several clusters'
                         disjoint column runs

Both are the same code path; basic mode is a 1-way merge.

The one property that makes remapping free: each shard's barcodes.tsv and the
global barcodes.tsv are both sorted, so local -> global is monotonic and column
indices stay ascending.  Nothing is re-sorted at any stage.

Note on output ordering: this writes each row's columns in ascending order.  The
matrices built by singlecell_tracking_to_sparse_matrix.py do NOT -- line 161
permutes a CSR's columns and nothing calls sort_indices() before tocoo():164, so
mmwrite emits them permuted.  Both are valid MatrixMarket and every common
reader accepts either, but it means output here is not byte-identical to the old
path even when the data is identical.  The acceptance criterion is same support
and same values, not identical bytes.
"""
import argparse, heapq, logging, os, shutil, subprocess, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sc_sparse_shard_io import (LEVELS, open_csr, read_lines, read_manifest,
                                write_lines)

logger = logging.getLogger(__name__)
WRITE_BATCH = 4_000_000


def _labelled(names, shard_idx):
    """(name, shard, local_index) for one shard.

    A generator expression would capture `shard_idx` from the enclosing
    comprehension's scope and every stream would report the LAST shard; only
    the outermost iterable of a genexp is bound eagerly.
    """
    return ((name, shard_idx, i) for i, name in enumerate(names))


def build_global_barcodes(shard_dirs):
    """k-way merge of the shards' sorted barcode lists.

    Returns the global sorted list and, per shard, an int32 local -> global map.
    """
    per_shard = [read_lines(os.path.join(d, "barcodes.tsv")) for d in shard_dirs]
    for d, names in zip(shard_dirs, per_shard):
        if any(names[i] > names[i + 1] for i in range(len(names) - 1)):
            raise SystemExit(f"Error, {d}/barcodes.tsv is not sorted; the merge "
                             "requires a monotonic local->global remap")
    maps = [np.full(len(names), -1, dtype=np.int32) for names in per_shard]
    streams = [_labelled(names, s) for s, names in enumerate(per_shard)]

    global_names, g = [], -1
    prev = None
    for name, s, i in heapq.merge(*streams):
        if name != prev:
            g += 1
            global_names.append(name)
            prev = name
        maps[s][i] = g
    logger.info("global barcodes: %d (from %d shards, %d local entries)",
                len(global_names), len(shard_dirs), sum(len(x) for x in per_shard))
    return global_names, maps


def build_global_features(shard_dirs, level, allow_shared):
    """Global sorted feature list plus, per global row, its (shard, local_row) sources."""
    per_shard = [read_lines(os.path.join(d, f"{level}.features.tsv")) for d in shard_dirs]
    streams = [_labelled(names, s) for s, names in enumerate(per_shard)]

    names_out, sources = [], []
    prev = None
    for name, s, i in heapq.merge(*streams):
        if name != prev:
            names_out.append(name)
            sources.append([(s, i)])
            prev = name
        else:
            sources[-1].append((s, i))
    widest = max((len(x) for x in sources), default=0)
    if not allow_shared and widest > 1:
        bad = names_out[max(range(len(sources)), key=lambda k: len(sources[k]))]
        raise SystemExit(
            f"Error, feature {bad!r} appears in {widest} shards at level {level}. "
            "Features must partition by contig; pass --shared-features if this run "
            "genuinely partitions cells instead (cluster-guided).")
    if "" in names_out:
        raise SystemExit(
            f"Error, an empty feature name is present at level {level}. Unassigned "
            "rows become \"\" and would appear in every shard, breaking the "
            "no-summing invariant this merge relies on.")
    logger.info("%s: %d global features, max %d shard(s) per feature",
                level, len(names_out), widest)
    return names_out, sources


def _flush(fh, buf_r, buf_c, buf_v):
    if not buf_r:
        return
    pd.DataFrame({"r": np.concatenate(buf_r),
                  "c": np.concatenate(buf_c),
                  "v": np.concatenate(buf_v)}).to_csv(
        fh, sep=" ", header=False, index=False, lineterminator="\n")
    buf_r.clear(); buf_c.clear(); buf_v.clear()


def merge_level(shard_dirs, level, bc_maps, n_barcodes, out_mtx, allow_shared):
    feat_names, sources = build_global_features(shard_dirs, level, allow_shared)

    handles, total_nnz = [], 0
    for d in shard_dirs:
        n_rows, _, nnz, indptr, indices, data = open_csr(
            os.path.join(d, f"{level}.csr"))
        handles.append((indptr, indices, data))
        total_nnz += nnz

    with open(out_mtx, "w") as fh:
        fh.write("%%MatrixMarket matrix coordinate real general\n%\n")
        fh.write(f"{len(feat_names)} {n_barcodes} {total_nnz}\n")

        buf_r, buf_c, buf_v, buffered = [], [], [], 0
        for g_row, srcs in enumerate(sources):
            cols, vals = [], []
            for s, local_row in srcs:
                indptr, indices, data = handles[s]
                lo, hi = int(indptr[local_row]), int(indptr[local_row + 1])
                if hi == lo:
                    continue
                cols.append(bc_maps[s][indices[lo:hi]])
                vals.append(np.asarray(data[lo:hi]))
            if not cols:
                continue
            if len(cols) == 1:
                c, v = cols[0], vals[0]      # already ascending: remap is monotonic
            else:
                c = np.concatenate(cols); v = np.concatenate(vals)
                order = np.argsort(c, kind="stable")   # disjoint runs, merge them
                c, v = c[order], v[order]
            buf_r.append(np.full(len(c), g_row + 1, dtype=np.int64))
            buf_c.append(c.astype(np.int64) + 1)       # MatrixMarket is 1-based
            buf_v.append(v)
            buffered += len(c)
            if buffered >= WRITE_BATCH:
                _flush(fh, buf_r, buf_c, buf_v); buffered = 0
        _flush(fh, buf_r, buf_c, buf_v)

    logger.info("%s: wrote %s (%d x %d, nnz %d)", level, out_mtx,
                len(feat_names), n_barcodes, total_nnz)
    return feat_names, total_nnz


def merge_mapping(shard_dirs, out_path):
    seen, rows = set(), []
    for d in shard_dirs:
        for line in read_lines(os.path.join(d, "mapping.tsv"))[1:]:
            if line not in seen:
                seen.add(line)
                rows.append(line)
    rows.sort()
    write_lines(out_path, ["\t".join(("gene_id", "transcript_id",
                                      "transcript_splice_hash_code",
                                      "num_exons"))] + rows)
    logger.info("mapping: %d entries -> %s", len(rows), out_path)
    return len(rows)


def gzip_file(path, level=1):
    pigz = shutil.which("pigz")
    if pigz:
        subprocess.run([pigz, f"-{level}", "-f", path], check=True)
    else:
        import gzip as _gz
        with open(path, "rb") as i, _gz.open(path + ".gz", "wb", compresslevel=level) as o:
            shutil.copyfileobj(i, o)
        os.remove(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--shard_dirs", nargs="+", required=True)
    ap.add_argument("--output_prefix", required=True)
    ap.add_argument("--gzip_level", type=int, default=1)
    ap.add_argument("--shared-features", action="store_true",
                    help="allow one feature to come from several shards "
                         "(cluster-guided, where cells partition instead)")
    ap.add_argument("--no_gzip", action="store_true", help="leave outputs plain (tests)")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, stream=sys.stderr,
                        format="%(asctime)s - %(levelname)s - %(message)s")

    dirs = list(args.shard_dirs)
    for d in dirs:
        read_manifest(d)                      # fails loudly on a malformed shard
    barcodes, bc_maps = build_global_barcodes(dirs)

    merge_mapping(dirs, f"{args.output_prefix}.gene_transcript_splicehashcode.tsv")

    for level in LEVELS:
        outdir = f"{args.output_prefix}^{level}-sparseM"
        os.makedirs(outdir, exist_ok=True)
        feats, _ = merge_level(dirs, level, bc_maps, len(barcodes),
                               os.path.join(outdir, "matrix.mtx"),
                               args.shared_features)
        write_lines(os.path.join(outdir, "features.tsv"), feats)
        write_lines(os.path.join(outdir, "barcodes.tsv"), barcodes)
        if not args.no_gzip:
            for fn in ("matrix.mtx", "features.tsv", "barcodes.tsv"):
                gzip_file(os.path.join(outdir, fn), args.gzip_level)


if __name__ == "__main__":
    main()
