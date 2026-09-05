#!/usr/bin/env python3
"""Build one shard's single-cell sparse matrices from that shard's quant.tracking.

Runs inside LRAA_runner_task, on the shard's own tracking file, so the work that
build_sc_sparse_matrices does once over the whole library instead happens 25 (or
clusters x 25) times in parallel.  The merge is then streaming I/O -- see
merge_shard_sparse_matrices.py.

Only the isoform level is accumulated.  gene and splice_pattern are exact row
aggregations of it, because transcript_id is a key of the mapping table and both
gene_id and the splice hash are functions of it (verified on run A: 327,187
transcripts, 0 mapping to two genes or two hashes).  That is a third of the
memory and a third of the finalize work.

No CSV parser is used: the tracking file is machine-generated TSV with no
quoting, and a plain split() loop measured 4.1x faster than pandas' python
engine and 1.7x faster than its C engine, with no exposure to the segfaults that
made the C engine unusable here.
"""
import argparse, gzip, json, logging, os, sys
from sys import intern

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sc_sparse_shard_io import (LEVELS, write_csr_with_shape, write_lines,
                                write_manifest)

logger = logging.getLogger(__name__)

NEEDED = ["gene_id", "transcript_id", "transcript_splice_hash_code",
          "num_exons", "read_name", "frac_assigned"]


def _index_dtype(n):
    """int32 while it is safe; the column axis is barcodes, which can grow."""
    return np.int32 if n < np.iinfo(np.int32).max else np.int64


class ShardSparseBuilder:
    """Accumulate one shard's isoform-level triples, then emit all three levels."""

    def __init__(self, shard):
        self.shard = shard
        self.tx2i, self.tx_labels = {}, []
        self.tx_gene, self.tx_hash, self.tx_nexon = [], [], []
        self.bc2i, self.bc_labels = {}, []
        self.rows, self.cols, self.vals = [], [], []
        self.n_rows_in = 0

    def add_chunk(self, gene, tx, hsh, nexon, barcode, frac):
        """All six arguments are equal-length sequences of one chunk's fields."""
        self.n_rows_in += len(tx)

        tx_arr = np.asarray(tx, dtype=object)
        tx_codes, tx_uniq = pd.factorize(tx_arr)
        tx_lut = np.empty(len(tx_uniq), dtype=np.int64)
        for k, name in enumerate(tx_uniq):
            idx = self.tx2i.get(name)
            if idx is None:
                idx = self.tx2i[name] = len(self.tx_labels)
                self.tx_labels.append(name)
                # mapping row captured once per transcript, not once per input row
                first = int(np.flatnonzero(tx_codes == k)[0])
                self.tx_gene.append(gene[first])
                self.tx_hash.append(hsh[first])
                self.tx_nexon.append(int(nexon[first]))
            tx_lut[k] = idx
        feat_idx = tx_lut[tx_codes]

        bc_arr = np.asarray(barcode, dtype=object)
        bc_codes, bc_uniq = pd.factorize(bc_arr)
        bc_lut = np.empty(len(bc_uniq), dtype=np.int64)
        for k, name in enumerate(bc_uniq):
            idx = self.bc2i.get(name)
            if idx is None:
                idx = self.bc2i[name] = len(self.bc_labels)
                self.bc_labels.append(name)
            bc_lut[k] = idx
        bc_idx = bc_lut[bc_codes]

        grouped = (pd.DataFrame({"f": feat_idx, "b": bc_idx,
                                 "v": np.asarray(frac, dtype=np.float32)})
                   .groupby(["f", "b"], observed=True)["v"].sum())
        if grouped.empty:
            return
        self.rows.append(grouped.index.get_level_values(0).to_numpy(dtype=np.int32))
        self.cols.append(grouped.index.get_level_values(1).to_numpy(dtype=np.int32))
        self.vals.append(grouped.to_numpy(dtype=np.float32))

    # ---------------------------------------------------------------- emit ---

    def _isoform_csr(self, bc_perm_inv):
        n_tx, n_bc = len(self.tx_labels), len(self.bc_labels)
        if self.rows:
            r = np.concatenate(self.rows); self.rows = None
            c = np.concatenate(self.cols); self.cols = None
            v = np.concatenate(self.vals); self.vals = None
        else:
            r = np.empty(0, np.int32); c = np.empty(0, np.int32)
            v = np.empty(0, np.float32)
        c = bc_perm_inv[c].astype(_index_dtype(n_bc), copy=False)
        # tocsr() sums duplicates itself; an explicit sum_duplicates() here would
        # cost a full lexsort over every non-zero and change nothing.
        return sparse.coo_matrix((v, (r, c)), shape=(n_tx, n_bc)).tocsr()

    @staticmethod
    def _regroup_rows(csr, new_row_of_row, n_new_rows):
        """Sum CSR rows into coarser rows.  Used to derive gene and splice."""
        counts = np.diff(csr.indptr)
        rows = np.repeat(new_row_of_row, counts)
        return sparse.coo_matrix((csr.data, (rows, csr.indices)),
                                 shape=(n_new_rows, csr.shape[1])).tocsr()

    @staticmethod
    def _sort_rows(csr, labels):
        order = np.argsort(np.asarray(labels, dtype=object), kind="stable")
        return csr[order], [labels[i] for i in order]

    def write(self, outdir):
        os.makedirs(outdir, exist_ok=True)

        # Sort this shard's barcodes explicitly.  The tracking file happens to be
        # read_name-ordered, but nothing states that as a contract; an unsorted
        # barcodes.tsv would make the merge's local->global remap non-monotonic
        # and produce a wrong matrix with nothing about it to look wrong.
        bc_arr = np.asarray(self.bc_labels, dtype=object)
        bc_order = np.argsort(bc_arr, kind="stable")
        bc_perm_inv = np.empty(len(bc_order), dtype=np.int64)
        bc_perm_inv[bc_order] = np.arange(len(bc_order))
        sorted_barcodes = [self.bc_labels[i] for i in bc_order]
        write_lines(os.path.join(outdir, "barcodes.tsv"), sorted_barcodes)

        iso = self._isoform_csr(bc_perm_inv)

        gene_names = sorted(set(self.tx_gene))
        gene_of = {g: i for i, g in enumerate(gene_names)}
        hash_names = sorted(set(self.tx_hash))
        hash_of = {h: i for i, h in enumerate(hash_names)}
        gene_row = np.fromiter((gene_of[g] for g in self.tx_gene),
                               dtype=np.int64, count=len(self.tx_gene))
        hash_row = np.fromiter((hash_of[h] for h in self.tx_hash),
                               dtype=np.int64, count=len(self.tx_hash))

        info = {}
        for level in LEVELS:
            if level == "isoform":
                m, labels = self._sort_rows(iso, self.tx_labels)
            elif level == "gene":
                m = self._regroup_rows(iso, gene_row, len(gene_names))
                labels = gene_names           # already sorted
            else:
                m = self._regroup_rows(iso, hash_row, len(hash_names))
                labels = hash_names           # already sorted
            m.sort_indices()
            write_lines(os.path.join(outdir, f"{level}.features.tsv"), labels)
            write_csr_with_shape(os.path.join(outdir, f"{level}.csr"), m)
            info[level] = {"n_features": int(m.shape[0]), "nnz": int(m.nnz)}
            logger.info("%s %s: %d x %d, nnz %d", self.shard, level,
                        m.shape[0], m.shape[1], m.nnz)
            del m

        write_lines(os.path.join(outdir, "mapping.tsv"),
                    ["\t".join(("gene_id", "transcript_id",
                                "transcript_splice_hash_code", "num_exons"))] +
                    ["\t".join((g, t, h, str(n))) for g, t, h, n in
                     sorted(zip(self.tx_gene, self.tx_labels,
                                self.tx_hash, self.tx_nexon),
                            key=lambda x: (x[0], x[1], x[2]))])
        write_manifest(outdir, self.shard, info, len(sorted_barcodes))
        return info


def stream_tracking(path, builder, chunksize=1_000_000):
    """Pure-python TSV reader; column positions resolved from the header line."""
    opener = gzip.open if path.endswith(".gz") else open
    fh = opener(path, "rt")
    try:
        while True:
            line = fh.readline()
            if not line:
                raise SystemExit(f"Error, {path} has no header line")
            if not line.startswith("#"):
                break
        hdr = line.rstrip("\n").split("\t")
        missing = [c for c in NEEDED if c not in hdr]
        if missing:
            raise SystemExit(f"Error, {path} lacks column(s): {missing}")
        ig, it, ih, ine, irn, ifr = (hdr.index(c) for c in NEEDED)

        while True:
            g = []; t = []; h = []; nx = []; b = []; f = []
            ga, ta, ha, na, ba, fa = (g.append, t.append, h.append,
                                      nx.append, b.append, f.append)
            for _ in range(chunksize):
                l = fh.readline()
                if not l:
                    break
                p = l.split("\t")
                # interning the low-cardinality ids makes every downstream dict
                # lookup a pointer compare; measured 1.29 s/Mrow off the loop
                ga(intern(p[ig])); ta(intern(p[it])); ha(intern(p[ih]))
                na(p[ine]); ba(intern(p[irn].partition("^")[0])); fa(p[ifr])
            if not t:
                break
            builder.add_chunk(g, t, h,
                              np.array(nx, dtype=np.int64), b,
                              np.array(f, dtype=np.float32))
            if len(t) < chunksize:
                break
    finally:
        fh.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tracking", required=True, help="this shard's quant.tracking[.gz]")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--shard", required=True, help="shard name, e.g. cluster_3/chr7")
    ap.add_argument("--chunksize", type=int, default=1_000_000)
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, stream=sys.stderr,
                        format="%(asctime)s - %(levelname)s - %(message)s")
    b = ShardSparseBuilder(args.shard)
    stream_tracking(args.tracking, b, args.chunksize)
    logger.info("%s: read %d rows, %d transcripts, %d barcodes",
                args.shard, b.n_rows_in, len(b.tx_labels), len(b.bc_labels))
    b.write(args.outdir)
    logger.info("%s: wrote %s", args.shard, args.outdir)


if __name__ == "__main__":
    main()
