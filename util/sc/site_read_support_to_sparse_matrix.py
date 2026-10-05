#!/usr/bin/env python3

"""Per-cell read support for LRAA's TSS and PolyA sites, as sparse matrices.

LRAA reports each TSS / PolyA site with one support value for the whole run (the
splice graph's sum of XW normalization weights). A test of whether site usage
differs between cells or cell clusters needs that support per cell. This utility
recounts it from the alignments with LRAA's own rules for which read ends may
support a site, so the counts are LRAA's site support broken down by cell barcode:

  - reads LRAA would keep (Util_funcs.quant_discard_reason: mapping quality,
    percent identity, secondary / supplementary / duplicate, long introns,
    polyA-terminal segments, ...), on their transcribed strand (minimap2 ts tag);
  - the read's ends taken as LRAA takes them (Pretty_alignment): the 5' end of the
    alignment is its TSS end and the 3' end its PolyA end, after LRAA strips a polyA
    tail or untemplated 5' G's from the soft clips;
  - an end supports a site only if its remaining soft clip is within
    max_soft_clip_at_TSS / max_soft_clip_at_PolyA (0 by default), and lies within
    int(max_dist_between_alt_{TSS,polyA}_sites / 2) nt of the site, the tolerance
    LRAA uses to call a read end the same site (sites of one run are more than the
    full distance apart, so an end can match at most one site; the nearest is taken).

Counts are reads (one per read end), not normalization weights, unless --weighted:
a bam thinned by LRAA's coverage normalization carries the XW weight, and summing it
estimates the unthinned count, as LRAA's own support does. Reads without a cell
barcode (config cell_barcode_tag, CB) are not counted.

Site beds are those LRAA writes (SplicePatternCollapse.write_site_bed: chrom, start,
end, name, score, strand, support, ...; the PolyA bed adds pas, pas_offset,
internal_priming), or the integrated beds of integrate_TSS_PolyA_sites.py (an extra
`source` column). Every site is counted, including PolyA sites flagged as internally
primed; the flag is carried into the site table for the caller to act on.

Outputs, for each site kind given (KIND = TSS | PolyA):
  <prefix>.<KIND>-sparseM/{matrix.mtx,features.tsv,barcodes.tsv}.gz
      sites x cells, in the layout of LRAA's other sparse matrices (feature = site
      name, e.g. TSS:chr1:29359:-). Sites with no supporting read are included.
  <prefix>.<KIND>.sites.tsv
      per site: bed fields, reads counted, cells with reads, and the bed's own support.
  <prefix>.<KIND>.cluster_counts.tsv   (with --cell_clusters)
      sites x clusters pseudobulk counts.
  <prefix>.site_read_support.summary.tsv
      reads seen / discarded by reason; per kind, ends at a site, ends rejected for
      soft clipping, ends at no site.
"""

import argparse
import bisect
import collections
import csv
import gzip
import logging
import os
import sys
from multiprocessing import Pool

import numpy as np
import pysam
from scipy import sparse
from scipy.io import mmwrite

sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../../pylib"]))

import LRAA_Globals  # noqa: E402
import Util_funcs  # noqa: E402
from Pretty_alignment import Pretty_alignment  # noqa: E402

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s", level=logging.INFO)

KINDS = ("TSS", "PolyA")
_DIST_KEY = {"TSS": "max_dist_between_alt_TSS_sites", "PolyA": "max_dist_between_alt_polyA_sites"}
_CLIP_KEY = {"TSS": "max_soft_clip_at_TSS", "PolyA": "max_soft_clip_at_PolyA"}

CHUNK_LEN = 5_000_000

# set in the parent before the pool forks
G = {}


def main():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bam", required=True, help="aligned reads with cell-barcode tags (indexed)")
    parser.add_argument("--TSS_bed", default=None, help="LRAA TSS site bed")
    parser.add_argument("--PolyA_bed", default=None, help="LRAA PolyA site bed")
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--cell_clusters", default=None,
                        help="optional: cell_barcode <tab> cluster (a header line is skipped); adds "
                             "sites x clusters pseudobulk counts. Cells are not restricted to these.")
    parser.add_argument("--HiFi", action="store_true",
                        help="apply LRAA's HiFi read-identity floor (min_per_id {}), as the LRAA run "
                             "given --HiFi did".format(LRAA_Globals.HIFI_MIN_PER_ID))
    parser.add_argument("--min_per_id", type=float, default=None,
                        help="override the percent-identity floor (default: LRAA config, or the HiFi one)")
    parser.add_argument("--min_mapping_quality", type=int, default=None,
                        help="override the mapping-quality floor (default: LRAA config min_mapping_quality)")
    parser.add_argument("--weighted", action="store_true",
                        help="sum XW normalization weights instead of counting reads")
    parser.add_argument("--CPU", type=int, default=4)
    args = parser.parse_args()

    beds = {k: b for k, b in (("TSS", args.TSS_bed), ("PolyA", args.PolyA_bed)) if b}
    if not beds:
        sys.exit("Error, give --TSS_bed and/or --PolyA_bed")

    min_per_id = args.min_per_id
    if min_per_id is None:
        min_per_id = LRAA_Globals.HIFI_MIN_PER_ID if args.HiFi else LRAA_Globals.config["min_per_id"]
    min_mapq = args.min_mapping_quality
    if min_mapq is None:
        min_mapq = int(LRAA_Globals.config["min_mapping_quality"])

    sites = {k: read_site_bed(b, k) for k, b in beds.items()}
    for k, s in sites.items():
        logger.info("%s: %d sites from %s", k, len(s), beds[k])

    G.update(bam=args.bam, site_index={k: index_sites(s) for k, s in sites.items()},
             tolerance={k: int(LRAA_Globals.config[_DIST_KEY[k]] / 2) for k in KINDS},
             max_clip={k: LRAA_Globals.config[_CLIP_KEY[k]] for k in KINDS},
             min_per_id=min_per_id, min_mapq=min_mapq, weighted=args.weighted,
             cb_tag=LRAA_Globals.config["cell_barcode_tag"])
    logger.info("read filter: min_per_id %s, min_mapping_quality %s; end-to-site tolerance %s; "
                "max soft clip %s; counting %s", min_per_id, min_mapq, G["tolerance"], G["max_clip"],
                "XW weights" if args.weighted else "reads")

    with pysam.AlignmentFile(args.bam) as bam:
        chunks = [(c, lo, min(lo + CHUNK_LEN, bam.get_reference_length(c)))
                  for c in bam.references for lo in range(0, bam.get_reference_length(c), CHUNK_LEN)]

    barcodes, cb_col = [], {}
    triplets = {k: ([], [], []) for k in sites}
    summary = collections.Counter()
    with Pool(args.CPU) as pool:
        for i, (chunk_barcodes, counts, summ) in enumerate(pool.imap_unordered(count_chunk, chunks)):
            remap = []
            for cb in chunk_barcodes:
                if cb not in cb_col:
                    cb_col[cb] = len(barcodes)
                    barcodes.append(cb)
                remap.append(cb_col[cb])
            remap = np.array(remap, dtype=np.int64)
            for k, (r, c, v) in counts.items():
                if len(r):
                    triplets[k][0].append(r)
                    triplets[k][1].append(remap[c])
                    triplets[k][2].append(v)
            summary.update(summ)
            if (i + 1) % 100 == 0:
                logger.info("%d / %d chunks, %d reads used", i + 1, len(chunks), summary["reads_used"])

    clusters = parse_cell_clusters(args.cell_clusters) if args.cell_clusters else None

    for k, s in sites.items():
        r, c, v = (np.concatenate(x) if x else np.zeros(0) for x in triplets[k])
        m = sparse.coo_matrix((v, (r.astype(np.int64), c.astype(np.int64))),
                              shape=(len(s), len(barcodes))).tocsr()  # sums duplicate entries
        if not args.weighted:
            m = m.astype(np.int64)
        write_sparseM(f"{args.output_prefix}.{k}-sparseM", m, [x["name"] for x in s], barcodes)
        write_site_table(f"{args.output_prefix}.{k}.sites.tsv", s, m)
        if clusters is not None:
            write_cluster_counts(f"{args.output_prefix}.{k}.cluster_counts.tsv", s, m, barcodes, clusters)
        logger.info("%s: %s read ends at %d sites in %d cells", k, f"{m.sum():,.0f}", (m.getnnz(axis=1) > 0).sum(),
                    (m.getnnz(axis=0) > 0).sum())

    with open(f"{args.output_prefix}.site_read_support.summary.tsv", "wt") as ofh:
        print("count\tvalue", file=ofh)
        for key, val in sorted(summary.items()):
            print(f"{key}\t{val}", file=ofh)
            logger.info("%s: %s", key, val)


def read_site_bed(path, kind):
    """sites in file order: chrom, pos (bed end, the 1-based site), name, strand, support, and
    for PolyA the internal-priming flag; '#' lines (provenance, column header) are skipped"""
    sites = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                sys.exit(f"{path}: expected an LRAA site bed (>= 7 columns), got: {line[:200]}")
            site = {"chrom": f[0], "pos": int(f[2]), "name": f[3], "strand": f[5], "support": f[6],
                    "internal_priming": f[11] if kind == "PolyA" and len(f) > 11 else "",
                    "source": f[-1] if f[-1] in ("cluster_guided", "basic") else ""}
            if site["strand"] not in ("+", "-"):
                sys.exit(f"{path}: site {site['name']} has strand {site['strand']!r}")
            sites.append(site)
    return sites


def index_sites(sites):
    """(chrom, strand) -> (sorted positions, row numbers)"""
    idx = collections.defaultdict(list)
    for i, s in enumerate(sites):
        idx[(s["chrom"], s["strand"])].append((s["pos"], i))
    return {k: (np.array([p for p, _ in v]), np.array([i for _, i in v]))
            for k, v in ((k, sorted(v)) for k, v in idx.items())}


def nearest_site(index, pos, tolerance):
    positions, rows = index
    i = bisect.bisect_left(positions, pos)
    best, best_d = None, None
    for j in (i - 1, i):
        if 0 <= j < len(positions):
            d = abs(int(positions[j]) - pos)
            if d <= tolerance and (best_d is None or d < best_d):
                best, best_d = rows[j], d
    return best


def count_chunk(chunk):
    contig, lo, hi = chunk
    site_index, tol, max_clip = G["site_index"], G["tolerance"], G["max_clip"]
    cb_tag, weighted = G["cb_tag"], G["weighted"]
    counts = {k: collections.Counter() for k in site_index}
    cb_col, barcodes = {}, []
    summ = collections.Counter()

    with pysam.AlignmentFile(G["bam"]) as bam:
        for read in bam.fetch(contig, lo, hi):
            # each read once: in the chunk holding its alignment start
            if read.reference_start < lo:
                continue
            summ["reads_seen"] += 1
            reason = Util_funcs.quant_discard_reason(read, None, min_mapping_quality=G["min_mapq"],
                                                     min_per_id=G["min_per_id"])
            if reason is not None:
                summ[f"reads_discarded:{reason}"] += 1
                continue
            if not read.has_tag(cb_tag):
                summ["reads_without_cell_barcode"] += 1
                continue
            summ["reads_used"] += 1

            pa = Pretty_alignment.get_pretty_alignment(read)
            strand = pa.get_strand()
            lend, rend = pa.get_alignment_span()
            plus = strand == "+"
            ends = {"TSS": (lend if plus else rend, pa.left_soft_clipping if plus else pa.right_soft_clipping),
                    "PolyA": (rend if plus else lend, pa.right_soft_clipping if plus else pa.left_soft_clipping)}
            col = None
            for kind, index in site_index.items():
                pos, clip = ends[kind]
                if clip > max_clip[kind]:
                    summ[f"{kind}_ends_soft_clipped"] += 1
                    continue
                idx = index.get((contig, strand))
                row = nearest_site(idx, pos, tol[kind]) if idx is not None else None
                if row is None:
                    summ[f"{kind}_ends_at_no_site"] += 1
                    continue
                summ[f"{kind}_ends_at_site"] += 1
                if col is None:
                    cb = read.get_tag(cb_tag)
                    col = cb_col.get(cb)
                    if col is None:
                        col = cb_col[cb] = len(barcodes)
                        barcodes.append(cb)
                counts[kind][(row, col)] += pa.get_normalization_weight() if weighted else 1

    out = {}
    for kind, cnt in counts.items():
        keys = list(cnt)
        out[kind] = (np.fromiter((k[0] for k in keys), dtype=np.int64, count=len(keys)),
                     np.fromiter((k[1] for k in keys), dtype=np.int64, count=len(keys)),
                     np.fromiter((cnt[k] for k in keys), dtype=float, count=len(keys)))
    return barcodes, out, summ


def write_sparseM(outdir, m, features, barcodes):
    """LRAA's sparse-matrix layout (as singlecell_tracking_to_sparse_matrix.py): matrix.mtx,
    features.tsv, barcodes.tsv, gzipped"""
    os.makedirs(outdir, exist_ok=True)
    mtx = os.path.join(outdir, "matrix.mtx")
    mmwrite(mtx, m)
    for fn, items in (("features.tsv", features), ("barcodes.tsv", barcodes)):
        with open(os.path.join(outdir, fn), "wt") as ofh:
            for x in items:
                print(x, file=ofh)
    for fn in ("matrix.mtx", "features.tsv", "barcodes.tsv"):
        path = os.path.join(outdir, fn)
        with open(path, "rb") as fin, gzip.open(path + ".gz", "wb") as fout:
            fout.writelines(fin)
        os.remove(path)


def write_site_table(path, sites, m):
    reads = np.asarray(m.sum(axis=1)).ravel()
    cells = m.getnnz(axis=1)
    with open(path, "wt") as ofh:
        w = csv.writer(ofh, delimiter="\t", lineterminator="\n")
        w.writerow(["site", "chrom", "pos", "strand", "source", "internal_priming", "bed_support",
                    "reads", "cells"])
        for s, r, c in zip(sites, reads, cells):
            w.writerow([s["name"], s["chrom"], s["pos"], s["strand"], s["source"], s["internal_priming"],
                        s["support"], round(float(r), 3), int(c)])


def parse_cell_clusters(path):
    clusters = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 2 or f[1].strip().lower() in ("cluster", "seurat_clusters", "cell_cluster"):
            continue
        clusters[f[0]] = f[1].strip()
    return clusters


def write_cluster_counts(path, sites, m, barcodes, clusters):
    labels = [clusters.get(cb) for cb in barcodes]
    names = sorted({x for x in labels if x is not None}, key=lambda x: (not x.lstrip("-").isdigit(),
                                                                         int(x) if x.lstrip("-").isdigit() else 0, x))
    col = {n: i for i, n in enumerate(names)}
    keep = [i for i, x in enumerate(labels) if x is not None]
    agg = sparse.csr_matrix((np.ones(len(keep)), (keep, [col[labels[i]] for i in keep])),
                            shape=(len(barcodes), len(names)))
    pb = (m @ agg).toarray()
    with open(path, "wt") as ofh:
        print("\t".join(["site"] + names), file=ofh)
        for s, row in zip(sites, pb):
            print("\t".join([s["name"]] + [f"{x:g}" for x in row]), file=ofh)


if __name__ == "__main__":
    main()
