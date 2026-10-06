#!/usr/bin/env python3

"""Turn site_read_support_to_sparse_matrix.py's per-cell TSS / PolyA read support into
the count files site_usage_dexseq.R reads, in place of count_site_read_ends.py.

The site support is then LRAA's own: LRAA's sites (the integrated beds), and an end
counted only as LRAA counts site support (no residual soft clip at the end, within
int(max_dist_between_alt_*_sites / 2) of the site). count_site_read_ends.py instead
recounts every read end within a window of its own.

Rows follow the site table from prep_site_table.py (run with no merging, so its site
ids are the beds' site names); a site with no supporting read is a zero row. Columns
are the clustered cells; reads from other cells are dropped.

Outputs, under --output_prefix (as count_site_read_ends.py writes them):
  .site_counts.mtx.gz   sites x cells MatrixMarket (integer), rows in site-table order
  .barcodes.tsv.gz      the cell barcodes, in column order, with their cluster
  .sites.tsv.gz         the site ids, in row order
  .summary.tsv          the support utility's read / end tallies
  .ends_at_no_site.tsv  header only: the support utility does not tally ends at no site per gene
"""

import argparse
import csv
import gzip
import logging
import os
import shutil
import sys

import numpy as np
from scipy import sparse
from scipy.io import mmread

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s", level=logging.INFO)

KINDS = ("TSS", "PolyA")


def main():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sites", required=True, help="site table from prep_site_table.py (no merging)")
    parser.add_argument("--support_prefix", required=True,
                        help="--output_prefix given to site_read_support_to_sparse_matrix.py")
    parser.add_argument("--cell_clusters", required=True,
                        help="cell_barcode <tab> cluster, one per line; a header line is skipped")
    parser.add_argument("--output_prefix", required=True)
    args = parser.parse_args()

    sites = list(csv.DictReader(open(args.sites), delimiter="\t"))
    row_of = {s["site_id"]: i for i, s in enumerate(sites)}
    merged = sum(1 for s in sites if int(s["n_merged"]) > 1)
    if merged:
        logger.warning("%d sites in %s merge several bed sites; only their representative's support is used",
                       merged, args.sites)

    barcodes, cluster_of = parse_cell_clusters(args.cell_clusters)
    col_of = {cb: j for j, cb in enumerate(barcodes)}

    blocks = []
    n_missing = 0
    for kind in KINDS:
        d = f"{args.support_prefix}.{kind}-sparseM"
        m = mmread(gzip.open(os.path.join(d, "matrix.mtx.gz"))).tocoo()
        features = [l.rstrip("\n") for l in gzip.open(os.path.join(d, "features.tsv.gz"), "rt")]
        cells = [l.rstrip("\n") for l in gzip.open(os.path.join(d, "barcodes.tsv.gz"), "rt")]
        frow = np.array([row_of.get(f, -1) for f in features])
        fcol = np.array([col_of.get(c, -1) for c in cells])
        n_missing += int((frow < 0).sum())
        r, c, v = frow[m.row], fcol[m.col], m.data
        keep = (r >= 0) & (c >= 0)
        blocks.append((r[keep], c[keep], v[keep]))
        logger.info("%s: %d sites, %d cells in the support matrix; %s of %s read ends kept (clustered cells, tabled sites)",
                    kind, len(features), len(cells), f"{v[keep].sum():,.0f}", f"{v.sum():,.0f}")
    if n_missing:
        logger.warning("%d support-matrix sites are not in the site table", n_missing)

    r, c, v = (np.concatenate(x) for x in zip(*blocks))
    counts = sparse.coo_matrix((v.astype(np.int64), (r, c)), shape=(len(sites), len(barcodes))).tocsr().tocoo()

    write_mtx(args.output_prefix + ".site_counts.mtx.gz", counts)
    with gzip.open(args.output_prefix + ".barcodes.tsv.gz", "wt") as ofh:
        for cb in barcodes:
            print(f"{cb}\t{cluster_of[cb]}", file=ofh)
    with gzip.open(args.output_prefix + ".sites.tsv.gz", "wt") as ofh:
        for s in sites:
            print(s["site_id"], file=ofh)
    shutil.copyfile(f"{args.support_prefix}.site_read_support.summary.tsv", args.output_prefix + ".summary.tsv")
    with open(args.output_prefix + ".ends_at_no_site.tsv", "wt") as ofh:
        print("kind\tgene_key\tcluster\tends_at_no_site", file=ofh)


def parse_cell_clusters(filename):
    barcodes, cluster_of = [], {}
    for line in open(filename):
        f = line.rstrip("\n").split("\t")
        if len(f) < 2 or not f[1].strip().lstrip("-").isdigit():
            continue
        if f[0] not in cluster_of:
            barcodes.append(f[0])
        cluster_of[f[0]] = "Cluster_" + f[1].strip()
    return barcodes, cluster_of


def write_mtx(filename, m):
    order = np.lexsort((m.row, m.col))
    with gzip.open(filename, "wt", compresslevel=4) as ofh:
        print("%%MatrixMarket matrix coordinate integer general", file=ofh)
        print(f"{m.shape[0]} {m.shape[1]} {m.nnz}", file=ofh)
        for r, c, v in zip(m.row[order] + 1, m.col[order] + 1, m.data[order]):
            ofh.write(f"{r} {c} {v}\n")


if __name__ == "__main__":
    main()
