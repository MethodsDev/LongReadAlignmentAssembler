#!/usr/bin/env python3

"""Count, per cell, the reads whose 5' end lies at each TSS site and whose 3' end
lies at each PolyA site.

One pass over the BAM, parallel over contigs. Each read's 5' end is assigned to
the nearest TSS site of the read's transcript strand within that site's window,
and its 3' end likewise to a PolyA site (windows and site spans from
prep_site_table.py). Reads are kept on their transcript strand: the alignment
orientation, flipped when minimap2's ts tag is "-"; reads without the tag are
taken as oriented to the transcript.

Read ends lying at no site are not counted toward any site but are tallied per
gene (the one gene symbol whose span covers the end on that strand) and cell
cluster, so a report can show what share of a gene's read ends the sites leave
out (internally primed A-runs, 5'-truncated reads and the like).

Outputs, under --output_prefix:
  .site_counts.mtx.gz   sites x cells MatrixMarket (integer), rows in site-table order
  .barcodes.tsv.gz      the cell barcodes, in column order, with their cluster
  .sites.tsv.gz         the site ids, in row order
  .ends_at_no_site.tsv  kind, gene_key, cluster, read ends at no site
  .summary.tsv          reads seen / used, ends assigned to a site or not, per kind
"""

import argparse
import bisect
import collections
import csv
import gzip
import logging
from multiprocessing import Pool

import numpy as np
import pysam

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s", level=logging.INFO)

KINDS = ("TSS", "PolyA")

# set in the parent before the pool forks
G = {}


def main():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sites", required=True, help="site table from prep_site_table.py")
    parser.add_argument("--gene_spans", required=True, help="gene spans from prep_site_table.py")
    parser.add_argument("--bam", required=True, help="aligned reads with CB cell-barcode tags (indexed)")
    parser.add_argument("--cell_clusters", required=True,
                        help="cell_barcode <tab> cluster, one per line; a header line is skipped; "
                             "reads from other cells are ignored")
    parser.add_argument("--CPU", type=int, default=8)
    parser.add_argument("--output_prefix", required=True)
    args = parser.parse_args()

    barcodes, cell_cluster = parse_cell_clusters(args.cell_clusters)
    sites = list(csv.DictReader(open(args.sites), delimiter="\t"))
    spans = list(csv.DictReader(open(args.gene_spans), delimiter="\t"))
    logger.info("%d cells, %d sites, %d gene spans", len(barcodes), len(sites), len(spans))

    G["bam"] = args.bam
    G["cb_col"] = {cb: i for i, cb in enumerate(barcodes)}
    G["site_index"] = index_sites(sites)
    G["span_index"] = index_spans(spans)
    G["cell_cluster"] = [cell_cluster[cb] for cb in barcodes]

    with pysam.AlignmentFile(args.bam) as bam:
        contigs = sorted(bam.references, key=lambda c: -bam.get_reference_length(c))

    rows, cols, vals = [], [], []
    nosite = collections.Counter()
    summary = collections.Counter()
    with Pool(args.CPU) as pool:
        for contig, r, c, v, ns, summ in pool.imap_unordered(count_contig, contigs):
            rows.append(r)
            cols.append(c)
            vals.append(v)
            nosite.update(ns)
            summary.update(summ)
            if summ["reads_used"]:
                logger.info("%s: %d reads used", contig, summ["reads_used"])

    rows, cols, vals = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)
    write_mtx(args.output_prefix + ".site_counts.mtx.gz", rows, cols, vals, len(sites), len(barcodes))

    with gzip.open(args.output_prefix + ".barcodes.tsv.gz", "wt") as ofh:
        for cb in barcodes:
            print(f"{cb}\t{cell_cluster[cb]}", file=ofh)
    with gzip.open(args.output_prefix + ".sites.tsv.gz", "wt") as ofh:
        for s in sites:
            print(s["site_id"], file=ofh)

    with open(args.output_prefix + ".ends_at_no_site.tsv", "wt") as ofh:
        print("kind\tgene_key\tcluster\tends_at_no_site", file=ofh)
        for (kind, gk, cl), n in sorted(nosite.items()):
            print(f"{kind}\t{gk}\t{cl}\t{n}", file=ofh)

    with open(args.output_prefix + ".summary.tsv", "wt") as ofh:
        print("count\tvalue", file=ofh)
        for k, v in sorted(summary.items()):
            print(f"{k}\t{v}", file=ofh)
            logger.info("%s: %d", k, v)


def parse_cell_clusters(filename):
    barcodes, cell_cluster = [], {}
    for line in open(filename):
        f = line.rstrip("\n").split("\t")
        if len(f) < 2 or not f[1].strip().lstrip("-").isdigit():
            continue  # header
        cb = f[0]
        if cb not in cell_cluster:
            barcodes.append(cb)
        cell_cluster[cb] = "Cluster_" + f[1].strip()
    return barcodes, cell_cluster


def index_sites(sites):
    """(kind, contig, strand) -> sorted span starts, span ends, windows, row numbers"""
    idx = collections.defaultdict(list)
    for i, s in enumerate(sites):
        idx[(s["kind"], s["chrom"], s["strand"])].append((int(s["span_lo"]), int(s["span_hi"]), int(s["window"]), i))
    out = {}
    for k, v in idx.items():
        v.sort()
        out[k] = tuple(np.array(x) for x in zip(*v))
    return out


def index_spans(spans):
    idx = collections.defaultdict(list)
    for s in spans:
        idx[(s["chrom"], s["strand"])].append((int(s["start"]), int(s["end"]), s["gene_key"]))
    out = {}
    for k, v in idx.items():
        v.sort()
        out[k] = ([x[0] for x in v], v, max(x[1] - x[0] for x in v))
    return out


def nearest_site(index, pos):
    """row of the site whose span is nearest pos, if pos is within that site's window"""
    lo, hi, win, row = index
    i = bisect.bisect_right(lo, pos)
    best, best_d = None, None
    # the candidates: the last span starting at or before pos, and the first starting after it
    for j in (i - 1, i):
        if 0 <= j < len(lo):
            d = 0 if lo[j] <= pos <= hi[j] else min(abs(pos - lo[j]), abs(pos - hi[j]))
            if d <= win[j] and (best_d is None or d < best_d):
                best, best_d = row[j], d
    return best


def covering_gene(span_index, contig, strand, pos):
    """the one gene_key whose span covers pos, or None when none or several do"""
    if (contig, strand) not in span_index:
        return None
    starts, v, longest = span_index[(contig, strand)]
    j = bisect.bisect_right(starts, pos) - 1
    hit = None
    while j >= 0 and v[j][0] >= pos - longest:
        if v[j][1] >= pos:
            if hit is not None and hit != v[j][2]:
                return None
            hit = v[j][2]
        j -= 1
    return hit


def count_contig(contig):
    cb_col, site_index, span_index, cell_cluster = G["cb_col"], G["site_index"], G["span_index"], G["cell_cluster"]
    counts = collections.Counter()
    nosite = collections.Counter()
    summ = collections.Counter()
    empty = None

    with pysam.AlignmentFile(G["bam"]) as bam:
        for read in bam.fetch(contig):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or not read.has_tag("CB"):
                continue
            col = cb_col.get(read.get_tag("CB"))
            if col is None:
                continue
            strand = "-" if read.is_reverse else "+"
            if read.has_tag("ts") and read.get_tag("ts") == "-":
                strand = "+" if strand == "-" else "-"
            summ["reads_used"] += 1
            if strand == "+":
                ends = (("TSS", read.reference_start + 1), ("PolyA", read.reference_end))
            else:
                ends = (("TSS", read.reference_end), ("PolyA", read.reference_start + 1))
            for kind, pos in ends:
                index = site_index.get((kind, contig, strand), empty)
                row = nearest_site(index, pos) if index is not None else None
                if row is not None:
                    counts[(row, col)] += 1
                    summ[f"{kind}_ends_at_site"] += 1
                else:
                    gk = covering_gene(span_index, contig, strand, pos)
                    if gk is not None:
                        nosite[(kind, gk, cell_cluster[col])] += 1
                        summ[f"{kind}_ends_at_no_site_in_gene"] += 1
                    else:
                        summ[f"{kind}_ends_at_no_site_no_gene"] += 1

    keys = list(counts.keys())
    r = np.fromiter((k[0] for k in keys), dtype=np.int64, count=len(keys))
    c = np.fromiter((k[1] for k in keys), dtype=np.int64, count=len(keys))
    v = np.fromiter((counts[k] for k in keys), dtype=np.int64, count=len(keys))
    return contig, r, c, v, nosite, summ


def write_mtx(filename, rows, cols, vals, n_rows, n_cols):
    order = np.lexsort((rows, cols))
    with gzip.open(filename, "wt", compresslevel=4) as ofh:
        print("%%MatrixMarket matrix coordinate integer general", file=ofh)
        print(f"{n_rows} {n_cols} {len(vals)}", file=ofh)
        for r, c, v in zip(rows[order] + 1, cols[order] + 1, vals[order]):
            ofh.write(f"{r} {c} {v}\n")


if __name__ == "__main__":
    main()
