#!/usr/bin/env python3

"""Pick the showcase site-switch events, for build_site_event_read_tracks.py, ahead of
knitting the site-usage notebook.

Same rule as the notebook's showcase tables (site_usage_funcs.R best_event_per_gene and
the notebook's showcase()): high-confidence events of one splicing group (alternative
terminal usage, or alternative splicing = terminal-exon or internal), both clusters of
at least --min_cluster_cells cells and at least --min_gene_reads of the gene's read ends
in each, the best event per gene (highest |delta usage|; ties by gene_key), ranked by
|delta usage|; the top --n_terminal_usage / --n_alt_splicing per site kind.

Output tsv: tag (<gene>.<kind>.<terminal_usage|alt_splicing>), gene_symbol, kind,
gained_site, lost_site, cluster_A, cluster_B, abs_delta.
"""

import argparse
import collections
import csv

import pandas as pd

GROUPS = {"terminal_usage": lambda c: c == "alt_terminal_usage",
          "alt_splicing": lambda c: str(c).startswith("alt_splicing")}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dexseq_prefix", required=True, help="<prefix> of <prefix>.<KIND>.events.tsv")
    parser.add_argument("--splicing", required=True, help="classify_site_pairs_by_splicing.py output")
    parser.add_argument("--cell_clusters", required=True, help="cell_barcode <tab> cluster (header skipped)")
    parser.add_argument("--min_cluster_cells", type=int, default=200)
    parser.add_argument("--min_gene_reads", type=int, default=50)
    parser.add_argument("--n_terminal_usage", type=int, default=6)
    parser.add_argument("--n_alt_splicing", type=int, default=15)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    sizes = collections.Counter()
    for line in open(args.cell_clusters):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[1].strip().lstrip("-").isdigit():
            sizes["Cluster_" + f[1].strip()] += 1

    split = pd.read_csv(args.splicing, sep="\t")
    rows = []
    for kind in ("TSS", "PolyA"):
        ev = pd.read_csv(f"{args.dexseq_prefix}.{kind}.events.tsv", sep="\t", keep_default_na=False)
        ev = ev.merge(split[split.kind == kind][["gene_key", "gained_site", "lost_site", "splicing_class"]],
                      on=["gene_key", "gained_site", "lost_site"], how="left")
        ev = ev[ev.high_confidence.astype(str) == "True"]
        ev = ev[(ev.cluster_A.map(sizes) >= args.min_cluster_cells) & (ev.cluster_B.map(sizes) >= args.min_cluster_cells)]
        ev = ev[(ev.gene_reads_A >= args.min_gene_reads) & (ev.gene_reads_B >= args.min_gene_reads)]
        for tag, in_group in GROUPS.items():
            g = ev[ev.splicing_class.map(in_group)].copy()
            g["reciprocal"] = g.switch_class == "reciprocal"
            best = (g.sort_values(["reciprocal", "abs_delta"], ascending=False, kind="mergesort")
                    .groupby("gene_key", sort=True).head(1)
                    .sort_values("gene_key", kind="mergesort")
                    .sort_values("abs_delta", ascending=False, kind="mergesort"))
            n = args.n_terminal_usage if tag == "terminal_usage" else args.n_alt_splicing
            for r in best.head(n).itertuples():
                rows.append({"tag": f"{r.gene_symbol}.{kind}.{tag}", "gene_symbol": r.gene_symbol, "kind": kind,
                             "gained_site": r.gained_site, "lost_site": r.lost_site,
                             "cluster_A": r.cluster_A, "cluster_B": r.cluster_B, "abs_delta": round(r.abs_delta, 4)})

    with open(args.output, "wt") as ofh:
        w = csv.DictWriter(ofh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
