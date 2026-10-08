#!/usr/bin/env python3

"""Pick the showcase site-switch events, for build_site_event_read_tracks.py, ahead of
knitting the site-usage notebook.

The notebook reads its showcase events from this file. Candidates: high-confidence
events of one splicing group (alternative terminal usage, or alternative splicing =
terminal-exon or internal), both clusters of at least --min_cluster_cells cells and at
least --min_gene_reads of the gene's read ends in each.

Only dominant switches are showcased, judged on the read ends: the gained site must be
the gene's most-used site of its kind in cluster_B, and the lost site the most-used in
cluster_A (<prefix>.<KIND>.cluster_usage.tsv.gz). A site share can also shift while
another site leads in both clusters; such events are real but not showcased.

Dominance is judged on sites, not on isoform quantifications: isoforms sharing their
introns and differing only at a terminus all fit the same reads, so the quantification
spreads reads among them whatever the reads' ends -- EMP3's upstream-TSS isoform stays
the top isoform in T cells, where only 7% of the gene's reads start at its TSS.

Per gene the best dominant event (reciprocal first, then highest |delta usage|; ties by
cluster pair), genes ranked by |delta usage|; the top --n_terminal_usage /
--n_alt_splicing per site kind. For each site the isoform to draw: of those carrying it
(site table transcript_ids; the gene's own, i.e. gtf ids prefixed with its symbol, in
any LRAA component) with >= --min_uniq_FSM unique full-splice-match reads (all of them
if none has that many), the one with the most reads assigned in the cluster favouring
the site.

Output tsv: tag (<gene>.<kind>.<terminal_usage|alt_splicing>), gene_symbol, kind,
gained_site, lost_site, cluster_A, cluster_B, abs_delta, splicing_class, gained_tx,
lost_tx (drawn by build_site_event_read_tracks.py), and the two sites' shares of the
gene's read ends in their clusters (gained_usage_B, lost_usage_A).
"""

import argparse
import collections
import csv
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from build_site_event_read_tracks import parse_FSM  # noqa: E402

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
    parser.add_argument("--sites", required=True, help="site table from prep_site_table.py (transcript_ids per site)")
    parser.add_argument("--gtf", required=True, help="LRAA gtf with SYMBOL^ transcript ids")
    parser.add_argument("--cluster_quant_tar", required=True, help="tar.gz of per-cluster quant.expr")
    parser.add_argument("--min_uniq_FSM", type=float, default=3,
                        help="unique FSM reads an isoform must hold to be drawn (summed over clusters)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    sizes = collections.Counter()
    for line in open(args.cell_clusters):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[1].strip().lstrip("-").isdigit():
            sizes["Cluster_" + f[1].strip()] += 1

    fsm, assigned, by_cluster = parse_FSM(args.cluster_quant_tar)
    gene_isoforms = multi_exon_isoforms_by_gene(args.gtf)
    carriers = {r["site_id"]: set(r["transcript_ids"].split(","))
                for r in csv.DictReader(open(args.sites), delimiter="\t")}

    def draw_isoform(gene_key, site, cl):
        own = set(gene_isoforms.get(gene_key, []))
        tids = [t for t in carriers.get(site, ()) if t in own] or [t for t in carriers.get(site, ()) if t]
        ok = [t for t in tids if fsm.get(t, 0) >= args.min_uniq_FSM] or tids
        reads = by_cluster.get(cl.replace("Cluster_", ""), {})
        return max(ok, key=lambda t: (reads.get(t, 0), assigned.get(t, 0), fsm.get(t, 0), t)) if ok else None

    split = pd.read_csv(args.splicing, sep="\t")
    # a site is dominant in a cluster if no site of the gene is used more (ties all count)
    top_sites = collections.defaultdict(set)
    for kind in ("TSS", "PolyA"):
        cu = pd.read_csv(f"{args.dexseq_prefix}.{kind}.cluster_usage.tsv.gz", sep="\t")
        cu = cu[cu.usage >= cu.groupby(["gene_key", "cluster"]).usage.transform("max")]
        for r in cu.itertuples():
            top_sites[(kind, r.gene_key, r.cluster)].add(r.site_id)
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
            g = g.sort_values(["gene_key", "cluster_A", "cluster_B"], kind="mergesort")
            g = g.sort_values(["reciprocal", "abs_delta"], ascending=False, kind="mergesort")
            best = {}
            for r in g.itertuples():
                if r.gene_key in best:
                    continue
                if (r.gained_site in top_sites[(kind, r.gene_key, r.cluster_B)]
                        and r.lost_site in top_sites[(kind, r.gene_key, r.cluster_A)]):
                    best[r.gene_key] = r
            n = args.n_terminal_usage if tag == "terminal_usage" else args.n_alt_splicing
            ranked = sorted(best.values(), key=lambda r: (-r.abs_delta, r.gene_key))
            for r in ranked[:n]:
                # a symbol can name genes at two loci (gene_key carries the locus): keep tags unique
                stem = r.gene_symbol
                if any(x["gene_symbol"] == r.gene_symbol and x["kind"] == kind and x["tag"].endswith(f".{tag}")
                       for x in rows):
                    _, chrom, strand = r.gene_key.split("|")
                    stem = f"{r.gene_symbol}_{chrom}{'plus' if strand == '+' else 'minus'}"
                rows.append({"tag": f"{stem}.{kind}.{tag}", "gene_symbol": r.gene_symbol, "kind": kind,
                             "gained_site": r.gained_site, "lost_site": r.lost_site,
                             "cluster_A": r.cluster_A, "cluster_B": r.cluster_B, "abs_delta": round(r.abs_delta, 4),
                             "splicing_class": r.splicing_class,
                             "gained_tx": draw_isoform(r.gene_key, r.gained_site, r.cluster_B),
                             "lost_tx": draw_isoform(r.gene_key, r.lost_site, r.cluster_A),
                             "gained_usage_B": round(r.gained_usage_B, 3), "lost_usage_A": round(r.lost_usage_A, 3)})
            print(f"{kind} {tag}: {g.gene_key.nunique()} genes with candidate events, "
                  f"{len(best)} with a dominant switch, {min(n, len(best))} showcased", file=sys.stderr)

    with open(args.output, "wt") as ofh:
        w = csv.DictWriter(ofh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def multi_exon_isoforms_by_gene(gtf):
    """gene_key (SYMBOL|chrom|strand) -> the multi-exon transcripts (bare ids) carrying that
    symbol in their gtf transcript id"""
    tid_re = re.compile(r'transcript_id "([^"]+)"')
    n_exons, key = collections.Counter(), {}
    for line in open(gtf):
        f = line.split("\t", 9)
        if len(f) < 9 or f[2] != "exon":
            continue
        full = tid_re.search(f[8]).group(1)
        if "^" not in full:
            continue
        symbol, bare = full.split("^", 1)
        n_exons[bare] += 1
        key[bare] = f"{symbol}|{f[0]}|{f[6]}"
    out = collections.defaultdict(list)
    for t, k in key.items():
        if n_exons[t] > 1:
            out[k].append(t)
    return out


if __name__ == "__main__":
    main()
