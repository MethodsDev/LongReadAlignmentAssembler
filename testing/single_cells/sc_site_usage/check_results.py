#!/usr/bin/env python3

"""Check the site-usage test outputs against switches known from the full PBMC analysis.

Asserts what is robust at this scale (seven genes): which genes are tested and stable,
which pairwise switches are found and between which clusters, their size (to within
--tol), which site gains, and how each switch is classified (alternative terminal usage
vs a switch coupled with alternative splicing). It does not assert p-values or the
expression-based switch class: with seven genes, DEXSeq's dispersion trend and the
per-cluster site-read totals behind the switch class differ from the genome-wide run.
"""

import argparse
import sys

import pandas as pd

FAILURES = []


def check(cond, msg):
    print(("ok    " if cond else "FAIL  ") + msg)
    if not cond:
        FAILURES.append(msg)
    return cond


# (kind, gene, cluster pair, gained site, |delta| in the full analysis, splicing class, event type)
EXPECTED_EVENTS = [
    ("TSS", "SELENOH", ("Cluster_7", "Cluster_11"), "TSS:chr11:57741491:+", 0.70, "alt_terminal_usage", "tandem_TSS"),
    ("TSS", "EMP3", ("Cluster_0", "Cluster_13"), "TSS:chr19:48325357:+", 0.68, "alt_terminal_usage", "tandem_TSS"),
    ("TSS", "CIAO2A", ("Cluster_2", "Cluster_1"), "TSS:chr15:64093795:-", 0.75, "alt_terminal_usage", "tandem_TSS"),
    ("TSS", "AIF1", ("Cluster_7", "Cluster_1"), "TSS:chr6:31615242:+", 0.78, "alt_splicing:terminal_exon", "alt_first_exon"),
    ("PolyA", "ELOVL5", ("Cluster_6", "Cluster_1"), "PolyA:chr6:53267405:-", 0.45, "alt_splicing:terminal_exon",
     "intronic_PolyA"),
]

STABLE = {"TSS": ["SELENOH", "EMP3", "CIAO2A", "AIF1"], "PolyA": ["ELOVL5", "POLR2K", "CMPK1", "AIF1"]}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--prefix", default="test")
    p.add_argument("--tol", type=float, default=0.05, help="allowed difference in |delta usage|")
    a = p.parse_args()
    dex = f"{a.prefix}.dexseq"

    # per-cell counts
    summ = dict(l.rstrip("\n").split("\t") for l in open(f"{a.prefix}.summary.tsv") if not l.startswith("count"))
    check(int(summ.get("TSS_ends_at_site", 0)) > 10000 and int(summ.get("PolyA_ends_at_site", 0)) > 5000,
          f"read ends counted at sites (TSS {summ.get('TSS_ends_at_site')}, PolyA {summ.get('PolyA_ends_at_site')})")

    split = pd.read_csv(f"{dex}.site_pairs.splicing.tsv", sep="\t")
    for kind, genes in STABLE.items():
        g = pd.read_csv(f"{dex}.{kind}.genes.tsv", sep="\t")
        for gene in genes:
            r = g[g.gene_symbol == gene]
            check(len(r) == 1 and bool(r.stable.iloc[0]),
                  f"{kind} {gene}: tested and stable ({int(r.n_seeds_significant.iloc[0]) if len(r) else 0}/5 seeds)")

    events = {k: pd.read_csv(f"{dex}.{k}.events.tsv", sep="\t", keep_default_na=False) for k in ("TSS", "PolyA")}
    for kind, gene, pair, gained, delta, sclass, etype in EXPECTED_EVENTS:
        e = events[kind]
        e = e[(e.gene_symbol == gene) & (e.cluster_A == pair[0]) & (e.cluster_B == pair[1])]
        label = f"{kind} {gene} {pair[0]}->{pair[1]}"
        if not check(len(e) == 1, f"{label}: switch event found"):
            continue
        e = e.iloc[0]
        check(e.gained_site == gained, f"{label}: gained site {e.gained_site} (expected {gained})")
        check(abs(float(e.abs_delta) - delta) <= a.tol,
              f"{label}: |delta usage| {float(e.abs_delta):.2f} (full analysis {delta:.2f})")
        check(e.event_type == etype, f"{label}: event type {e.event_type} (expected {etype})")
        s = split[(split.kind == kind) & (split.gene_key == e.gene_key) & (split.gained_site == e.gained_site)
                  & (split.lost_site == e.lost_site)]
        check(len(s) == 1 and s.splicing_class.iloc[0] == sclass,
              f"{label}: splicing class {s.splicing_class.iloc[0] if len(s) else None} (expected {sclass})")

    # showcase read tracks
    m = pd.read_csv("read_tracks/manifest.tsv", sep="\t")
    check(len(m) >= 4, f"read-track data built for {len(m)} showcase events")
    for r in m.itertuples():
        reads = pd.read_csv(f"read_tracks/{r.tag}.reads.tsv", sep="\t")
        ends = pd.read_csv(f"read_tracks/{r.tag}.ends.tsv", sep="\t")
        check(reads.read_name.nunique() > 0 and ends.reads.sum() > 0,
              f"read tracks {r.tag}: {reads.read_name.nunique()} reads drawn, {ends.reads.sum()} read ends counted")

    print(f"\n{len(FAILURES)} check(s) failed" if FAILURES else "\nall checks passed")
    sys.exit(1 if FAILURES else 0)


if __name__ == "__main__":
    main()
