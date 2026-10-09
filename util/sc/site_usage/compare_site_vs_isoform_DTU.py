#!/usr/bin/env python3

"""Compare site-level switch calls (site_usage_dexseq.R + annotate_site_usage_events.py)
with isoform-level DTU calls (sc_pseudobulk_test_isoform_DiffUsage.py, chi-square on
cluster pseudobulks) under shared rules, per gene x cluster pair x site kind, and say
why each call made by only one method is missed by the other.

Shared rules (applied to both methods' calls; "harmonized" calls):
  - the switch is between two distinct annotated LRAA sites of the kind that the site
    analysis can test (site table, competing). Site level: by construction. Isoform
    level: the two isoforms of a significant DTU row (the first transcript of each side)
    carry different sites of the kind (the models' own TSS / PolyA claims, site table
    transcript_ids); a pair differing at both ends counts for both kinds;
  - PolyA: both sites have evidence of being real polyadenylation sites (a PAS hexamer;
    a PolyASite 2.0 cluster within --site_tolerance used by >= --min_polyasite_frac of
    its samples; or <= --max_downstream_A A's in the 20 genomic bases past it). Site
    level: no polyA_site_unsupported flag (same rule, annotate_site_usage_events.py);
  - effect size >= --min_delta: site level |delta usage| (by construction); isoform
    level |delta pi|;
  - depth: site level >= 20 gene read ends per cluster, isoform level >= 25 gene reads
    per cluster (each method's own floor; close enough not to matter).

Each harmonized call made by only one method gets the reason the other misses it:

  site-level only -- isoform level:
    no DTU row        the chi-square never tested the gene x pair; with --isoform_DTU_skip_reasons
                      (gene x pair -> skip_reason, from sc_pseudobulk_test_isoform_DiffUsage.py
                      --save_annotated_isoform_details) the gate that stopped it: fewer than 25
                      gene reads in a cluster, a single isoform, its top isoforms' reciprocal
                      delta pi < 0.1, an isoform < 25 reads, isoforms in < 5% of a cluster's
                      cells, or the gene absent from its table
    not significant   tested, adjusted p > its FDR (0.001)
    small shift       significant at two sites of this kind but |delta pi| < --min_delta
    site filtered     significant at two sites of this kind, but a site isn't testable
                      or fails the PolyA evidence rule
    other isoforms    significant, but the isoform pair tested doesn't differ at two
                      annotated sites of this kind (its top isoforms differ elsewhere)
  isoform-level only -- site level:
    event filtered    an event exists but fails the PolyA evidence rule
    gene not tested   fewer than 2 of the gene's sites pass the usage filter
    sites not tested  the gene is tested, but not both of these sites
    gene not stable   gene q >= 0.05 in more than 1 of the 5 pseudo-replicate dealings
    not significant   pairwise padj >= --fdr at every site of the gene for this pair
    small shift       significant, but |delta usage| < --min_delta
    low depth         < --min_gene_reads gene read ends in a cluster

Output tsv, one row per gene x cluster pair x kind called by either method (raw calls
included, with their harmonized status): gene, gene_key, cluster_A, cluster_B, kind,
site_call, site_harmonized, site_high_confidence, site_group, site_sites, isoform_call,
isoform_harmonized, isoform_group, isoform_sites, membership (of the harmonized calls;
empty when neither call is harmonized), same_site_pair, reason.
"""

import argparse
import collections
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from annotate_site_usage_events import load_polyasite, polyasite_frac, downstream_A  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dexseq_prefix", required=True,
                   help="<prefix> of <prefix>.<KIND>.{events,genes}.tsv, .pairwise.tsv.gz, .seed1.sites.tsv.gz")
    p.add_argument("--splicing", required=True, help="classify_site_pairs_by_splicing.py output")
    p.add_argument("--isoform_DTU", required=True, help="sc_pseudobulk_test_isoform_DiffUsage.py diff_iso.tsv")
    p.add_argument("--sites", required=True, help="site table from prep_site_table.py")
    p.add_argument("--isoform_DTU_skip_reasons", default=None,
                   help="tsv gene_symbol, cluster_A, cluster_B, skip_reason: why the chi-square skipped a gene x pair")
    p.add_argument("--polyasite_atlas", default=None, help="PolyASite 2.0 atlas clusters bed(.gz)")
    p.add_argument("--genome_fa", default=None, help="genome fasta, for the A count past PolyA sites")
    p.add_argument("--min_polyasite_frac", type=float, default=0.10)
    p.add_argument("--max_downstream_A", type=int, default=7)
    p.add_argument("--site_tolerance", type=int, default=25)
    p.add_argument("--min_delta", type=float, default=0.2)
    p.add_argument("--fdr", type=float, default=0.05, help="site-level pairwise padj")
    p.add_argument("--min_gene_reads", type=int, default=20, help="site-level gene read ends per cluster")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    sites = pd.read_csv(args.sites, sep="\t", low_memory=False, keep_default_na=False).set_index("site_id")
    tx_site = collections.defaultdict(dict)  # transcript -> kind -> site_id
    for sid, tids, kind in zip(sites.index, sites.transcript_ids, sites.kind):
        for t in str(tids).split(","):
            if t:
                tx_site[t][kind] = sid
    competing = set(sites.index[sites.competing.astype(str) == "True"])

    atlas = load_polyasite(args.polyasite_atlas) if args.polyasite_atlas else None
    genome = None
    if args.genome_fa:
        import pysam
        genome = pysam.FastaFile(args.genome_fa)
    pa_ok = {}

    def polyA_supported(sid):
        if sid not in pa_ok:
            _, chrom, pos, strand = sid.split(":")
            pas = str(sites.at[sid, "pas"])
            ok = pas not in ("", "none", "NA", "nan")
            if not ok and atlas is not None:
                ok = polyasite_frac(atlas, chrom, strand, int(pos), args.site_tolerance) >= args.min_polyasite_frac
            if not ok and genome is not None:
                ok = downstream_A(genome, chrom, int(pos), strand == "+") <= args.max_downstream_A
            pa_ok[sid] = ok
        return pa_ok[sid]

    # ---- site level
    split = pd.read_csv(args.splicing, sep="\t")
    ctx = dict(genes_tested={}, genes_stable={}, tested_sites={}, pw={}, dtu_any=set(), dtu_sig=set(), skip=None)
    if args.isoform_DTU_skip_reasons:
        sk = pd.read_csv(args.isoform_DTU_skip_reasons, sep="\t", keep_default_na=False)
        ctx["skip"] = {(g, *sorted((a, b))): r for g, a, b, r in zip(sk.gene_symbol, sk.cluster_A, sk.cluster_B, sk.skip_reason)}
    rows = []
    for kind in ("TSS", "PolyA"):
        ev = pd.read_csv(f"{args.dexseq_prefix}.{kind}.events.tsv", sep="\t", keep_default_na=False)
        ev = ev.merge(split[split.kind == kind][["gene_key", "gained_site", "lost_site", "splicing_class"]],
                      on=["gene_key", "gained_site", "lost_site"], how="left")
        for r in ev.itertuples():
            a, b = sorted((r.cluster_A, r.cluster_B))
            rows.append(dict(gene=r.gene_symbol, gene_key=r.gene_key, cluster_A=a, cluster_B=b, kind=kind,
                             site_harmonized="polyA_site_unsupported" not in str(r.flags),
                             site_high_confidence=str(r.high_confidence) == "True",
                             site_group=site_group(r.splicing_class),
                             site_sites=frozenset((r.gained_site, r.lost_site))))
        g = pd.read_csv(f"{args.dexseq_prefix}.{kind}.genes.tsv", sep="\t")
        ctx["genes_tested"][kind] = set(g.gene_key)
        ctx["genes_stable"][kind] = set(g[g.stable.astype(str) == "True"].gene_key)
        ctx["tested_sites"][kind] = set(pd.read_csv(f"{args.dexseq_prefix}.{kind}.seed1.sites.tsv.gz", sep="\t").site_id)
        pw = pd.read_csv(f"{args.dexseq_prefix}.{kind}.pairwise.tsv.gz", sep="\t")
        pw["a"] = pw[["cluster_A", "cluster_B"]].min(axis=1)
        pw["b"] = pw[["cluster_A", "cluster_B"]].max(axis=1)
        pw["sig"] = pw.padj < args.fdr
        pw["sig_big"] = pw.sig & (pw.delta_usage.abs() >= args.min_delta)
        pw["deep"] = (pw.gene_reads_A >= args.min_gene_reads) & (pw.gene_reads_B >= args.min_gene_reads)
        agg = pw.groupby(["gene_key", "a", "b"]).agg(sig=("sig", "any"), sig_big=("sig_big", "any"), deep=("deep", "all"))
        for (gk, a, b), r in zip(agg.index, agg.itertuples(index=False)):
            ctx["pw"][(kind, gk, a, b)] = r
    site = pd.DataFrame(rows)
    # a gene x pair x kind can switch several site pairs: harmonized first, then high confidence
    site = (site.sort_values(["site_harmonized", "site_high_confidence"], ascending=False, kind="mergesort")
            .drop_duplicates(["gene_key", "cluster_A", "cluster_B", "kind"]))
    site["site_call"] = True

    # ---- isoform level
    d = pd.read_csv(args.isoform_DTU, sep="\t", low_memory=False)
    rows = []
    for r in d.itertuples():
        dom, alt = str(r.dominant_transcript_ids).split(","), str(r.alternate_transcript_ids).split(",")
        syms = {i.split("^", 1)[0] for i in dom + alt if "^" in i}
        if len(syms) != 1:
            continue  # isoforms of two genes (read-through component)
        gene = next(iter(syms))
        a, b = sorted((r.cluster_A, r.cluster_B))
        ctx["dtu_any"].add((gene, a, b))
        if str(r.significant) != "True":
            continue
        ctx["dtu_sig"].add((gene, a, b))
        t1, t2 = dom[0].split("^", 1)[-1], alt[0].split("^", 1)[-1]
        group = ("terminal usage" if r.dominant_splice_hashcodes == r.alternate_splice_hashcodes
                 else "alternative splicing")
        for kind in ("TSS", "PolyA"):
            s1, s2 = tx_site[t1].get(kind), tx_site[t2].get(kind)
            if s1 is None or s2 is None or s1 == s2:
                continue
            big = abs(r.delta_pi) >= args.min_delta
            sites_ok = s1 in competing and s2 in competing and (
                kind == "TSS" or (polyA_supported(s1) and polyA_supported(s2)))
            rows.append(dict(gene=gene, gene_key_iso=sites.at[s1, "gene_key"] or sites.at[s2, "gene_key"],
                             cluster_A=a, cluster_B=b, kind=kind, isoform_group=group,
                             isoform_harmonized=big and sites_ok, isoform_big=big, isoform_sites_ok=sites_ok,
                             isoform_sites=frozenset((s1, s2))))
    iso = pd.DataFrame(rows)
    iso = (iso.sort_values(["isoform_harmonized", "isoform_big", "isoform_sites_ok"], ascending=False, kind="mergesort")
           .drop_duplicates(["gene", "cluster_A", "cluster_B", "kind"]))
    iso["isoform_call"] = True

    m = site.merge(iso, on=["gene", "cluster_A", "cluster_B", "kind"], how="outer")
    for c in ("site_call", "isoform_call", "site_harmonized", "site_high_confidence", "isoform_harmonized",
              "isoform_big"):
        m[c] = m[c].eq(True)
    m["gene_key"] = m.gene_key.fillna(m.gene_key_iso)
    sh, ih = m.site_call & m.site_harmonized, m.isoform_call & m.isoform_harmonized
    m["membership"] = ["both" if s and i else "site only" if s else "isoform only" if i else ""
                       for s, i in zip(sh, ih)]
    m["same_site_pair"] = [(s == i) if isinstance(s, frozenset) and isinstance(i, frozenset) else ""
                           for s, i in zip(m.site_sites, m.isoform_sites)]
    m["reason"] = [reason(r, ctx) for r in m.itertuples()]
    for c in ("site_sites", "isoform_sites"):
        m[c] = [",".join(sorted(x)) if isinstance(x, frozenset) else "" for x in m[c]]
    m = m[["gene", "gene_key", "cluster_A", "cluster_B", "kind", "site_call", "site_harmonized", "site_high_confidence",
           "site_group", "site_sites", "isoform_call", "isoform_harmonized", "isoform_group", "isoform_sites",
           "membership", "same_site_pair", "reason"]]
    m.to_csv(args.output, sep="\t", index=False)
    h = m[m.membership != ""]
    print(h.groupby(["kind", "membership"]).size().to_string())
    print(h[h.membership != "both"].groupby(["kind", "membership", "reason"]).size().to_string())


def reason(r, ctx):
    if r.membership == "site only":
        key = (r.gene, r.cluster_A, r.cluster_B)
        if r.isoform_call:  # significant at two annotated sites of this kind, but not harmonized
            return ("isoform: small shift (|delta pi| < 0.2)" if not r.isoform_big
                    else "isoform: a site untestable or failing PolyA evidence")
        if key not in ctx["dtu_any"]:
            if ctx["skip"] is None:
                return "isoform: no DTU row (not tested / gates)"
            return "isoform: " + SKIP_REASONS.get(ctx["skip"].get(key, "absent"), "not tested (other gate)")
        if key not in ctx["dtu_sig"]:
            return "isoform: tested, not significant"
        return "isoform: significant, other isoforms (not these site kinds)"
    if r.membership == "isoform only":
        kind, gk = r.kind, r.gene_key
        if r.site_call:
            return "site: event fails PolyA evidence"
        if gk not in ctx["genes_tested"][kind]:
            return "site: gene not tested (< 2 sites pass filter)"
        if not all(s in ctx["tested_sites"][kind] for s in r.isoform_sites):
            return "site: these sites not tested"
        if gk not in ctx["genes_stable"][kind]:
            return "site: gene not stable across seeds"
        b = ctx["pw"].get((kind, gk, r.cluster_A, r.cluster_B))
        if b is None or not b.sig:
            return "site: pair not significant"
        if not b.sig_big:
            return "site: small shift (|delta usage| < 0.2)"
        if not b.deep:
            return "site: low depth (< 20 gene reads in a cluster)"
        return "site: other"
    return ""


SKIP_REASONS = {
    "insufficient_gene_reads": "not tested: < 25 gene reads in a cluster",
    "zero_total_counts": "not tested: gene absent from a cluster",
    "single_isoform": "not tested: a single isoform expressed",
    "delta_pi_fail": "not tested: top isoforms' reciprocal delta pi < 0.1",
    "dominant_isoform_low_reads": "not tested: an isoform with < 25 reads",
    "alternate_isoform_low_reads": "not tested: an isoform with < 25 reads",
    "cell_fraction_fail": "not tested: isoforms in < 5% of a cluster's cells",
    "absent": "not tested: gene not in its table (symbol / unspliced only)",
}


def site_group(c):
    c = str(c)
    if c == "alt_terminal_usage":
        return "terminal usage"
    if c.startswith("alt_splicing"):
        return "alternative splicing"
    return "unresolved"


if __name__ == "__main__":
    main()
