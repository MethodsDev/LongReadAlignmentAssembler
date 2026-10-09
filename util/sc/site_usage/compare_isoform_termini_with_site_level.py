#!/usr/bin/env python3

"""IsoformTermini (isoforms of one splice pattern differing in TSS or PolyA site, EM counts)
against the site-level analysis (read ends counted at each site; site_usage_dexseq.R +
annotate_site_usage_events.py + classify_site_pairs_by_splicing.py outputs in --site_dir).

Per site kind, each gene with a high-confidence IsoformTermini event at that end (TSS:
differing_end TSS or both; PolyA: PolyA or both) is put in the first class that applies:

  site_HC_same_sites    a high-confidence site-level event between the same two sites
                        (either direction, any cluster pair)
  site_event_same_sites a site-level event (not high confidence) between the same two sites
  site_event_other_sites site-level events for the gene, none between these sites
  site_stable_no_event  the gene is stable at site level but has no event
  site_not_stable       tested at site level, not stable (q < 0.05 in < 4 of 5 seeds)
  site_not_tested       not tested at site level (< 2 of its sites pass the read filters)

and the reverse: site-level high-confidence genes whose best event is alternative terminal
usage (tandem sites: the case IsoformTermini addresses) and whether IsoformTermini has a
high-confidence event at that end for the gene.

With --bam (and --gtf, --features, --cluster_quant_tar, --cell_clusters, --genome_fa, and
--HiFi / --rdna_mask_bed as given to the site-usage counting), each
gene's best IsoformTermini event is also checked against the alignments
(alt_termini_read_check.py): the share of the reads at the differing end (carrying the
intron next to it, or unspliced in that terminal exon) that end at the gained isoform's
terminus rather than the lost one's, in each cluster. read_ends_switch: that share rises
by >= 0.2 from cluster_A to cluster_B. The isoform used per feature is the one with the
most reads assigned in the cluster favouring it.

Output: <output_prefix>.<kind>.genes.tsv per kind and a summary on stdout.
"""

import argparse
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "diff_iso_usage"))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--isoform_termini_events", required=True)
    p.add_argument("--site_dir", required=True, help="site-usage analysis dir (PBMCs.dexseq.<kind>.* files)")
    p.add_argument("--prefix", default="PBMCs")
    p.add_argument("--kinds", default="TSS,PolyA")
    p.add_argument("--output_prefix", required=True)
    p.add_argument("--bam")
    p.add_argument("--gtf")
    p.add_argument("--features")
    p.add_argument("--cluster_quant_tar")
    p.add_argument("--cell_clusters")
    p.add_argument("--genome_fa")
    p.add_argument("--HiFi", action="store_true", help="LRAA's HiFi read-identity floor for the BAM read check")
    p.add_argument("--rdna_mask_bed", default=None, help="LRAA's rDNA mask bed: masked reads are not counted")
    args = p.parse_args()
    checker = ReadEndChecker(args) if args.bam else None

    it = pd.read_csv(args.isoform_termini_events, sep="\t", low_memory=False)
    split = pd.read_csv(f"{args.site_dir}/{args.prefix}.dexseq.site_pairs.splicing.tsv", sep="\t")
    for kind in args.kinds.split(","):
        dex = f"{args.site_dir}/{args.prefix}.dexseq.{kind}"
        se = pd.read_csv(f"{dex}.events.tsv", sep="\t", low_memory=False)
        sg = pd.read_csv(f"{dex}.genes.tsv", sep="\t")
        sk = split[split.kind == kind][["gene_key", "gained_site", "lost_site", "splicing_class"]] \
            .drop_duplicates(["gene_key", "gained_site", "lost_site"])
        se = se.merge(sk, how="left", on=["gene_key", "gained_site", "lost_site"])
        se["pair"] = [frozenset((a, b)) for a, b in zip(se.gained_site, se.lost_site)]

        ends = (kind, "both")
        hc = it[(it.high_confidence == True) & it.differing_end.isin(ends)].copy()  # noqa: E712
        hc["pair"] = [frozenset((a, b)) for a, b in zip(hc[f"gained_{kind}_site"], hc[f"lost_{kind}_site"])]
        rows = []
        # chromosome by chromosome: the read check holds one contig's sequence at a time
        by_gene = sorted(hc.groupby("gene_symbol"), key=lambda gx: (gx[1][f"gained_{kind}_site"].iloc[0].split(":")[1], gx[0]))
        for g, x in by_gene:
            pairs = set(x.pair)
            s = se[se.gene_symbol == g]
            genes = sg[sg.gene_symbol == g]
            if len(s[(s.high_confidence == True) & s.pair.isin(pairs)]):  # noqa: E712
                c = "site_HC_same_sites"
            elif len(s[s.pair.isin(pairs)]):
                c = "site_event_same_sites"
            elif len(s):
                c = "site_event_other_sites"
            elif genes.stable.any():
                c = "site_stable_no_event"
            elif len(genes):
                c = "site_not_stable"
            else:
                c = "site_not_tested"
            best = x.sort_values("abs_delta", ascending=False).iloc[0]
            rc = checker.check(best, kind) if checker else {}
            rows.append({**rc, "gene_symbol": g, "class": c, "n_IT_HC_events": len(x), "IT_best_abs_delta": round(best.abs_delta, 3),
                         "IT_best_clusters": f"{best.cluster_A}>{best.cluster_B}",
                         "IT_sites": f"{best[f'gained_{kind}_site']} / {best[f'lost_{kind}_site']}",
                         "site_n_events": len(s), "site_n_HC_events": int((s.high_confidence == True).sum()),  # noqa: E712
                         "site_best_abs_delta": round(s.abs_delta.max(), 3) if len(s) else None,
                         "site_gene_n_seeds_significant": genes.n_seeds_significant.max() if len(genes) else None})
        out = pd.DataFrame(rows).sort_values("gene_symbol").reset_index(drop=True)
        out.to_csv(f"{args.output_prefix}.{kind}.genes.tsv", sep="\t", index=False)
        print(f"== {kind}: {len(out)} genes with a high-confidence IsoformTermini event at the {kind}")
        if checker:
            t = out.groupby("class").agg(genes=("gene_symbol", "size"),
                                         read_ends_switch=("read_ends_switch", "sum"),
                                         median_IT_delta=("IT_best_abs_delta", "median"),
                                         median_read_end_delta=("read_end_delta", "median"))
            print(t.sort_values("genes", ascending=False).to_string())
        else:
            print(out["class"].value_counts().to_string())

        # reverse: site-level HC tandem-site genes
        shc = se[se.high_confidence == True]  # noqa: E712
        tandem = set(shc[shc.splicing_class == "alt_terminal_usage"].gene_symbol)
        it_genes = set(hc.gene_symbol)
        it_any = set(it[it.differing_end.isin(ends)].gene_symbol)
        print(f"site-level HC genes with a tandem ({kind}) event: {len(tandem)}; with an IsoformTermini HC event "
              f"at the {kind}: {len(tandem & it_genes)}; with any IsoformTermini event there: {len(tandem & it_any)}")
        alt = set(shc[shc.splicing_class.astype(str).str.startswith("alt_splicing")].gene_symbol) - tandem
        print(f"site-level HC genes whose {kind} switch comes with alternative splicing only: {len(alt)} "
              f"(outside IsoformTermini's scope: different splice patterns); of these with an IsoformTermini HC event: "
              f"{len(alt & it_genes)}")


class ReadEndChecker:
    def __init__(self, args):
        import collections
        import pysam
        from build_site_event_read_tracks import parse_FSM, parse_gtf
        from alt_termini_read_check import check_pair, make_read_filter
        self.check_pair = check_pair
        self.read_filter = make_read_filter(HiFi=args.HiFi, rdna_mask_bed=args.rdna_mask_bed)
        self.parse_gtf = parse_gtf
        self.gtf = args.gtf
        self.feats = pd.read_csv(args.features, sep="\t", low_memory=False).set_index("site_id")
        _, self.assigned, self.by_cluster = parse_FSM(args.cluster_quant_tar)
        self.cluster_of = {}
        for line in open(args.cell_clusters):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[1].strip().lstrip("-").isdigit():
                self.cluster_of[f[0]] = "Cluster_" + f[1].strip()
        self.bam = pysam.AlignmentFile(args.bam)
        self.genome = pysam.FastaFile(args.genome_fa)
        self.exon_cache = {}

    def pick(self, feature, cluster):
        tids = [t for t in str(self.feats.loc[feature, "transcript_ids"]).split(",") if t]
        in_cl = self.by_cluster.get(cluster.replace("Cluster_", ""), {})
        return max(tids, key=lambda t: (in_cl.get(t, 0), self.assigned.get(t, 0), t))

    def check(self, ev, kind):
        g_tx, l_tx = self.pick(ev.gained_feature, ev.cluster_B), self.pick(ev.lost_feature, ev.cluster_A)
        need = {g_tx, l_tx} - set(self.exon_cache)
        if need:
            ex, gid = self.parse_gtf(self.gtf, need)
            for t in need:
                self.exon_cache[t] = (gid[t], ex[t])
        (gfull, gex), (lfull, lex) = self.exon_cache[g_tx], self.exon_cache[l_tx]
        res = self.check_pair({"alt_terminus": kind, "dominant_transcript_ids": gfull, "alternate_transcript_ids": lfull,
                               "cluster_A": ev.cluster_A, "cluster_B": ev.cluster_B},
                              {gfull: sorted(gex["exons"]), lfull: sorted(lex["exons"])},
                              {gfull: gex["strand"], lfull: lex["strand"]}, {gfull: gex["chrom"], lfull: lex["chrom"]},
                              self.bam, self.genome, self.cluster_of, 50, read_filter=self.read_filter)
        sa, sb = res["read_dom_share_A"], res["read_dom_share_B"]
        delta = round(sb - sa, 3) if sa != "" and sb != "" else None
        return {"read_end_gained_share_A": sa, "read_end_gained_share_B": sb, "read_end_delta": delta,
                "read_ends_A": res["reads_dom_A"] + res["reads_alt_A"], "read_ends_B": res["reads_dom_B"] + res["reads_alt_B"],
                "read_frac_elsewhere": res["read_frac_elsewhere"],
                "read_ends_switch": delta is not None and delta >= 0.2}


if __name__ == "__main__":
    main()
