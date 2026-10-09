#!/usr/bin/env python3

"""Turn site_usage_dexseq.R's pairwise contrasts on splice patterns (kind SplicePattern) or
on terminal-site combinations within a splice pattern (kind IsoformTermini), as made by
feature_counts_from_sparse_matrix.py, into switching events and annotate them.

Events are called as annotate_site_usage_events.py calls them (its call_events): a
stable group x cluster pair with a feature at pairwise padj < --fdr and |delta usage| >=
--min_delta, and >= --min_group_reads of the group's counts in each cluster; oriented so
the significant feature with the largest |delta| gains from cluster_A to cluster_B (the
gained feature); the lost feature is the other one whose usage falls the most.

Each event is annotated with:
  gene_symbol, the two features, their usages, counts and the delta.
  switch_class   from the two features' counts per million (all of the matrix's counts in
                 the cluster, --library; +1), 1.5-fold, as in annotate_site_usage_events.py:
                 reciprocal, concordant, one site changes, neither.
  FSM support    <side>_uniq_FSM: unique FSM reads summed over the cluster quantifications
                 and over the feature's member isoforms.
  SplicePattern  each pattern's exon count and main TSS / PolyA (the site with most support
                 in the splice-pattern gtf); TSS_differs / PolyA_differs (main sites more than
                 --site_tolerance apart; NA when a pattern has none at that end).
  IsoformTermini each combination's TSS and PolyA site ids, differing_end (TSS, PolyA, both),
                 the two sites' separation (nt) at each end, and for a PolyA change the
                 PolyA site evidence of the two sites as in annotate_site_usage_events.py: PAS hexamer, or a PolyASite atlas cluster within
                 --site_tolerance with >= --min_polyasite_frac of its samples, or <=
                 --max_downstream_A A's in the 20 genomic bases past the site.
  flags          polyA_site_unsupported (IsoformTermini PolyA change: a site without that
                 evidence); differing_end_unannotated (IsoformTermini: one of the two has no
                 annotated site at an end where they differ, e.g. a TSS site vs 'none').
  high_confidence  reciprocal, >= --min_FSM unique FSM reads for both features, no flags.
"""

import argparse
import logging
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from annotate_site_usage_events import (call_events, parse_FSM, load_polyasite, polyasite_frac,  # noqa: E402
                                        downstream_A, add_switch_class)

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s", level=logging.INFO)


def main():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dexseq_prefix", required=True, help="--output_prefix given to site_usage_dexseq.R")
    parser.add_argument("--kind", required=True, choices=["SplicePattern", "IsoformTermini"])
    parser.add_argument("--features", required=True, help="feature table from feature_counts_from_sparse_matrix.py")
    parser.add_argument("--library", required=True,
                        help="<prefix>.cluster_library.tsv from feature_counts_from_sparse_matrix.py (CPM library)")
    parser.add_argument("--cluster_quant_tar", required=True,
                        help="tar.gz of per-cluster LRAA quant.expr files (uniq_FSM_reads)")
    parser.add_argument("--sites", default=None, help="site table from prep_site_table.py (PAS of PolyA sites)")
    parser.add_argument("--polyasite_atlas", default=None, help="PolyASite 2.0 atlas clusters bed(.gz)")
    parser.add_argument("--genome_fa", default=None, help="genome fasta, for the A-content past PolyA sites")
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--min_delta", type=float, default=0.2)
    parser.add_argument("--min_group_reads", type=int, default=20)
    parser.add_argument("--min_FSM", type=int, default=5)
    parser.add_argument("--min_polyasite_frac", type=float, default=0.10)
    parser.add_argument("--max_downstream_A", type=int, default=7)
    parser.add_argument("--site_tolerance", type=int, default=25)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    kind = args.kind
    pref = f"{args.dexseq_prefix}.{kind}"
    pw = pd.read_csv(f"{pref}.pairwise.tsv.gz", sep="\t")
    genes = pd.read_csv(f"{pref}.genes.tsv", sep="\t")
    conf = pd.read_csv(f"{pref}.sites.tsv", sep="\t").set_index("site_id")["n_seeds_confirmed"]
    usage = pd.read_csv(f"{pref}.cluster_usage.tsv.gz", sep="\t")
    feat = pd.read_csv(args.features, sep="\t", keep_default_na=False, low_memory=False).set_index("site_id")

    events = call_events(pw, args.fdr, args.min_delta, args.min_group_reads)
    if events.empty:
        logger.info("%s: no events", kind)
        pd.DataFrame(columns=["gene_key", "gene_symbol", "cluster_A", "cluster_B", "gained_feature", "lost_feature",
                              "abs_delta", "high_confidence", "switch_class", "flags"]).to_csv(args.output, sep="\t", index=False)
        return
    events = events.merge(genes[["gene_key", "n_seeds_significant", "median_q"]], on="gene_key")
    add_switch_class(events, usage, args.library)

    fsm = parse_FSM(args.cluster_quant_tar)
    events["gene_symbol"] = feat.loc[events.gained_site, "gene_symbol"].values
    for side in ("gained", "lost"):
        f = feat.loc[events[f"{side}_site"]]
        tids = [t.split(",") if t else [] for t in f.transcript_ids]
        events[f"{side}_n_isoforms"] = [len(t) for t in tids]
        events[f"{side}_uniq_FSM"] = [int(sum(fsm.get(x, 0) for x in t)) for t in tids]
        events[f"{side}_n_seeds_confirmed"] = conf.reindex(events[f"{side}_site"]).values
        events[f"{side}_transcript_ids"] = f.transcript_ids.values
        if kind == "SplicePattern":
            events[f"{side}_num_exons"] = f.num_exons.values
            events[f"{side}_main_TSS"] = f.main_TSS.values
            events[f"{side}_main_PolyA"] = f.main_PolyA.values
        else:
            events[f"{side}_TSS_site"] = f.TSS_site.values
            events[f"{side}_PolyA_site"] = f.PolyA_site.values

    flags = [[] for _ in range(len(events))]
    if kind == "SplicePattern":
        for end in ("TSS", "PolyA"):
            events[f"{end}_differs"] = [differs(a, b, args.site_tolerance)
                                        for a, b in zip(events[f"gained_main_{end}"], events[f"lost_main_{end}"])]
        events["exon_count_differs"] = events.gained_num_exons != events.lost_num_exons
    else:
        tss_d = events.gained_TSS_site != events.lost_TSS_site
        pa_d = events.gained_PolyA_site != events.lost_PolyA_site
        events["differing_end"] = np.select([tss_d & pa_d, tss_d, pa_d], ["both", "TSS", "PolyA"], "none")
        for end in ("TSS", "PolyA"):
            events[f"{end}_separation"] = [separation(a, b) for a, b in zip(events[f"gained_{end}_site"],
                                                                            events[f"lost_{end}_site"])]
        for i, r in enumerate(events.itertuples()):
            for end, d in (("TSS", tss_d.iloc[i]), ("PolyA", pa_d.iloc[i])):
                if d and "none" in (getattr(r, f"gained_{end}_site"), getattr(r, f"lost_{end}_site")):
                    flags[i].append("differing_end_unannotated")
                    break
        add_polyA_evidence(events, pa_d, flags, args)

    events["flags"] = [",".join(f) for f in flags]
    events["both_FSM"] = (events.gained_uniq_FSM >= args.min_FSM) & (events.lost_uniq_FSM >= args.min_FSM)
    events["high_confidence"] = (events.switch_class == "reciprocal") & events.both_FSM & (events["flags"] == "")

    events = events.rename(columns={"gained_site": "gained_feature", "lost_site": "lost_feature",
                                    "n_sites": "n_features_tested", "n_sites_sig": "n_features_sig",
                                    "gene_reads_A": "group_reads_A", "gene_reads_B": "group_reads_B"})
    first = ["gene_key", "gene_symbol", "cluster_A", "cluster_B", "gained_feature", "lost_feature", "abs_delta",
             "high_confidence", "switch_class", "flags", "gained_uniq_FSM", "lost_uniq_FSM"]
    events = events[first + [c for c in events.columns if c not in first]]
    events = events.sort_values(["high_confidence", "abs_delta"], ascending=[False, False])
    events.to_csv(args.output, sep="\t", index=False)
    hc = events[events.high_confidence]
    logger.info("%s: %d events in %d groups (%d genes); %d high-confidence events in %d groups (%d genes)", kind,
                len(events), events.gene_key.nunique(), events.gene_symbol.nunique(),
                len(hc), hc.gene_key.nunique(), hc.gene_symbol.nunique())


def separation(a, b):
    """nt between two site ids (chrom:pos in the id); NA when either is 'none'"""
    if a == "none" or b == "none":
        return "NA"
    return abs(int(a.split(":")[2]) - int(b.split(":")[2]))


def differs(a, b, tol):
    if a in ("", None) or b in ("", None) or (isinstance(a, float) and np.isnan(a)) or (isinstance(b, float) and np.isnan(b)):
        return "NA"
    return str(abs(int(float(a)) - int(float(b))) > tol)


def add_polyA_evidence(events, pa_differs, flags, args):
    """for events whose PolyA site differs: each annotated PolyA site's evidence, by the rule of
    annotate_site_usage_events.py (a PAS, PolyASite >= --min_polyasite_frac, or <= --max_downstream_A A's)"""
    pas = {}
    if args.sites:
        s = pd.read_csv(args.sites, sep="\t", keep_default_na=False, low_memory=False, usecols=["site_id", "pas"])
        pas = dict(zip(s.site_id, s.pas))
    else:
        logger.warning("no --sites: PolyA sites' PAS not known")
    atlas = load_polyasite(args.polyasite_atlas) if args.polyasite_atlas else None
    genome = None
    if args.genome_fa:
        import pysam
        genome = pysam.FastaFile(args.genome_fa)
    for side in ("gained", "lost"):
        ev, ok = [], []
        for d, sid in zip(pa_differs, events[f"{side}_PolyA_site"]):
            if not d or sid == "none":
                ev.append("")
                ok.append(np.nan)
                continue
            _, chrom, pos, strand = sid.split(":")
            pos = int(pos)
            p = str(pas.get(sid, ""))
            has_pas = p not in ("", "none", "NA", "nan")
            frac = polyasite_frac(atlas, chrom, strand, pos, args.site_tolerance) if atlas else None
            n_a = downstream_A(genome, chrom, pos, strand == "+") if genome else None
            ok.append(has_pas or (frac is not None and frac >= args.min_polyasite_frac)
                      or (n_a is not None and n_a <= args.max_downstream_A))
            ev.append(";".join([f"PAS={p if has_pas else 'none'}", f"PolyASite={'NA' if frac is None else round(frac, 2)}",
                                f"A={'NA' if n_a is None else n_a}"]))
        events[f"{side}_PolyA_site_evidence"] = ev
        events[f"{side}_PolyA_site_supported"] = ok
    for i, (g, l) in enumerate(zip(events.gained_PolyA_site_supported, events.lost_PolyA_site_supported)):
        if g == False or l == False:  # noqa: E712 (nan: not assessed)
            flags[i].append("polyA_site_unsupported")


if __name__ == "__main__":
    main()
