#!/usr/bin/env python3

"""Pick showcase events of a feature-usage analysis (annotate_feature_usage_events.py:
SplicePattern or IsoformTermini) and write them in the event format of
build_site_event_read_tracks.py, which draws them.

The test runs on LRAA's EM-assigned counts; a showcase must also be borne out by reads
that can be placed without the EM:

  1. candidates: high-confidence events, both clusters of >= --min_cluster_cells cells
     and >= --min_group_reads of the group's counts in each, that are dominant switches
     on the EM usage: the gained feature is the group's most-used feature in cluster_B,
     the lost feature the most-used in cluster_A. Per gene the best event (reciprocal
     first, then |delta usage|); genes ranked by |delta usage|, the top --n_candidates
     checked (plus any --anchor_genes).
  2. the isoform drawn per feature: of the feature's isoforms with >= --min_uniq_FSM
     unique full-splice-match (FSM) reads (summed over the cluster quantifications), the
     one with the most reads assigned in the cluster favouring the feature (cluster_B for
     the gained feature, cluster_A for the lost one). The two drawn isoforms must
     overlap on the genome.
  3. FSM support, from the tracking file: a read is a full-splice match of a splice
     pattern when it is FSM to an isoform of the pattern and every isoform it is
     assigned to carries that pattern. SplicePattern: the gained pattern has
     >= --min_pattern_FSM such reads in cluster_B and the lost one in cluster_A, and the
     gained pattern's share of the pair's FSM reads rises from cluster_A to cluster_B
     by >= --min_read_delta, crossing one half (the FSM reads switch too).
     IsoformTermini: both features are one splice pattern, so FSM reads cannot tell them
     apart; step 4 decides.
  4. read ends. IsoformTermini: alt_termini_read_check.py's count, at the end that
     differs, of the reads carrying the intron next to it (and unspliced reads in that
     terminal exon) ending at the gained vs the lost isoform's terminus: the gained
     terminus' share must rise by >= --min_read_delta from cluster_A to cluster_B and
     cross one half. Both: for each drawn isoform, the share of its FSM reads (in the two
     clusters) whose 5' end lies at a TSS, and whose 3' end at a PolyA site, of the two
     features (an annotated site carried by one of their isoforms, within its counting
     window, or the drawn model's own end +- --end_tolerance) is reported
     (read_TSS_agree / read_PolyA_agree);
     SplicePattern needs >= --min_end_agree at an end the two patterns' models differ
     at (where the pair's distinguishing terminus is), IsoformTermini at the differing end.

Showcases: passing candidates, anchor genes first (IsoformTermini: then the best
--min_PolyA PolyA events; SplicePattern: then the best --min_same_termini switches whose
drawn isoforms share both termini), then by |delta usage|, at most
--max_per_pair per cluster pair, --n in all. For build_site_event_read_tracks.py each
event gets a site kind and a site per side: the end that differs (TSS when both differ
or neither does; then both sides share the site and only full-splice-match reads are
drawn), the sites those of the site table whose counting windows hold the drawn models'
ends. Genes matching --exclude_genes_regex (immunoglobulin / T-cell receptor segments by
default) are never showcased.

Outputs: --output (showcase events: tag, gene_symbol, kind, gained_site, lost_site,
cluster_A, cluster_B, gained_tx, lost_tx and the checks) and --output.candidates.tsv
(every checked candidate with its checks and why it failed).
"""

import argparse
import collections
import csv
import gzip
import logging
import os
import sys

import pandas as pd
import pysam

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from build_site_event_read_tracks import parse_FSM, parse_gtf, blocks  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "diff_iso_usage"))
from alt_termini_read_check import check_pair, make_read_filter, transcript_strand  # noqa: E402

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s", level=logging.INFO)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dexseq_prefix", required=True, help="<prefix> of <prefix>.<kind>.events.tsv / .cluster_usage.tsv.gz")
    p.add_argument("--kind", required=True, choices=["SplicePattern", "IsoformTermini"])
    p.add_argument("--features", required=True, help="feature table (feature_counts_from_sparse_matrix.py)")
    p.add_argument("--sites", required=True, help="site table (prep_site_table.py)")
    p.add_argument("--gtf", required=True, help="LRAA gtf with SYMBOL^ transcript ids (whole exons)")
    p.add_argument("--cluster_quant_tar", required=True)
    p.add_argument("--tracking", required=True, help="merged per-cluster quant.tracking(.gz)")
    p.add_argument("--bam", required=True)
    p.add_argument("--genome_fa", required=True)
    p.add_argument("--cell_clusters", required=True)
    p.add_argument("--min_cluster_cells", type=int, default=200)
    p.add_argument("--min_group_reads", type=float, default=50)
    p.add_argument("--min_uniq_FSM", type=float, default=3)
    p.add_argument("--min_pattern_FSM", type=float, default=10,
                   help="SplicePattern: full-splice-match reads of the gained pattern in cluster_B and of the lost one "
                        "in cluster_A")
    p.add_argument("--min_read_delta", type=float, default=0.2)
    p.add_argument("--end_tolerance", type=int, default=50)
    p.add_argument("--min_end_agree", type=float, default=0.3)
    p.add_argument("--n_candidates", type=int, default=60)
    p.add_argument("--anchor_genes", default="", help="comma-separated genes checked and showcased first when they pass")
    p.add_argument("--n", type=int, default=10)
    p.add_argument("--HiFi", action="store_true", help="LRAA's HiFi read-identity floor for the BAM read checks")
    p.add_argument("--rdna_mask_bed", default=None, help="LRAA's rDNA mask bed: masked reads are not counted")
    p.add_argument("--max_per_pair", type=int, default=2)
    p.add_argument("--min_PolyA", type=int, default=2,
                   help="IsoformTermini: showcase at least this many PolyA events when that many pass (TSS events "
                        "otherwise fill the list on |delta usage|)")
    p.add_argument("--min_same_termini", type=int, default=2,
                   help="SplicePattern: showcase at least this many switches whose drawn isoforms share both "
                        "termini (internal splicing only) when that many pass")
    p.add_argument("--exclude_genes_regex", default=r"^(IG[HKL][VDJC]|TR[ABDG][VDJC])",
                   help="genes never showcased (default: immunoglobulin / T-cell receptor segments, whose "
                        "'splice patterns' follow V(D)J recombination)")
    p.add_argument("--output", required=True)
    args = p.parse_args()
    anchors = [g for g in args.anchor_genes.split(",") if g]

    cluster_of, n_cells = {}, collections.Counter()
    for line in open(args.cell_clusters):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[1].strip().lstrip("-").isdigit():
            cluster_of[f[0]] = "Cluster_" + f[1].strip()
            n_cells[cluster_of[f[0]]] += 1

    ev = pd.read_csv(f"{args.dexseq_prefix}.{args.kind}.events.tsv", sep="\t", low_memory=False)
    usage = pd.read_csv(f"{args.dexseq_prefix}.{args.kind}.cluster_usage.tsv.gz", sep="\t")
    feats = pd.read_csv(args.features, sep="\t", low_memory=False).set_index("site_id")

    # 1. candidates
    top = usage.sort_values("usage", ascending=False).drop_duplicates(["gene_key", "cluster"]) \
        .set_index(["gene_key", "cluster"])["site_id"]
    e = ev[(ev.high_confidence == True)  # noqa: E712
           & ev.cluster_A.map(n_cells).ge(args.min_cluster_cells) & ev.cluster_B.map(n_cells).ge(args.min_cluster_cells)
           & (ev.group_reads_A >= args.min_group_reads) & (ev.group_reads_B >= args.min_group_reads)].copy()
    e["dominant"] = [top.get((g, b)) == gf and top.get((g, a)) == lf
                     for g, a, b, gf, lf in zip(e.gene_key, e.cluster_A, e.cluster_B, e.gained_feature, e.lost_feature)]
    e = e[e.dominant & ~e.gene_symbol.astype(str).str.match(args.exclude_genes_regex)]
    e = e.assign(recip=e.switch_class.eq("reciprocal")) \
        .sort_values(["recip", "abs_delta"], ascending=False).drop_duplicates("gene_symbol") \
        .sort_values("abs_delta", ascending=False)
    cand = pd.concat([e[e.gene_symbol.isin(anchors)], e.head(args.n_candidates)]).drop_duplicates("gene_key")
    missing = [g for g in anchors if g not in set(cand.gene_symbol)]
    if missing:
        logger.warning("anchor genes without a dominant high-confidence event: %s", ", ".join(missing))
    logger.info("%d dominant high-confidence genes; %d candidates checked", len(e), len(cand))

    # 2. drawn isoforms
    fsm, assigned, by_cluster = parse_FSM(args.cluster_quant_tar)
    rows = []
    for _, r in cand.iterrows():
        d = r.to_dict()
        d["fail"] = []
        for side, cl in (("gained", r.cluster_B), ("lost", r.cluster_A)):
            tids = [t for t in str(feats.loc[r[f"{side}_feature"], "transcript_ids"]).split(",") if t]
            ok = [t for t in tids if fsm.get(t, 0) >= args.min_uniq_FSM]
            in_cl = by_cluster.get(cl.replace("Cluster_", ""), {})
            d[f"{side}_tx"] = max(ok, key=lambda t: (in_cl.get(t, 0), assigned.get(t, 0), t)) if ok else ""
            d[f"{side}_tx_uniq_FSM"] = fsm.get(d[f"{side}_tx"], 0) if ok else max([fsm.get(t, 0) for t in tids] or [0])
            if not ok:
                d["fail"].append(f"{side}_no_isoform_with_{args.min_uniq_FSM:g}_uniq_FSM")
        rows.append(d)

    want = {d[k] for d in rows for k in ("gained_tx", "lost_tx") if d[k]}
    exons, gtf_id = parse_gtf(args.gtf, want)
    for d in rows:
        if d["gained_tx"] and d["lost_tx"]:
            a, b = exons[d["gained_tx"]]["exons"], exons[d["lost_tx"]]["exons"]
            if max(x for _, x in a) < min(s for s, _ in b) or max(x for _, x in b) < min(s for s, _ in a):
                d["fail"].append("isoforms_do_not_overlap")

    # 3. FSM reads by splice pattern and cluster, and each drawn isoform's FSM read names
    comps = {t.rsplit(":", 1)[0].replace("t:", "g:", 1) for t in want}
    pattern_of = {}  # tracking hash of each drawn isoform
    read_hashes = collections.defaultdict(set)
    read_fsm_hash = {}
    read_cluster = {}
    fsm_reads_of_tx = collections.defaultdict(set)  # tx -> {(cluster, molecule)}
    opener = gzip.open if args.tracking.endswith(".gz") else open
    with opener(args.tracking, "rt") as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("gene_id"):
                continue
            g = line[:line.index("\t")]
            if g not in comps:
                continue
            f = line.rstrip("\n").split("\t")
            tx, h, name = f[1], f[2], f[5]
            read_hashes[name].add(h)
            if tx in want:
                pattern_of[tx] = h
            if f[9] == "1":
                read_fsm_hash.setdefault(name, h)
                if tx in want:
                    parts = name.split("^")
                    cl = cluster_of.get(parts[0])
                    if cl:
                        fsm_reads_of_tx[tx].add((cl, parts[-1]))
            if name not in read_cluster:
                read_cluster[name] = cluster_of.get(name.split("^")[0])
    fsm_by_pattern = collections.Counter()
    for name, h in read_fsm_hash.items():
        if read_hashes[name] == {h} and read_cluster.get(name):
            fsm_by_pattern[(h, read_cluster[name])] += 1
    logger.info("tracking: %d reads in the candidates' components", len(read_hashes))

    # annotated sites: counting windows, the ones carried by each isoform, and lookup by position
    site_rows = list(csv.DictReader(open(args.sites), delimiter="\t"))
    sites_of_tx = collections.defaultdict(list)  # tx -> [(kind, lo, hi)]
    sites_by_strand = collections.defaultdict(list)  # (kind, chrom, strand) -> [(pos, lo, hi, site_id)]
    for srow in site_rows:
        lo_w, hi_w = int(srow["span_lo"]) - int(srow["window"]), int(srow["span_hi"]) + int(srow["window"])
        sites_by_strand[(srow["kind"], srow["chrom"], srow["strand"])].append((int(srow["pos"]), lo_w, hi_w, srow["site_id"]))
        for t in srow["transcript_ids"].split(","):
            if t:
                sites_of_tx[t].append((srow["kind"], lo_w, hi_w))

    def site_at(kind, chrom, strand, pos):
        """the site whose counting window holds pos (the nearest one), else None"""
        hits = [x for x in sites_by_strand.get((kind, chrom, strand), []) if x[1] <= pos <= x[2]]
        return min(hits, key=lambda x: abs(x[0] - pos))[3] if hits else None

    def feature_windows(*features):
        w = {"TSS": [], "PolyA": []}
        for ft in features:
            for t in str(feats.loc[ft, "transcript_ids"]).split(","):
                for kind, lo_w, hi_w in sites_of_tx.get(t, []):
                    w[kind].append((lo_w, hi_w))
        return w

    bam = pysam.AlignmentFile(args.bam)
    genome = pysam.FastaFile(args.genome_fa)
    # the reads LRAA would use, on LRAA's transcribed strand (alt_termini_read_check)
    read_filter = make_read_filter(HiFi=args.HiFi, rdna_mask_bed=args.rdna_mask_bed)
    for d in rows:
        if not (d["gained_tx"] and d["lost_tx"]):
            continue
        A, B = d["cluster_A"], d["cluster_B"]
        g_tx, l_tx = d["gained_tx"], d["lost_tx"]
        # FSM reads of the two splice patterns in the two clusters
        gh, lh = pattern_of.get(g_tx), pattern_of.get(l_tx)
        for side, h in (("gained", gh), ("lost", lh)):
            for ab, cl in (("A", A), ("B", B)):
                d[f"{side}_pattern_FSM_{ab}"] = fsm_by_pattern.get((h, cl), 0)
        if args.kind == "SplicePattern":
            sh = {}
            for ab in "AB":
                n = d[f"gained_pattern_FSM_{ab}"] + d[f"lost_pattern_FSM_{ab}"]
                sh[ab] = d[f"gained_pattern_FSM_{ab}"] / n if n else None
                d[f"gained_FSM_share_{ab}"] = round(sh[ab], 3) if sh[ab] is not None else ""
            if d["gained_pattern_FSM_B"] < args.min_pattern_FSM or d["lost_pattern_FSM_A"] < args.min_pattern_FSM:
                d["fail"].append("few_pattern_FSM_reads")
            elif not (sh["A"] is not None and sh["B"] is not None and sh["B"] - sh["A"] >= args.min_read_delta
                      and sh["B"] > 0.5 > sh["A"]):
                d["fail"].append("FSM_reads_do_not_switch")

        # where each drawn isoform's FSM reads start and end: at one of the annotated sites
        # of the two features (or the drawn model's own end, +- end_tolerance)
        ex = {t: exons[t] for t in (g_tx, l_tx)}
        lo = min(min(s for s, _ in ex[t]["exons"]) for t in ex) - 200
        hi = max(max(x for _, x in ex[t]["exons"]) for t in ex) + 200
        names = {t: {m for c, m in fsm_reads_of_tx[t] if c in (A, B)} for t in ex}
        agree = {t: [0, 0, 0] for t in ex}  # n, 5' ok, 3' ok
        strand = ex[g_tx]["strand"]
        win = feature_windows(d["gained_feature"], d["lost_feature"])
        for t in ex:
            xs = ex[t]["exons"]
            tss = min(s for s, _ in xs) if strand == "+" else max(x for _, x in xs)
            pa = max(x for _, x in xs) if strand == "+" else min(s for s, _ in xs)
            win["TSS"].append((tss - args.end_tolerance, tss + args.end_tolerance))
            win["PolyA"].append((pa - args.end_tolerance, pa + args.end_tolerance))
        for read in bam.fetch(ex[g_tx]["chrom"], max(0, lo), hi):
            if not read_filter(read) or transcript_strand(read, genome) != strand:
                continue
            for t in ex:
                if read.query_name in names[t]:
                    five = read.reference_start + 1 if strand == "+" else read.reference_end
                    three = read.reference_end if strand == "+" else read.reference_start + 1
                    agree[t][0] += 1
                    agree[t][1] += any(a <= five <= b for a, b in win["TSS"])
                    agree[t][2] += any(a <= three <= b for a, b in win["PolyA"])
        for side, t in (("gained", g_tx), ("lost", l_tx)):
            n, a5, a3 = agree[t]
            d[f"{side}_FSM_reads_in_pair_clusters"] = n
            d[f"{side}_read_TSS_agree"] = round(a5 / n, 3) if n else ""
            d[f"{side}_read_PolyA_agree"] = round(a3 / n, 3) if n else ""

        # the end the event is about
        def end_of(t, which):
            xs = ex[t]["exons"]
            first = min(s for s, _ in xs) if strand == "+" else max(x for _, x in xs)
            last = max(x for _, x in xs) if strand == "+" else min(s for s, _ in xs)
            return first if which == "TSS" else last
        if args.kind == "IsoformTermini":
            de = d["differing_end"]
            ends = ["TSS", "PolyA"] if de == "both" else [de]
        else:
            ends = [w for w in ("TSS", "PolyA") if abs(end_of(g_tx, w) - end_of(l_tx, w)) > 25]
        d["checked_ends"] = ",".join(ends)
        for w in ends:
            if any(d[f"{s}_read_{w}_agree"] == "" or d[f"{s}_read_{w}_agree"] < args.min_end_agree
                   for s in ("gained", "lost")):
                d["fail"].append(f"FSM_reads_disagree_with_{w}")

        # read ends at the differing terminus (IsoformTermini)
        if args.kind == "IsoformTermini":
            w = max(ends, key=lambda w: abs(end_of(g_tx, w) - end_of(l_tx, w))) if len(ends) > 1 else ends[0]
            gfull, lfull = gtf_id[g_tx], gtf_id[l_tx]
            gex = {gfull: sorted(ex[g_tx]["exons"]), lfull: sorted(ex[l_tx]["exons"])}
            res = check_pair({"alt_terminus": w, "dominant_transcript_ids": gfull, "alternate_transcript_ids": lfull,
                              "cluster_A": A, "cluster_B": B}, gex,
                             {gfull: strand, lfull: strand}, {gfull: ex[g_tx]["chrom"], lfull: ex[g_tx]["chrom"]},
                             bam, genome, cluster_of, 50, read_filter=read_filter)
            d["read_end_kind"] = w
            for k in ("n_reads", "read_frac_at_dom", "read_frac_at_alt", "reads_dom_A", "reads_alt_A",
                      "read_dom_share_A", "reads_dom_B", "reads_alt_B", "read_dom_share_B"):
                d["readend_" + k.replace("dom", "gained").replace("alt", "lost")] = res[k]
            sa, sb = res["read_dom_share_A"], res["read_dom_share_B"]
            if sa == "" or sb == "" or not (sb - sa >= args.min_read_delta and sb > 0.5 > sa):
                d["fail"].append("read_ends_do_not_switch")

    # site kind and sites for the read tracks: the sites whose windows hold the drawn models' ends
    for d in rows:
        if d["fail"]:
            continue
        g_tx, l_tx = d["gained_tx"], d["lost_tx"]
        chrom, strand = exons[g_tx]["chrom"], exons[g_tx]["strand"]

        def model_site(t, kind):
            xs = exons[t]["exons"]
            first = min(s for s, _ in xs) if strand == "+" else max(x for _, x in xs)
            last = max(x for _, x in xs) if strand == "+" else min(s for s, _ in xs)
            return site_at(kind, chrom, strand, first if kind == "TSS" else last)

        if args.kind == "IsoformTermini":
            kind = d["read_end_kind"]
            gs, ls = d[f"gained_{kind}_site"], d[f"lost_{kind}_site"]
        else:
            ends = d["checked_ends"].split(",") if d["checked_ends"] else []
            kinds = [k for k in (ends or ["TSS", "PolyA"]) if model_site(g_tx, k) and model_site(l_tx, k)]
            if not kinds:
                d["fail"].append("drawn_isoforms_end_at_no_site")
                continue
            kind = kinds[0]
            gs, ls = model_site(g_tx, kind), model_site(l_tx, kind)
        d["track_kind"], d["gained_site"], d["lost_site"] = kind, gs, ls

    for d in rows:
        d["passed"] = not d["fail"]
        d["fail"] = ";".join(d["fail"])

    out = pd.DataFrame(rows)
    out["anchor"] = out.gene_symbol.isin(anchors)
    out.to_csv(args.output + ".candidates.tsv", sep="\t", index=False)

    ok = out[out.passed].sort_values(["anchor", "abs_delta"], ascending=False)
    if args.kind == "IsoformTermini" and len(ok):
        # anchors, then the best PolyA events, then the rest
        polyA = ok[~ok.anchor & ok.track_kind.eq("PolyA")].head(args.min_PolyA)
        ok = pd.concat([ok[ok.anchor], polyA, ok.drop(index=polyA.index)]).drop_duplicates("gene_key")
    if args.kind == "SplicePattern" and len(ok):
        same = ok[~ok.anchor & ok.checked_ends.fillna("").eq("")].head(args.min_same_termini)
        ok = pd.concat([ok[ok.anchor], same, ok.drop(index=same.index)]).drop_duplicates("gene_key")
    chosen, per_pair = [], collections.Counter()
    for _, r in ok.iterrows():
        pair = tuple(sorted((r.cluster_A, r.cluster_B)))
        if per_pair[pair] >= args.max_per_pair and not r.anchor:
            continue
        per_pair[pair] += 1
        chosen.append(r)
        if len(chosen) >= args.n:
            break
    sc = pd.DataFrame(chosen)
    sc.insert(0, "tag", [f"{g}.{args.kind}" for g in sc.gene_symbol])
    sc = sc.rename(columns={"kind": "_kind"}) if "kind" in sc.columns else sc
    sc.insert(2, "kind", sc.track_kind)
    lead = ["tag", "gene_symbol", "kind", "gained_site", "lost_site", "cluster_A", "cluster_B", "gained_tx", "lost_tx",
            "gained_feature", "lost_feature", "abs_delta", "switch_class"]
    sc = sc[lead + [c for c in sc.columns if c not in lead and c != "track_kind"]]
    sc.to_csv(args.output, sep="\t", index=False)
    logger.info("%d candidates pass; %d showcases: %s", len(ok), len(sc), ", ".join(sc.gene_symbol))


if __name__ == "__main__":
    main()
