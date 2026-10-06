#!/usr/bin/env python3

"""Turn site_usage_dexseq.R's pairwise contrasts into site-switching events and
annotate them.

An event is a stable gene x cluster pair with at least one site at pairwise
padj < --fdr and |delta usage| >= --min_delta, and at least --min_gene_reads of
the gene's read ends (at its tested sites) in each of the two clusters: with a
handful of reads in a cluster, usage jumps to 0 or 1 and the pairwise test,
running on dispersions fitted across all clusters, can still call it. It is oriented so the site with
the largest significant |delta| gains usage from cluster_A to cluster_B (the
gained site); the lost site is the one whose usage falls the most.

Each event is annotated with:
  event_type   PolyA: tandem_3UTR (both sites in the last exon of one isoform),
               intronic_PolyA (the upstream site lies in an intron of an isoform
               ending at the downstream one), alt_last_exon (otherwise).
               TSS: tandem_TSS (both in the first exon of one isoform),
               alt_first_exon (otherwise).
  FSM support  the most unique FSM reads among the isoforms carrying each site,
               summed over the cluster quantifications.
  switch_class from the two sites' read ends per million site-ending reads in
               each cluster (+1), 1.5-fold: reciprocal (gained site up, lost
               site down), concordant (both up or both down), one site changes,
               neither. Direction-independent, as for the isoform-level classes.
  flags        monoexonic (a site carried only by single-exon isoforms),
               downstream_TSS_no_FSM (the gained TSS lies downstream of the lost
               one and no isoform starting there has --min_FSM unique FSM reads:
               5'-truncated reads look like this), close_sites (separation <
               --min_separation), A_rich_downstream (PolyA: >= 12 A of the 20
               genomic bases past the gained site, an oligo-dT priming template).
  high_confidence  reciprocal, >= --min_FSM at both sites, no flags.
  isoform-level DTU on the same gene x cluster pair (alt-termini rows, where
               the two isoforms share a splice pattern), when --isoform_DTU is given.
"""

import argparse
import collections
import csv
import io
import logging
import re
import tarfile

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s", level=logging.INFO)

FOLD = np.log2(1.5)
A_RICH = 12


def main():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dexseq_prefix", required=True, help="--output_prefix given to site_usage_dexseq.R")
    parser.add_argument("--kind", required=True, choices=["TSS", "PolyA"])
    parser.add_argument("--sites", required=True, help="site table from prep_site_table.py")
    parser.add_argument("--gtf", required=True, help="LRAA gtf (SYMBOL^t:... transcript ids), for exon structures")
    parser.add_argument("--cluster_quant_tar", required=True,
                        help="tar.gz of per-cluster LRAA quant.expr files (uniq_FSM_reads)")
    parser.add_argument("--genome_fa", default=None, help="genome fasta, for the A-content past PolyA sites")
    parser.add_argument("--isoform_DTU", default=None, help="isoform-level DTU tsv (sc_pseudobulk_test_isoform_DiffUsage.py)")
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--min_delta", type=float, default=0.2)
    parser.add_argument("--min_gene_reads", type=int, default=20,
                        help="gene read ends (at tested sites) required in each cluster of an event")
    parser.add_argument("--min_FSM", type=int, default=5)
    parser.add_argument("--min_separation", type=int, default=30)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    kind = args.kind
    pref = f"{args.dexseq_prefix}.{kind}"
    pw = pd.read_csv(f"{pref}.pairwise.tsv.gz", sep="\t")
    genes = pd.read_csv(f"{pref}.genes.tsv", sep="\t")
    site_conf = pd.read_csv(f"{pref}.sites.tsv", sep="\t").set_index("site_id")["n_seeds_confirmed"]
    usage = pd.read_csv(f"{pref}.cluster_usage.tsv.gz", sep="\t")
    sites = pd.read_csv(args.sites, sep="\t", low_memory=False, keep_default_na=False).set_index("site_id")

    events = call_events(pw, args.fdr, args.min_delta, args.min_gene_reads)
    logger.info("%s: %d events in %d genes", kind, len(events), events.gene_key.nunique())
    if events.empty:
        events.to_csv(args.output, sep="\t", index=False)
        return

    events = events.merge(genes[["gene_key", "gene_symbol", "n_seeds_significant", "median_q"]], on="gene_key")

    fsm = parse_FSM(args.cluster_quant_tar)
    want_symbols = set(events.gene_symbol)
    exons = parse_exons(args.gtf, want_symbols)

    for side in ("gained", "lost"):
        sid = events[f"{side}_site"]
        s = sites.loc[sid]
        events[f"{side}_pos"] = s.pos.astype(int).values
        events[f"{side}_source"] = s.source.values
        events[f"{side}_n_merged"] = s.n_merged.values
        tids = [t.split(",") if t else [] for t in s.transcript_ids]
        events[f"{side}_n_isoforms"] = [len(t) for t in tids]
        events[f"{side}_max_isoform_uniq_FSM"] = [max([fsm.get(x, 0) for x in t], default=0) for t in tids]
        events[f"{side}_monoexonic_only"] = [bool(t) and all(len(exons.get(x, [[0, 0], [0, 0]])) == 1 for x in t)
                                             for t in tids]
        events[f"{side}_n_seeds_confirmed"] = site_conf.reindex(sid).values
        if kind == "PolyA":
            events[f"{side}_PAS"] = s.pas.values

    strand = events.gene_key.str.split("|").str[2]
    plus = strand == "+"
    events["separation"] = (events.gained_pos - events.lost_pos).abs()
    events["gained_is_downstream"] = np.where(plus, events.gained_pos > events.lost_pos,
                                              events.gained_pos < events.lost_pos)

    events["event_type"] = [event_type(kind, r, sites, exons) for r in events.itertuples()]

    if kind == "PolyA" and args.genome_fa:
        import pysam
        genome = pysam.FastaFile(args.genome_fa)
        chrom = events.gene_key.str.split("|").str[1]
        events["gained_downstream_A_of_20"] = [downstream_A(genome, c, p, s == "+")
                                               for c, p, s in zip(chrom, events.gained_pos, strand)]
        events["lost_downstream_A_of_20"] = [downstream_A(genome, c, p, s == "+")
                                             for c, p, s in zip(chrom, events.lost_pos, strand)]

    add_switch_class(events, usage)

    flags = []
    for r in events.to_dict("records"):
        f = []
        if r["gained_monoexonic_only"] or r["lost_monoexonic_only"]:
            f.append("monoexonic")
        if kind == "TSS" and r["gained_is_downstream"] and r["gained_max_isoform_uniq_FSM"] < args.min_FSM:
            f.append("downstream_TSS_no_FSM")
        if r["separation"] < args.min_separation:
            f.append("close_sites")
        if kind == "PolyA" and r.get("gained_downstream_A_of_20", 0) >= A_RICH:
            f.append("A_rich_downstream")
        flags.append(",".join(f))
    events["flags"] = flags

    events["both_sites_FSM"] = (events.gained_max_isoform_uniq_FSM >= args.min_FSM) & \
                               (events.lost_max_isoform_uniq_FSM >= args.min_FSM)
    events["high_confidence"] = (events.switch_class == "reciprocal") & events.both_sites_FSM & (events["flags"] == "")

    if args.isoform_DTU:
        add_isoform_DTU(events, args.isoform_DTU, kind)

    events = events.sort_values(["high_confidence", "abs_delta"], ascending=[False, False])
    events.to_csv(args.output, sep="\t", index=False)
    logger.info("%s: %d high-confidence events in %d genes", kind, events.high_confidence.sum(),
                events.loc[events.high_confidence, "gene_key"].nunique())


def call_events(pw, fdr, min_delta, min_gene_reads):
    pw = pw[(pw.gene_reads_A >= min_gene_reads) & (pw.gene_reads_B >= min_gene_reads)].copy()
    pw["sig"] = (pw.padj < fdr) & (pw.delta_usage.abs() >= min_delta)
    rows = []
    for (gk, ca, cb), g in pw.groupby(["gene_key", "cluster_A", "cluster_B"], sort=False):
        sig = g[g.sig]
        if sig.empty:
            continue
        top = sig.loc[sig.delta_usage.abs().idxmax()]
        flip = top.delta_usage < 0
        d = -g.delta_usage if flip else g.delta_usage
        gained = g.loc[d.idxmax()]
        lost = g.loc[d.idxmin()]
        A, B = (cb, ca) if flip else (ca, cb)
        sfx_A, sfx_B = ("_B", "_A") if flip else ("_A", "_B")
        rows.append({
            "gene_key": gk, "cluster_A": A, "cluster_B": B, "n_sites": len(g), "n_sites_sig": len(sig),
            "gained_site": gained.site_id, "lost_site": lost.site_id,
            "abs_delta": abs(top.delta_usage),
            "gained_usage_A": gained["usage" + sfx_A], "gained_usage_B": gained["usage" + sfx_B],
            "lost_usage_A": lost["usage" + sfx_A], "lost_usage_B": lost["usage" + sfx_B],
            "gained_reads_A": gained["reads" + sfx_A], "gained_reads_B": gained["reads" + sfx_B],
            "lost_reads_A": lost["reads" + sfx_A], "lost_reads_B": lost["reads" + sfx_B],
            "gene_reads_A": gained["gene_reads" + sfx_A], "gene_reads_B": gained["gene_reads" + sfx_B],
            "gained_padj": gained.padj, "lost_padj": lost.padj,
        })
    return pd.DataFrame(rows)


def parse_FSM(tar):
    fsm = collections.Counter()
    with tarfile.open(tar) as tf:
        for mem in tf.getmembers():
            if mem.name.endswith("quant.expr"):
                q = pd.read_csv(io.BytesIO(tf.extractfile(mem).read()), sep="\t", comment="#",
                                usecols=["transcript_id", "uniq_FSM_reads"])
                fsm.update(dict(zip(q.transcript_id, q.uniq_FSM_reads)))
    return fsm


def parse_exons(gtf, want_symbols):
    """transcript_id (without the symbol prefix) -> merged sorted exons, for the wanted genes"""
    exons = collections.defaultdict(list)
    tid_re = re.compile(r'transcript_id "([^"]+)"')
    for line in open(gtf):
        f = line.split("\t", 9)
        if len(f) < 9 or f[2] != "exon":
            continue
        tid = tid_re.search(f[8]).group(1)
        sym, _, bare = tid.partition("^")
        if not bare or sym not in want_symbols:
            continue
        exons[bare].append((int(f[3]), int(f[4])))
    merged = {}
    for tid, segs in exons.items():
        segs.sort()
        m = [list(segs[0])]
        for s, e in segs[1:]:
            if s <= m[-1][1] + 1:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        merged[tid] = m
    return merged


def event_type(kind, r, sites, exons):
    plus = r.gene_key.split("|")[2] == "+"
    tids = set()
    for sid in (r.gained_site, r.lost_site):
        t = sites.at[sid, "transcript_ids"]
        tids.update(t.split(",") if t else [])
    lo, hi = sorted((r.gained_pos, r.lost_pos))
    terminal_is_last = (kind == "PolyA") == plus  # the varying terminus is at the high coordinate end
    for t in tids:
        ex = exons.get(t)
        if not ex:
            continue
        term = ex[-1] if terminal_is_last else ex[0]
        if term[0] - 25 <= lo and hi <= term[1] + 25:
            return "tandem_3UTR" if kind == "PolyA" else "tandem_TSS"
    if kind == "PolyA":
        upstream = lo if plus else hi
        downstream = hi if plus else lo
        for t in tids:
            ex = exons.get(t)
            if not ex or len(ex) < 2:
                continue
            end = ex[-1][1] if plus else ex[0][0]
            if abs(end - downstream) > 25:
                continue
            if any(a[1] < upstream < b[0] for a, b in zip(ex, ex[1:])):
                return "intronic_PolyA"
        return "alt_last_exon"
    return "alt_first_exon"


def downstream_A(genome, contig, pos, plus):
    if plus:
        return genome.fetch(contig, pos, pos + 20).upper().count("A")
    return genome.fetch(contig, max(0, pos - 21), pos - 1).upper().count("T")


def add_switch_class(events, usage):
    lib = usage.groupby("cluster").reads.sum()
    reads = usage.set_index(["site_id", "cluster"]).reads

    def log2fc(site, a, b):
        cpm = lambda c: reads.get((site, c), 0) / lib[c] * 1e6
        return np.log2((cpm(b) + 1) / (cpm(a) + 1))

    events["gained_log2FC"] = [log2fc(s, a, b) for s, a, b in zip(events.gained_site, events.cluster_A, events.cluster_B)]
    events["lost_log2FC"] = [log2fc(s, a, b) for s, a, b in zip(events.lost_site, events.cluster_A, events.cluster_B)]
    g, l = events.gained_log2FC, events.lost_log2FC
    events["switch_class"] = np.select(
        [(g >= FOLD) & (l <= -FOLD),
         ((g >= FOLD) & (l >= FOLD)) | ((g <= -FOLD) & (l <= -FOLD)),
         (g.abs() >= FOLD) | (l.abs() >= FOLD)],
        ["reciprocal", "concordant", "one site changes"], "neither")


def add_isoform_DTU(events, filename, kind):
    d = pd.read_csv(filename, sep="\t", low_memory=False)
    d = d[d.significant.astype(str) == "True"]
    alt_termini = d[d.dominant_splice_hashcodes == d.alternate_splice_hashcodes]
    pairs_any = {(s, frozenset((a, b))) for s, a, b in zip(d.gene_symbol, d.cluster_A, d.cluster_B)}
    pairs_at = {(s, frozenset((a, b))) for s, a, b in zip(alt_termini.gene_symbol, alt_termini.cluster_A, alt_termini.cluster_B)}
    key = [(s, frozenset((a, b))) for s, a, b in zip(events.gene_symbol, events.cluster_A, events.cluster_B)]
    events["isoform_alt_termini_DTU_same_pair"] = [k in pairs_at for k in key]
    events["isoform_any_DTU_same_pair"] = [k in pairs_any for k in key]
    events["isoform_alt_termini_DTU_gene"] = events.gene_symbol.isin(set(alt_termini.gene_symbol))


if __name__ == "__main__":
    main()
