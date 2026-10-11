#!/usr/bin/env python3
"""
5' completeness of reads, and the untemplated-G (cap) signature at their 5' ends, per cell
cluster: what share of polyA-anchored molecules reach the TSS, as a function of distance
from the 3' end, and whether reads starting at a TSS carry the template-switch G's more
often than reads starting inside a transcript body. Measurement only: nothing here changes
how LRAA calls TSS sites.

Reads are taken as alt_termini_read_check.py / read_end_profile_by_cluster.py take them:
LRAA's read filters (Util_funcs.quant_discard_reason; --HiFi identity floor, --rdna_mask_bed),
LRAA's transcribed-strand rule, intron chains from the CIGAR. The 5' clip and whether LRAA
strips it as untemplated G's come from Pretty_alignment (max_untemplated_G_at_TSS,
min_proximal_untemplated_G_at_TSS), not reimplemented here.

1. Single-TSS genes (completeness curves, per-gene truncation).
   A gene (gene_key of the dataset's site table, `--sites`) qualifies if it has exactly one
   TSS site there, its symbol is a gene_name in `--gtf` on the same contig and strand, and
   a GENCODE transcript of the gene has the site in (or within 50 nt upstream of) its first
   exon. The main transcript is, among those, the one whose full intron chain most reads
   reproduce (ties: more reads fitting, then longer). Restricting to transcripts that start
   at the site matters: a fragment model starting downstream would make truncated reads
   look full length. Its CAGE
   status (`--cage_peaks`, FANTOM5): peaks on the gene's strand over the main transcript
   plus 100 nt upstream, merged within 100 nt into clusters; "single" if the top cluster
   holds >= 90% of the peak score and lies within 100 nt of the LRAA site, "multi" if not,
   "none" if there are no peaks.
   Reads fitting the main transcript (introns all among its introns) whose 3' end lies within
   50 nt of the observed 3' end (the most common read end, 10 nt bins) are polyA-anchored.
   For each, extent = distance along the transcript from its 5' end to the observed 3' end;
   reads starting within 50 nt of the TSS (or upstream of it) reach it (extent = gene length
   L = observed end - TSS).

2. Cap signature. Per read, the transcript-5' soft clip: "noclip" (none), "G" (clipped,
   and LRAA strips it all as untemplated G's), "other" (clip LRAA leaves). An untemplated G
   aligns when the genome has a G just upstream, so the genomic base upstream of the read's
   5' end is recorded (upstream_G). Classes of 5' end: "TSS" (reaches it), "near" (50-100 nt
   downstream), "body" (> 100 nt downstream).

3. Event sites (`--events`, high-confidence site-level TSS events): for each event's
   upstream (U) and downstream (D) site, reads whose 5' end is within 25 nt of the site
   (the site-counting tolerance), by cluster x clip class x upstream_G; and, from the reads
   starting at U, whether D lies in their aligned exonic blocks (so a truncated U molecule
   can start at D) and the median aligned distances U -> read 3' end and D -> read 3' end.

Outputs (<prefix>.):
  genes.tsv            one row per candidate gene: status, main transcript, L, CAGE status,
                       first-exon length, GC of the first 100 nt, reads
  extent.tsv.gz        gene_key, cluster, n_anchored, n_reach, then extent histogram rows
                       (25 nt bins) of the anchored reads that do not reach the TSS
  capsig.tsv.gz        gene_key, cluster, five_class, upstream_G, clip, reads
  event_pairs.tsv      per (gene_key, U site, D site): D_in_U_body share, d_U, d_D (median nt)
  event_sites.tsv.gz   site_id, cluster, clip, upstream_G, reads
"""

import argparse
import collections
import csv
import gzip
import logging
import multiprocessing
import os
import statistics
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "diff_iso_usage"))
sys.path.insert(0, HERE)
from alt_termini_read_check import (make_read_filter, parse_cell_clusters, read_geometry,  # noqa: E402
                                    transcript_strand, register_contig)
from read_end_profile_by_cluster import parse_ref_gtf, introns_of, chain_in  # noqa: E402

sys.path.insert(0, os.path.join(HERE, "../../pylib"))
import Util_funcs  # noqa: E402
from Pretty_alignment import Pretty_alignment  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

END_TOL = 50      # read 3' end within this of the observed end: polyA-anchored
REACH_TOL = 50    # read 5' end within this of the TSS: reaches it
SITE_TOL = 25     # site-counting tolerance (half of max_dist_between_alt_TSS_sites)
BIN = 25
COMP = str.maketrans("ACGTN", "TGCAN")

G = {}  # per-worker globals


def tx_coord(pos, exons, strand):
    """distance (nt) of genomic pos from the transcript's 5' end along its exons; None in an
    intron. Beyond the transcript's ends the genomic distance continues the scale (negative
    upstream of the 5' end, > length past the 3' end) rather than clamping, so reads whose
    ends run past the annotated model keep their true positions."""
    tx_len = sum(e - s + 1 for s, e in exons)
    off = 0
    for s, e in exons:
        if s <= pos <= e:
            d = off + (pos - s)
            return d if strand == "+" else tx_len - 1 - d
        off += e - s + 1
    if pos < exons[0][0]:
        return -(exons[0][0] - pos) if strand == "+" else tx_len - 1 + (exons[0][0] - pos)
    if pos > exons[-1][1]:
        return tx_len - 1 + (pos - exons[-1][1]) if strand == "+" else -(pos - exons[-1][1])
    return None


def init_worker(args):
    G["args"] = args
    G["keep"] = make_read_filter(HiFi=args.HiFi, rdna_mask_bed=args.rdna_mask_bed)
    G["cells"] = parse_cell_clusters(args.cell_clusters)
    G["genome"] = pysam.FastaFile(args.genome)
    G["bam"] = pysam.AlignmentFile(args.bam)


def read_records(contig, lo, hi, strand):
    """filtered reads on the strand overlapping [lo, hi] (1-based): (read, cluster, lend, rend, introns)"""
    keep, cells, genome = G["keep"], G["cells"], G["genome"]
    for r in G["bam"].fetch(contig, max(0, lo - 1), hi):
        if not r.has_tag("CB"):
            continue
        cl = cells.get(r.get_tag("CB"))
        if cl is None or not keep(r) or transcript_strand(r, genome) != strand:
            continue
        lend, rend, introns = read_geometry(r)
        yield r, cl, lend, rend, introns


def clip_class(r, strand):
    """noclip / G (5' clip that LRAA strips as untemplated G's) / other"""
    pa = Pretty_alignment.get_pretty_alignment(r)
    if strand == "+":
        orig, resid = pa.orig_left_soft_clipping, pa.left_soft_clipping
    else:
        orig, resid = pa.orig_right_soft_clipping, pa.right_soft_clipping
    if not orig:
        return "noclip"
    return "G" if resid == 0 else "other"


def upstream_is_G(contig, lend, rend, strand):
    """is the genomic base just 5' of the read's 5' end a G (on the transcript strand)"""
    seq = Util_funcs.contig_seq_for_strand_check(contig)
    if seq is None:
        return "NA"
    if strand == "+":
        i = lend - 2
        b = seq[i] if 0 <= i < len(seq) else "N"
    else:
        i = rend
        b = seq[i].translate(COMP) if 0 <= i < len(seq) else "N"
    return "1" if b == "G" else "0"


def blocks_of(lend, rend, introns):
    out, s = [], lend
    for a, b in introns:
        out.append((s, a - 1))
        s = b + 1
    out.append((s, rend))
    return out


def mode10(vals):
    c = collections.Counter(v // 10 for v in vals)
    return c.most_common(1)[0][0] * 10 + 5 if c else None


def gc_first(contig, exons, strand, start_tx, n=100):
    seq = Util_funcs.contig_seq_for_strand_check(contig)
    if seq is None:
        return ""
    s = "".join(seq[a - 1:b] for a, b in exons)
    if strand == "-":
        s = s[::-1].translate(COMP)
    w = s[start_tx:start_tx + n]
    return round((w.count("G") + w.count("C")) / len(w), 3) if w else ""


def cage_status(contig, strand, exons, site_pos, peaks):
    lo, hi = exons[0][0], exons[-1][1]
    if strand == "+":
        lo -= 100
    else:
        hi += 100
    ps = sorted((s, e, sc) for s, e, sc in peaks.get((contig, strand), []) if e >= lo and s <= hi)
    if not ps:
        return "none"
    clusters = []
    for s, e, sc in ps:
        if clusters and s - clusters[-1][1] <= 100:
            clusters[-1][1] = max(clusters[-1][1], e)
            clusters[-1][2] += sc
        else:
            clusters.append([s, e, sc])
    total = sum(c[2] for c in clusters)
    top = max(clusters, key=lambda c: c[2])
    near = top[0] - 100 <= site_pos <= top[1] + 100
    return "single" if top[2] >= 0.9 * total and near else "multi"


def do_contig(job):
    contig, genes, pairs, sites = job
    args = G["args"]
    register_contig(G["genome"], contig)
    gene_rows, ext_rows, cap = [], [], collections.Counter()
    for g in genes:
        recs, strand, site_pos = g["recs"], g["strand"], g["site_pos"]
        row = {"gene_key": g["gene_key"], "gene_symbol": g["symbol"], "status": "", "main_transcript": "",
               "L": "", "tss_tx": "", "end_tx": "", "tx_len": "", "first_exon_len": "", "gc100": "",
               "cage": "", "n_fit_reads": 0, "n_anchored": 0}
        lo = min(ex[0][0] for _, _, ex in recs.values())
        hi = max(ex[-1][1] for _, _, ex in recs.values())
        tx_introns = {t: introns_of(ex) for t, (_, _, ex) in recs.items()}
        reads = []
        for r, cl, lend, rend, introns in read_records(contig, lo, hi, strand):
            if rend < lo or lend > hi:
                continue
            fits = [t for t, ti in tx_introns.items() if chain_in(introns, ti)]
            if fits:
                reads.append((r, cl, lend, rend, introns, fits))
        if not reads:
            row["status"] = "no_reads"
            gene_rows.append(row)
            continue
        tx_len = {t: sum(e - s + 1 for s, e in ex) for t, (_, _, ex) in recs.items()}
        # candidates: transcripts whose first exon holds the TSS site (a model starting
        # downstream, e.g. a fragment, would make truncated reads look full length)
        def first_exon_len(t):
            e = recs[t][2]
            return (e[0][1] - e[0][0] + 1) if strand == "+" else (e[-1][1] - e[-1][0] + 1)
        cands = []
        for t in recs:
            c = tx_coord(site_pos, recs[t][2], strand)
            if c is not None and c <= first_exon_len(t) + REACH_TOL and \
                    (recs[t][2][0][0] - REACH_TOL <= site_pos if strand == "+" else site_pos <= recs[t][2][-1][1] + REACH_TOL):
                cands.append(t)
        if not cands:
            row["status"] = "no_transcript_starting_at_site"
            gene_rows.append(row)
            continue
        fsm = collections.Counter(t for _, _, _, _, introns, fits in reads if introns for t in fits
                                  if len(introns) == len(tx_introns[t]))
        subset = collections.Counter(t for *_, fits in reads for t in fits)
        main_t = max(cands, key=lambda t: (fsm[t], subset[t], tx_len[t]))
        ex = recs[main_t][2]
        row.update(main_transcript=main_t, tx_len=tx_len[main_t],
                   first_exon_len=(ex[0][1] - ex[0][0] + 1) if strand == "+" else (ex[-1][1] - ex[-1][0] + 1),
                   cage=cage_status(contig, strand, ex, site_pos, G["peaks"]))
        tss_tx = tx_coord(site_pos, ex, strand)
        if tss_tx is None:
            row["status"] = "site_not_exonic_in_main"
            gene_rows.append(row)
            continue
        mine = []
        for r, cl, lend, rend, introns, fits in reads:
            if main_t not in fits:
                continue
            five, three = (lend, rend) if strand == "+" else (rend, lend)
            f5, f3 = tx_coord(five, ex, strand), tx_coord(three, ex, strand)
            if f5 is None or f3 is None:
                continue
            mine.append((r, cl, lend, rend, f5, f3))
        row["n_fit_reads"] = len(mine)
        end_tx = mode10([f3 for *_, f3 in mine]) if mine else None
        if end_tx is None or end_tx - tss_tx < 200:
            row["status"] = "too_short_or_no_reads"
            gene_rows.append(row)
            continue
        L = end_tx - tss_tx
        row.update(tss_tx=tss_tx, end_tx=end_tx, L=L, gc100=gc_first(contig, ex, strand, tss_tx))
        per = collections.defaultdict(lambda: [0, 0, collections.Counter()])
        for r, cl, lend, rend, f5, f3 in mine:
            if abs(f3 - end_tx) > END_TOL:
                continue
            d5 = f5 - tss_tx
            five_class = "TSS" if d5 <= REACH_TOL else ("near" if d5 <= 100 else "body")
            rec = per[cl]
            rec[0] += 1
            if five_class == "TSS":
                rec[1] += 1
            else:
                rec[2][(end_tx - f5) // BIN * BIN] += 1
            cap[(g["gene_key"], cl, five_class, upstream_is_G(contig, lend, rend, strand),
                 clip_class(r, strand))] += 1
        row["n_anchored"] = sum(v[0] for v in per.values())
        row["status"] = "ok"
        gene_rows.append(row)
        for cl, (n, nr, hist) in per.items():
            ext_rows.append((g["gene_key"], cl, n, nr, ";".join(f"{b}:{k}" for b, k in sorted(hist.items()))))

    # event sites: reads starting at each site; D inside U-initiated reads' bodies
    site_counts = collections.Counter()
    for sid, pos, strand in sites:
        for r, cl, lend, rend, introns in read_records(contig, pos - SITE_TOL, pos + SITE_TOL, strand):
            five = lend if strand == "+" else rend
            if abs(five - pos) > SITE_TOL:
                continue
            site_counts[(sid, cl, clip_class(r, strand), upstream_is_G(contig, lend, rend, strand))] += 1
    pair_rows = []
    for gene_key, u_sid, u_pos, d_sid, d_pos, strand in pairs:
        n_u, n_cover, dus, dds = 0, 0, [], []
        for r, cl, lend, rend, introns in read_records(contig, u_pos - SITE_TOL, u_pos + SITE_TOL, strand):
            five = lend if strand == "+" else rend
            if abs(five - u_pos) > SITE_TOL:
                continue
            n_u += 1
            bl = blocks_of(lend, rend, introns)
            total = sum(e - s + 1 for s, e in bl)
            dus.append(total)
            if any(s - 10 <= d_pos <= e + 10 for s, e in bl):
                n_cover += 1
                if strand == "+":
                    dds.append(sum(max(0, e - max(s, d_pos) + 1) for s, e in bl))
                else:
                    dds.append(sum(max(0, min(e, d_pos) - s + 1) for s, e in bl))
        pair_rows.append({"gene_key": gene_key, "U_site": u_sid, "D_site": d_sid, "U_pos": u_pos, "D_pos": d_pos,
                          "U_start_reads": n_u, "D_in_U_body": round(n_cover / n_u, 3) if n_u else "",
                          "d_U": statistics.median(dus) if dus else "", "d_D": statistics.median(dds) if dds else ""})
    logger.info("%s: %d genes, %d event pairs, %d event sites", contig, len(genes), len(pairs), len(sites))
    return gene_rows, ext_rows, cap, pair_rows, site_counts


def load_peaks(path):
    peaks = collections.defaultdict(list)
    if not path:
        return peaks
    with gzip.open(path, "rt") if path.endswith(".gz") else open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 6:
                continue
            peaks[(f[0], f[5].strip())].append((int(f[1]) + 1, int(f[2]), float(f[4])))
    return peaks


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True, help="reference annotation with gene_name attributes")
    ap.add_argument("--sites", required=True, help="the dataset's site table (site-usage <prefix>.sites.tsv)")
    ap.add_argument("--cage_peaks", default=None, help="FANTOM5 CAGE peak bed (hg38_fair+new_CAGE_peaks_phase1and2.bed.gz)")
    ap.add_argument("--events", default=None, help="site-level TSS events (only high_confidence rows are used)")
    ap.add_argument("--cell_clusters", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--HiFi", action="store_true")
    ap.add_argument("--rdna_mask_bed", default=None)
    ap.add_argument("--CPU", type=int, default=4)
    ap.add_argument("--contigs", default=None, help="comma-separated contigs to restrict to (testing)")
    ap.add_argument("--output_prefix", required=True)
    args = ap.parse_args()

    tss = collections.defaultdict(list)
    site_pos = {}
    sym = {}
    with open(args.sites) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["kind"] != "TSS":
                continue
            tss[r["gene_key"]].append(r["site_id"])
            site_pos[r["site_id"]] = int(r["pos"])
            sym[r["gene_key"]] = r["gene_symbol"]
    single = {k: v[0] for k, v in tss.items() if len(v) == 1}
    txs = parse_ref_gtf(args.gtf, {sym[k] for k in single})
    by_contig = collections.defaultdict(lambda: [[], [], []])
    n_cand = 0
    for k, sid in single.items():
        s = sym[k]
        _, contig, strand = k.split("|")
        recs = {t: v for t, v in txs.get(s, {}).items() if v[0] == contig and v[1] == strand}
        if not recs:
            continue
        by_contig[contig][0].append({"gene_key": k, "symbol": s, "strand": strand, "recs": recs,
                                     "site_pos": site_pos[sid]})
        n_cand += 1
    logger.info("%d genes with one TSS site, %d of them in the gtf", len(single), n_cand)

    if args.events:
        seen_pairs, seen_sites = set(), set()
        with open(args.events) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["high_confidence"] != "True":
                    continue
                _, contig, strand = r["gene_key"].split("|")
                down = r["gained_is_downstream"] == "True"
                u_sid, d_sid = (r["lost_site"], r["gained_site"]) if down else (r["gained_site"], r["lost_site"])
                u_pos, d_pos = int(u_sid.split(":")[2]), int(d_sid.split(":")[2])
                key = (r["gene_key"], u_sid, d_sid)
                if key not in seen_pairs:
                    seen_pairs.add(key)
                    by_contig[contig][1].append((r["gene_key"], u_sid, u_pos, d_sid, d_pos, strand))
                for sid, pos in ((u_sid, u_pos), (d_sid, d_pos)):
                    if sid not in seen_sites:
                        seen_sites.add(sid)
                        by_contig[contig][2].append((sid, pos, strand))

    peaks = load_peaks(args.cage_peaks)
    G["peaks"] = peaks  # inherited by forked workers
    only = set(args.contigs.split(",")) if args.contigs else None
    jobs = sorted(((c, *v) for c, v in by_contig.items() if only is None or c in only), key=lambda j: -len(j[1]))
    gene_rows, ext_rows, cap, pair_rows, site_counts = [], [], collections.Counter(), [], collections.Counter()
    with multiprocessing.get_context("fork").Pool(args.CPU, initializer=init_worker, initargs=(args,)) as pool:
        for gr, er, cp, pr, sc in pool.imap_unordered(do_contig, jobs):
            gene_rows += gr
            ext_rows += er
            cap.update(cp)
            pair_rows += pr
            site_counts.update(sc)

    p = args.output_prefix
    cols = ["gene_key", "gene_symbol", "status", "main_transcript", "tx_len", "tss_tx", "end_tx", "L",
            "first_exon_len", "gc100", "cage", "n_fit_reads", "n_anchored"]
    with open(f"{p}.genes.tsv", "w") as out:
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(gene_rows)
    with gzip.open(f"{p}.extent.tsv.gz", "wt") as out:
        out.write("gene_key\tcluster\tn_anchored\tn_reach\textent_hist\n")
        for r in ext_rows:
            out.write("\t".join(map(str, r)) + "\n")
    with gzip.open(f"{p}.capsig.tsv.gz", "wt") as out:
        out.write("gene_key\tcluster\tfive_class\tupstream_G\tclip\treads\n")
        for k, v in sorted(cap.items()):
            out.write("\t".join(map(str, k)) + f"\t{v}\n")
    if args.events:
        with open(f"{p}.event_pairs.tsv", "w") as out:
            w = csv.DictWriter(out, fieldnames=["gene_key", "U_site", "D_site", "U_pos", "D_pos", "U_start_reads",
                                                "D_in_U_body", "d_U", "d_D"], delimiter="\t")
            w.writeheader()
            w.writerows(pair_rows)
        with gzip.open(f"{p}.event_sites.tsv.gz", "wt") as out:
            out.write("site_id\tcluster\tclip\tupstream_G\treads\n")
            for k, v in sorted(site_counts.items()):
                out.write("\t".join(map(str, k)) + f"\t{v}\n")
    st = collections.Counter(r["status"] for r in gene_rows)
    logger.info("genes by status: %s; CAGE status of ok genes: %s", dict(st),
                dict(collections.Counter(r["cage"] for r in gene_rows if r["status"] == "ok")))


if __name__ == "__main__":
    main()
