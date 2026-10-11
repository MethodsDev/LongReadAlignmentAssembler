#!/usr/bin/env python3
"""
Features of tandem TSS pairs: two TSSs of the same splice pattern (same intron chain), the
downstream one (D) inside the upstream one's (U) first exon. Here a read from a molecule
that started at U but lost its 5' end looks exactly like a read that started at D, so these
are the pairs where a degradation / truncation product can pass for an alternative TSS.
Measurement only: nothing here changes how LRAA calls TSS sites. The untemplated-G (cap)
signature is deliberately not a feature (it depends on the library chemistry).

Pairs (one row each, `<prefix>.pairs.tsv`; `set` says where a pair came from, a pair can
be in several):
  tandem     every pair of TSS-defined 5' ends (TSS "True") of transcripts sharing an intron
             chain in `--lraa_gtf` (ends within 25 nt merged), at least 30 nt apart.
             `cage_both` marks those with distinct FANTOM5 CAGE peaks at both ends.
  negative   read-start peaks inside the first exon of single-TSS genes (`--single_tss_genes`,
             the genes.tsv of tss_5prime_completeness.py; status ok, CAGE single or none, or
             `--negative_cage single` for single only): positions >= 50 nt downstream of
             the gene's TSS site and >= 10 nt before the first donor, in first exons >= 150
             nt, that pass LRAA's own TSS peak calling on the pooled read starts
             (Splice_graph.aggregate_sites_within_window, window
             max_dist_between_alt_TSS_sites, min min_alignments_define_TSS_site, then
             filter_non_peaky_positions), at most 3 per gene. False by construction: the
             gene has one TSS, so a read starting there is a truncated molecule.
  isoT_event site pairs of IsoformTermini events with differing_end TSS (`--isoT_events`).
  site_event site-level TSS events (`--site_events`) whose two sites are 5' ends of a
             shared splice pattern (found through `--sites` transcript_ids), any class.
  tandem_ref another dataset's CAGE-supported tandem pairs (`--ref_pairs`, its pairs.tsv),
             measured in this dataset's reads (cross-platform checks without enough of the
             dataset's own models).

Reads: LRAA's read filters and transcribed-strand rule (as alt_termini_read_check.py),
cell clusters from CB tags (`--cell_clusters`; without it every read is one group, "all"),
and, as LRAA counts a read 5' end at a TSS, only reads with no residual 5' soft clip after
LRAA's clip handling (Pretty_alignment; used here only to decide eligibility). A read
counts for a pair if its first intron (in transcript orientation) is the pattern's first
intron (+/- 2 nt): it comes from the pattern and starts in the first exon.

Positions are offsets from U along the transcript (D_off = D's offset; the first donor at
donor_off). Read-start windows: U +/- 25, D +/- 25; gap (25, D_off - 25); downstream
(D_off + 25, min(D_off + 125, donor_off - 5)].

Features (pooled over cells; per-cluster counts in `<prefix>.pair_clusters.tsv.gz`):
  bg_log2      log2((D_win + 1) / (expected + 1)): expected = start density in the gap and
               downstream windows x 51 (D's excess over the truncation background around it)
  D_sharp      share of the D-window starts within +/- 2 nt of their mode (U_sharp likewise)
  mode_repro   share of clusters (>= 8 D-window starts) whose mode is within 3 nt of the
               pooled mode
  d_U          median aligned exonic length of the reads starting at U (U's distance from
               the 3' end); d_D = d_U - D_off
  sequence (transcript strand, at D and U): inr (YR at -1/+1), tata (best match of
               TATAWAWR in -35..-20, of 8), gc250 and cpg_oe250 (+/- 250 nt), gc40up and
               mfe40up (the 40 nt 5' of the site, ViennaRNA), mfe60 (-50..+10), dinuc (-1,+1)
  cage_dist / reftss_dist: nt to the nearest FANTOM5 CAGE peak / refTSS on the strand

Outputs (<prefix>.): pairs.tsv, pair_clusters.tsv.gz, event_map.tsv (source, event key,
pair_id), summary.tsv (median aligned read length over a genome-wide sample).
"""

import argparse
import bisect
import collections
import csv
import gzip
import logging
import multiprocessing
import os
import random
import statistics
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "diff_iso_usage"))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "../../pylib"))
from alt_termini_read_check import (make_read_filter, parse_cell_clusters, read_geometry,  # noqa: E402
                                    transcript_strand, register_contig)
import Util_funcs  # noqa: E402
import LRAA_Globals  # noqa: E402
from Pretty_alignment import Pretty_alignment  # noqa: E402
from Splice_graph import aggregate_sites_within_window, filter_non_peaky_positions  # noqa: E402

try:
    import RNA  # ViennaRNA
except ImportError:
    RNA = None

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

WIN = 25           # site window (half of max_dist_between_alt_TSS_sites)
INTRON_TOL = 2
UPSTREAM_PAD = 100
COMP = str.maketrans("ACGTN", "TGCAN")
G = {}


# ---------------------------------------------------------------- inputs

def gtf_attr(a, key):
    k = key + ' "'
    return a.split(k, 1)[1].split('"', 1)[0] if k in a else None


def parse_lraa_gtf(gtf):
    """transcript_id -> dict(contig, strand, exons, tss_flag, symbol)"""
    tx = {}
    for line in open(gtf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        tid = gtf_attr(f[8], "transcript_id")
        if tid is None:
            continue
        if f[2] == "transcript":
            d = tx.setdefault(tid, {"exons": []})
            d.update(contig=f[0], strand=f[6], tss=gtf_attr(f[8], "TSS") == "True",
                     symbol=tid.split("^", 1)[0] if "^" in tid else "")
        elif f[2] == "exon":
            d = tx.setdefault(tid, {"exons": []})
            d["exons"].append((int(f[3]), int(f[4])))
            d.setdefault("contig", f[0])
            d.setdefault("strand", f[6])
    for d in tx.values():
        d["exons"].sort()
        d.setdefault("tss", False)
        d.setdefault("symbol", "")
    return tx


def chain_of(exons):
    return tuple((exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1))


def five_end(exons, strand):
    return exons[0][0] if strand == "+" else exons[-1][1]


def first_intron_donor(exons, strand):
    ch = chain_of(exons)
    if strand == "+":
        return ch[0], exons[0][1]
    return ch[-1], exons[-1][0]


def load_intervals(bed, has_header=False):
    """(contig, strand) -> sorted [(start1, end, id, score)]"""
    idx = collections.defaultdict(list)
    op = gzip.open if bed.endswith(".gz") else open
    with op(bed, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or not f[1].isdigit():
                continue
            try:
                sc = float(f[4])
            except ValueError:
                sc = 0.0
            idx[(f[0], f[5])].append((int(f[1]) + 1, int(f[2]), f[3], sc))
    for k in idx:
        idx[k].sort()
    starts = {k: [x[0] for x in v] for k, v in idx.items()}
    return idx, starts


def nearest(intervals, contig, strand, pos):
    """(distance, id) to the nearest interval on the strand"""
    idx, starts = intervals
    v = idx.get((contig, strand))
    if not v:
        return None, ""
    i = bisect.bisect_right(starts[(contig, strand)], pos)
    best = (None, "")
    for j in range(max(0, i - 50), min(len(v), i + 3)):
        s, e, pid, _ = v[j]
        d = 0 if s <= pos <= e else min(abs(pos - s), abs(pos - e))
        if best[0] is None or d < best[0]:
            best = (d, pid)
    return best


# ---------------------------------------------------------------- pair building

def build_pairs(args, cage):
    """groups: (contig, strand, first_intron, donor) -> dict(pairs=[...], negatives=[...])"""
    tx = parse_lraa_gtf(args.lraa_gtf)
    pats = collections.defaultdict(list)    # (contig, strand, chain) -> [tid]
    for tid, d in tx.items():
        if len(d["exons"]) < 2 or "contig" not in d:
            continue
        pats[(d["contig"], d["strand"], chain_of(d["exons"]))].append(tid)
    tid_pat = {}
    for k, tids in pats.items():
        for t in tids:
            tid_pat[t] = k
    by_suffix = {t.split("^", 1)[-1]: t for t in tx}

    pairs = {}      # (contig, strand, U, D) -> row
    event_map = []

    def add(contig, strand, a, b, pat, sets, symbol):
        if abs(a - b) < 30:
            return None
        U, D = (min(a, b), max(a, b)) if strand == "+" else (max(a, b), min(a, b))
        intron, donor = first_intron_donor(tx[pats[pat][0]]["exons"], strand)
        donor_off = (donor - U) if strand == "+" else (U - donor)
        D_off = abs(D - U)
        if D_off >= donor_off:    # D not inside the first exon
            return None
        key = (contig, strand, U, D)
        row = pairs.get(key)
        if row is None:
            row = pairs[key] = {"contig": contig, "strand": strand, "U": U, "D": D, "intron": intron,
                                "donor": donor, "sets": set(), "gene_symbol": symbol}
        row["sets"] |= set(sets)
        if symbol and not row["gene_symbol"]:
            row["gene_symbol"] = symbol
        return key

    # tandem: TSS-defined 5' ends sharing a chain
    for pat, tids in pats.items():
        contig, strand, _ = pat
        ends = sorted({five_end(tx[t]["exons"], strand) for t in tids if tx[t]["tss"]})
        merged = []
        for e in ends:
            if merged and e - merged[-1] <= WIN:
                continue
            merged.append(e)
        syms = collections.Counter(tx[t]["symbol"] for t in tids if tx[t]["symbol"])
        sym = syms.most_common(1)[0][0] if syms else ""
        for i in range(len(merged)):
            for j in range(i + 1, len(merged)):
                add(contig, strand, merged[i], merged[j], pat, ["tandem"], sym)

    def pat_of_tids(tids):
        out = []
        for t in tids:
            full = tx.get(t) and t or by_suffix.get(t)
            if full and full in tid_pat:
                out.append(tid_pat[full])
        return out

    def tx_five_end(tids):
        ends = collections.Counter()
        cs = None
        for t in tids:
            full = tx.get(t) and t or by_suffix.get(t)
            if full and "contig" in tx[full]:
                d = tx[full]
                ends[five_end(d["exons"], d["strand"])] += 1
                cs = (d["contig"], d["strand"])
        if not ends:
            return None
        return cs[0], ends.most_common(1)[0][0], cs[1]

    def site_pos(site_id):
        p = (site_id or "").split(":")
        if len(p) != 4 or not p[2].isdigit():
            return None
        return p[1], int(p[2]), p[3]

    if args.isoT_events:
        with open(args.isoT_events) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("differing_end") != "TSS":
                    continue
                ekey = "|".join([r["gene_key"], r["cluster_A"], r["cluster_B"], r["gained_feature"], r["lost_feature"]])
                ga, lo_ = r["gained_transcript_ids"].split(","), r["lost_transcript_ids"].split(",")
                a, b = site_pos(r["gained_TSS_site"]), site_pos(r["lost_TSS_site"])
                # a feature without an annotated TSS site: its transcripts' most common 5' end
                a = a or tx_five_end(ga)
                b = b or tx_five_end(lo_)
                if a is None or b is None:
                    event_map.append(("isoT", ekey, "TSS_site_missing"))
                    continue
                c1, p1, s1 = a
                p2 = b[1]
                pa = set(pat_of_tids(ga))
                pb = set(pat_of_tids(lo_))
                shared = pa & pb
                if not shared:
                    event_map.append(("isoT", ekey, "no_shared_pattern_in_gtf"))
                    continue
                key = add(c1, s1, p1, p2, next(iter(shared)), ["isoT_event"], r["gene_symbol"])
                event_map.append(("isoT", ekey, key if key else "D_not_in_first_exon"))

    if args.site_events and args.sites:
        site_tids = {}
        with open(args.sites) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["kind"] == "TSS":
                    site_tids[r["site_id"]] = [t for t in (r.get("transcript_ids") or "").split(",") if t]
        with open(args.site_events) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                ekey = "|".join([r["gene_key"], r["cluster_A"], r["cluster_B"], r["gained_site"], r["lost_site"]])
                pa = set(pat_of_tids(site_tids.get(r["gained_site"], [])))
                pb = set(pat_of_tids(site_tids.get(r["lost_site"], [])))
                shared = pa & pb
                if not shared:
                    event_map.append(("site", ekey, "no_shared_pattern"))
                    continue
                a, b = site_pos(r["gained_site"]), site_pos(r["lost_site"])
                if a is None or b is None:
                    event_map.append(("site", ekey, "TSS_site_missing"))
                    continue
                c1, p1, s1 = a
                p2 = b[1]
                key = add(c1, s1, p1, p2, next(iter(shared)), ["site_event"], r["gene_symbol"])
                event_map.append(("site", ekey, key if key else "D_not_in_first_exon"))

    # another dataset's CAGE-supported tandem pairs, measured in these reads ("tandem_ref")
    if args.ref_pairs:
        n0 = len(pairs)
        with open(args.ref_pairs) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if "tandem" not in r["sets"] or r["cage_both"] != "1" or "negative" in r["sets"]:
                    continue
                key = (r["contig"], r["strand"], int(r["U"]), int(r["D"]))
                row = pairs.get(key)
                if row is None:
                    row = pairs[key] = {"contig": key[0], "strand": key[1], "U": key[2], "D": key[3],
                                        "intron": (int(r["intron_lo"]), int(r["intron_hi"])), "donor": int(r["donor"]),
                                        "sets": set(), "gene_symbol": r["gene_symbol"]}
                row["sets"].add("tandem_ref")
        logger.info("%d pairs added from %s", len(pairs) - n0, args.ref_pairs)

    # negatives: genes whose peaks are found later from the reads
    negs = []
    if args.single_tss_genes:
        ok_cage = {"single"} if args.negative_cage == "single" else {"single", "none"}
        want = {}
        with open(args.single_tss_genes) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["status"] == "ok" and r["cage"] in ok_cage:
                    want[r["main_transcript"]] = r
        site_of = {}
        if args.negative_sites:
            with open(args.negative_sites) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r["kind"] == "TSS":
                        site_of[r["gene_key"]] = int(r["pos"])
        # a genes table may carry the TSS position itself (site_pos), e.g. one pooled from datasets
        for r in want.values():
            if r.get("site_pos"):
                site_of[r["gene_key"]] = int(r["site_pos"])
        ex = collections.defaultdict(list)
        meta = {}
        for line in open(args.ref_gtf):
            f = line.split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            t = gtf_attr(f[8], "transcript_id")
            if t in want:
                ex[t].append((int(f[3]), int(f[4])))
                meta[t] = (f[0], f[6])
        for t, r in want.items():
            if t not in ex or r["gene_key"] not in site_of:
                continue
            exons = sorted(ex[t])
            if len(exons) < 2:
                continue
            contig, strand = meta[t]
            intron, donor = first_intron_donor(exons, strand)
            U = site_of[r["gene_key"]]
            donor_off = (donor - U) if strand == "+" else (U - donor)
            if donor_off < 150:
                continue
            negs.append({"contig": contig, "strand": strand, "U": U, "intron": intron, "donor": donor,
                         "gene_symbol": r["gene_symbol"], "gene_key": r["gene_key"]})
    logger.info("%d pairs from the gtf / events, %d single-TSS genes to search for negatives", len(pairs), len(negs))
    return pairs, negs, event_map


# ---------------------------------------------------------------- reads

def init_worker(args):
    G["args"] = args
    G["keep"] = make_read_filter(HiFi=args.HiFi, rdna_mask_bed=args.rdna_mask_bed)
    G["cells"] = parse_cell_clusters(args.cell_clusters) if args.cell_clusters else None
    G["genome"] = pysam.FastaFile(args.genome)
    G["bam"] = pysam.AlignmentFile(args.bam)


def eligible(r, strand):
    """no residual 5' soft clip once LRAA has handled the clip (LRAA's TSS counting rule)"""
    ct = r.cigartuples
    clipped = (ct[0][0] == 4) if strand == "+" else (ct[-1][0] == 4)
    if not clipped:
        return True
    pa = Pretty_alignment.get_pretty_alignment(r)
    resid = pa.left_soft_clipping if strand == "+" else pa.right_soft_clipping
    return resid == 0


def group_starts(contig, strand, intron, lo, hi):
    """read 5' starts in [lo, hi] of reads whose transcript-first intron is `intron`:
    cluster -> Counter(pos), and pos -> [exonic aligned length]"""
    keep, cells, genome = G["keep"], G["cells"], G["genome"]
    starts = collections.defaultdict(collections.Counter)
    lens = collections.defaultdict(list)
    for r in G["bam"].fetch(contig, max(0, lo - 1), hi):
        if cells is not None:
            if not r.has_tag("CB"):
                continue
            cl = cells.get(r.get_tag("CB"))
            if cl is None:
                continue
        else:
            cl = "all"
        if not keep(r) or transcript_strand(r, genome) != strand:
            continue
        lend, rend, introns = read_geometry(r)
        if not introns:
            continue
        fi = introns[0] if strand == "+" else introns[-1]
        if abs(fi[0] - intron[0]) > INTRON_TOL or abs(fi[1] - intron[1]) > INTRON_TOL:
            continue
        p5 = lend if strand == "+" else rend
        if p5 < lo or p5 > hi or not eligible(r, strand):
            continue
        starts[cl][p5] += 1
        lens[p5].append((rend - lend + 1) - sum(b - a + 1 for a, b in introns))
    return starts, lens


def oseq(contig, strand, pos, a, b):
    """transcript-strand sequence from offset a to b (inclusive) relative to pos (+1 = pos)"""
    seq = Util_funcs.contig_seq_for_strand_check(contig)
    if seq is None:
        return ""
    if strand == "+":
        lo, hi = pos + a, pos + b
        s = seq[max(0, lo - 1):max(0, hi)]
    else:
        lo, hi = pos - b, pos - a
        s = seq[max(0, lo - 1):max(0, hi)][::-1].translate(COMP)
    return s.upper()


def seq_features(contig, strand, pos):
    m1 = oseq(contig, strand, pos, -1, -1)
    p1 = oseq(contig, strand, pos, 0, 0)
    tata_win = oseq(contig, strand, pos, -35, -20)
    motif = ["T", "A", "T", "A", "AT", "A", "AT", "AG"]
    tata = 0
    for i in range(0, max(0, len(tata_win) - 7)):
        tata = max(tata, sum(tata_win[i + k] in motif[k] for k in range(8)))
    w = oseq(contig, strand, pos, -250, 250)
    gc = (w.count("G") + w.count("C")) / len(w) if w else None
    cpg = w.count("CG")
    c, g = w.count("C"), w.count("G")
    cpg_oe = cpg * len(w) / (c * g) if c and g else None
    up = oseq(contig, strand, pos, -40, -1)
    gc40 = (up.count("G") + up.count("C")) / len(up) if up else None
    mid = oseq(contig, strand, pos, -50, 10)
    mfe40 = RNA.fold(up.replace("T", "U"))[1] if RNA and up else None
    mfe60 = RNA.fold(mid.replace("T", "U"))[1] if RNA and mid else None
    return {"inr": int(m1 in "CT" and p1 in "AG" and m1 != "" and p1 != ""), "tata": tata,
            "gc250": gc, "cpg_oe250": cpg_oe, "gc40up": gc40, "mfe40up": mfe40, "mfe60": mfe60,
            "dinuc": m1 + p1}


def pair_features(contig, strand, U, D, donor, starts, lens):
    off = (lambda p: p - U) if strand == "+" else (lambda p: U - p)
    D_off = off(D)
    donor_off = off(donor)
    gap_lo, gap_hi = WIN + 1, D_off - WIN - 1                    # inclusive
    dn_lo, dn_hi = D_off + WIN + 1, min(D_off + WIN + 100, donor_off - 5)
    gap_len = max(0, gap_hi - gap_lo + 1)
    dn_len = max(0, dn_hi - dn_lo + 1)
    per = []
    pooled_D = collections.Counter()
    pooled_U = collections.Counter()
    tot = dict(U_win=0, D_win=0, gap_n=0, down_n=0)
    for cl, cnt in starts.items():
        u = d = gp = dn = 0
        dpos = collections.Counter()
        for p, n in cnt.items():
            o = off(p)
            if abs(o) <= WIN:
                u += n
                pooled_U[o] += n
            if abs(o - D_off) <= WIN:
                d += n
                dpos[o] += n
            if gap_lo <= o <= gap_hi:
                gp += n
            if dn_lo <= o <= dn_hi:
                dn += n
        pooled_D.update(dpos)
        mode = dpos.most_common(1)[0][0] if dpos else None
        per.append((cl, u, d, gp, dn, "" if mode is None else mode - D_off))
        tot["U_win"] += u
        tot["D_win"] += d
        tot["gap_n"] += gp
        tot["down_n"] += dn

    def sharp(c):
        n = sum(c.values())
        if n < 5:
            return None, None
        m = c.most_common(1)[0][0]
        return sum(v for k, v in c.items() if abs(k - m) <= 2) / n, m

    D_sharp, D_mode = sharp(pooled_D)
    U_sharp, _ = sharp(pooled_U)
    bg_len = gap_len + dn_len
    if bg_len >= 20:
        expected = (tot["gap_n"] + tot["down_n"]) / bg_len * (2 * WIN + 1)
        import math
        bg_log2 = math.log2((tot["D_win"] + 1) / (expected + 1))
    else:
        expected, bg_log2 = None, None
    repro, n_cl = None, 0
    if D_mode is not None:
        ok = 0
        for cl, cnt in starts.items():
            c = collections.Counter({off(p): n for p, n in cnt.items() if abs(off(p) - D_off) <= WIN})
            if sum(c.values()) >= 8:
                n_cl += 1
                if abs(c.most_common(1)[0][0] - D_mode) <= 3:
                    ok += 1
        repro = ok / n_cl if n_cl >= 2 else None
    ulens = [x for p, v in lens.items() if abs(off(p)) <= WIN for x in v]
    d_U = statistics.median(ulens) if len(ulens) >= 5 else None
    row = dict(tot, D_off=D_off, donor_off=donor_off, gap_len=gap_len, down_len=dn_len,
               expected_bg=expected, bg_log2=bg_log2, D_sharp=D_sharp, U_sharp=U_sharp,
               D_mode_off=(None if D_mode is None else D_mode - D_off), mode_repro=repro,
               n_clusters_mode=n_cl, d_U=d_U, n_U_len=len(ulens),
               d_D=(None if d_U is None else d_U - D_off))
    return row, per


def do_contig(job):
    contig, groups = job
    args = G["args"]
    register_contig(G["genome"], contig)
    out_rows, out_per = [], []
    for g in groups:
        strand, intron, donor = g["strand"], g["intron"], g["donor"]
        us = [p["U"] for p in g["pairs"]] + ([g["neg"]["U"]] if g["neg"] else [])
        if strand == "+":
            lo, hi = min(us) - UPSTREAM_PAD, donor
        else:
            lo, hi = donor, max(us) + UPSTREAM_PAD
        starts, lens = group_starts(contig, strand, intron, lo, hi)
        pairs = list(g["pairs"])
        if g["neg"]:
            ng = g["neg"]
            pooled = collections.Counter()
            for c in starts.values():
                pooled.update(c)
            if pooled:
                grouped = aggregate_sites_within_window(pooled, LRAA_Globals.config["max_dist_between_alt_TSS_sites"],
                                                        LRAA_Globals.config["min_alignments_define_TSS_site"])
                grouped = filter_non_peaky_positions(grouped, pooled, contig, strand)
                U = ng["U"]
                donor_off = (donor - U) if strand == "+" else (U - donor)
                cands = []
                for pos, cnt in grouped:
                    o = (pos - U) if strand == "+" else (U - pos)
                    if 50 <= o <= donor_off - 10:
                        cands.append((cnt, pos))
                for cnt, pos in sorted(cands, reverse=True)[:3]:
                    pairs.append({"contig": contig, "strand": strand, "U": U, "D": pos, "intron": intron,
                                  "donor": donor, "sets": {"negative"}, "gene_symbol": ng["gene_symbol"],
                                  "gene_key": ng["gene_key"]})
        for p in pairs:
            row, per = pair_features(contig, strand, p["U"], p["D"], donor, starts, lens)
            pid = f"{contig}:{strand}:{p['U']}:{p['D']}"
            row.update(pair_id=pid, contig=contig, strand=strand, U=p["U"], D=p["D"], donor=donor,
                       intron_lo=intron[0], intron_hi=intron[1],
                       gene_symbol=p.get("gene_symbol", ""), sets=",".join(sorted(p["sets"])))
            for site, pos in (("D", p["D"]), ("U", p["U"])):
                for k, v in seq_features(contig, strand, pos).items():
                    row[f"{k}_{site}"] = v
                cd, cid = nearest(G["cage"], contig, strand, pos)
                row[f"cage_dist_{site}"], row[f"cage_id_{site}"] = cd, cid
                if G.get("reftss"):
                    rd, _ = nearest(G["reftss"], contig, strand, pos)
                    row[f"reftss_dist_{site}"] = rd
            out_rows.append(row)
            for cl, u, d, gp, dn, m in per:
                out_per.append((pid, cl, u, d, gp, dn, m))
    return out_rows, out_per


def sample_read_length(args, n=100000, seed=1):
    random.seed(seed)
    bam = pysam.AlignmentFile(args.bam)
    keep = make_read_filter(HiFi=args.HiFi, rdna_mask_bed=args.rdna_mask_bed)
    contigs = [(c, l) for c, l in zip(bam.references, bam.lengths)
               if c.startswith("chr") and "_" not in c and c != "chrM"]
    total = sum(l for _, l in contigs)
    cum, acc = [], 0
    for _, l in contigs:
        acc += l
        cum.append(acc)
    lengths = []
    tries = 0
    while len(lengths) < n and tries < 20000:
        tries += 1
        x = random.randrange(total)
        i = bisect.bisect_right(cum, x)
        c = contigs[i][0]
        start = x - (cum[i - 1] if i else 0)
        k = 0
        for r in bam.fetch(c, start, start + 200000):
            if r.reference_start < start or not keep(r):
                continue
            lengths.append(r.query_alignment_length)
            k += 1
            if k >= 500:
                break
    return statistics.median(lengths) if lengths else None, len(lengths)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--lraa_gtf", required=True, help="LRAA gtf: transcripts with TSS attributes")
    ap.add_argument("--cage_peaks", required=True, help="FANTOM5 CAGE peak bed")
    ap.add_argument("--reftss", default=None, help="refTSS bed (optional)")
    ap.add_argument("--cell_clusters", default=None)
    ap.add_argument("--HiFi", action="store_true")
    ap.add_argument("--rdna_mask_bed", default=None)
    ap.add_argument("--isoT_events", default=None, help="IsoformTermini events tsv")
    ap.add_argument("--site_events", default=None, help="site-level TSS events tsv")
    ap.add_argument("--sites", default=None, help="site table with transcript_ids (for --site_events)")
    ap.add_argument("--single_tss_genes", default=None, help="genes.tsv of tss_5prime_completeness.py")
    ap.add_argument("--negative_sites", default=None,
                    help="site table giving those genes' TSS position (or a site_pos column in --single_tss_genes)")
    ap.add_argument("--ref_gtf", default=None, help="reference gtf with the genes' main transcripts")
    ap.add_argument("--ref_pairs", default=None,
                    help="pairs.tsv of another dataset: its CAGE-supported tandem pairs are measured here too (set tandem_ref)")
    ap.add_argument("--negative_cage", choices=["single", "single_or_none"], default="single_or_none")
    ap.add_argument("--length_sample", type=int, default=100000)
    ap.add_argument("--CPU", type=int, default=4)
    ap.add_argument("--contigs", default=None)
    ap.add_argument("--output_prefix", required=True)
    args = ap.parse_args()

    cage = load_intervals(args.cage_peaks)
    pairs, negs, event_map = build_pairs(args, cage)
    groups = {}

    def gkey(contig, strand, intron, donor):
        return (contig, strand, intron, donor)
    for key, p in pairs.items():
        k = gkey(p["contig"], p["strand"], p["intron"], p["donor"])
        groups.setdefault(k, {"strand": p["strand"], "intron": p["intron"], "donor": p["donor"], "pairs": [], "neg": None})
        groups[k]["pairs"].append(p)
    for n in negs:
        k = gkey(n["contig"], n["strand"], n["intron"], n["donor"]) + ("neg",)
        groups[k] = {"strand": n["strand"], "intron": n["intron"], "donor": n["donor"], "pairs": [], "neg": n}
    by_contig = collections.defaultdict(list)
    for k, g in groups.items():
        by_contig[k[0]].append(g)
    contigs = sorted(by_contig, key=lambda c: -len(by_contig[c]))
    if args.contigs:
        want = set(args.contigs.split(","))
        contigs = [c for c in contigs if c in want]

    G["cage"] = cage
    G["reftss"] = load_intervals(args.reftss) if args.reftss else None
    rows, per = [], []
    with multiprocessing.Pool(args.CPU, initializer=init_worker, initargs=(args,)) as pool:
        for i, (r, p) in enumerate(pool.imap_unordered(do_contig, [(c, by_contig[c]) for c in contigs])):
            rows.extend(r)
            per.extend(p)
            logger.info("%d / %d contigs, %d pairs", i + 1, len(contigs), len(rows))

    # one row per pair: a negative can coincide with a tandem pair; merge their sets
    merged = {}
    for r in rows:
        m = merged.get(r["pair_id"])
        if m is None:
            merged[r["pair_id"]] = r
        else:
            m["sets"] = ",".join(sorted(set(m["sets"].split(",")) | set(r["sets"].split(","))))
    rows = list(merged.values())
    for r in rows:
        r["cage_both"] = int(r.get("cage_dist_U") is not None and r.get("cage_dist_D") is not None and
                             r["cage_dist_U"] <= 50 and r["cage_dist_D"] <= 50 and r["cage_id_U"] != r["cage_id_D"])
    cols = ["pair_id", "sets", "gene_symbol", "contig", "strand", "U", "D", "donor", "intron_lo", "intron_hi",
            "D_off", "donor_off",
            "U_win", "D_win", "gap_n", "gap_len", "down_n", "down_len", "expected_bg", "bg_log2",
            "D_sharp", "U_sharp", "D_mode_off", "mode_repro", "n_clusters_mode", "d_U", "d_D", "n_U_len",
            "cage_both"]
    extra = sorted({k for r in rows for k in r} - set(cols))
    with open(f"{args.output_prefix}.pairs.tsv", "wt") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(cols + extra)
        for r in rows:
            w.writerow(["" if r.get(c) is None else (round(r[c], 4) if isinstance(r[c], float) else r[c])
                        for c in cols + extra])
    seen = set()
    with gzip.open(f"{args.output_prefix}.pair_clusters.tsv.gz", "wt") as out:
        out.write("pair_id\tcluster\tU_win\tD_win\tgap_n\tdown_n\tD_mode_off\n")
        for t in per:
            if (t[0], t[1]) in seen:
                continue
            seen.add((t[0], t[1]))
            out.write("\t".join(map(str, t)) + "\n")
    with open(f"{args.output_prefix}.event_map.tsv", "wt") as out:
        out.write("source\tevent_key\tpair_id\n")
        for src, ek, key in event_map:
            pid = f"{key[0]}:{key[1]}:{key[2]}:{key[3]}" if isinstance(key, tuple) else key
            out.write(f"{src}\t{ek}\t{pid}\n")
    med, n = sample_read_length(args, args.length_sample) if args.length_sample else (None, 0)
    with open(f"{args.output_prefix}.summary.tsv", "wt") as out:
        out.write(f"median_aligned_read_length\t{med}\nreads_sampled\t{n}\npairs\t{len(rows)}\n")
    logger.info("done: %d pairs, median read length %s", len(rows), med)


if __name__ == "__main__":
    main()
