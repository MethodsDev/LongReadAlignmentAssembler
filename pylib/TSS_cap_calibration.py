#!/usr/bin/env python3

"""Library-level calibration of the TSS cap signature, before any splice graph.

A read whose 5' end is a capped RNA end carries one more untemplated G than the
library's uncapped reads. Whether a library shows that signature, and its baseline G
run, are properties of the library. Deciding them per splice graph (each contig strand
of each chunk, each cluster's reads) let a chunk dominated by one deep uncapped 5'
truncation pile-up, or a cluster with few reads, switch the cap rules off locally
(v0.46.0: FTD nuclei chr10:10-20 Mb plus strand, VIM; the smallest clusters).

This pass reads the bam once, one worker per contig, and keeps the strands apart. Each
worker returns raw counts: per contig strand, the 5' G runs of spliced informative read
starts away from every TSS peak, and the G runs at each peak. The parent pools them into
one baseline and one decision, which reach every graph build as TSS_cap_signature
("on"/"off") and TSS_cap_G_baseline (an int): Splice_graph then estimates nothing.

Read-level rules are Splice_graph's and Pretty_alignment's: the transcribed strand
(Util_funcs.transcribed_strand), the 5' G run read from the clip before the
untemplated-G strip, a read defines a TSS only if what is left of its 5' clip is within
max_soft_clip_at_TSS (five_prime_clip), and a read whose genomic base just 5' of its start is G carries no
cap information. Peaks are called as a graph build calls them
(aggregate_sites_within_window, filter_non_peaky_positions); only spliced reads enter
the statistics, as in _calibrate_TSS_cap_signature.
"""

import sys, os
import bisect
import json
import logging
import multiprocessing
from collections import defaultdict

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import LRAA_Globals
import Util_funcs

logger = logging.getLogger(__name__)

MAX_G_RUN = 8  # as Splice_graph's TSS_cap_counter


def five_prime_clip(clip_seq, tx_forward):
    """(G run, residual clip length) of a read's 5' soft clip, as Pretty_alignment.

    clip_seq is the clip as stored (genomic frame): the LEFT clip of a + strand read,
    the RIGHT clip of a - strand read, whose G run reads as C's at its head. The G run
    is the run touching the alignment, recorded before stripping. The clip is stripped
    when it is nothing but a short G run (<= max_untemplated_G_at_TSS) or when its
    proximal G run is >= min_proximal_untemplated_G_at_TSS (ONT: adapter + primer + G
    run). test_TSS_cap_calibration pins the agreement with Pretty_alignment.
    """
    cfg = LRAA_Globals.config
    seq = clip_seq.upper()
    clip = len(seq)
    if tx_forward:
        G_run = clip - len(seq.rstrip("G"))
        base = "G"
    else:
        G_run = clip - len(seq.lstrip("C"))
        base = "C"
    max_G = cfg["max_untemplated_G_at_TSS"]
    if max_G > 0 and 0 < clip <= max_G and set(seq) == {base}:
        return G_run, 0
    min_prox = cfg["min_proximal_untemplated_G_at_TSS"]
    if min_prox > 0 and clip > 0 and G_run >= min_prox:
        return G_run, 0
    return G_run, clip


def _five_prime_start(read, strand, seq):
    """(TSS position, G run, spliced) of a read, or None if it defines no TSS.

    G run is None when the read is not informative (genomic G just 5' of the start).
    """
    cig = read.cigartuples
    if not cig or read.query_sequence is None:
        return None
    qseq = read.query_sequence
    forward = strand == "+"
    if forward:
        clip = cig[0][1] if cig[0][0] == 4 else 0
        clip_seq = qseq[:clip]
        pos = read.reference_start + 1
    else:
        clip = cig[-1][1] if cig[-1][0] == 4 else 0
        clip_seq = qseq[len(qseq) - clip :] if clip else ""
        pos = read.reference_end
    G_run, residual = five_prime_clip(clip_seq, forward)
    if residual > LRAA_Globals.config["max_soft_clip_at_TSS"]:
        return None

    spliced = any(op == 3 for op, _ in cig)
    if forward:
        i = pos - 2
        genomic_G = 0 <= i < len(seq) and seq[i] in "Gg"
    else:
        i = pos
        genomic_G = 0 <= i < len(seq) and seq[i] in "Cc"
    return pos, (None if genomic_G else min(G_run, MAX_G_RUN)), spliced


def scan_contig(bam_file, genome_fa, contig, config):
    """One contig, both strands: off-peak G-run counts and per-peak G-run counts.

    Returns {strand: {"off_peak": {run: n}, "peaks": [(pos, {run: n}), ...],
    "n_starts": n}}.
    """
    # imported here so a worker process loads the graph module only when it scans
    from Splice_graph import aggregate_sites_within_window, filter_non_peaky_positions

    LRAA_Globals.config.update(config)
    cfg = LRAA_Globals.config
    seq = Util_funcs.retrieve_contig_seq_from_fasta_file(contig, genome_fa)
    Util_funcs.register_contig_seq_for_strand_check(contig, seq)

    starts = {"+": defaultdict(int), "-": defaultdict(int)}
    cap = {"+": defaultdict(lambda: defaultdict(int)), "-": defaultdict(lambda: defaultdict(int))}
    min_mapq = cfg["min_mapping_quality"]
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(contig):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            strand = Util_funcs.transcribed_strand(read)
            s = _five_prime_start(read, strand, seq)
            if s is None:
                continue
            pos, G_run, spliced = s
            starts[strand][pos] += 1
            if G_run is not None and spliced:
                cap[strand][pos][G_run] += 1

    window = cfg["max_dist_between_alt_TSS_sites"]
    result = {}
    for strand in ("+", "-"):
        counter = starts[strand]
        grouped = aggregate_sites_within_window(
            counter, window, cfg["min_alignments_define_TSS_site"]
        )
        peaks = sorted(
            p for p, _ in filter_non_peaky_positions(grouped, counter, contig, strand)
        )
        positions = sorted(cap[strand])
        off_peak = defaultdict(int)
        for pos in positions:
            k = bisect.bisect_left(peaks, pos - window)
            if not (k < len(peaks) and peaks[k] <= pos + window):
                for run, n in cap[strand][pos].items():
                    off_peak[run] += n
        peak_runs = []
        for site in peaks:
            lo = bisect.bisect_left(positions, site - window)
            hi = bisect.bisect_right(positions, site + window)
            runs = defaultdict(int)
            for pos in positions[lo:hi]:
                for run, n in cap[strand][pos].items():
                    runs[run] += n
            if runs:
                peak_runs.append((site, dict(runs)))
        result[strand] = {
            "off_peak": dict(off_peak),
            "peaks": peak_runs,
            "n_starts": sum(counter.values()),
        }
    return result


def _scan_contig_star(task):
    contig = task[2]
    return contig, scan_contig(*task)


def _split_two_modes(fracs):
    """Otsu split of sorted capped fractions: the threshold maximising the
    between-group variance (grid of 0.01). Returns the threshold, or None."""
    n = len(fracs)
    best, best_t = -1.0, None
    total = sum(fracs)
    for i in range(1, 100):
        t = i / 100.0
        k = bisect.bisect_left(fracs, t)
        if k == 0 or k == n:
            continue
        s_lo = sum(fracs[:k]) if k < n / 2 else total - sum(fracs[k:])
        m0, m1 = s_lo / k, (total - s_lo) / (n - k)
        v = (k / n) * ((n - k) / n) * (m1 - m0) ** 2
        if v > best:
            best, best_t = v, t
    return best_t


HIST_BIN = 0.02


def decide(per_contig, config=None):
    """Pool the contigs' counts into the library's baseline and decision.

    Baseline: the most common 5' G run among spliced informative starts away from
    every peak; capped = more than baseline. The decision reads the per-peak capped
    fractions (peaks with TSS_cap_prepass_min_peak_reads or more informative spliced
    starts). With the signature they are bimodal -- an uncapped mode near the
    transcript-body rate (5' truncation and degradation pile-ups) and a capped mode
    (promoters): MEASURED, medians uncapped / capped: FTD nuclei 1.3 / 48%, PBMC
    1.4 / 65%, BT474 and K562 PacBio 1.6 / 77%, MCF7 and A549 ONT 4.5-6 / 39-43%;
    97% of capped-mode peaks at a FANTOM5 CAGE peak, 14% of uncapped-mode ones.
    Without it (SIRVs, simulated reads, trimmed MORF reads) there is one mode at 0.
    The fractions are split into two groups where they separate best (Otsu); the
    signature is present when the capped group holds at least
    TSS_cap_prepass_min_capped_mode_frac of the peaks and
    TSS_cap_prepass_min_capped_mode_peaks peaks, and its median exceeds the uncapped
    group's by TSS_cap_prepass_min_mode_separation.
    """
    cfg = dict(LRAA_Globals.config)
    if config:
        cfg.update(config)
    off_peak = defaultdict(int)
    for strands in per_contig.values():
        for d in strands.values():
            for run, n in d["off_peak"].items():
                off_peak[int(run)] += n
    n_off = sum(off_peak.values())
    stats = dict(n_contigs=len(per_contig), n_off_peak=n_off)
    if n_off == 0:
        stats.update(decision="undetermined", reason="no informative off-peak read starts")
        return None, False, stats
    baseline = cfg.get("TSS_cap_G_baseline", "auto")
    if baseline in ("auto", None):
        baseline = max(off_peak.items(), key=lambda x: (x[1], -x[0]))[0]
    baseline = int(baseline)
    frac_off = sum(n for run, n in off_peak.items() if run > baseline) / n_off

    min_reads = cfg["TSS_cap_prepass_min_peak_reads"]
    fracs = []
    for strands in per_contig.values():
        for d in strands.values():
            for _site, runs in d["peaks"]:
                n = sum(runs.values())
                if n >= min_reads:
                    fracs.append(
                        sum(w for run, w in runs.items() if int(run) > baseline) / n
                    )
    fracs.sort()
    nbins = int(round(1 / HIST_BIN))
    hist = [0] * nbins
    for x in fracs:
        hist[min(int(x / HIST_BIN), nbins - 1)] += 1
    stats.update(
        baseline_G_run=baseline,
        G_run_off_peak={str(k): off_peak[k] for k in sorted(off_peak)},
        capped_off_peak=round(frac_off, 4),
        peak_min_informative_reads=min_reads,
        peaks_judged=len(fracs),
        histogram_bin=HIST_BIN,
        histogram=hist,
        per_site_threshold=cfg["min_TSS_cap_frac"],
    )
    if len(fracs) < cfg["TSS_cap_prepass_min_peaks"]:
        stats.update(decision="undetermined", reason="too few peaks to judge")
        return baseline, False, stats

    def median(v):
        return v[(len(v) - 1) // 2] if v else None

    split = _split_two_modes(fracs)
    if split is None:
        lo, hi = fracs, []
    else:
        k = bisect.bisect_left(fracs, split)
        lo, hi = fracs[:k], fracs[k:]
    m_lo, m_hi = median(lo), median(hi)
    hi_frac = len(hi) / len(fracs)
    separation = (m_hi - m_lo) if (lo and hi) else 0.0
    # the valley: the emptiest histogram bin between the two medians
    valley = None
    if lo and hi:
        a_, b_ = int(m_lo / HIST_BIN), min(int(m_hi / HIST_BIN), nbins - 1)
        if b_ > a_ + 1:
            v = min(range(a_ + 1, b_), key=lambda i: (hist[i], i))
            valley = round((v + 0.5) * HIST_BIN, 3)
    active = (
        hi_frac >= cfg["TSS_cap_prepass_min_capped_mode_frac"]
        and len(hi) >= cfg["TSS_cap_prepass_min_capped_mode_peaks"]
        and separation >= cfg["TSS_cap_prepass_min_mode_separation"]
    )
    # too few capped-mode peaks with an otherwise clear signature is thin data, not
    # absence: undetermined, so each graph build decides for itself (a 10 Mb test
    # chunk of FTD reads: 29 capped-mode peaks, 56% of those judged)
    thin = (
        not active
        and hi_frac >= cfg["TSS_cap_prepass_min_capped_mode_frac"]
        and separation >= cfg["TSS_cap_prepass_min_mode_separation"]
    )
    reasons = []
    if not active:
        if hi_frac < cfg["TSS_cap_prepass_min_capped_mode_frac"]:
            reasons.append("capped mode {:.1%} of peaks".format(hi_frac))
        if len(hi) < cfg["TSS_cap_prepass_min_capped_mode_peaks"]:
            reasons.append("capped mode {} peaks".format(len(hi)))
        if separation < cfg["TSS_cap_prepass_min_mode_separation"]:
            reasons.append("modes {:.2f} apart".format(separation))
    stats.update(
        mode_split=split,
        uncapped_mode_peaks=len(lo),
        uncapped_mode_median=None if m_lo is None else round(m_lo, 4),
        capped_mode_peaks=len(hi),
        capped_mode_frac=round(hi_frac, 4),
        capped_mode_median=None if m_hi is None else round(m_hi, 4),
        mode_separation=round(separation, 4),
        valley=valley,
        decision="on" if active else ("undetermined" if thin else "off"),
    )
    if reasons:
        stats["reason"] = "; ".join(reasons)
    return baseline, active, stats


SCAN_COLUMNS = ["kind", "contig", "strand", "pos"] + ["G%d" % g for g in range(MAX_G_RUN + 1)]


def write_scan_table(path, per_contig, min_peak_reads=0):
    """The scan's counts as a tsv (gzipped if path ends .gz): one 'off_peak' row per
    contig strand (pos '.') and one 'peak' row per TSS peak, each with its spliced
    informative read starts by 5' G run (G0..G8). decide() needs nothing else, so a
    decision can be re-made, or tested, from this table without the bam. With
    min_peak_reads, peaks with fewer informative reads are left out (they are never
    judged at that floor or above)."""
    import gzip

    op = gzip.open if path.endswith(".gz") else open
    with op(path, "wt") as ofh:
        ofh.write("\t".join(SCAN_COLUMNS) + "\n")
        for contig in sorted(per_contig):
            for strand in ("+", "-"):
                d = per_contig[contig].get(strand)
                if d is None:
                    continue
                off = {int(k): v for k, v in d["off_peak"].items()}
                ofh.write("\t".join(["off_peak", contig, strand, "."] + [str(off.get(g, 0)) for g in range(MAX_G_RUN + 1)]) + "\n")
                for site, runs in d["peaks"]:
                    runs = {int(k): v for k, v in runs.items()}
                    if sum(runs.values()) < min_peak_reads:
                        continue
                    ofh.write("\t".join(["peak", contig, strand, str(site)] + [str(runs.get(g, 0)) for g in range(MAX_G_RUN + 1)]) + "\n")


def read_scan_table(path):
    """write_scan_table's tsv back into scan()'s structure."""
    import gzip

    op = gzip.open if path.endswith(".gz") else open
    per_contig = {}
    with op(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header != SCAN_COLUMNS:
            raise ValueError("{}: not a TSS cap scan table (header {})".format(path, header[:5]))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            kind, contig, strand, pos = f[:4]
            runs = {g: int(n) for g, n in enumerate(f[4:]) if int(n) > 0}
            d = per_contig.setdefault(contig, {}).setdefault(
                strand, {"off_peak": {}, "peaks": [], "n_starts": None}
            )
            if kind == "off_peak":
                for g, n in runs.items():
                    d["off_peak"][g] = d["off_peak"].get(g, 0) + n
            elif kind == "peak":
                d["peaks"].append((int(pos), runs))
            else:
                raise ValueError("{}: unknown row kind {!r}".format(path, kind))
    return per_contig


def write_report(prefix, overrides, stats):
    """<prefix>.TSS_cap_calibration.{json,hist.tsv,svg}: the decision and config keys,
    the per-peak capped-fraction histogram, and its plot (plain SVG, no plotting
    library needed)."""
    out = {}
    j = prefix + ".TSS_cap_calibration.json"
    with open(j, "wt") as ofh:
        json.dump(dict(config_update=overrides, stats=stats), ofh, indent=2, sort_keys=True)
    out["json"] = j
    hist = stats.get("histogram")
    if hist is None:
        return out
    w = stats["histogram_bin"]
    t = prefix + ".TSS_cap_calibration.hist.tsv"
    with open(t, "wt") as ofh:
        ofh.write("capped_frac_lo\tcapped_frac_hi\tpeaks\n")
        for i, n in enumerate(hist):
            ofh.write("{:.2f}\t{:.2f}\t{}\n".format(i * w, (i + 1) * w, n))
    out["tsv"] = t
    out["svg"] = _write_svg(prefix + ".TSS_cap_calibration.svg", stats)
    return out


def _write_svg(path, stats):
    import math

    hist, w = stats["histogram"], stats["histogram_bin"]
    W, H, L, R, T, B = 640, 360, 60, 20, 50, 50
    pw, ph = W - L - R, H - T - B
    ymax = math.log10(max(hist) + 1) if hist and max(hist) > 0 else 1
    X = lambda f: L + f * pw
    Y = lambda n: T + ph - (math.log10(n + 1) / ymax) * ph
    e = []
    e.append('<rect x="0" y="0" width="{}" height="{}" fill="white"/>'.format(W, H))
    bw = pw / len(hist)
    for i, n in enumerate(hist):
        if n:
            y = Y(n)
            e.append(
                '<rect x="{:.1f}" y="{:.1f}" width="{:.1f}" height="{:.1f}" fill="#888"/>'.format(
                    L + i * bw, y, bw - 0.5, T + ph - y
                )
            )
    def vline(f, colour, dash, label, row):
        x = X(f)
        e.append(
            '<line x1="{0:.1f}" y1="{1}" x2="{0:.1f}" y2="{2}" stroke="{3}" stroke-dasharray="{4}"/>'.format(
                x, T, T + ph, colour, dash
            )
        )
        e.append(
            '<text x="{:.1f}" y="{}" font-size="10" fill="{}">{}</text>'.format(
                x + 3, T + 12 + 12 * row, colour, label
            )
        )
    if stats.get("per_site_threshold") is not None:
        vline(stats["per_site_threshold"], "#c00", "4,3", "per-site threshold", 0)
    if stats.get("valley") is not None:
        vline(stats["valley"], "#06c", "2,2", "valley", 1)
    # axes
    e.append('<line x1="{0}" y1="{1}" x2="{2}" y2="{1}" stroke="black"/>'.format(L, T + ph, L + pw))
    e.append('<line x1="{0}" y1="{1}" x2="{0}" y2="{2}" stroke="black"/>'.format(L, T, T + ph))
    for f in (0, 0.2, 0.4, 0.6, 0.8, 1.0):
        e.append('<text x="{:.1f}" y="{}" font-size="10" text-anchor="middle">{:.1f}</text>'.format(X(f), T + ph + 14, f))
    for k in range(0, int(ymax) + 1):
        e.append('<text x="{}" y="{:.1f}" font-size="10" text-anchor="end">{}</text>'.format(L - 4, Y(10 ** k - 1) + 3, 10 ** k))
    e.append('<text x="{}" y="{}" font-size="11" text-anchor="middle">capped fraction of informative spliced read starts at the peak</text>'.format(L + pw / 2, H - 12))
    e.append('<text x="14" y="{}" font-size="11" transform="rotate(-90 14 {})" text-anchor="middle">TSS peaks (log)</text>'.format(T + ph / 2, T + ph / 2))
    title = "TSS cap signature: {} (baseline {} G; {} peaks with >= {} reads)".format(
        stats.get("decision"), stats.get("baseline_G_run"), stats.get("peaks_judged"), stats.get("peak_min_informative_reads")
    )
    parts = []
    if stats.get("capped_mode_median") is not None:
        parts.append(
            "uncapped mode median {:.1%}, capped mode median {:.1%} ({:.0%} of peaks)".format(
                stats["uncapped_mode_median"], stats["capped_mode_median"], stats["capped_mode_frac"]
            )
        )
    if stats.get("reason"):
        parts.append(stats["reason"])
    sub = "; ".join(parts)
    e.append('<text x="{}" y="18" font-size="13" font-weight="bold">{}</text>'.format(L, title))
    e.append('<text x="{}" y="34" font-size="11">{}</text>'.format(L, sub))
    with open(path, "wt") as ofh:
        ofh.write('<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" font-family="sans-serif">\n'.format(W, H))
        ofh.write("\n".join(e))
        ofh.write("\n</svg>\n")
    return path


def scan(bam_file, genome_fa, cpu=1, contigs=None):
    """Scan the bam, one worker per contig: {contig: scan_contig's result}."""
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        if not bam.has_index():
            raise RuntimeError(
                "TSS cap calibration needs an indexed bam: {}".format(bam_file)
            )
        mapped = {s.contig: s.mapped for s in bam.get_index_statistics()}
        names = list(bam.references)
    if contigs:
        names = [c for c in names if c in set(contigs)]
    names = [c for c in names if mapped.get(c, 0) > 0]
    # reads on a contig the genome lacks (BT474/K562 E0 bams: chrEBV) carry no
    # sequence to check: skipped, as the graph builds never visit them either
    with pysam.FastaFile(genome_fa) as fa:
        in_genome = set(fa.references)
    missing = [c for c in names if c not in in_genome]
    if missing:
        logger.warning(
            "TSS cap calibration: skipping %d bam contig(s) absent from the genome: %s",
            len(missing),
            ", ".join(missing[:10]),
        )
        names = [c for c in names if c in in_genome]
    # largest first, so the longest scans start early
    names.sort(key=lambda c: -mapped.get(c, 0))
    config = dict(LRAA_Globals.config)
    tasks = [(bam_file, genome_fa, c, config) for c in names]
    per_contig = {}
    if cpu > 1 and len(tasks) > 1:
        ctx = multiprocessing.get_context("spawn")
        with ctx.Pool(min(cpu, len(tasks))) as pool:
            for contig, res in pool.imap_unordered(_scan_contig_star, tasks):
                per_contig[contig] = res
    else:
        for t in tasks:
            contig, res = _scan_contig_star(t)
            per_contig[contig] = res
    return per_contig


def calibrate(bam_file, genome_fa, cpu=1, contigs=None):
    """Scan the bam (one worker per contig) and decide. Returns
    (config_overrides, stats): the overrides set TSS_cap_signature and
    TSS_cap_G_baseline for every later graph build."""
    per_contig = scan(bam_file, genome_fa, cpu=cpu, contigs=contigs)
    return overrides_from(*decide(per_contig))


def overrides_from(baseline, active, stats):
    """decide()'s result as LRAA config keys, and the stats. Undetermined (too little
    data): no keys, so TSS_cap_signature stays "auto"."""
    overrides = {}
    if stats["decision"] in ("on", "off"):
        overrides["TSS_cap_signature"] = stats["decision"]
        if baseline is not None:
            overrides["TSS_cap_G_baseline"] = baseline
    return overrides, stats


def main(argv=None):
    """Standalone calibration, for a workflow's pre-run step (util/calibrate_TSS_cap_signature.py)."""
    import argparse

    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--bam", help="coordinate-sorted, indexed bam")
    p.add_argument("--genome", help="genome fasta (indexed)")
    p.add_argument(
        "--scan_table",
        help="decide from a saved scan (<prefix>.TSS_cap_calibration.scan.tsv.gz, "
        "written by every bam run) instead of scanning a bam",
    )
    p.add_argument("--cpu", type=int, default=1, help="contigs scanned in parallel")
    p.add_argument(
        "--config_update",
        default=None,
        help="LRAA config overrides (file or inline JSON) applied before scanning, as "
        "LRAA would; with --merge, also carried into the output",
    )
    p.add_argument(
        "--output_prefix",
        required=True,
        help="writes <prefix>.TSS_cap_calibration.{json,hist.tsv,svg,scan.tsv.gz} and "
        "<prefix>.TSS_cap_calibration.config.json, for LRAA --config_update: the "
        "decision's keys (none if undetermined) and TSS_cap_prepass false",
    )
    p.add_argument(
        "--merge",
        action="store_true",
        help="the .config.json holds --config_update's keys plus the calibration's",
    )
    p.add_argument(
        "--contigs", default=None, help="comma-separated contigs to scan (default all)"
    )
    args = p.parse_args(argv)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s"
    )

    user = {}
    if args.config_update:
        v = args.config_update.strip()
        user = json.loads(v) if v.startswith("{") else json.load(open(v))
    unknown = [k for k in user if k not in LRAA_Globals.config]
    if unknown:
        logger.warning("ignoring unknown config keys: %s", ", ".join(unknown))
    LRAA_Globals.config.update({k: v for k, v in user.items() if k not in unknown})

    if args.scan_table:
        if args.bam or args.genome:
            p.error("--scan_table replaces --bam/--genome")
        per_contig = read_scan_table(args.scan_table)
    else:
        if not (args.bam and args.genome):
            p.error("--bam and --genome are required (or --scan_table)")
        per_contig = scan(
            args.bam,
            args.genome,
            cpu=args.cpu,
            contigs=args.contigs.split(",") if args.contigs else None,
        )
        write_scan_table(args.output_prefix + ".TSS_cap_calibration.scan.tsv.gz", per_contig)
    overrides, stats = overrides_from(*decide(per_contig))
    log_summary(overrides, stats)
    write_report(args.output_prefix, overrides, stats)
    out = dict(user) if args.merge else {}
    out.update(overrides)
    # the runs handed this file must not calibrate again on their own subset of the
    # reads (a workflow's shards, a cell cluster): an undetermined library then falls
    # back to the per-graph decision, as in LRAA itself
    out["TSS_cap_prepass"] = False
    with open(args.output_prefix + ".TSS_cap_calibration.config.json", "wt") as ofh:
        json.dump(out, ofh, indent=2, sort_keys=True)


def log_summary(overrides, stats):
    logger.info(
        "TSS cap calibration: %s (config %s) -- baseline G run %s; %s peaks judged; uncapped mode "
        "median %s, capped mode median %s (%s of peaks)%s",
        stats.get("decision"),
        overrides or "unchanged: auto",
        stats.get("baseline_G_run"),
        stats.get("peaks_judged"),
        stats.get("uncapped_mode_median"),
        stats.get("capped_mode_median"),
        stats.get("capped_mode_frac"),
        "; " + stats["reason"] if stats.get("reason") else "",
    )


if __name__ == "__main__":
    main()
