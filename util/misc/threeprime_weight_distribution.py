#!/usr/bin/env python3
"""Measure what the 3'-agreement read weighting actually does, per ambiguity level.

Measurement instrumentation for evaluating config['weight_reads_by_3prime_agreement'].
Reads a LRAA quant.tracking(.gz) produced with the weighting ON and reports, bucketed by
N = the number of transcripts a read is compatible with, how much relative spread the
weights carry and how much assigned read mass sits at each N.

The weighting is

    w_i = 1 - (d_i + P) / (sum_j d_j + N*P)      renormalized to sum to 1

so sum_j(d_j) grows with N while any one (d_i + P) does not. The prediction under test is
that the RELATIVE spread between a read's weights collapses as N grows, i.e. the weighting
discriminates only where ambiguity is mild and goes uniform exactly where it is worst.

Emits a TSV to stdout. Nothing here changes LRAA behaviour.
"""

import argparse
import gzip
import sys
from collections import defaultdict


BUCKETS = [(1, 1), (2, 2), (3, 3), (4, 5), (6, 10), (11, 10**9)]
BUCKET_LABELS = ["1", "2", "3", "4-5", "6-10", ">10"]

# The tracking column is printed at 6 decimal places, so weights that differ only in the
# 7th place are indistinguishable in this file and are treated as equal rather than as
# evidence of discrimination.
SPREAD_EPS = 1e-5


def bucket_of(n):
    for i, (lo, hi) in enumerate(BUCKETS):
        if lo <= n <= hi:
            return i
    raise AssertionError(f"unbucketable N={n}")


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def load_read_groups(path):
    """read_name -> (list of weights, list of frac_assigned).

    A read lives in exactly one multipath, so its rows are exactly its compatible
    transcript set. Grouping on the read name rather than the mp id keeps the unit of
    analysis the thing the question is about: one read, and the weights EM will use to
    apportion it.
    """
    weights = defaultdict(list)
    fracs = defaultdict(list)
    mp_of_read = {}
    multi_mp_reads = 0
    with opener(path) as fh:
        header_seen = False
        for line in fh:
            if line.startswith("#"):
                continue
            if not header_seen:
                header_seen = True
                if line.startswith("gene_id\t"):
                    continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            mp_id, readname, frac, w = f[4], f[5], f[6], f[7]
            prev = mp_of_read.setdefault(readname, mp_id)
            if prev != mp_id:
                multi_mp_reads += 1
            weights[readname].append(float(w))
            fracs[readname].append(float(frac))
    return weights, fracs, multi_mp_reads


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tracking", help="LRAA *.quant.tracking or .gz, weighting ON")
    ap.add_argument("--label", default="", help="sample label for the output rows")
    args = ap.parse_args()

    weights, fracs, multi_mp_reads = load_read_groups(args.tracking)
    if not weights:
        sys.exit(f"Error, no tracking rows parsed from {args.tracking}")

    n_b = [0] * len(BUCKETS)
    mass_b = [0.0] * len(BUCKETS)
    ratios_b = [[] for _ in BUCKETS]
    cvs_b = [[] for _ in BUCKETS]
    wmax_b = [[] for _ in BUCKETS]
    wmin_b = [[] for _ in BUCKETS]
    nonuniform_reads_b = [0] * len(BUCKETS)
    nonuniform_mass_b = [0.0] * len(BUCKETS)
    strong_reads_b = [0] * len(BUCKETS)
    strong_mass_b = [0.0] * len(BUCKETS)

    for readname, ws in weights.items():
        n = len(ws)
        b = bucket_of(n)
        n_b[b] += 1
        mass = sum(fracs[readname])
        mass_b[b] += mass
        lo, hi = min(ws), max(ws)
        mean = sum(ws) / n
        # max/min is undefined against a zero floor; report it only where it is defined
        # and count the degenerate case separately rather than silently dropping it.
        if lo > 0:
            ratios_b[b].append(hi / lo)
        if mean > 0:
            var = sum((w - mean) ** 2 for w in ws) / n
            cvs_b[b].append((var ** 0.5) / mean)
        wmax_b[b].append(hi)
        wmin_b[b].append(lo)
        # A read whose compatible transcripts are all equidistant from its 3' end gets
        # identical weights, so the weighting is a no-op for it no matter what N is.
        # Counting those separately distinguishes "flattened by large N" from "never had
        # anything to say", which are different defects with different fixes.
        if n > 1 and hi > lo * (1.0 + SPREAD_EPS):
            nonuniform_reads_b[b] += 1
            nonuniform_mass_b[b] += mass
            if hi > lo * 1.10:
                strong_reads_b[b] += 1
                strong_mass_b[b] += mass

    total_reads = sum(n_b)
    total_mass = sum(mass_b)

    def med(v):
        if not v:
            return float("nan")
        s = sorted(v)
        m = len(s) // 2
        return s[m] if len(s) % 2 else 0.5 * (s[m - 1] + s[m])

    def p90(v):
        if not v:
            return float("nan")
        s = sorted(v)
        return s[min(len(s) - 1, int(0.9 * len(s)))]

    def pct(num, den):
        return f"{100.0 * num / den:.3f}" if den else "nan"

    cols = ["sample", "N_bucket", "n_reads", "pct_reads", "assigned_mass",
            "pct_mass", "median_wmax_over_wmin", "p90_wmax_over_wmin",
            "median_CV_of_weights", "median_wmax", "median_wmin",
            "pct_reads_nonuniform", "pct_of_total_mass_nonuniform",
            "pct_reads_spread_gt10pct", "pct_of_total_mass_spread_gt10pct"]
    print("\t".join(cols))
    for i, lab in enumerate(BUCKET_LABELS):
        print("\t".join([
            args.label, lab, f"{n_b[i]}",
            pct(n_b[i], total_reads),
            f"{mass_b[i]:.2f}",
            pct(mass_b[i], total_mass),
            f"{med(ratios_b[i]):.4f}", f"{p90(ratios_b[i]):.4f}",
            f"{med(cvs_b[i]):.5f}", f"{med(wmax_b[i]):.5f}", f"{med(wmin_b[i]):.5f}",
            pct(nonuniform_reads_b[i], n_b[i]),
            pct(nonuniform_mass_b[i], total_mass),
            pct(strong_reads_b[i], n_b[i]),
            pct(strong_mass_b[i], total_mass),
        ]))
    print("\t".join([
        args.label, "ALL", f"{total_reads}", "100.000", f"{total_mass:.2f}", "100.000",
        "", "", "", "", "",
        pct(sum(nonuniform_reads_b), total_reads),
        pct(sum(nonuniform_mass_b), total_mass),
        pct(sum(strong_reads_b), total_reads),
        pct(sum(strong_mass_b), total_mass),
    ]))
    print(f"# total_reads={total_reads} total_assigned_mass={total_mass:.2f} "
          f"reads_spanning_multiple_mps={multi_mp_reads}", file=sys.stderr)


if __name__ == "__main__":
    main()
