#!/usr/bin/env python3
"""How much expression actually moves between two LRAA quant.expr files.

Measurement instrumentation. A benchmark metric answers "did the answer get better";
this answers the prior question "did the answer change at all", against no truth set and
with no scorer in the loop. Useful when a config is suspected of being inert: if total
variation between two arms is a fraction of a percent, no accuracy metric computed from
them can be reporting a real effect either.

Compares the read-count column the benchmark scorer consumes (0-indexed col 3,
all_reads), renormalized per file to sum to 1, and reports the total variation distance

    TVD = 0.5 * sum_t |p_t - q_t|

which is exactly the fraction of total expression mass that changed hands.
"""

import argparse
import sys


def load(path, id_col=1, val_col=3):
    vals = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(id_col, val_col) or f[id_col] == "transcript_id":
                continue
            try:
                vals[f[id_col]] = float(f[val_col])
            except ValueError:
                continue
    total = sum(vals.values())
    if total > 0:
        vals = {k: v / total for k, v in vals.items()}
    return vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--label", default="")
    ap.add_argument("--header", action="store_true")
    args = ap.parse_args()

    A, B = load(args.a), load(args.b)
    keys = set(A) | set(B)
    if not keys:
        sys.exit(f"Error, no rows parsed from {args.a} / {args.b}")

    diffs = [(k, A.get(k, 0.0), B.get(k, 0.0)) for k in keys]
    tvd = 0.5 * sum(abs(a - b) for _, a, b in diffs)
    n_moved = sum(1 for _, a, b in diffs if a != b)
    # Largest single reallocation, as a share of the library.
    k_max, a_max, b_max = max(diffs, key=lambda t: abs(t[1] - t[2]))
    # Sign flips between zero and nonzero are what turn into false negatives/positives.
    a_zero_b_not = sum(1 for _, a, b in diffs if a == 0 and b > 0)
    b_zero_a_not = sum(1 for _, a, b in diffs if b == 0 and a > 0)

    if args.header:
        print("label\tn_transcripts\tTVD_frac_mass_moved\tn_changed\tmax_single_shift"
              "\tmax_shift_transcript\tzero_in_A_only\tzero_in_B_only")
    print(f"{args.label}\t{len(keys)}\t{tvd:.6e}\t{n_moved}\t{abs(a_max - b_max):.6e}"
          f"\t{k_max}\t{a_zero_b_not}\t{b_zero_a_not}")


if __name__ == "__main__":
    main()
