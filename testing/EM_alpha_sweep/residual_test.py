#!/usr/bin/env python3
"""Test alpha's mechanism on the slice it actually controls.

alpha does not move the whole library. Reads compatible with exactly one
transcript land there whatever theta says, and they are 88-92% of assigned mass.
What alpha interpolates is the split of the remaining AMBIGUOUS mass, between

    likelihood split   the allocation EM reaches on its own at alpha=0
    profile split      allocation proportional to each transcript's ambiguous
                       support, which is where theta goes as alpha -> infinity

So the claim to test is not "is the ambiguity profile close to the truth" -- it
is not, it is farther from truth than a flat prior on every sample -- but "which
of those two splits is closer to the ambiguous mass truth actually allots". That
is what should predict the sign of the alpha effect, and it is computable
without a sweep.

The target, in COUNT space throughout so no effective-length normalization can
creep into it:

    u_t         unique read mass on t          (quant.expr uniq_reads)
    N           total assigned mass            (sum of all_reads)
    T_t         truth mass on t, rescaled to N
    residual_t  T_t - u_t, the ambiguous mass truth says belongs to t

residual_t CAN BE NEGATIVE, where unique reads already exceed truth's whole
allotment for that transcript. Those are not clipped. A negative residual is not
something alpha can fix -- it means the unique-read assignment itself disagrees
with truth -- so their count and total is reported as a diagnostic of the corpus
rather than folded into the comparison.

Distance is L1 over total ambiguous mass, i.e. the fraction of ambiguous mass a
split misallocates. It is well defined with negative targets, which the
symmetric-relative-difference form used for MARD is not.
"""

import argparse
import csv
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep_config as C  # noqa: E402
from ambiguity_profile import PROBE, load_truth, profile_from_tracking  # noqa: E402


def read_expr(path):
    """(uniq_reads, all_reads) per transcript from a quant.expr."""
    uniq, allr = {}, {}
    with open(path) as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
            uniq[r["transcript_id"]] = float(r["uniq_reads"])
            allr[r["transcript_id"]] = float(r["all_reads"])
    return uniq, allr


def expr_path(sample, arm):
    return os.path.join(
        C.WORK, sample, arm, f"{sample}.LRAA.{arm}.LRAA.quant-only.quant.expr"
    )


def misallocation(split, target, total):
    """Fraction of ambiguous mass placed on the wrong transcript."""
    return float(np.abs(split - target).sum() / (2.0 * total)) if total else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples")
    ap.add_argument("--out", default=os.path.join(C.BASE, "results", "residual_test.tsv"))
    args = ap.parse_args()

    samples = (
        [C.SAMPLES_BY_NAME[s] for s in args.samples.split(",")]
        if args.samples else list(C.SAMPLES)
    )

    rows = []
    for smp in samples:
        s = smp["sample"]
        trk = os.path.join(PROBE, s, "prof.LRAA.quant-only.quant.tracking.gz")
        p0 = expr_path(s, C.arm_name(0.0, True))
        if not (os.path.exists(trk) and os.path.exists(p0)):
            continue
        amb, _n_reads, _n_amb = profile_from_tracking(trk)
        uniq, all0 = read_expr(p0)
        truth = load_truth(smp["truth_quant"])

        keys = sorted(set(uniq) | set(truth) | set(amb))
        u = np.array([uniq.get(k, 0.0) for k in keys])
        a = np.array([amb.get(k, 0.0) for k in keys])
        a0 = np.array([all0.get(k, 0.0) for k in keys])
        t = np.array([truth.get(k, 0.0) for k in keys])
        N = a0.sum()
        t = t / t.sum() * N

        residual = t - u
        A = N - u.sum()                       # total ambiguous mass, by construction
        split_like = a0 - u                   # what EM does with it at alpha=0
        split_prof = a / a.sum() * A if a.sum() else np.zeros_like(a)

        # Where a limit run exists, the measured alpha->infinity split, which is
        # the profile split only up to EM's own redistribution at that theta.
        limit_arm = C.arm_name(C.LIMIT_ALPHA, True)
        split_limit = None
        if os.path.exists(expr_path(s, limit_arm)):
            _uL, allL = read_expr(expr_path(s, limit_arm))
            split_limit = np.array([allL.get(k, 0.0) for k in keys]) - u

        neg = residual < 0
        r = dict(
            sample=s, corpus=smp["corpus"], platform=smp["platform"],
            n_tx=len(keys), total_assigned=N, ambiguous_mass=A,
            frac_ambiguous=A / N,
            n_residual_negative=int(neg.sum()),
            mass_residual_negative=float(-residual[neg].sum()),
            frac_amb_mass_negative=float(-residual[neg].sum() / A) if A else float("nan"),
            misalloc_likelihood=misallocation(split_like, residual, A),
            misalloc_profile=misallocation(split_prof, residual, A),
            spearman_likelihood=float(spearmanr(split_like, residual).statistic),
            spearman_profile=float(spearmanr(split_prof, residual).statistic),
        )
        if split_limit is not None:
            r["misalloc_limit_run"] = misallocation(split_limit, residual, A)
            r["spearman_limit_run"] = float(spearmanr(split_limit, residual).statistic)
        r["profile_beats_likelihood"] = r["misalloc_profile"] < r["misalloc_likelihood"]
        rows.append(r)

    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    pd.set_option("display.width", 250, "display.max_columns", 40)
    print(out.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
