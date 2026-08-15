#!/usr/bin/env python3
"""Does a sample's sensitivity to EM_alpha scale with its ambiguous read mass?

The prediction under test. alpha enters the M-step as
`base_alpha * ambiguous_read_counts`, so it can only move mass that some read
left ambiguous. A library with twice the ambiguous mass should therefore move
about twice as far when alpha is swept, and a library with none should not move
at all. If that fails, alpha is doing something other than what the mechanism
says, and no fitted optimum is trustworthy.

The covariate is measured, not modelled. Per-sample ambiguous-mass fractions
come from Quant3Prime's read-level extraction of LRAA's own quant.tracking
(v0.20.0, quant-only, cpu_budget 3, alpha 0.01, 3' weighting on; mass == read
count, unnormalized bams, 1.0 per read). Two columns are supplied:

  all_ambig    every read compatible with >1 transcript. This is the
               mechanistically faithful regressor: EM.py:290 counts every
               multipath with len>1 into ambiguous_read_counts regardless of
               weight, so this is exactly alpha's reach.
  irreducible  the subset whose candidate transcripts share a 3' end exactly
               and differ by internal splicing. Such a read cannot span the
               junction that separates its candidates, so no per-read evidence
               of any design can resolve it and alpha decides it unopposed.
               The interpretively interesting column, not the faithful one.

Three constraints on what a positive result licenses, all of them real:

  * The two regressors are COLLINEAR here -- BT474 is the maximum on both
    within E2, UHRR is the maximum on both within E1 -- so a correlation
    cannot attribute effect size to the irreducible fraction specifically.
  * n=4 per E-stratum. Hitting the predicted maximum by chance is p~0.25 per
    stratum, ~0.06 for both. Reportable, not established.
  * Truth flatness is a confound and it differs between E1 and E2 but not
    within, so correlations are computed WITHIN E-level. The pooled figure is
    printed for completeness and should not be read as controlling for it.

Sensitivity is measured over a common alpha set across the samples being
compared, because a spread taken over different grids is not a comparable
number.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep_config as C  # noqa: E402

# Quant3Prime, per-sample, % of total library mass. morf2_pacbio was still
# running when this table was handed over; it is absent rather than imputed.
AMBIG = {
    #  sample                irreducible  all_ambig
    "CL_BT474_E1_sirv":     (3.529, 6.276),
    "CL_HG002_E1_sirv":     (4.004, 7.207),
    "CL_K562_E1_sirv":      (3.780, 6.808),
    "CL_UHRR_E1_sirv":      (5.890, 9.811),
    "CL_BT474_E2_sirv":     (7.580, 13.439),
    "CL_HG002_E2_sirv":     (3.808, 8.135),
    "CL_K562_E2_sirv":      (3.823, 7.786),
    "CL_UHRR_E2_sirv":      (3.755, 7.578),
    "morf2_ont":            (8.754, 11.353),
}

PREDICTED_MAX = {"E1": "CL_UHRR_E1_sirv", "E2": "CL_BT474_E2_sirv"}


# Central difference for d MARD / d log10(alpha) at the default.  These two
# alphas bracket 0.01 by exactly one decade each way.
DERIV_LO, DERIV_HI, DERIV_AT = 0.003, 0.03, 0.01


def sensitivity(df, alphas, threep=1):
    """Per-sample alpha sensitivity over a fixed alpha set.

    Three measures, because they are not interchangeable:

    spread   max - min MARD across the set.  Intuitive, but it is NOT a clean
             sensitivity measure when different samples put their optimum at
             opposite ends of the grid: for a sample optimal at alpha=0 it
             measures how bad the top of the range gets, and for a sample
             optimal at the top it measures how much headroom remains.  Those
             are different quantities wearing one name, which is why the
             spread-versus-mass relationship comes out superlinear in one
             stratum and flat in the other.
    signed   MARD(max alpha) - MARD(min alpha).  Direction as well as size.
    dmard    |d MARD / d log10 alpha| at the shipped default, by central
             difference over 0.003 and 0.03.  Well defined for every sample
             wherever its optimum sits, including off-grid, and it is the
             textbook definition of local sensitivity.  This is the measure to
             regress ambiguous mass against.
    """
    rows = []
    for smp, g in df[df.threep == threep].groupby("sample"):
        g = g[g.alpha.isin(alphas)]
        if len(g) != len(alphas):
            continue
        m = g.set_index("alpha")["mard"]
        have_deriv = DERIV_LO in m.index and DERIV_HI in m.index
        rows.append(dict(
            sample=smp, corpus=g.corpus.iloc[0], platform=g.platform.iloc[0],
            elevel=("E1" if "_E1_" in smp else "E2" if "_E2_" in smp else "-"),
            mard_spread=float(m.max() - m.min()),
            mard_signed=float(m.loc[max(alphas)] - m.loc[min(alphas)]),
            dmard_dlogalpha=(
                abs(float(m.loc[DERIV_HI] - m.loc[DERIV_LO]))
                / (np.log10(DERIV_HI) - np.log10(DERIV_LO))
                if have_deriv else np.nan
            ),
            best_alpha=float(m.idxmin()),
        ))
    out = pd.DataFrame(rows)
    out["irreducible"] = out["sample"].map(lambda s: AMBIG.get(s, (np.nan,) * 2)[0])
    out["all_ambig"] = out["sample"].map(lambda s: AMBIG.get(s, (np.nan,) * 2)[1])
    return out


def correlate(sub, xcol, ycol):
    s = sub.dropna(subset=[xcol, ycol])
    if len(s) < 3:
        return dict(n=len(s), pearson_r=np.nan, pearson_p=np.nan,
                    spearman_r=np.nan, spearman_p=np.nan)
    pr = pearsonr(s[xcol], s[ycol])
    sr = spearmanr(s[xcol], s[ycol])
    return dict(n=len(s), pearson_r=pr.statistic, pearson_p=pr.pvalue,
                spearman_r=sr.statistic, spearman_p=sr.pvalue)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics", default=os.path.join(C.BASE, "results", "grid_metrics.tsv"))
    ap.add_argument("--outdir", default=os.path.join(C.BASE, "results"))
    ap.add_argument("--threep", type=int, default=1)
    args = ap.parse_args()

    df = pd.read_csv(args.metrics, sep="\t")
    for c in ("alpha", "mard"):
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Common alpha set: the widest set present on every sample that has a
    # covariate. Comparing spreads taken over different grids would compare
    # grid coverage, not sensitivity.
    have = df[df["sample"].isin(AMBIG)]
    counts = have[have.threep == args.threep].groupby("alpha")["sample"].nunique()
    nsamp = have[have.threep == args.threep]["sample"].nunique()
    alphas = sorted(counts[counts == nsamp].index)
    if not alphas:
        print("no alpha point is present on every covariate sample yet")
        return 0

    sens = sensitivity(have, alphas, args.threep)
    sens = sens.sort_values(["elevel", "sample"])
    os.makedirs(args.outdir, exist_ok=True)
    sens.to_csv(os.path.join(args.outdir, "ambiguity_sensitivity.tsv"), sep="\t", index=False)

    pd.set_option("display.width", 200, "display.max_columns", 40)
    print(f"alpha set used ({len(alphas)} points): {alphas}")
    print(f"3' weighting: {'on' if args.threep else 'off'}\n")
    print(sens.to_string(index=False, float_format=lambda x: f"{x:.6g}"))

    rows = []
    # irreducible first: it is the PRE-REGISTERED, theoretically motivated
    # regressor.  all_ambig is reported second as the mechanistically faithful
    # but merely correlated alternative.  dmard_dlogalpha first for the same
    # reason -- it is the measure that is comparable across samples whose
    # optima sit at opposite ends of the grid.
    for ycol in ("dmard_dlogalpha", "mard_spread"):
        for xcol in ("irreducible", "all_ambig"):
            for label, sub in (
                ("E1 only", sens[sens.elevel == "E1"]),
                ("E2 only", sens[sens.elevel == "E2"]),
                ("SIRV pooled (flatness NOT controlled)", sens[sens.corpus == "SIRV"]),
                ("all samples (flatness NOT controlled)", sens),
            ):
                r = dict(sensitivity=ycol, regressor=xcol, stratum=label)
                r.update(correlate(sub, xcol, ycol))
                rows.append(r)
    corr = pd.DataFrame(rows)
    corr.to_csv(os.path.join(args.outdir, "ambiguity_correlation.tsv"), sep="\t", index=False)
    print("\n===== correlation of alpha sensitivity with ambiguous mass =====")
    print(corr.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    print("\n===== pre-registered predicted maxima =====")
    for e, want in PREDICTED_MAX.items():
        sub = sens[sens.elevel == e]
        if sub.empty:
            print(f"  {e}: not yet measured")
            continue
        got = sub.loc[sub.mard_spread.idxmax(), "sample"]
        n = len(sub)
        print(f"  {e}: predicted {want}, observed {got} -> "
              f"{'HIT' if got == want else 'MISS'} (n={n}, chance 1/{n})")
    sirv = sens[sens.corpus == "SIRV"]
    ont = sens[sens["sample"] == "morf2_ont"]
    if not ont.empty and not sirv.empty:
        v = float(ont.mard_spread.iloc[0])
        mx = float(sirv.mard_spread.max())
        print(f"  morf2_ont (8.754% irreducible) spread {v:.6g} vs max SIRV {mx:.6g} -> "
              f"{'HIT' if v > mx else 'MISS'} ({len(sirv)} SIRV samples compared)")
    else:
        print("  morf2_ont: not yet measured")
    return 0


if __name__ == "__main__":
    sys.exit(main())
