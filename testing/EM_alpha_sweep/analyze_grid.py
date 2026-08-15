#!/usr/bin/env python3
"""Paired analysis of the EM_alpha x 3'-weighting grid.

Everything here is PAIRED BY SAMPLE.  LRAA quant is deterministic, so there is
no replicate noise to average over, but sample-to-sample variation in MARD is
two orders of magnitude larger than the effects being chased -- a mean MARD
compared across samples between two configs would be dominated by which samples
happen to be in the corpus.  So every comparison is a within-sample delta
against one reference arm, and the headline number is how many samples moved
the same way, not the mean.

Emits, under <BASE>/results/:
  paired_deltas.tsv      per (sample, arm) delta vs the reference arm
  response_curves.tsv    per (stratum, 3p, alpha) paired summary + sign test
  per_sample_optima.tsv  per (sample, 3p) argmin alpha for each metric
  interaction.tsv        alpha optimum with 3p on vs off, per sample
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep_config as C  # noqa: E402

# Lower is better for these; higher is better for the rest.
LOWER_BETTER = {"mard", "rmse", "nrmse_mean_truth", "n_pred_zero", "n_pred_lt_1pct", "n_pred_lt_10pct"}
METRICS = ["mard", "spearman", "pearson_log1p", "ccc_log1p", "nrmse_mean_truth",
           "n_pred_zero", "n_pred_lt_1pct", "n_pred_lt_10pct"]


def load(path):
    df = pd.read_csv(path, sep="\t")
    for c in METRICS + ["alpha"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def balanced(df):
    """Drop (alpha, 3p) cells not measured on every sample IN THIS FRAME.

    Applied per stratum, not once corpus-wide.  MORF was cut to five alpha
    points for wall-clock and SIRV was extended to eleven, so a single global
    balance would throw away the extension everywhere -- including from the
    SIRV stratum that paid for it.  Balancing inside each stratum keeps every
    paired comparison over a complete set of samples while letting different
    strata cover different alpha ranges.
    """
    per_alpha = df.groupby(["alpha", "threep"])["sample"].nunique()
    full = df["sample"].nunique()
    keep = set(per_alpha[per_alpha == full].index)
    return df[[(a, w) in keep for a, w in zip(df["alpha"], df["threep"])]].copy()


def paired_deltas(df, ref_alpha, ref_threep):
    ref = (
        df[(df.alpha == ref_alpha) & (df.threep == ref_threep)]
        .set_index("sample")[METRICS]
        .add_suffix("_ref")
    )
    out = df.join(ref, on="sample")
    for m in METRICS:
        out[f"d_{m}"] = out[m] - out[f"{m}_ref"]
    return out.drop(columns=[f"{m}_ref" for m in METRICS])


def sign_summary(vals, lower_better):
    """(n_better, n_worse, n_tie, sign-test p) for a vector of deltas vs reference."""
    v = np.asarray(vals, dtype=float)
    v = v[~np.isnan(v)]
    better = int((v < 0).sum()) if lower_better else int((v > 0).sum())
    worse = int((v > 0).sum()) if lower_better else int((v < 0).sum())
    ties = int((v == 0).sum())
    n = better + worse
    p = binomtest(better, n, 0.5).pvalue if n else float("nan")
    return better, worse, ties, p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics", default=os.path.join(C.BASE, "results", "grid_metrics.tsv"))
    ap.add_argument("--outdir", default=os.path.join(C.BASE, "results"))
    ap.add_argument("--ref_alpha", type=float, default=C.DEFAULT_ALPHA)
    ap.add_argument("--ref_threep", type=int, default=1)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    raw = load(args.metrics)

    # SIRV E1 and E2 are split out because they are not two samples of one
    # thing.  Within an E-level the truth vector is bit-identical across all
    # four cell lines, so the eight SIRV libraries carry exactly two distinct
    # truths -- E1 with an 8-fold dynamic range (CV 0.68) and E2 with 128-fold
    # (CV 1.37).  Pooling them hides a disagreement that turns out to be the
    # main result.
    def sub(mask):
        return balanced(raw[mask].copy())

    strata = {
        "pooled": balanced(raw.copy()),
        "HiFi": sub(raw.platform == "HiFi"),
        "ONT": sub(raw.platform == "ONT"),
        "SIRV": sub(raw.corpus == "SIRV"),
        "SIRV_E1": sub(raw["sample"].str.contains("_E1_")),
        "SIRV_E2": sub(raw["sample"].str.contains("_E2_")),
        "MORF": sub(raw.corpus == "MORF"),
    }
    df = strata["pooled"]

    # ---- paired deltas vs the reference arm ------------------------------
    dd = paired_deltas(raw, args.ref_alpha, args.ref_threep)
    keep = ["sample", "corpus", "platform", "arm", "alpha", "threep"] + METRICS + [f"d_{m}" for m in METRICS]
    dd[keep].sort_values(["corpus", "sample", "threep", "alpha"], ascending=[True, True, False, True]).to_csv(
        os.path.join(args.outdir, "paired_deltas.tsv"), sep="\t", index=False
    )

    # ---- response curves --------------------------------------------------
    rows = []
    for sname, sdf in strata.items():
        if sdf.empty:
            continue
        sdd = paired_deltas(sdf, args.ref_alpha, args.ref_threep)
        nsamp = sdf["sample"].nunique()
        for w in sorted(sdf.threep.unique(), reverse=True):
            for a in sorted(sdf.alpha.unique()):
                cell = sdd[(sdd.alpha == a) & (sdd.threep == w)]
                if cell.empty:
                    continue
                r = dict(stratum=sname, n_samples=nsamp, threep=int(w), alpha=a)
                for m in METRICS:
                    r[f"mean_{m}"] = cell[m].mean()
                    r[f"meand_{m}"] = cell[f"d_{m}"].mean()
                    b, wo, t, p = sign_summary(cell[f"d_{m}"], m in LOWER_BETTER)
                    r[f"better_{m}"] = b
                    r[f"worse_{m}"] = wo
                    r[f"tie_{m}"] = t
                    r[f"p_{m}"] = p
                # Wilcoxon signed-rank on MARD deltas, when there is anything to rank.
                v = cell["d_mard"].dropna().values
                r["wilcoxon_p_mard"] = (
                    wilcoxon(v).pvalue if len(v) >= 5 and np.any(v != 0) else float("nan")
                )
                rows.append(r)
    curves = pd.DataFrame(rows)
    curves.to_csv(os.path.join(args.outdir, "response_curves.tsv"), sep="\t", index=False)

    # ---- per-sample optima -------------------------------------------------
    rows = []
    for (smp, w), g in raw.groupby(["sample", "threep"]):
        r = dict(sample=smp, corpus=g.corpus.iloc[0], platform=g.platform.iloc[0], threep=int(w))
        for m in METRICS:
            gg = g.dropna(subset=[m])
            if gg.empty or gg[m].nunique() == 1:
                r[f"best_alpha_{m}"] = float("nan")
                continue
            idx = gg[m].idxmin() if m in LOWER_BETTER else gg[m].idxmax()
            r[f"best_alpha_{m}"] = gg.loc[idx, "alpha"]
            r[f"best_{m}"] = gg.loc[idx, m]
        rows.append(r)
    optima = pd.DataFrame(rows).sort_values(["corpus", "sample", "threep"], ascending=[True, True, False])
    optima.to_csv(os.path.join(args.outdir, "per_sample_optima.tsv"), sep="\t", index=False)

    # ---- alpha x 3p interaction --------------------------------------------
    rows = []
    for smp, g in raw.groupby("sample"):
        on, off = g[g.threep == 1], g[g.threep == 0]
        if on.empty or off.empty:
            continue
        a_on = on.loc[on.mard.idxmin(), "alpha"]
        a_off = off.loc[off.mard.idxmin(), "alpha"]
        # Interaction contrast: does moving alpha from default to its optimum
        # buy the same amount with 3p on as with 3p off?
        d_on = on.mard.min() - float(on[on.alpha == args.ref_alpha].mard.iloc[0])
        d_off = off.mard.min() - float(off[off.alpha == args.ref_alpha].mard.iloc[0])
        rows.append(dict(
            sample=smp, corpus=g.corpus.iloc[0], platform=g.platform.iloc[0],
            best_alpha_3p_on=a_on, best_alpha_3p_off=a_off,
            mard_best_on=on.mard.min(), mard_best_off=off.mard.min(),
            gain_from_alpha_on=d_on, gain_from_alpha_off=d_off,
            interaction=d_off - d_on,
            mard_3p_effect_at_default=(
                float(off[off.alpha == args.ref_alpha].mard.iloc[0])
                - float(on[on.alpha == args.ref_alpha].mard.iloc[0])
            ),
        ))
    inter = pd.DataFrame(rows)
    inter.to_csv(os.path.join(args.outdir, "interaction.tsv"), sep="\t", index=False)

    # ---- console summary ---------------------------------------------------
    pd.set_option("display.width", 200, "display.max_columns", 60)
    for sname in ("pooled", "HiFi", "ONT", "SIRV_E1", "SIRV_E2", "MORF"):
        c = curves[curves.stratum == sname]
        if c.empty:
            continue
        print(f"\n===== {sname}  (n={c.n_samples.iloc[0]} samples) =====")
        show = c[["threep", "alpha", "mean_mard", "meand_mard", "better_mard", "worse_mard",
                  "p_mard", "mean_spearman", "mean_n_pred_lt_1pct"]]
        print(show.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\n===== per-sample MARD optima =====")
    print(optima[["sample", "platform", "threep", "best_alpha_mard", "best_mard",
                  "best_alpha_spearman"]].to_string(index=False))
    print("\n===== alpha x 3p interaction =====")
    print(inter.to_string(index=False, float_format=lambda x: f"{x:.6g}"))


if __name__ == "__main__":
    main()
