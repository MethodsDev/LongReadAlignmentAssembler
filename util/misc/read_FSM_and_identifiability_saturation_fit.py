#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute saturation curves for:
  (1) Identifiable (single isoform in 'matching_isoforms')
  (2) FSM subset (sqanti_cat == 'FSM' & single isoform)

Then fit saturation models to each curve and estimate:
  - Predicted Vmax
    - Reads to reach 80–95% of Vmax (step 5; x80_pred, x85_pred, x90_pred, x95_pred)
Also report observed:
  - Observed max (last y)
    - Reads to reach 80–95% of observed max (x80_obs, x85_obs, x90_obs, x95_obs)

Outputs:
  - TSV (thinned curve samples): read_index, cum_unique_identifiable, cum_unique_FSM
  - TSV (fit summary): one row per curve with fit + observed stats
  - PNG overlay: data curves, fitted curves, 80% lines

Author: you + ChatGPT
"""

import argparse
import logging
import os
import time
from typing import Optional, Dict, Tuple, Set

import numpy as np
import pandas as pd


# --------------- IO ---------------


def read_columns(
    path: str,
    iso_col: str = "matching_isoforms",
    cat_col: str = "sqanti_cat",
    nrows: Optional[int] = None,
) -> pd.DataFrame:
    engine = "pyarrow"
    usecols = [iso_col, cat_col]
    try:
        pd.read_csv(path, sep="\t", usecols=usecols, nrows=0, engine=engine)
    except Exception:
        engine = "c"

    df = pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        dtype={iso_col: "string", cat_col: "string"},
        engine=engine,
        nrows=nrows,
        compression="infer",
    ).copy()

    df[iso_col] = df[iso_col].str.strip().mask(lambda s: s == "", other=pd.NA)
    df[cat_col] = df[cat_col].str.strip().mask(lambda s: s == "", other=pd.NA)
    return df


def auto_paths(
    in_path: str, out_tsv: Optional[str], out_png: Optional[str], out_fit: Optional[str]
):
    base = os.path.basename(in_path)
    root = base
    for ext in (".tsv.gz", ".txt.gz", ".csv.gz", ".tsv", ".txt", ".csv"):
        if root.endswith(ext):
            root = root[: -len(ext)]
            break
    if out_tsv is None:
        out_tsv = f"{root}.saturation_identifiable_vs_FSM.thin.tsv.gz"
    if out_png is None:
        out_png = f"{root}.saturation_identifiable_vs_FSM.fit.png"
    if out_fit is None:
        out_fit = f"{root}.saturation_identifiable_vs_FSM.fit_summary.tsv"
    return out_tsv, out_png, out_fit


# --------------- Curves ---------------


def cumulative_unique_curve(ids: pd.Series, eligible_mask: np.ndarray) -> np.ndarray:
    n = len(ids)
    y_full = np.zeros(n, dtype=np.int64)
    if n == 0 or not eligible_mask.any():
        return y_full
    idx = np.nonzero(eligible_mask)[0]
    ids_sub = ids.iloc[idx]
    is_new = ~ids_sub.duplicated(keep="first")
    y_full[idx] = is_new.cumsum().to_numpy(dtype=np.int64)
    return np.maximum.accumulate(y_full)


def build_curves(
    df: pd.DataFrame,
    iso_col: str,
    cat_col: str,
    fsm_value: str,
    allowed_isos: Optional[Set[str]] = None,
) -> Dict[str, np.ndarray]:
    iso = df[iso_col]
    cat = df[cat_col]

    mask_ident = iso.notna() & ~iso.str.contains(",", na=False)
    if allowed_isos is not None:
        # Only consider identifiable reads whose isoform passes the RPM threshold
        mask_ident = mask_ident & iso.isin(allowed_isos)
    mask_fsm = mask_ident & (cat == fsm_value)

    y_ident = cumulative_unique_curve(iso, mask_ident.to_numpy())
    y_fsm = cumulative_unique_curve(iso, mask_fsm.to_numpy())

    return {
        "ident": y_ident,
        "fsm": y_fsm,
        "mask_counts": (int(mask_ident.sum()), int(mask_fsm.sum())),
    }


# --------------- Expression thresholding (RPM with fractional counting) ---------------


def compute_fractional_iso_counts(df: pd.DataFrame, iso_col: str) -> pd.Series:
    """
    Compute fractional read support per isoform.
    For each row with k assigned transcripts, each transcript gets +1/k.
    Returns a Series indexed by isoform ID with float counts.
    """
    s = df[iso_col].dropna()
    if s.empty:
        return pd.Series(dtype=np.float64)
    splits = s.str.split(",", regex=False)
    lens = splits.str.len().astype("int64")
    exploded = splits.explode()
    # Strip whitespace from tokens and drop empty tokens
    exploded = exploded.str.strip()
    exploded = exploded[exploded.ne("")]
    if exploded.empty:
        return pd.Series(dtype=np.float64)
    weights_per_row = 1.0 / lens
    weights = weights_per_row.loc[exploded.index]
    counts = weights.groupby(exploded).sum()
    # Ensure float64 dtype
    return counts.astype(np.float64)


def allowed_isoforms_from_rpm(
    counts: pd.Series, total_reads: int, min_rpm: float
) -> Optional[Set[str]]:
    """
    Given fractional counts per isoform and total reads, return isoforms with
    counts >= min_rpm * total_reads / 1e6.
    """
    if min_rpm is None or min_rpm <= 0 or total_reads <= 0 or counts.empty:
        # No filtering to apply
        return None
    thr = (min_rpm * float(total_reads)) / 1_000_000.0
    passed = counts[counts >= thr]
    return set(passed.index.astype(str))


# --------------- Fitting ---------------


def _exp_model(x, vmax, k):
    # y = vmax * (1 - exp(-k x))
    return vmax * (1.0 - np.exp(-k * x))


def _mm_model(x, vmax, K):
    # y = vmax * x / (K + x)
    return vmax * (x / (K + x))


def _safe_curve_fit(model, x, y, p0):
    from scipy.optimize import curve_fit

    try:
        popt, pcov = curve_fit(model, x, y, p0=p0, maxfev=20000)
        yhat = model(x, *popt)
        rss = float(np.sum((y - yhat) ** 2))
        return True, popt, yhat, rss
    except Exception as e:
        logging.debug("curve_fit failed: %s", e)
        return False, None, None, np.inf


def _aic(n, rss, k_params):
    # Gaussian log-likelihood approximation
    return n * np.log(rss / n if rss > 0 else 1e-12) + 2 * k_params


def fit_best_model(x: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    """
    Try exponential and Michaelis–Menten models; pick by AIC.
    Returns dict with model name, params, rss, aic, vmax, x80_pred.
    """
    n = len(x)
    if n < 5 or np.nanmax(y) <= 0:
        return {
            "model": "none",
            "rss": float("nan"),
            "aic": float("nan"),
            "vmax": float(np.nan),
            "param1": float("nan"),
            "x80_pred": float(np.nan),
        }

    # Initial guesses
    y_max = float(np.nanmax(y))
    # For exp: k ~ initial slope / vmax; rough slope from first few points
    idx_front = max(5, n // 200)  # up to ~0.5%
    slope_guess = max(
        1e-9,
        (y[min(n - 1, idx_front)] - y[0]) / (x[min(n - 1, idx_front)] - x[0] + 1e-9),
    )
    k0 = min(1.0, max(1e-9, slope_guess / max(1.0, y_max)))
    p0_exp = (max(y_max, 1.0), k0)

    # For MM: K ~ x at half of ymax
    half = 0.5 * y_max
    # Find closest index where y crosses half
    if np.any(y >= half):
        j = int(np.argmax(y >= half))
        K0 = max(1.0, x[j])
    else:
        K0 = np.median(x)
    p0_mm = (max(y_max, 1.0), K0)

    # Fit both
    ok_e, pexp, yhat_e, rss_e = _safe_curve_fit(_exp_model, x, y, p0_exp)
    ok_m, pmm, yhat_m, rss_m = _safe_curve_fit(_mm_model, x, y, p0_mm)

    # Score
    aic_e = _aic(n, rss_e, 2) if ok_e else np.inf
    aic_m = _aic(n, rss_m, 2) if ok_m else np.inf

    if aic_e < aic_m:
        vmax, k = pexp if ok_e else (np.nan, np.nan)
        x80 = (np.log(5.0) / k) if ok_e and k > 0 else np.nan  # ln(5) ~ 1.609
        return {
            "model": "exp",
            "rss": rss_e,
            "aic": aic_e,
            "vmax": float(vmax),
            "param1": float(k),
            "x80_pred": float(x80),
        }
    elif aic_m < np.inf:
        vmax, K = pmm
        x80 = 4.0 * K
        return {
            "model": "mm",
            "rss": rss_m,
            "aic": aic_m,
            "vmax": float(vmax),
            "param1": float(K),
            "x80_pred": float(x80),
        }
    else:
        return {
            "model": "none",
            "rss": float("nan"),
            "aic": float("nan"),
            "vmax": float("nan"),
            "param1": float("nan"),
            "x80_pred": float("nan"),
        }


def observed_x80(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """
    Return (x80_obs, y80_obs) where y reaches 80% of observed max (first crossing).
    """
    if len(y) == 0:
        return float("nan"), float("nan")
    y_obs_max = float(np.nanmax(y))
    thr = 0.8 * y_obs_max
    idx = np.argmax(y >= thr)
    if y[idx] < thr:
        return float("nan"), float("nan")
    return float(x[idx]), float(thr)


def observed_x_at_fraction(x: np.ndarray, y: np.ndarray, frac: float) -> float:
    """
    First x where y >= frac * observed_max. Returns NaN if never reached.
    """
    if len(y) == 0 or not (0.0 < frac < 1.0):
        return float("nan")
    y_obs_max = float(np.nanmax(y))
    thr = frac * y_obs_max
    idx = int(np.argmax(y >= thr))
    if y[idx] < thr:
        return float("nan")
    return float(x[idx])


def predicted_x_at_fraction(fit: Dict[str, float], frac: float) -> float:
    """
    Compute predicted reads to reach a given fraction of Vmax based on chosen model.
    """
    model = fit.get("model", "none")
    if not (0.0 < frac < 1.0):
        return float("nan")
    if model == "exp":
        k = fit.get("param1", float("nan"))
        if k and k > 0:
            return float(-np.log(1.0 - frac) / k)
        return float("nan")
    if model == "mm":
        K = fit.get("param1", float("nan"))
        if np.isfinite(K):
            return float((frac * K) / (1.0 - frac))
        return float("nan")
    return float("nan")


# --------------- Main ---------------


def main():
    p = argparse.ArgumentParser(
        description="Saturation analysis for Identifiable vs FSM (fit + 80% reads)."
    )
    p.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input TSV(.gz) with 'matching_isoforms' and 'sqanti_cat'.",
    )
    p.add_argument(
        "--iso-col", default="matching_isoforms", help="Column with isoform IDs."
    )
    p.add_argument(
        "--cat-col", default="sqanti_cat", help="Column with SQANTI category."
    )
    p.add_argument(
        "--fsm-value",
        default="FSM",
        help="Value in sqanti_cat indicating full splice match.",
    )
    p.add_argument("--seed", type=int, default=42, help="Shuffle seed.")
    p.add_argument("--no-shuffle", action="store_true", help="Disable shuffling.")
    p.add_argument(
        "--thin",
        type=int,
        default=10_000,
        help="Target points kept for curve TSV/plot.",
    )
    p.add_argument(
        "--limit-rows", type=int, default=0, help="Process only first N rows (0=all)."
    )
    p.add_argument(
        "-o", "--out-tsv", default=None, help="Output TSV for thinned curves."
    )
    p.add_argument("--png-out", default=None, help="Output PNG plot.")
    p.add_argument("--fit-out", default=None, help="Output TSV summarizing fits.")
    p.add_argument(
        "--min-rpm",
        type=float,
        default=0.0,
        help=
        "Minimum expression threshold in reads-per-million total reads."
        " Isoforms must meet this threshold (with fractional counting) to be included",
    )
    p.add_argument("-q", "--quiet", action="store_true", help="Reduce logging.")
    args = p.parse_args()

    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
    )

    nrows = args.limit_rows if args.limit_rows and args.limit_rows > 0 else None

    t0 = time.perf_counter()
    df = read_columns(args.input, args.iso_col, args.cat_col, nrows=nrows)
    n = len(df)
    logging.info("Loaded %d rows.", n)

    if not args.no_shuffle and n > 0:
        rng = np.random.default_rng(args.seed)
        df = df.take(rng.permutation(n)).reset_index(drop=True)
        logging.info("Shuffled with seed=%d", args.seed)
    else:
        df = df.reset_index(drop=True)

    # Fractional counting pass and RPM filtering
    frac_counts = compute_fractional_iso_counts(df, args.iso_col)
    allowed_isos = allowed_isoforms_from_rpm(frac_counts, n, args.min_rpm)
    n_allowed_iso = (
        len(allowed_isos)
        if allowed_isos is not None
        else int(len(frac_counts))
    )
    if args.min_rpm and args.min_rpm > 0:
        logging.info(
            "RPM filter: min_rpm=%.6g; total_reads=%d; allowed_isoforms=%d",
            args.min_rpm,
            n,
            n_allowed_iso,
        )

    curves = build_curves(
        df, args.iso_col, args.cat_col, args.fsm_value, allowed_isos=allowed_isos
    )
    y_ident, y_fsm = curves["ident"], curves["fsm"]
    n_ident, n_fsm = curves["mask_counts"]

    # Build x
    x = np.arange(1, n + 1, dtype=np.float64)

    # Fit on a thinned grid for speed/robustness, but evaluate stats on full vectors
    # (We fit using the thinned subset but compute observed x80 on the full curve.)
    fit_step = max(1, n // max(1, args.thin))
    xf = x[::fit_step]
    yi_f = y_ident[::fit_step].astype(np.float64)
    yf_f = y_fsm[::fit_step].astype(np.float64)

    fit_ident = fit_best_model(xf, yi_f)
    fit_fsm = fit_best_model(xf, yf_f)

    # Observed stats (full curves)
    x80i_obs, y80i_obs = observed_x80(x, y_ident)
    x80f_obs, y80f_obs = observed_x80(x, y_fsm)

    # Additional quantiles (80 to 95 step 5)
    fracs = [0.80, 0.85, 0.90, 0.95]
    obs_ident = {int(f * 100): observed_x_at_fraction(x, y_ident, f) for f in fracs}
    obs_fsm = {int(f * 100): observed_x_at_fraction(x, y_fsm, f) for f in fracs}
    pred_ident = {int(f * 100): predicted_x_at_fraction(fit_ident, f) for f in fracs}
    pred_fsm = {int(f * 100): predicted_x_at_fraction(fit_fsm, f) for f in fracs}

    # Compose fit summary TSV
    out_tsv, out_png, out_fit = auto_paths(
        args.input, args.out_tsv, args.png_out, args.fit_out
    )
    fit_df = pd.DataFrame(
        [
            {
                "curve": "Identifiable",
                "eligible_read_count": n_ident,
                "observed_max": int(y_ident[-1]) if n else 0,
                "observed_x80_reads": x80i_obs,
                # Additional observed quantiles
                "observed_x85_reads": obs_ident.get(85, float("nan")),
                "observed_x90_reads": obs_ident.get(90, float("nan")),
                "observed_x95_reads": obs_ident.get(95, float("nan")),
                "model": fit_ident["model"],
                "vmax_pred": fit_ident["vmax"],
                "param1": fit_ident["param1"],  # k for exp; K for MM
                "rss": fit_ident["rss"],
                "aic": fit_ident["aic"],
                "pred_x80_reads": fit_ident["x80_pred"],
                # Additional predicted quantiles
                "pred_x85_reads": pred_ident.get(85, float("nan")),
                "pred_x90_reads": pred_ident.get(90, float("nan")),
                "pred_x95_reads": pred_ident.get(95, float("nan")),
            },
            {
                "curve": "FSM",
                "eligible_read_count": n_fsm,
                "observed_max": int(y_fsm[-1]) if n else 0,
                "observed_x80_reads": x80f_obs,
                "observed_x85_reads": obs_fsm.get(85, float("nan")),
                "observed_x90_reads": obs_fsm.get(90, float("nan")),
                "observed_x95_reads": obs_fsm.get(95, float("nan")),
                "model": fit_fsm["model"],
                "vmax_pred": fit_fsm["vmax"],
                "param1": fit_fsm["param1"],
                "rss": fit_fsm["rss"],
                "aic": fit_fsm["aic"],
                "pred_x80_reads": fit_fsm["x80_pred"],
                "pred_x85_reads": pred_fsm.get(85, float("nan")),
                "pred_x90_reads": pred_fsm.get(90, float("nan")),
                "pred_x95_reads": pred_fsm.get(95, float("nan")),
            },
        ]
    )
    # Add metadata columns (same values for both rows)
    rpm_thr_count = (args.min_rpm * float(n)) / 1_000_000.0 if n else 0.0
    fit_df.insert(1, "total_reads", n)
    fit_df.insert(2, "min_rpm", args.min_rpm)
    fit_df.insert(3, "rpm_count_threshold", rpm_thr_count)
    fit_df.insert(4, "n_isoforms_above_threshold", n_allowed_iso)
    fit_df.to_csv(out_fit, sep="\t", index=False)
    logging.info("Wrote fit summary: %s", out_fit)

    # Thinned curve TSV for plotting elsewhere
    step_curve = max(1, n // max(1, args.thin))
    xt = np.arange(1, n + 1, step_curve, dtype=np.int64)
    yi_t = y_ident[::step_curve]
    yf_t = y_fsm[::step_curve]
    thin_df = pd.DataFrame(
        {
            "read_index": xt,
            "cum_unique_identifiable": yi_t,
            "cum_unique_FSM": yf_t,
        }
    )
    thin_df.to_csv(out_tsv, sep="\t", index=False, compression="infer")
    logging.info("Wrote thinned curves: %s  (%d rows)", out_tsv, len(thin_df))

    # Plot overlay with fitted curves + 80% lines and subtle extra quantiles (85/90/95)
    import matplotlib.pyplot as plt

    plt.figure(figsize=(8.6, 5.4), dpi=120)
    # Empirical (thinned)
    plt.plot(xt, yi_t, label="Identifiable (empirical)", linewidth=1.6)
    plt.plot(xt, yf_t, label="FSM (empirical)", linewidth=1.6)

    # Fitted curves (evaluated on xt for display)
    # Extra fractions beyond 80% for subtle guides
    fracs_extra = [0.85, 0.90, 0.95]

    if fit_ident["model"] == "exp":
        vmax, k = fit_ident["vmax"], fit_ident["param1"]
        yi_fit = vmax * (1.0 - np.exp(-k * xt))
        y80 = 0.8 * vmax
        x80 = fit_ident["x80_pred"]
        plt.plot(
            xt,
            yi_fit,
            linestyle="--",
            linewidth=1.2,
            label=f"Identifiable fit ({fit_ident['model']})",
        )
        plt.axhline(y80, linestyle=":", linewidth=1.0)
        if np.isfinite(x80):
            plt.axvline(x80, linestyle=":", linewidth=1.0)
        # Subtle extra quantile guides
        for f in fracs_extra:
            y_f = f * vmax
            x_f = predicted_x_at_fraction(fit_ident, f)
            plt.axhline(y_f, linestyle=":", linewidth=0.8, alpha=0.25, color="gray")
            if np.isfinite(x_f):
                plt.axvline(x_f, linestyle=":", linewidth=0.8, alpha=0.25, color="gray")
    elif fit_ident["model"] == "mm":
        vmax, K = fit_ident["vmax"], fit_ident["param1"]
        yi_fit = vmax * (xt / (K + xt))
        y80 = 0.8 * vmax
        x80 = fit_ident["x80_pred"]
        plt.plot(
            xt,
            yi_fit,
            linestyle="--",
            linewidth=1.2,
            label=f"Identifiable fit ({fit_ident['model']})",
        )
        plt.axhline(y80, linestyle=":", linewidth=1.0)
        if np.isfinite(x80):
            plt.axvline(x80, linestyle=":", linewidth=1.0)
        for f in fracs_extra:
            y_f = f * vmax
            x_f = predicted_x_at_fraction(fit_ident, f)
            plt.axhline(y_f, linestyle=":", linewidth=0.8, alpha=0.25, color="gray")
            if np.isfinite(x_f):
                plt.axvline(x_f, linestyle=":", linewidth=0.8, alpha=0.25, color="gray")

    if fit_fsm["model"] == "exp":
        vmax, k = fit_fsm["vmax"], fit_fsm["param1"]
        yf_fit = vmax * (1.0 - np.exp(-k * xt))
        y80 = 0.8 * vmax
        x80 = fit_fsm["x80_pred"]
        plt.plot(
            xt,
            yf_fit,
            linestyle="--",
            linewidth=1.2,
            label=f"FSM fit ({fit_fsm['model']})",
        )
        plt.axhline(y80, linestyle=":", linewidth=1.0)
        if np.isfinite(x80):
            plt.axvline(x80, linestyle=":", linewidth=1.0)
        for f in fracs_extra:
            y_f = f * vmax
            x_f = predicted_x_at_fraction(fit_fsm, f)
            plt.axhline(y_f, linestyle=":", linewidth=0.8, alpha=0.25, color="gray")
            if np.isfinite(x_f):
                plt.axvline(x_f, linestyle=":", linewidth=0.8, alpha=0.25, color="gray")
    elif fit_fsm["model"] == "mm":
        vmax, K = fit_fsm["vmax"], fit_fsm["param1"]
        yf_fit = vmax * (xt / (K + xt))
        y80 = 0.8 * vmax
        x80 = fit_fsm["x80_pred"]
        plt.plot(
            xt,
            yf_fit,
            linestyle="--",
            linewidth=1.2,
            label=f"FSM fit ({fit_fsm['model']})",
        )
        plt.axhline(y80, linestyle=":", linewidth=1.0)
        if np.isfinite(x80):
            plt.axvline(x80, linestyle=":", linewidth=1.0)
        for f in fracs_extra:
            y_f = f * vmax
            x_f = predicted_x_at_fraction(fit_fsm, f)
            plt.axhline(y_f, linestyle=":", linewidth=0.8, alpha=0.25, color="gray")
            if np.isfinite(x_f):
                plt.axvline(x_f, linestyle=":", linewidth=0.8, alpha=0.25, color="gray")

    plt.xlabel("Number of reads")
    plt.ylabel("Unique transcript discovered")
    plt.title("Saturation: Identifiable vs FSM (fit + 80% thresholds)")
    plt.legend(frameon=False, ncol=2)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()
    logging.info("Wrote plot: %s  (elapsed: %.2fs)", out_png, time.perf_counter() - t0)


if __name__ == "__main__":
    main()
