#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute saturation curves for:
  (1) Identifiable (single isoform in 'matching_isoforms')
  (2) FSM subset (sqanti_cat == 'FSM' & single isoform)

Then fit saturation models to each curve and estimate:
  - Predicted Vmax
  - Reads to reach 80–95% of Vmax (x80/x85/x90/x95, from the chosen model)
Also report observed:
  - Observed max (last y)
  - Reads to reach 80–95% of observed max (x80/x85/x90/x95)

Outputs:
  - TSV (thinned curve samples): read_index, cum_unique_identifiable, cum_unique_FSM
  - TSV (fit summary): one row per curve with fit + observed stats
  - PNG overlay: data curves, fitted curves, 80% line (with optional 85/90/95 guides)

Updates in this version:
  * Predicted thresholds are computed consistently via predicted_x_at_fraction()
  * Optional --fit-max-reads / --fit-max-frac to constrain fit range
  * Optional --model {auto,exp,mm} to force a model across runs
"""

import argparse
import logging
import os
import tempfile
import shutil
import time
import gzip
import heapq
from typing import Optional, Dict, Tuple, Set, List

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
    Fractional read support per isoform.
    For each row with k assigned transcripts, each transcript gets +1/k.
    """
    s = df[iso_col].dropna()
    if s.empty:
        return pd.Series(dtype=np.float64)
    splits = s.str.split(",", regex=False)
    lens = splits.str.len().astype("int64")
    exploded = splits.explode().str.strip()
    exploded = exploded[exploded.ne("")]
    if exploded.empty:
        return pd.Series(dtype=np.float64)
    weights = (1.0 / lens).loc[exploded.index]
    counts = weights.groupby(exploded).sum()
    return counts.astype(np.float64)


def allowed_isoforms_from_rpm(
    counts: pd.Series, total_reads: int, min_rpm: float
) -> Optional[Set[str]]:
    if min_rpm is None or min_rpm <= 0 or total_reads <= 0 or counts.empty:
        return None
    thr = (min_rpm * float(total_reads)) / 1_000_000.0
    passed = counts[counts >= thr]
    return set(passed.index.astype(str))


# --------------- Streaming utilities ---------------


def stream_fractional_counts(
    path: str,
    iso_col: str,
    chunksize: int = 5_000_000,
) -> Tuple[pd.Series, int]:
    total = 0
    counts_sum: Optional[pd.Series] = None
    usecols = [iso_col]
    read_iter = pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        dtype={iso_col: "string"},
        chunksize=chunksize,
        engine="c",
        compression="infer",
    )
    for chunk in read_iter:
        total += len(chunk)
        s = (
            chunk[iso_col]
            .astype("string")
            .str.strip()
            .mask(lambda x: x == "", other=pd.NA)
            .dropna()
        )
        if s.empty:
            continue
        splits = s.str.split(",", regex=False)
        lens = splits.str.len().astype("int64")
        exploded = splits.explode().str.strip()
        exploded = exploded[exploded.ne("")]
        if exploded.empty:
            continue
        weights = (1.0 / lens).loc[exploded.index]
        counts_chunk = weights.groupby(exploded).sum()
        if counts_sum is None:
            counts_sum = counts_chunk.astype(np.float64)
        else:
            counts_sum = counts_sum.add(counts_chunk, fill_value=0.0)
    if counts_sum is None:
        counts_sum = pd.Series(dtype=np.float64)
    else:
        counts_sum = counts_sum.astype(np.float64)
    return counts_sum, int(total)


def stream_curves_thin_and_fit(
    path: str,
    iso_col: str,
    cat_col: str,
    fsm_value: str,
    allowed_isos: Optional[Set[str]],
    total_reads: int,
    thin_target: int,
    chunksize: int = 5_000_000,
) -> Tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    int,
    int,
    int,
    int,
]:
    """
    Stream to build thinned curves and produce the (x,y) arrays used for fitting.
    Returns (xt, yi_t, yf_t, xf, yi_f, yf_f, n_ident_reads, n_fsm_reads, y_ident, y_fsm)
    """
    step_curve = max(1, total_reads // max(1, thin_target))
    xt = np.arange(1, total_reads + 1, step_curve, dtype=np.int64)
    xf = xt.astype(np.float64)

    yi_t = np.zeros_like(xt, dtype=np.int64)
    yf_t = np.zeros_like(xt, dtype=np.int64)
    yi_f = np.zeros_like(xf, dtype=np.float64)
    yf_f = np.zeros_like(xf, dtype=np.float64)

    seen_ident: Set[str] = set()
    seen_fsm: Set[str] = set()
    y_ident = 0
    y_fsm = 0
    n_ident_reads = 0
    n_fsm_reads = 0

    next_curve_idx = 0
    next_fit_idx = 0

    usecols = [iso_col, cat_col]
    read_iter = pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        dtype={iso_col: "string", cat_col: "string"},
        chunksize=chunksize,
        engine="c",
        compression="infer",
    )

    i_global = 0
    for chunk in read_iter:
        iso_s = chunk[iso_col].astype("string").str.strip()
        cat_s = chunk[cat_col].astype("string").str.strip()
        for iso_val, cat_val in zip(iso_s, cat_s):
            i_global += 1
            if pd.isna(iso_val) or ("," in str(iso_val)):
                pass
            else:
                iso_token = str(iso_val)
                if allowed_isos is None or iso_token in allowed_isos:
                    n_ident_reads += 1
                    if iso_token not in seen_ident:
                        seen_ident.add(iso_token)
                        y_ident += 1
                    if (not pd.isna(cat_val)) and (str(cat_val) == fsm_value):
                        n_fsm_reads += 1
                        if iso_token not in seen_fsm:
                            seen_fsm.add(iso_token)
                            y_fsm += 1
            while next_curve_idx < len(xt) and i_global >= int(xt[next_curve_idx]):
                yi_t[next_curve_idx] = y_ident
                yf_t[next_curve_idx] = y_fsm
                next_curve_idx += 1
            while next_fit_idx < len(xf) and i_global >= int(xf[next_fit_idx]):
                yi_f[next_fit_idx] = float(y_ident)
                yf_f[next_fit_idx] = float(y_fsm)
                next_fit_idx += 1

    return (
        xt,
        yi_t,
        yf_t,
        xf,
        yi_f,
        yf_f,
        n_ident_reads,
        n_fsm_reads,
        int(y_ident),
        int(y_fsm),
    )


def stream_observed_quantiles(
    path: str,
    iso_col: str,
    cat_col: str,
    fsm_value: str,
    allowed_isos: Optional[Set[str]],
    thresholds_ident: Dict[int, int],
    thresholds_fsm: Dict[int, int],
    chunksize: int = 5_000_000,
) -> Tuple[Dict[int, float], Dict[int, float]]:
    """
    Return dicts percent -> read index (float) where y first crosses frac*observed_max.
    """
    pending_i = {k: True for k in thresholds_ident}
    pending_f = {k: True for k in thresholds_fsm}
    out_i: Dict[int, float] = {k: float("nan") for k in thresholds_ident}
    out_f: Dict[int, float] = {k: float("nan") for k in thresholds_fsm}

    seen_ident: Set[str] = set()
    seen_fsm: Set[str] = set()
    y_ident = 0
    y_fsm = 0

    i_global = 0
    usecols = [iso_col, cat_col]
    read_iter = pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        dtype={iso_col: "string", cat_col: "string"},
        chunksize=chunksize,
        engine="c",
        compression="infer",
    )
    for chunk in read_iter:
        iso_s = chunk[iso_col].astype("string").str.strip()
        cat_s = chunk[cat_col].astype("string").str.strip()
        for iso_val, cat_val in zip(iso_s, cat_s):
            i_global += 1
            if not pd.isna(iso_val) and "," not in str(iso_val):
                iso_token = str(iso_val)
                if allowed_isos is None or iso_token in allowed_isos:
                    if iso_token not in seen_ident:
                        seen_ident.add(iso_token)
                        y_ident += 1
                        for p, thr in thresholds_ident.items():
                            if pending_i.get(p) and y_ident >= thr:
                                out_i[p] = float(i_global)
                                pending_i[p] = False
                    if (not pd.isna(cat_val)) and (str(cat_val) == fsm_value):
                        if iso_token not in seen_fsm:
                            seen_fsm.add(iso_token)
                            y_fsm += 1
                            for p, thr in thresholds_fsm.items():
                                if pending_f.get(p) and y_fsm >= thr:
                                    out_f[p] = float(i_global)
                                    pending_f[p] = False
        if not any(pending_i.values()) and not any(pending_f.values()):
            break
    return out_i, out_f


# --------------- External full shuffle (optional) ---------------


def _open_read_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def _open_write_text(path: str, compresslevel: Optional[int] = None):
    if path.endswith(".gz"):
        if compresslevel is not None:
            return gzip.open(path, "wt", compresslevel=compresslevel)
        return gzip.open(path, "wt")
    return open(path, "wt")


def _write_sorted_run(
    batch_lines: List[str], rng: np.random.Generator, tmp_dir: str, run_paths: List[str]
):
    n = len(batch_lines)
    if n == 0:
        return
    keys = rng.random(n)
    order = np.argsort(keys, kind="mergesort")
    run_path = os.path.join(tmp_dir, f"run_{len(run_paths)}.run")
    with open(run_path, "wt") as f:
        for idx in order:
            key_str = f"{keys[idx]:0.17f}"
            f.write(key_str)
            f.write("\t")
            f.write(batch_lines[idx])
    run_paths.append(run_path)


def _kway_merge_runs(
    run_paths: List[str],
    out_path: str,
    append: bool = False,
    compresslevel: Optional[int] = None,
):
    mode = "at" if append else "wt"
    if out_path.endswith(".gz"):
        fout = (
            gzip.open(out_path, mode, compresslevel=compresslevel)
            if compresslevel is not None
            else gzip.open(out_path, mode)
        )
    else:
        fout = open(out_path, mode)
    files = [open(p, "rt") for p in run_paths]
    try:
        heap = []
        for i, f in enumerate(files):
            line = f.readline()
            if not line:
                continue
            key, row = line.split("\t", 1)
            heap.append((key, i, row))
        heapq.heapify(heap)
        while heap:
            key, i, row = heapq.heappop(heap)
            if out_path.endswith(".run"):
                fout.write(key)
                fout.write("\t")
                fout.write(row)
            else:
                fout.write(row)
            nxt = files[i].readline()
            if nxt:
                k2, r2 = nxt.split("\t", 1)
                heapq.heappush(heap, (k2, i, r2))
    finally:
        for f in files:
            try:
                f.close()
            except Exception:
                pass
        try:
            fout.close()
        except Exception:
            pass


def _kway_merge_runs_to_handle(run_paths: List[str], fout):
    files = [open(p, "rt") for p in run_paths]
    try:
        heap = []
        for i, f in enumerate(files):
            line = f.readline()
            if not line:
                continue
            key, row = line.split("\t", 1)
            heap.append((key, i, row))
        heapq.heapify(heap)
        while heap:
            key, i, row = heapq.heappop(heap)
            fout.write(row)
            nxt = files[i].readline()
            if nxt:
                k2, r2 = nxt.split("\t", 1)
                heapq.heappush(heap, (k2, i, r2))
    finally:
        for f in files:
            try:
                f.close()
            except Exception:
                pass


def external_shuffle_to_file(
    in_path: str,
    out_path: str,
    seed: int,
    lines_per_run: int = 200_000,
    fan_in: int = 64,
    tmp_dir: Optional[str] = None,
    gzip_level: Optional[int] = None,
) -> str:
    logging.info(
        "Shuffle params: seed=%d, lines_per_run=%d, fan_in=%d",
        seed,
        lines_per_run,
        fan_in,
    )
    if tmp_dir is None:
        cwd = os.getcwd()
        tmp_dir = os.path.join(cwd, f"shuffle_runs_{int(time.time())}")
        os.makedirs(tmp_dir, exist_ok=True)
        auto_cleanup = True
    else:
        os.makedirs(tmp_dir, exist_ok=True)
        auto_cleanup = False

    rng = np.random.default_rng(seed)
    run_paths: List[str] = []
    total_lines = 0

    with _open_read_text(in_path) as fin:
        header = fin.readline()
        if header == "":
            with _open_write_text(out_path) as fout:
                pass
            return out_path
        if not header.endswith("\n"):
            header = header + "\n"

        batch_lines: List[str] = []
        for line in fin:
            batch_lines.append(line)
            total_lines += 1
            if len(batch_lines) >= lines_per_run:
                _write_sorted_run(batch_lines, rng, tmp_dir, run_paths)
                batch_lines.clear()
        if batch_lines:
            _write_sorted_run(batch_lines, rng, tmp_dir, run_paths)
            batch_lines.clear()

    cur_runs = run_paths
    level = 0
    while len(cur_runs) > fan_in:
        next_runs: List[str] = []
        for i in range(0, len(cur_runs), fan_in):
            group = cur_runs[i : i + fan_in]
            out_run = os.path.join(tmp_dir, f"merge_level{level}_group{i//fan_in}.run")
            _kway_merge_runs(group, out_run)
            next_runs.append(out_run)
            for rp in group:
                try:
                    os.remove(rp)
                except OSError:
                    pass
        cur_runs = next_runs
        level += 1

    with _open_write_text(out_path, compresslevel=gzip_level) as fout:
        fout.write(header)
        _kway_merge_runs_to_handle(cur_runs, fout)

    for rp in cur_runs:
        try:
            os.remove(rp)
        except OSError:
            pass
    if auto_cleanup:
        try:
            shutil.rmtree(tmp_dir)
        except Exception as e:
            logging.warning("Failed to clean temp dir %s: %s", tmp_dir, e)
    return out_path


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
    return n * np.log(rss / n if rss > 0 else 1e-12) + 2 * k_params


def fit_best_model(
    x: np.ndarray, y: np.ndarray, force: str = "auto"
) -> Dict[str, float]:
    """
    Fit exponential and Michaelis–Menten; pick by AIC unless 'force' is 'exp' or 'mm'.
    Returns: {model, rss, aic, vmax, param1}
      param1 = k (exp) or K (mm)
    """
    n = len(x)
    if n < 5 or np.nanmax(y) <= 0:
        return {
            "model": "none",
            "rss": float("nan"),
            "aic": float("nan"),
            "vmax": float("nan"),
            "param1": float("nan"),
        }

    y_max = float(np.nanmax(y))
    idx_front = max(5, n // 200)
    slope_guess = max(
        1e-9,
        (y[min(n - 1, idx_front)] - y[0]) / (x[min(n - 1, idx_front)] - x[0] + 1e-9),
    )
    k0 = min(1.0, max(1e-9, slope_guess / max(1.0, y_max)))
    p0_exp = (max(y_max, 1.0), k0)

    half = 0.5 * y_max
    if np.any(y >= half):
        j = int(np.argmax(y >= half))
        K0 = max(1.0, x[j])
    else:
        K0 = np.median(x)
    p0_mm = (max(y_max, 1.0), K0)

    ok_e, pexp, yhat_e, rss_e = _safe_curve_fit(_exp_model, x, y, p0_exp)
    ok_m, pmm, yhat_m, rss_m = _safe_curve_fit(_mm_model, x, y, p0_mm)

    aic_e = _aic(n, rss_e, 2) if ok_e else np.inf
    aic_m = _aic(n, rss_m, 2) if ok_m else np.inf

    if force == "exp" and ok_e:
        vmax, k = pexp
        return {
            "model": "exp",
            "rss": rss_e,
            "aic": aic_e,
            "vmax": float(vmax),
            "param1": float(k),
        }
    if force == "mm" and ok_m:
        vmax, K = pmm
        return {
            "model": "mm",
            "rss": rss_m,
            "aic": aic_m,
            "vmax": float(vmax),
            "param1": float(K),
        }

    if aic_e < aic_m:
        vmax, k = pexp if ok_e else (np.nan, np.nan)
        return {
            "model": "exp",
            "rss": rss_e,
            "aic": aic_e,
            "vmax": float(vmax),
            "param1": float(k),
        }
    elif aic_m < np.inf:
        vmax, K = pmm
        return {
            "model": "mm",
            "rss": rss_m,
            "aic": aic_m,
            "vmax": float(vmax),
            "param1": float(K),
        }
    else:
        return {
            "model": "none",
            "rss": float("nan"),
            "aic": float("nan"),
            "vmax": float("nan"),
            "param1": float("nan"),
        }


def predicted_x_at_fraction(fit: Dict[str, float], frac: float) -> float:
    """
    Compute predicted reads to reach a given fraction of Vmax for the chosen model.
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
        description="Saturation analysis for Identifiable vs FSM (fit + thresholds)."
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
        help="Minimum reads-per-million threshold (fractional counting).",
    )
    p.add_argument(
        "--chunksize",
        type=int,
        default=5_000_000,
        help="Rows per chunk when streaming.",
    )
    p.add_argument(
        "--shuffle-buffer", type=int, default=2_000_000, help="(unused in this path)"
    )
    p.add_argument(
        "--no-plot", action="store_true", help="Skip generating the PNG plot."
    )
    p.add_argument(
        "--no-curves", action="store_true", help="Skip writing the thinned curves TSV."
    )
    p.add_argument(
        "--shuffle-out",
        default=None,
        help="Path to write the shuffled file (default: auto).",
    )
    p.add_argument(
        "--skip-shuffle",
        action="store_true",
        help="Skip external shuffling; use input as-is.",
    )
    p.add_argument(
        "--shuffle-gzip-level",
        type=int,
        default=3,
        help="Gzip level for shuffled file.",
    )
    p.add_argument(
        "--shuffle-plain",
        action="store_true",
        help="Write shuffled output uncompressed (.tsv).",
    )
    p.add_argument(
        "--shuffle-run-lines",
        type=int,
        default=1_000_000,
        help="Approx. lines per run.",
    )
    p.add_argument("--shuffle-fan-in", type=int, default=256, help="Run merge fan-in.")
    p.add_argument("--shuffle-tmpdir", default=None, help="Temp dir for shuffle runs.")
    # NEW knobs for comparability
    p.add_argument(
        "--fit-max-reads",
        type=float,
        default=None,
        help="Only use reads <= this value for model fitting (still plots full curve).",
    )
    p.add_argument(
        "--fit-max-frac",
        type=float,
        default=None,
        help="Alternatively, use the first frac (0..1) of reads for fitting.",
    )
    p.add_argument(
        "--model",
        choices=["auto", "exp", "mm"],
        default="auto",
        help="Force model selection across runs; default auto by AIC.",
    )
    p.add_argument("-q", "--quiet", action="store_true", help="Reduce logging.")
    args = p.parse_args()

    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
    )

    t0 = time.perf_counter()
    in_path = args.input

    # External shuffle (unless skipped)
    if not args.skip_shuffle:
        base = os.path.basename(in_path)
        root = base
        for ext in (".tsv.gz", ".txt.gz", ".csv.gz", ".tsv", ".txt", ".csv"):
            if root.endswith(ext):
                root = root[: -len(ext)]
                break
        if args.shuffle_plain:
            shuffle_out = args.shuffle_out or f"{root}.shuffled.tsv"
            gzip_level = None
        else:
            shuffle_out = args.shuffle_out or f"{root}.shuffled.tsv.gz"
            gzip_level = (
                int(args.shuffle_gzip_level)
                if args.shuffle_gzip_level is not None
                else 3
            )
        logging.info("External shuffle starting -> %s", shuffle_out)
        external_shuffle_to_file(
            in_path,
            shuffle_out,
            seed=args.seed,
            lines_per_run=args.shuffle_run_lines,
            fan_in=args.shuffle_fan_in,
            tmp_dir=args.shuffle_tmpdir,
            gzip_level=gzip_level,
        )
        in_path = shuffle_out
        logging.info("External shuffle complete: %s", in_path)
    else:
        logging.info("Skipping external shuffle; using input as-is: %s", in_path)

    # Pass 1: fractional counts + total reads
    frac_counts, n = stream_fractional_counts(
        in_path, args.iso_col, chunksize=args.chunksize
    )
    allowed_isos = allowed_isoforms_from_rpm(frac_counts, n, args.min_rpm)
    n_allowed_iso = (
        len(allowed_isos) if allowed_isos is not None else int(len(frac_counts))
    )
    if args.min_rpm and args.min_rpm > 0:
        logging.info(
            "RPM filter: min_rpm=%.6g; total_reads=%d; allowed_isoforms=%d",
            args.min_rpm,
            n,
            n_allowed_iso,
        )

    # Pass 2: thinned curves + fit arrays
    xt, yi_t, yf_t, xf, yi_f, yf_f, n_ident, n_fsm, obs_max_ident, obs_max_fsm = (
        stream_curves_thin_and_fit(
            in_path,
            args.iso_col,
            args.cat_col,
            args.fsm_value,
            allowed_isos,
            n,
            args.thin,
            chunksize=args.chunksize,
        )
    )

    # Optionally restrict the *fit* domain for stability/comparability
    fit_xmax = None
    if args.fit_max_reads is not None:
        fit_xmax = float(args.fit_max_reads)
    elif args.fit_max_frac is not None:
        fit_xmax = float(args.fit_max_frac) * float(n)
    if fit_xmax is not None and np.isfinite(fit_xmax):
        m = xf <= fit_xmax
        if np.count_nonzero(m) >= 5:
            xf_use, yi_use, yf_use = xf[m], yi_f[m], yf_f[m]
            fit_xmax_used = float(np.max(xf_use))
        else:
            logging.warning("Requested fit range too small; using full domain.")
            xf_use, yi_use, yf_use = xf, yi_f, yf_f
            fit_xmax_used = float(np.max(xf))
    else:
        xf_use, yi_use, yf_use = xf, yi_f, yf_f
        fit_xmax_used = float(np.max(xf))

    # Fit (optionally force model)
    fit_ident = fit_best_model(xf_use, yi_use, force=args.model)
    fit_fsm = fit_best_model(xf_use, yf_use, force=args.model)

    # Pass 3: exact observed thresholds
    fracs = [0.80, 0.85, 0.90, 0.95]
    thr_i = {int(f * 100): int(np.ceil(f * obs_max_ident)) for f in fracs}
    thr_f = {int(f * 100): int(np.ceil(f * obs_max_fsm)) for f in fracs}
    obs_ident, obs_fsm = stream_observed_quantiles(
        in_path,
        args.iso_col,
        args.cat_col,
        args.fsm_value,
        allowed_isos,
        thr_i,
        thr_f,
        chunksize=args.chunksize,
    )

    # Predicted thresholds (ALL via one routine)
    pred_ident = {int(f * 100): predicted_x_at_fraction(fit_ident, f) for f in fracs}
    pred_fsm = {int(f * 100): predicted_x_at_fraction(fit_fsm, f) for f in fracs}

    # Compose outputs
    out_tsv, out_png, out_fit = auto_paths(
        args.input, args.out_tsv, args.png_out, args.fit_out
    )
    rpm_thr_count = (args.min_rpm * float(n)) / 1_000_000.0 if n else 0.0

    fit_df = pd.DataFrame(
        [
            {
                "curve": "Identifiable",
                "total_reads": n,
                "min_rpm": args.min_rpm,
                "rpm_count_threshold": rpm_thr_count,
                "n_isoforms_above_threshold": n_allowed_iso,
                "eligible_read_count": n_ident,
                "observed_max": int(obs_max_ident),
                "observed_x80_reads": obs_ident.get(80, float("nan")),
                "observed_x85_reads": obs_ident.get(85, float("nan")),
                "observed_x90_reads": obs_ident.get(90, float("nan")),
                "observed_x95_reads": obs_ident.get(95, float("nan")),
                "model": fit_ident["model"],
                "vmax_pred": fit_ident["vmax"],
                "param1": fit_ident["param1"],  # k (exp) or K (mm)
                "rss": fit_ident["rss"],
                "aic": fit_ident["aic"],
                "fit_max_reads_used": fit_xmax_used,
                "pred_x80_reads": pred_ident.get(80, float("nan")),
                "pred_x85_reads": pred_ident.get(85, float("nan")),
                "pred_x90_reads": pred_ident.get(90, float("nan")),
                "pred_x95_reads": pred_ident.get(95, float("nan")),
            },
            {
                "curve": "FSM",
                "total_reads": n,
                "min_rpm": args.min_rpm,
                "rpm_count_threshold": rpm_thr_count,
                "n_isoforms_above_threshold": n_allowed_iso,
                "eligible_read_count": n_fsm,
                "observed_max": int(obs_max_fsm),
                "observed_x80_reads": obs_fsm.get(80, float("nan")),
                "observed_x85_reads": obs_fsm.get(85, float("nan")),
                "observed_x90_reads": obs_fsm.get(90, float("nan")),
                "observed_x95_reads": obs_fsm.get(95, float("nan")),
                "model": fit_fsm["model"],
                "vmax_pred": fit_fsm["vmax"],
                "param1": fit_fsm["param1"],
                "rss": fit_fsm["rss"],
                "aic": fit_fsm["aic"],
                "fit_max_reads_used": fit_xmax_used,
                "pred_x80_reads": pred_fsm.get(80, float("nan")),
                "pred_x85_reads": pred_fsm.get(85, float("nan")),
                "pred_x90_reads": pred_fsm.get(90, float("nan")),
                "pred_x95_reads": pred_fsm.get(95, float("nan")),
            },
        ]
    )
    fit_df.to_csv(out_fit, sep="\t", index=False)
    logging.info("Wrote fit summary: %s", out_fit)

    if not args.no_curves:
        thin_df = pd.DataFrame(
            {"read_index": xt, "cum_unique_identifiable": yi_t, "cum_unique_FSM": yf_t}
        )
        thin_df.to_csv(out_tsv, sep="\t", index=False, compression="infer")
        logging.info("Wrote thinned curves: %s  (%d rows)", out_tsv, len(thin_df))

    if not args.no_plot:
        try:
            import matplotlib.pyplot as plt
        except Exception as e:
            logging.warning("matplotlib not available (%s); skipping plot.", e)
        else:
            plt.figure(figsize=(8.6, 5.4), dpi=120)
            # Empirical (thinned)
            plt.plot(xt, yi_t, label="Identifiable (empirical)", linewidth=1.6)
            plt.plot(xt, yf_t, label="FSM (empirical)", linewidth=1.6)

            fracs_extra = [0.85, 0.90, 0.95]

            # Identifiable fit + guides
            if fit_ident["model"] == "exp":
                vmax, k = fit_ident["vmax"], fit_ident["param1"]
                yi_fit = vmax * (1.0 - np.exp(-k * xt))
            elif fit_ident["model"] == "mm":
                vmax, K = fit_ident["vmax"], fit_ident["param1"]
                yi_fit = vmax * (xt / (K + xt))
            else:
                yi_fit = None

            if yi_fit is not None:
                plt.plot(
                    xt,
                    yi_fit,
                    linestyle="--",
                    linewidth=1.2,
                    label=f"Identifiable fit ({fit_ident['model']})",
                )
                y80 = 0.8 * fit_ident["vmax"]
                x80 = predicted_x_at_fraction(fit_ident, 0.80)
                plt.axhline(y80, linestyle=":", linewidth=1.0)
                if np.isfinite(x80):
                    plt.axvline(x80, linestyle=":", linewidth=1.0)
                for f in fracs_extra:
                    y_f = f * fit_ident["vmax"]
                    x_f = predicted_x_at_fraction(fit_ident, f)
                    plt.axhline(
                        y_f, linestyle=":", linewidth=0.8, alpha=0.25, color="gray"
                    )
                    if np.isfinite(x_f):
                        plt.axvline(
                            x_f, linestyle=":", linewidth=0.8, alpha=0.25, color="gray"
                        )

            # FSM fit + guides
            if fit_fsm["model"] == "exp":
                vmax, k = fit_fsm["vmax"], fit_fsm["param1"]
                yf_fit = vmax * (1.0 - np.exp(-k * xt))
            elif fit_fsm["model"] == "mm":
                vmax, K = fit_fsm["vmax"], fit_fsm["param1"]
                yf_fit = vmax * (xt / (K + xt))
            else:
                yf_fit = None

            if yf_fit is not None:
                plt.plot(
                    xt,
                    yf_fit,
                    linestyle="--",
                    linewidth=1.2,
                    label=f"FSM fit ({fit_fsm['model']})",
                )
                y80 = 0.8 * fit_fsm["vmax"]
                x80 = predicted_x_at_fraction(fit_fsm, 0.80)
                plt.axhline(y80, linestyle=":", linewidth=1.0)
                if np.isfinite(x80):
                    plt.axvline(x80, linestyle=":", linewidth=1.0)
                for f in fracs_extra:
                    y_f = f * fit_fsm["vmax"]
                    x_f = predicted_x_at_fraction(fit_fsm, f)
                    plt.axhline(
                        y_f, linestyle=":", linewidth=0.8, alpha=0.25, color="gray"
                    )
                    if np.isfinite(x_f):
                        plt.axvline(
                            x_f, linestyle=":", linewidth=0.8, alpha=0.25, color="gray"
                        )

            plt.xlabel("Number of reads")
            plt.ylabel("Unique transcript discovered")
            plt.title("Saturation: Identifiable vs FSM (fit + thresholds)")
            plt.legend(frameon=False, ncol=2)
            plt.tight_layout()
            plt.savefig(out_png)
            plt.close()
            logging.info(
                "Wrote plot: %s  (elapsed: %.2fs)", out_png, time.perf_counter() - t0
            )


if __name__ == "__main__":
    main()
