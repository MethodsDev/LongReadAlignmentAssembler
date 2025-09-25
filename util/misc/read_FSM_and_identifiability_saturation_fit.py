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


# --------------- Streaming utilities for low-memory operation ---------------


def stream_fractional_counts(
    path: str,
    iso_col: str,
    chunksize: int = 5_000_000,
) -> Tuple[pd.Series, int]:
    """
    Stream the file to compute fractional counts per isoform and total reads.
    Returns (counts_series, total_rows).
    """
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
        s = chunk[iso_col].astype("string").str.strip().mask(lambda x: x == "", other=pd.NA).dropna()
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
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, int, int, int]:
    """
    Stream to build thinned curves and perform fit on the same thinned grid.
    Returns (xt, yi_t, yf_t, xf, yi_f, yf_f, n_ident_reads, n_fsm_reads) where:
      - xt: x positions for thinned curve TSV
      - yi_t, yf_t: y values (int) for Identifiable and FSM at xt
      - xf: x positions for fitting
      - yi_f, yf_f: y values (float) for fitting
      - n_ident_reads, n_fsm_reads: counts of eligible reads per curve
    """
    # Determine grids
    step_curve = max(1, total_reads // max(1, thin_target))
    xt = np.arange(1, total_reads + 1, step_curve, dtype=np.int64)
    fit_step = max(1, total_reads // max(1, thin_target))
    xf = xt.astype(np.float64)

    yi_t = np.zeros_like(xt, dtype=np.int64)
    yf_t = np.zeros_like(xt, dtype=np.int64)
    yi_f = np.zeros_like(xf, dtype=np.float64)
    yf_f = np.zeros_like(xf, dtype=np.float64)

    # Streaming counters
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
            # Identifiable if single isoform token (no comma) and not NA/empty
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
            # Record thinned points
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
    Stream to find first read index crossing each absolute threshold.
    thresholds_* map percent (e.g., 80) -> y threshold count (absolute).
    Returns dicts mapping percent -> read index (float) or NaN if never crossed.
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


# --------------- External disk-based full randomization (random-key sort) ---------------


def _open_read_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def _open_write_text(path: str):
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "wt")


def _write_sorted_run(batch_lines: List[str], rng: np.random.Generator, tmp_dir: str, run_paths: List[str]):
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


def _kway_merge_runs(run_paths: List[str], out_path: str, append: bool = False):
    mode = "at" if append else "wt"
    # Use gzip when writing to a gzip file to avoid corrupting output
    fout = gzip.open(out_path, mode) if out_path.endswith(".gz") else open(out_path, mode)
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
        try:
            fout.close()
        except Exception:
            pass


def external_shuffle_to_file(
    in_path: str,
    out_path: str,
    seed: int,
    lines_per_run: int = 200_000,
    fan_in: int = 64,
    tmp_dir: Optional[str] = None,
) -> str:
    if tmp_dir is None:
        tdir_obj = tempfile.TemporaryDirectory(prefix="shuffle_runs_")
        tmp_dir = tdir_obj.name
        auto_cleanup = True
    else:
        os.makedirs(tmp_dir, exist_ok=True)
        tdir_obj = None
        auto_cleanup = False

    rng = np.random.default_rng(seed)

    run_paths: List[str] = []
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

    with _open_write_text(out_path) as fout:
        fout.write(header)
    _kway_merge_runs(cur_runs, out_path, append=True)
    for rp in cur_runs:
        try:
            os.remove(rp)
        except OSError:
            pass

    if auto_cleanup and tdir_obj is not None:
        tdir_obj.cleanup()
    return out_path


def randomized_row_iter(
    path: str,
    iso_col: str,
    cat_col: str,
    chunksize: int,
    buffer_max: int,
    seed: int,
):
    """
    Yield (iso, cat) rows in approximately random order using a bounded shuffle buffer.
    Memory usage is O(buffer_max) rows.
    """
    rng = np.random.default_rng(seed)
    buffer = []  # list of tuples (iso, cat)
    usecols = [iso_col, cat_col]
    reader = pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        dtype={iso_col: "string", cat_col: "string"},
        chunksize=chunksize,
        engine="c",
        compression="infer",
    )
    # Helper to fill buffer up to capacity
    def extend_buffer_from_chunk(df_chunk: pd.DataFrame):
        nonlocal buffer
        if df_chunk is None or df_chunk.empty:
            return
        iso_s = df_chunk[iso_col].astype("string").str.strip()
        cat_s = df_chunk[cat_col].astype("string").str.strip()
        # Append tuples; avoid creating large intermediate DataFrame
        for iso_val, cat_val in zip(iso_s, cat_s):
            buffer.append((iso_val, cat_val))

    chunk_iter = iter(reader)
    next_chunk = None

    # Prime buffer
    while len(buffer) < buffer_max:
        try:
            next_chunk = next(chunk_iter)
        except StopIteration:
            next_chunk = None
            break
        extend_buffer_from_chunk(next_chunk)
        next_chunk = None

    # Main loop: pop random items, refill from stream
    while buffer:
        # Refill if below capacity and chunks remain
        while len(buffer) < buffer_max:
            try:
                next_chunk = next(chunk_iter)
            except StopIteration:
                next_chunk = None
                break
            extend_buffer_from_chunk(next_chunk)
            next_chunk = None
        # Randomly pop one
        j = int(rng.integers(len(buffer)))
        iso_val, cat_val = buffer[j]
        buffer[j] = buffer[-1]
        buffer.pop()
        yield (iso_val, cat_val)
    # Drain remaining chunks (if any) — occurs when initial file smaller than buffer
    for chunk in chunk_iter:
        extend_buffer_from_chunk(chunk)
        while buffer:
            j = int(rng.integers(len(buffer)))
            iso_val, cat_val = buffer[j]
            buffer[j] = buffer[-1]
            buffer.pop()
            yield (iso_val, cat_val)


def stream_curves_thin_and_fit_from_rows(
    row_iter,
    fsm_value: str,
    allowed_isos: Optional[Set[str]],
    total_reads: int,
    thin_target: int,
):
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

    i_global = 0
    for iso_val, cat_val in row_iter:
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

    return xt, yi_t, yf_t, xf, yi_f, yf_f, n_ident_reads, n_fsm_reads, int(y_ident), int(y_fsm)


def stream_observed_quantiles_from_rows(
    row_iter,
    fsm_value: str,
    allowed_isos: Optional[Set[str]],
    thresholds_ident: Dict[int, int],
    thresholds_fsm: Dict[int, int],
):
    pending_i = {k: True for k in thresholds_ident}
    pending_f = {k: True for k in thresholds_fsm}
    out_i: Dict[int, float] = {k: float("nan") for k in thresholds_ident}
    out_f: Dict[int, float] = {k: float("nan") for k in thresholds_fsm}

    seen_ident: Set[str] = set()
    seen_fsm: Set[str] = set()
    y_ident = 0
    y_fsm = 0

    i_global = 0
    for iso_val, cat_val in row_iter:
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
    p.add_argument(
        "--chunksize",
        type=int,
        default=5_000_000,
        help="Rows per chunk when streaming (memory/perf trade-off).",
    )
    p.add_argument(
        "--shuffle-buffer",
        type=int,
        default=2_000_000,
        help="Approximate shuffle buffer size for streaming mode (O(buffer) memory).",
    )
    p.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip generating the PNG plot (saves memory and dependencies).",
    )
    p.add_argument(
        "--no-curves",
        action="store_true",
        help="Skip writing the thinned curves TSV (saves I/O).",
    )
    p.add_argument(
        "--shuffle-out",
        default=None,
        help="Path to write the shuffled file (default: auto next to input). Use .gz to compress.",
    )
    p.add_argument(
        "--shuffle-run-lines",
        type=int,
        default=200_000,
        help="Approx. number of lines per run (controls memory/temp file size in external shuffle).",
    )
    p.add_argument(
        "--shuffle-fan-in",
        type=int,
        default=64,
        help="Run merge fan-in (controls number of runs merged at once).",
    )
    p.add_argument(
        "--shuffle-tmpdir",
        default=None,
        help="Temporary directory for external shuffle runs (default: system temp).",
    )
    p.add_argument("-q", "--quiet", action="store_true", help="Reduce logging.")
    args = p.parse_args()

    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
    )

    t0 = time.perf_counter()
    in_path = args.input
    # Always do an external full shuffle pre-step
    base = os.path.basename(in_path)
    root = base
    for ext in (".tsv.gz", ".txt.gz", ".csv.gz", ".tsv", ".txt", ".csv"):
        if root.endswith(ext):
            root = root[: -len(ext)]
            break
    shuffle_out = args.shuffle_out or f"{root}.shuffled.tsv.gz"
    logging.info("External shuffle starting -> %s", shuffle_out)
    external_shuffle_to_file(
        in_path,
        shuffle_out,
        seed=args.seed,
        lines_per_run=args.shuffle_run_lines,
        fan_in=args.shuffle_fan_in,
        tmp_dir=args.shuffle_tmpdir,
    )
    in_path = shuffle_out
    logging.info("External shuffle complete: %s", in_path)

    # Always run in streaming mode
    # Pass 1: fractional counts and total reads
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

    # Pass 2: thinned curves and fit (sequential since input is externally shuffled)
    xt, yi_t, yf_t, xf, yi_f, yf_f, n_ident, n_fsm, obs_max_ident, obs_max_fsm = stream_curves_thin_and_fit(
        in_path,
        args.iso_col,
        args.cat_col,
        args.fsm_value,
        allowed_isos,
        n,
        args.thin,
        chunksize=args.chunksize,
    )
    fit_ident = fit_best_model(xf, yi_f)
    fit_fsm = fit_best_model(xf, yf_f)

    # Compute observed quantiles (exact) with pass 3
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
    # Extract x80 from dicts
    x80i_obs = obs_ident.get(80, float("nan"))
    x80f_obs = obs_fsm.get(80, float("nan"))
    pred_ident = {int(f * 100): predicted_x_at_fraction(fit_ident, f) for f in fracs}
    pred_fsm = {int(f * 100): predicted_x_at_fraction(fit_fsm, f) for f in fracs}
    # Prepare arrays for output already computed: xt, yi_t, yf_t
    y_ident = None  # not used further
    y_fsm = None

    # Compose fit summary TSV
    out_tsv, out_png, out_fit = auto_paths(args.input, args.out_tsv, args.png_out, args.fit_out)
    # Observed maxima per curve for summary
    observed_max_ident_val = int(obs_max_ident)
    observed_max_fsm_val = int(obs_max_fsm)

    fit_df = pd.DataFrame(
        [
            {
                "curve": "Identifiable",
                "eligible_read_count": n_ident,
                "observed_max": observed_max_ident_val,
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
                "observed_max": observed_max_fsm_val,
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

    if not args.no_curves:
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
    if not args.no_plot:
        try:
            import matplotlib.pyplot as plt  # type: ignore
        except Exception as e:
            logging.warning("matplotlib not available (%s); skipping plot.", e)
        else:
            plt.figure(figsize=(8.6, 5.4), dpi=120)
            # Empirical (thinned)
            plt.plot(xt, yi_t, label="Identifiable (empirical)", linewidth=1.6)
            plt.plot(xt, yf_t, label="FSM (empirical)", linewidth=1.6)

            # Fitted curves (evaluated on xt for display)
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
            plt.title("Saturation: Identifiable vs FSM (fit + thresholds)")
            plt.legend(frameon=False, ncol=2)
            plt.tight_layout()
            plt.savefig(out_png)
            plt.close()
            logging.info("Wrote plot: %s  (elapsed: %.2fs)", out_png, time.perf_counter() - t0)


if __name__ == "__main__":
    main()
