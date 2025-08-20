#!/usr/bin/env python3

# written by chatgpt-5, mod by bhaas

"""
Compute expression breadth/specificity metrics per transcript across clusters.

Input
-----
A table with rows = transcripts (or genes) and columns = clusters.
Values should be nonnegative pseudobulk expression (raw counts or CPM/TPM).
An ID column may be included (default: first column is ID).

Output
------
TSV with metrics for each transcript. Note: "evenness" was removed as it is
equivalent to normalized entropy.

Usage
-----
python specificity_metrics.py \
    --input pseudobulk.tsv \
    --output metrics.tsv
"""

import argparse
import math
import numpy as np
import pandas as pd
import sys

# ------------------------
# Metric functions
# ------------------------


def shannon_entropy(p: np.ndarray) -> float:
    """Shannon entropy H = -sum p log p (natural log).

    Interpretation:
    - Low H (~0): expression concentrated in one/few clusters (specific).
    - High H (up to log K): expression spread uniformly across clusters (broad).
    """
    p = p[(p > 0) & np.isfinite(p)]
    if p.size == 0:
        return np.nan
    return float(-np.sum(p * np.log(p)))


def normalized_entropy(p: np.ndarray) -> float:
    """H / log(K), ranges in [0,1] for K>1.

    Interpretation:
    - 0: highly specific (mass on a single cluster).
    - 1: perfectly uniform across clusters (maximal breadth).
    """
    K = p.size
    if K <= 1:
        return np.nan
    H = shannon_entropy(p)
    if not np.isfinite(H):
        return np.nan
    return float(H / math.log(K))


def gini(x: np.ndarray) -> float:
    """Gini coefficient over raw (nonnegative) expression x.

    Interpretation:
    - 0: perfectly equal expression across clusters.
    - 1: all expression in one cluster (max inequality/specificity).
    """
    x = x[np.isfinite(x)]
    x = x[x >= 0]
    n = x.size
    if n == 0:
        return np.nan
    s = x.sum()
    if s == 0:
        return np.nan
    xs = np.sort(x)
    idx = np.arange(1, n + 1)
    return float((2.0 * np.sum(idx * xs)) / (n * s) - (n + 1.0) / n)


def tau_metric(x: np.ndarray) -> float:
    """Tau (Yanai et al. 2005). 0=uniform, 1=single-cluster specific.

    Interpretation:
    - 0: uniform (broadly expressed).
    - 1: expressed only in the max-expressing cluster (high specificity).
    """
    x = x[np.isfinite(x)]
    K = x.size
    if K <= 1:
        return np.nan
    m = x.max(initial=0.0)
    if m <= 0:
        return np.nan
    return float(np.sum(1.0 - (x / m)) / (K - 1))


def kl_divergence_to_uniform(p: np.ndarray) -> float:
    """KL divergence to the uniform distribution.

    Interpretation:
    - 0: identical to uniform (broad expression).
    - Higher: deviates from uniform (more specific/skewed).
    """
    K = p.size
    if K == 0:
        return np.nan
    H = shannon_entropy(p)
    if not np.isfinite(H):
        return np.nan
    return float(-H + math.log(K))


def simpson_index(p: np.ndarray) -> float:
    """Simpson’s index: sum p_i^2.

    Interpretation:
    - Low (near 1/K): uniform distribution (broad expression).
    - High (up to 1): dominated by a single cluster (specific).
    """
    if p.size == 0:
        return np.nan
    return float(np.sum(np.square(p)))


def as_proportions(x: np.ndarray, pseudocount: float = 0.0):
    """Convert to proportions after adding pseudocount."""
    if pseudocount < 0:
        raise ValueError("pseudocount must be >= 0")
    x = np.asarray(x, dtype=float)
    if np.any(x < 0):
        raise ValueError("expression values must be nonnegative")
    x = x + pseudocount
    s = x.sum()
    if s == 0:
        return np.full_like(x, np.nan), s
    return x / s, s


def compute_metrics_for_row(
    x: np.ndarray, cluster_names: np.ndarray, pseudocount: float
):
    p, row_sum = as_proportions(x, pseudocount=pseudocount)

    if not np.isfinite(p).all():
        H = H_norm = KL = D = np.nan
    else:
        H = shannon_entropy(p)
        H_norm = normalized_entropy(p)
        KL = kl_divergence_to_uniform(p)
        D = simpson_index(p)

    G = gini(x)
    tau = tau_metric(x)

    if np.all(~np.isfinite(p)) or np.all(np.isnan(p)):
        max_frac = np.nan
        max_cluster = np.nan
    else:
        max_idx = np.nanargmax(p)
        max_frac = float(p[max_idx])
        max_cluster = cluster_names[max_idx]

    return {
        "row_sum": row_sum,
        "entropy": H,
        "entropy_normalized": H_norm,
        "kl_to_uniform": KL,
        "simpson": D,
        "gini": G,
        "tau": tau,
        "max_fraction": max_frac,
        "max_cluster": max_cluster,
    }


# ------------------------
# CPM normalization
# ------------------------


def cpm_normalize(df_num: pd.DataFrame, per_million: float = 1e6) -> pd.DataFrame:
    """Column-wise CPM normalization."""
    col_sums = df_num.sum(axis=0)
    scale = pd.Series(0.0, index=df_num.columns)
    pos = col_sums > 0
    scale.loc[pos] = per_million / col_sums.loc[pos]
    return df_num.multiply(scale, axis=1)


# ------------------------
# Main
# ------------------------


def main():
    ap = argparse.ArgumentParser(
        description="Compute breadth/specificity metrics per transcript."
    )
    ap.add_argument(
        "--input",
        required=True,
        help="Input table (TSV/CSV). Rows=transcripts, Cols=clusters.",
    )
    ap.add_argument("--output", required=True, help="Output TSV with metrics.")
    ap.add_argument(
        "--pseudocount",
        type=float,
        default=0.0,
        help="Pseudocount added before converting to proportions (default: 0.0).",
    )
    args = ap.parse_args()

    # Detect delimiter
    sep = "\t"

    df = pd.read_csv(args.input, sep=sep)
    # Preserve the original first column (e.g., transcript_id)
    id_col_name = df.columns[0]
    id_values = df.iloc[:, 0].copy()

    # Keep only numeric cluster columns
    cluster_cols = list(df.columns[1:])

    # Always CPM normalize cluster columns
    df[cluster_cols] = cpm_normalize(df[cluster_cols], per_million=1e6)

    X = df[cluster_cols].to_numpy(dtype=float)
    cluster_names = np.array(cluster_cols, dtype=object)

    results = []
    for ridx in range(X.shape[0]):
        x = X[ridx, :]
        metrics = compute_metrics_for_row(
            x, cluster_names, pseudocount=args.pseudocount
        )
        results.append(metrics)

    out = pd.DataFrame(results, index=df.index)
    cols = [
        "row_sum",
        "entropy",
        "entropy_normalized",
        "kl_to_uniform",
        "simpson",
        "gini",
        "tau",
        "max_fraction",
        "max_cluster",
    ]
    out = out[cols]
    # Insert original ID column as the first column and write without index
    out.insert(0, id_col_name, id_values.values)
    out.to_csv(args.output, sep="\t", index=False)
    print(
        f"Wrote metrics for {out.shape[0]} transcripts across {len(cluster_cols)} clusters to {args.output}"
    )


if __name__ == "__main__":
    main()
