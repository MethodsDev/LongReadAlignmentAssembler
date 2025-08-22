#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
import pandas as pd


# Optional progress bar support
try:
    from tqdm.auto import tqdm
except Exception:  # fallback no-op if tqdm not installed

    def tqdm(x=None, total=None, desc=None, disable=False):
        return x if x is not None else range(0)


import numpy as np

REQUIRED_COLS = [
    "transcript_id",
    "p_val",
    "avg_log2FC",
    "pct.1",
    "pct.2",
    "p_val_adj",
    "cluster1",
    "cluster2",
    "gene_id",
    "comparison",
    "isoform_fraction_cluster1",
    "isoform_fraction_cluster2",
    "isoform_fraction_diff",
]


def parse_args():
    ap = argparse.ArgumentParser(
        description="Screen for differential transcript usage (isoform switching) using DE + isoform fractions, and summarize dominant isoforms per cluster."
    )
    ap.add_argument(
        "-i", "--input", required=True, help="Input TSV with the required columns."
    )
    ap.add_argument(
        "-o",
        "--out_prefix",
        required=True,
        help="Prefix for output files (writes *_rows.tsv, *_groups.tsv, *_dominant.tsv).",
    )
    ap.add_argument(
        "--min_log2fc",
        type=float,
        default=0.25,
        help="Minimum |avg_log2FC| to consider an isoform 'changed' (default: 0.25).",
    )
    ap.add_argument(
        "--min_frac_diff",
        type=float,
        default=0.05,
        help="Minimum |isoform_fraction_diff| to consider meaningful (default: 0.05).",
    )
    ap.add_argument(
        "--allow-sign-inconsistency",
        action="store_true",
        help="By default, require sign(avg_log2FC) matches sign(isoform_fraction_diff). "
        "Set this flag to DISABLE that requirement (less strict).",
    )
    ap.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress bars and periodic status messages.",
    )
    ap.add_argument(
        "--progress-interval",
        type=int,
        default=1000,
        help="If tqdm is unavailable, print a status line every N groups (default: 1000).",
    )
    return ap.parse_args()


def check_columns(df):
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: Input is missing required columns: {missing}")


def annotate_rows(df, min_log2fc, min_frac_diff, require_sign_consistency=True):
    df = df.copy()
    df["direction"] = np.where(
        df["avg_log2FC"] > min_log2fc,
        "up",
        np.where(df["avg_log2FC"] < -min_log2fc, "down", "neutral"),
    )
    df["passes_fc"] = df["avg_log2FC"].abs() >= min_log2fc
    df["passes_frac"] = df["isoform_fraction_diff"].abs() >= min_frac_diff

    # Sign consistency: up => frac_diff > 0; down => frac_diff < 0
    sign_consistent = np.full(len(df), True, dtype=bool)
    up_mask = df["direction"] == "up"
    down_mask = df["direction"] == "down"
    sign_consistent[up_mask] = df.loc[up_mask, "isoform_fraction_diff"] > 0
    sign_consistent[down_mask] = df.loc[down_mask, "isoform_fraction_diff"] < 0
    df["sign_consistent"] = sign_consistent

    if require_sign_consistency:
        df["qualifies"] = (
            df["passes_fc"]
            & df["passes_frac"]
            & df["sign_consistent"]
            & (df["direction"] != "neutral")
        )
    else:
        df["qualifies"] = (
            df["passes_fc"] & df["passes_frac"] & (df["direction"] != "neutral")
        )

    df["abs_log2fc"] = df["avg_log2FC"].abs()
    df["abs_frac_diff"] = df["isoform_fraction_diff"].abs()
    return df


def summarize_groups(df, show_progress=True, fallback_interval=1000):
    group_cols = ["gene_id", "comparison"]

    # Estimate number of groups up-front for better progress bars
    n_groups = df[group_cols].drop_duplicates().shape[0]
    results = []

    iterator = df.groupby(group_cols, sort=False)
    if hasattr(tqdm, "__call__"):
        iterator = tqdm(
            iterator,
            total=n_groups,
            desc="Summarizing groups",
            disable=not show_progress,
        )

    for (gene_id, comparison), g in iterator:
        n = len(g)
        n_up = (g["direction"] == "up").sum()
        n_down = (g["direction"] == "down").sum()
        n_up_qual = ((g["direction"] == "up") & g["qualifies"]).sum()
        n_down_qual = ((g["direction"] == "down") & g["qualifies"]).sum()

        switching = (n_up_qual >= 1) and (n_down_qual >= 1)

        g_sorted = g.sort_values(
            ["abs_frac_diff", "abs_log2fc"], ascending=[False, False]
        )
        top_up = g_sorted[g_sorted["direction"] == "up"].head(3)
        top_down = g_sorted[g_sorted["direction"] == "down"].head(3)

        results.append(
            {
                "gene_id": gene_id,
                "comparison": comparison,
                "n_transcripts": int(n),
                "n_up": int(n_up),
                "n_down": int(n_down),
                "n_up_qualifying": int(n_up_qual),
                "n_down_qualifying": int(n_down_qual),
                "switching_candidate": bool(switching),
                "top_up_transcripts": ",".join(top_up["transcript_id"].astype(str)),
                "top_down_transcripts": ",".join(top_down["transcript_id"].astype(str)),
                "max_abs_frac_diff": float(g["abs_frac_diff"].max()),
                "max_abs_log2fc": float(g["abs_log2fc"].max()),
            }
        )

        # Fallback periodic logging if tqdm isn't available
        if (
            not hasattr(tqdm, "__call__")
            and show_progress
            and (len(results) % fallback_interval == 0)
        ):
            done = len(results)
            print(
                f"[progress] summarize_groups: processed {done}/{n_groups} groups...",
                file=sys.stderr,
            )

    summary = pd.DataFrame(results)
    summary = summary.sort_values(
        ["switching_candidate", "max_abs_frac_diff", "max_abs_log2fc"],
        ascending=[False, False, False],
    )
    return summary


def dominant_summary(df, show_progress=True, fallback_interval=1000):
    """
    For each (gene_id, comparison) group, identify the dominant isoform in each cluster
    using isoform_fraction_cluster1/2. Also report coverage diagnostics
    (sum of fractions in the input per group).
    """
    group_cols = ["gene_id", "comparison"]

    def pick_dominant(g, frac_col):
        # Highest fraction -> dominant. In case of ties, break by higher abs_frac_diff then abs_log2fc.
        order_cols = [frac_col, "abs_frac_diff", "abs_log2fc"]
        g_sorted = g.sort_values(order_cols, ascending=[False, False, False])
        return g_sorted.iloc[0]

    # Estimate number of groups up-front
    n_groups = df[group_cols].drop_duplicates().shape[0]
    records = []

    iterator = df.groupby(group_cols, sort=False)
    if hasattr(tqdm, "__call__"):
        iterator = tqdm(
            iterator,
            total=n_groups,
            desc="Dominant isoform summary",
            disable=not show_progress,
        )

    for (gene_id, comparison), g in iterator:
        sum1 = float(g["isoform_fraction_cluster1"].sum())
        sum2 = float(g["isoform_fraction_cluster2"].sum())

        dom1 = pick_dominant(g, "isoform_fraction_cluster1")
        dom2 = pick_dominant(g, "isoform_fraction_cluster2")

        records.append(
            {
                "gene_id": gene_id,
                "comparison": comparison,
                "dom_isoform_cluster1": dom1["transcript_id"],
                "dom1_fraction": float(dom1["isoform_fraction_cluster1"]),
                "dom1_log2fc": float(dom1["avg_log2FC"]),
                "dom1_direction": dom1["direction"],
                "dom_isoform_cluster2": dom2["transcript_id"],
                "dom2_fraction": float(dom2["isoform_fraction_cluster2"]),
                "dom2_log2fc": float(dom2["avg_log2FC"]),
                "dom2_direction": dom2["direction"],
                "switched": dom1["transcript_id"] != dom2["transcript_id"],
                "sum_frac_cluster1": sum1,
                "sum_frac_cluster2": sum2,
            }
        )

        if (
            not hasattr(tqdm, "__call__")
            and show_progress
            and (len(records) % fallback_interval == 0)
        ):
            done = len(records)
            print(
                f"[progress] dominant_summary: processed {done}/{n_groups} groups...",
                file=sys.stderr,
            )

    dom_df = pd.DataFrame.from_records(records)
    dom_df["max_dom_fraction"] = dom_df[["dom1_fraction", "dom2_fraction"]].max(axis=1)
    dom_df = dom_df.sort_values(
        ["switched", "max_dom_fraction"], ascending=[False, False]
    ).drop(columns=["max_dom_fraction"])
    return dom_df


def main():
    args = parse_args()
    in_path = Path(args.input)
    if not in_path.exists():
        sys.exit(f"ERROR: Input file not found: {in_path}")

    df = pd.read_csv(in_path, sep="	")
    if not args.no_progress:
        print(
            f"Loaded input with {df.shape[0]:,} rows. Annotating rows...",
            file=sys.stderr,
        )
    check_columns(df)

    require_sign_consistency = not args.allow_sign_inconsistency

    df_ann = annotate_rows(
        df,
        min_log2fc=args.min_log2fc,
        min_frac_diff=args.min_frac_diff,
        require_sign_consistency=require_sign_consistency,
    )

    group_summary = summarize_groups(
        df_ann,
        show_progress=(not args.no_progress),
        fallback_interval=args.progress_interval,
    )
    dom_summary = dominant_summary(
        df_ann,
        show_progress=(not args.no_progress),
        fallback_interval=args.progress_interval,
    )

    rows_out = f"{args.out_prefix}_rows.tsv"
    groups_out = f"{args.out_prefix}_groups.tsv"
    dom_out = f"{args.out_prefix}_dominant.tsv"

    df_ann.to_csv(rows_out, sep="\t", index=False)
    group_summary.to_csv(groups_out, sep="\t", index=False)
    dom_summary.to_csv(dom_out, sep="\t", index=False)

    n_groups = group_summary.shape[0]
    n_switching = int(group_summary["switching_candidate"].sum())
    n_switched_dom = int(dom_summary["switched"].sum())

    print(f"Wrote per-row annotations: {rows_out}")
    print(f"Wrote per-group summary:   {groups_out}")
    print(f"Wrote dominant summary:    {dom_out}")
    print(
        f"Groups analyzed: {n_groups}; switching candidates: {n_switching}; dominant switched: {n_switched_dom}"
    )


if __name__ == "__main__":
    main()
