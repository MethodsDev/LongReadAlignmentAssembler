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
        help="Prefix for output file (writes only *_dominant.tsv).",
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
        default=0.1,
        help="Minimum |isoform_fraction_diff| to consider meaningful (default: 0.1).",
    )
    ap.add_argument(
        "--min_cells_fraction",
        type=float,
        default=0.05,
        help="Minimum fraction of cells (pct.*) required in the relevant cluster for a dominant isoform to count toward both_dom_qualify (default: 0.05).",
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


def annotate_rows(df, min_log2fc, min_frac_diff):
    df = df.copy()
    # Expression direction (change in expression level)
    df["expr_direction"] = np.where(
        df["avg_log2FC"] > min_log2fc,
        "up",
        np.where(df["avg_log2FC"] < -min_log2fc, "down", "neutral"),
    )
    df["passes_fc"] = df["avg_log2FC"].abs() >= min_log2fc
    df["passes_frac"] = df["isoform_fraction_diff"].abs() >= min_frac_diff
    # Qualification now ignores sign consistency between fold-change direction and fraction difference.
    df["qualifies"] = (
        df["passes_fc"] & df["passes_frac"] & (df["expr_direction"] != "neutral")
    )

    df["abs_log2fc"] = df["avg_log2FC"].abs()
    df["abs_frac_diff"] = df["isoform_fraction_diff"].abs()
    # Isoform fraction direction: increase or decrease in isoform fraction
    df["isoform_fraction_direction"] = np.where(
        df["isoform_fraction_diff"] > 0,
        "up",
        np.where(df["isoform_fraction_diff"] < 0, "down", "neutral"),
    )
    return df


def dominant_summary(df, min_cells_fraction, show_progress=True, fallback_interval=1000):
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

        # Determine cell fraction pass criteria for dominant isoforms in their respective clusters
        dom1_cells_frac_cluster1 = float(dom1["pct.1"])  # cluster1 fraction for domIso1
        dom2_cells_frac_cluster2 = float(dom2["pct.2"])  # cluster2 fraction for domIso2
        dom1_cells_pass = dom1_cells_frac_cluster1 >= min_cells_fraction
        dom2_cells_pass = dom2_cells_frac_cluster2 >= min_cells_fraction

        records.append(
            {
                "gene_id": gene_id,
                "comparison": comparison,
                "domIso1_transcript_id": dom1["transcript_id"],
                "domIso1_fraction_cluster1": float(dom1["isoform_fraction_cluster1"]),
                "domIso1_fraction_cluster2": float(dom1["isoform_fraction_cluster2"]),
                "domIso1_isoform_fraction_diff": float(dom1["isoform_fraction_diff"]),
                "domIso1_isoform_fraction_direction": dom1["isoform_fraction_direction"],
                "domIso1_log2fc": float(dom1["avg_log2FC"]),
                "domIso1_expr_direction": dom1["expr_direction"],
                "domIso1_cells_fraction1": dom1_cells_frac_cluster1,
                "domIso1_cells_fraction2": float(dom1["pct.2"]),
                "domIso1_cells_fraction_pass": dom1_cells_pass,
                "domIso1_p_val": float(dom1["p_val"]),
                "domIso1_p_val_adj": float(dom1["p_val_adj"]),
                "domIso1_qualifies": bool(dom1["qualifies"]),
                "domIso2_transcript_id": dom2["transcript_id"],
                "domIso2_fraction_cluster1": float(dom2["isoform_fraction_cluster1"]),
                "domIso2_fraction_cluster2": float(dom2["isoform_fraction_cluster2"]),
                "domIso2_isoform_fraction_diff": float(dom2["isoform_fraction_diff"]),
                "domIso2_isoform_fraction_direction": dom2["isoform_fraction_direction"],
                "domIso2_log2fc": float(dom2["avg_log2FC"]),
                "domIso2_expr_direction": dom2["expr_direction"],
                "domIso2_cells_fraction1": float(dom2["pct.1"]),
                "domIso2_cells_fraction2": dom2_cells_frac_cluster2,
                "domIso2_cells_fraction_pass": dom2_cells_pass,
                "domIso2_p_val": float(dom2["p_val"]),
                "domIso2_p_val_adj": float(dom2["p_val_adj"]),
                "domIso2_qualifies": bool(dom2["qualifies"]),
                "switched": dom1["transcript_id"] != dom2["transcript_id"],
                "sum_frac_cluster1": sum1,
                "sum_frac_cluster2": sum2,
                # PASS if both dominant isoforms qualify and meet cell fraction thresholds
                "PASS": bool(
                    dom1["qualifies"]
                    and dom2["qualifies"]
                    and dom1_cells_pass
                    and dom2_cells_pass
                    and (
                        (dom1["isoform_fraction_direction"] == "up" and dom2["isoform_fraction_direction"] == "down")
                        or (dom1["isoform_fraction_direction"] == "down" and dom2["isoform_fraction_direction"] == "up")
                    )
                ),
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

    df_ann = annotate_rows(
        df,
        min_log2fc=args.min_log2fc,
        min_frac_diff=args.min_frac_diff,
    )

    dom_summary = dominant_summary(
        df_ann,
        min_cells_fraction=args.min_cells_fraction,
        show_progress=(not args.no_progress),
        fallback_interval=args.progress_interval,
    )

    # Round numeric columns to three decimal places for output consistency
    def round_numeric(df, ndigits=3):
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        df[numeric_cols] = df[numeric_cols].round(ndigits)
        return df

    dom_summary = round_numeric(dom_summary)
    dom_out = f"{args.out_prefix}_dominant.tsv"
    both_out = f"{args.out_prefix}_dominant.PASS.tsv"
    # Only writing dominant summary as final output
    dom_summary.to_csv(dom_out, sep="\t", index=False)
    dom_summary[dom_summary["PASS"]].to_csv(both_out, sep="\t", index=False)

    n_switched_dom = int(dom_summary["switched"].sum())
    n_both_qual = int((dom_summary["switched"] & dom_summary["PASS"]).sum())
    print(f"Wrote dominant summary: {dom_out}")
    print(f"Wrote PASS subset: {both_out}")
    print(
        f"Dominant isoform switches: {n_switched_dom}; switched comparisons with both dominant isoforms qualifying: {n_both_qual}"
    )

    # Done


if __name__ == "__main__":
    main()
