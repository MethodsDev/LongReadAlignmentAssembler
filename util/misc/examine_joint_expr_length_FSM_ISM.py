#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
examine_joint_expr_length_FSM.v3.py

Improved and fully stable version:
  • Correct TPM bin ordering (low → high)
  • Fixed categorical fillna issue
  • Added logging and timing
  • Silenced pandas deprecation warnings
  • Clean FSM/ISM side-by-side heatmaps

Author: Adapted from Rmd by ChatGPT
"""

import gzip
import time
import logging
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from argparse import ArgumentParser

# ---------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------
logging.basicConfig(
    format="%(asctime)s [%(levelname)s] %(message)s",
    level=logging.INFO,
    datefmt="%H:%M:%S",
)


# ---------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------
def read_tsv_maybe_gz(path: str, usecols=None) -> pd.DataFrame:
    """Read TSV (optionally gzipped)."""
    start = time.time()
    logging.info(f"Reading file: {path}")
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            df = pd.read_csv(fh, sep="\t", usecols=usecols)
    else:
        df = pd.read_csv(path, sep="\t", usecols=usecols)
    logging.info(f"Loaded {len(df):,} rows in {time.time() - start:.1f}s")
    return df


def bin_data(df: pd.DataFrame):
    """Apply fixed bins for TPM and transcript length."""
    tpm_breaks = [0, 0.25, 1, 4, 16, 64, 256, np.inf]
    tpm_labels = ["0+", "0.25+", "1+", "4+", "16+", "64+", "256+"]
    length_breaks = [0, 500, 1000, 2000, 3000, 5000, np.inf]
    length_labels = ["0bp+", "500bp+", "1kbp+", "2kbp+", "3kbp+", "5kbp+"]

    df = df.copy()
    df["tpm_bin"] = pd.cut(
        df["TPM"], bins=tpm_breaks, labels=tpm_labels, include_lowest=True
    )
    df["length_bin"] = pd.cut(
        df["transcript_length"],
        bins=length_breaks,
        labels=length_labels,
        include_lowest=True,
    )

    df_total_ref = (
        df.groupby(["length_bin", "tpm_bin"], dropna=False, observed=True)
        .size()
        .reset_index(name="total_n")
    )
    return df, df_total_ref


def plot_fraction_heatmap_fixed(df_all, df_total_ref, sample_frac=1.0, seed=1234):
    """Simulate sequencing depth by sampling transcripts and plot FSM/ISM dominance."""
    logging.info(f"→ Sampling {sample_frac*100:.0f}% of transcripts ...")
    np.random.seed(seed)
    df_sub = df_all.sample(frac=sample_frac, random_state=seed)

    # Collapse FSM > ISM dominance
    df_sub_collapsed = (
        df_sub.groupby("transcript_id", as_index=False, observed=True)
        .agg(
            {
                "transcript_length": "first",
                "TPM": "first",
                "sqanti_cat": lambda x: (
                    "FSM"
                    if "FSM" in x.values
                    else ("ISM" if "ISM" in x.values else np.nan)
                ),
            }
        )
        .dropna(subset=["sqanti_cat"])
    )

    # Re-bin
    df_binned, _ = bin_data(df_sub_collapsed)

    # Count FSM/ISM
    df_counts = (
        df_binned.groupby(
            ["length_bin", "tpm_bin"], dropna=False, observed=True, group_keys=False
        )
        .apply(
            lambda g: pd.Series(
                {
                    "n_FSM": (g["sqanti_cat"] == "FSM").sum(),
                    "n_ISM": (g["sqanti_cat"] == "ISM").sum(),
                }
            ),
            include_groups=False,
        )
        .reset_index()
    )

    # Merge with totals
    df_frac = pd.merge(
        df_total_ref, df_counts, on=["length_bin", "tpm_bin"], how="left"
    )

    # Only fill numeric columns
    numeric_cols = df_frac.select_dtypes(include=[np.number]).columns
    df_frac[numeric_cols] = df_frac[numeric_cols].fillna(0)

    df_frac["frac_FSM"] = df_frac["n_FSM"] / df_frac["total_n"]
    df_frac["frac_ISM"] = df_frac["n_ISM"] / df_frac["total_n"]

    df_frac_melt = df_frac.melt(
        id_vars=["length_bin", "tpm_bin"],
        value_vars=["frac_FSM", "frac_ISM"],
        var_name="sqanti_cat",
        value_name="fraction",
    )
    df_frac_melt["sqanti_cat"] = df_frac_melt["sqanti_cat"].str.replace("frac_", "")

    # Ensure consistent ordering
    tpm_order = ["0+", "0.25+", "1+", "4+", "16+", "64+", "256+"]
    length_order = ["0bp+", "500bp+", "1kbp+", "2kbp+", "3kbp+", "5kbp+"]

    # Plot (two side-by-side panels)
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharex=True, sharey=True)
    for ax, cat in zip(axes, ["FSM", "ISM"]):
        sub = df_frac_melt[df_frac_melt["sqanti_cat"] == cat]
        pivot = sub.pivot(
            index="tpm_bin", columns="length_bin", values="fraction"
        ).reindex(index=tpm_order, columns=length_order)
        sns.heatmap(
            pivot,
            ax=ax,
            cmap="viridis",
            vmin=0,
            vmax=1,
            cbar=True,
            cbar_kws={"label": "Fraction in bin"},
        )
        ax.invert_yaxis()  # <-- reverse y-axis: low → high
        ax.set_title(cat)
        ax.set_xlabel("cDNA length bin")
        ax.set_ylabel("TPM bin")

    plt.suptitle(
        f"FSM > ISM dominance (sample = {int(sample_frac*100)}%)",
        fontsize=14,
        weight="bold",
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    return fig


# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------
def main():
    parser = ArgumentParser(
        description="Examine joint expression-length FSM dominance patterns."
    )
    parser.add_argument(
        "--iso_cats", required=True, help="Input .iso_cats.tsv or .tsv.gz file"
    )
    parser.add_argument(
        "--transcript_lengths",
        required=True,
        help="Transcript length file (FASTA seqlens)",
    )
    parser.add_argument("--expr", required=True, help="Expression file with TPM values")
    parser.add_argument("--prefix", default="sample", help="Output file prefix")
    parser.add_argument("--outdir", default="plots", help="Output directory")
    args = parser.parse_args()

    start_total = time.time()
    Path(args.outdir).mkdir(exist_ok=True)

    # --- Load isoform categories ---
    cols = [
        "feature_name",
        "sqanti_cat",
        "read_length",
        "alignment_length",
        "num_exon_segments",
        "matching_isoforms",
    ]
    iso_cats_df = read_tsv_maybe_gz(args.iso_cats, usecols=cols)
    iso_cats_df = iso_cats_df[iso_cats_df["sqanti_cat"].isin(["FSM", "ISM"])]
    iso_cats_df["unique_match"] = ~iso_cats_df["matching_isoforms"].str.contains(
        ",", na=False
    )
    iso_cats_df = iso_cats_df[iso_cats_df["unique_match"]].copy()
    iso_cats_df = iso_cats_df.rename(
        columns={"matching_isoforms": "transcript_id"}
    ).drop(columns=["read_length", "alignment_length"])
    logging.info(f"After FSM/ISM + unique filtering: {len(iso_cats_df):,} reads")

    # --- Load transcript lengths ---
    transcript_lengths = pd.read_csv(
        args.transcript_lengths,
        sep="\t",
        header=None,
        names=["transcript_id", "transcript_length"],
    )
    iso_cats_df = pd.merge(
        iso_cats_df, transcript_lengths, on="transcript_id", how="left"
    ).dropna(subset=["transcript_length"])
    logging.info(f"After joining transcript lengths: {len(iso_cats_df):,}")

    # --- Load expression ---
    transcript_expr_info = pd.read_csv(args.expr, sep="\t")
    transcript_expr_info = transcript_expr_info[["gene_id", "transcript_id", "TPM"]]
    iso_cats_df = pd.merge(
        iso_cats_df, transcript_expr_info, on="transcript_id", how="left"
    )
    df = iso_cats_df[iso_cats_df["TPM"] >= 0.1].copy()
    df_all = df[
        ["transcript_id", "transcript_length", "sqanti_cat", "TPM"]
    ].drop_duplicates()
    logging.info(f"Filtered to {len(df_all):,} transcripts with TPM >= 0.1")

    # --- Bin ---
    logging.info("Computing global bin reference ...")
    df_binned_all, df_total_ref = bin_data(df_all)
    logging.info(f"Total bins: {len(df_total_ref):,}")

    # --- Plot each sample fraction ---
    sample_fracs = [0.10, 0.25, 0.50, 0.75, 1.00]
    for frac in sample_fracs:
        fig = plot_fraction_heatmap_fixed(df_all, df_total_ref, sample_frac=frac)
        outpath = (
            Path(args.outdir) / f"{args.prefix}.FSM_ISM_heatmap_{int(frac*100)}pct.png"
        )
        plt.savefig(outpath, dpi=200)
        plt.close(fig)
        logging.info(f"Saved {outpath}")

    logging.info(f"✅ Done in {(time.time() - start_total)/60:.1f} min")


if __name__ == "__main__":
    main()
