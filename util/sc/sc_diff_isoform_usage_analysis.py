#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests


def differential_isoform_tests(df):

    ## method written by chatgpt based on methods description in:
    ## https://www.nature.com/articles/s41467-020-20343-5
    ## A spatially resolved brain region- and cell type-specific isoform atlas of the postnatal mouse brain
    ## Jogelkar et al., Nature Communications volume 12, Article number: 463 (2021)

    results = []
    grouped = df.groupby("gene_id")

    for gene_id, group in grouped:
        if group[["count_A", "count_B"]].sum().sum() < 25:
            continue  # Skip genes with insufficient depth

        # Sort isoforms by abundance
        group = group.sort_values(by=["count_A", "count_B"], ascending=False)

        # Construct isoform × category matrix (max 11 × 2)
        matrix = group.iloc[:10][["count_A", "count_B"]].values

        if len(group) > 10:
            other_counts = (
                group.iloc[10:][["count_A", "count_B"]].sum().values.reshape(1, -1)
            )
            matrix = np.vstack([matrix, other_counts])

        # Perform chi-squared test
        chi2, pvalue, _, _ = chi2_contingency(matrix)

        # Calculate ΔΠ (percent isoform changes)
        total_counts_A = group["count_A"].sum()
        total_counts_B = group["count_B"].sum()

        if total_counts_A == 0 or total_counts_B == 0:
            continue  # Avoid division errors

        pi_A = group["count_A"] / total_counts_A
        pi_B = group["count_B"] / total_counts_B

        delta_pi = (pi_B - pi_A).abs().nlargest(2).sum()

        results.append([gene_id, pvalue, delta_pi])

    # Multiple testing correction
    if results:
        results_df = pd.DataFrame(results, columns=["gene_id", "pvalue", "delta_pi"])
        results_df["adjusted_pvalue"] = multipletests(
            results_df["pvalue"], method="fdr_bh"
        )[1]

        # Identify significant differential splicing events
        results_df["significant"] = (results_df["adjusted_pvalue"] <= 0.05) & (
            results_df["delta_pi"] > 0.1
        )
        return results_df

    return pd.DataFrame(
        columns=["gene_id", "pvalue", "delta_pi", "adjusted_pvalue", "significant"]
    )
