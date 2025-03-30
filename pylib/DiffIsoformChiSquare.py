#!/usr/bin/env python3

import os, sys, re
import logging
import argparse
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
import glob

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def differential_isoform_tests(
    df, min_reads_per_gene=25, min_delta_pi=0.1, top_isoforms_each=5
):

    ## Requirements:
    ## df had columns 'gene_id', 'count_A', and 'count_B' with transcripts as separate rows.

    ## method initially written by chatgpt based on methods description in:
    ## https://www.nature.com/articles/s41467-020-20343-5
    ## A spatially resolved brain region- and cell type-specific isoform atlas of the postnatal mouse brain
    ## Jogelkar et al., Nature Communications volume 12, Article number: 463 (2021)
    ## modified by bhaas

    results = []
    grouped = df.groupby("gene_id")

    for gene_id, group in grouped:

        if len(group) < 2:
            continue

        if group["count_A"].sum() == 0 or group["count_B"].sum() == 0:
            continue

        if group[["count_A", "count_B"]].sum().sum() < min_reads_per_gene:
            continue  # Skip genes with insufficient depth

        # print(group)

        top_countA = group.nlargest(top_isoforms_each, "count_A")

        top_countB = group.nlargest(top_isoforms_each, "count_B")

        group = pd.concat([top_countA, top_countB]).drop_duplicates()

        # Sort isoforms by abundance
        group = group.sort_values(by=["count_A", "count_B"], ascending=False)

        # Calculate ΔΠ (percent isoform changes)
        total_counts_A = group["count_A"].sum()
        total_counts_B = group["count_B"].sum()

        if total_counts_A == 0 or total_counts_B == 0:
            continue  # Avoid division errors

        # Construct isoform × category matrix (max 11 × 2)
        matrix = group.iloc[:10][["count_A", "count_B"]].values

        # if len(group) > 10:
        #    other_counts = (
        #        group.iloc[10:][["count_A", "count_B"]].sum().values.reshape(1, -1)
        #    )
        #    matrix = np.vstack([matrix, other_counts])

        # print(matrix)

        pi_A = group["count_A"] / total_counts_A
        pi_B = group["count_B"] / total_counts_B

        delta_pi = pi_B - pi_A

        # Separate positive and negative changes
        positive_changes = delta_pi[delta_pi > 0].nlargest(2).sum()
        negative_changes = delta_pi[delta_pi < 0].nsmallest(2).sum()

        if not (
            abs(positive_changes) > min_delta_pi or abs(negative_changes) > min_delta_pi
        ):
            continue

        # Determine the dominant direction
        if abs(positive_changes) > abs(negative_changes):
            dominant_delta_pi = positive_changes
        else:
            dominant_delta_pi = negative_changes

        # Perform chi-squared test
        chi2, pvalue, _, _ = chi2_contingency(matrix)

        results.append([gene_id, pvalue, dominant_delta_pi])

    # Multiple testing correction
    if results:
        results_df = pd.DataFrame(results, columns=["gene_id", "pvalue", "delta_pi"])
        results_df["adj_pvalue"] = multipletests(results_df["pvalue"], method="fdr_bh")[
            1
        ]

        # Identify significant differential splicing events
        results_df["significant"] = (results_df["adj_pvalue"] <= 0.001) & (
            results_df["delta_pi"].abs() > 0.1
        )
        return results_df

    return None
    # return pd.DataFrame(
    #    columns=["gene_id", "pvalue", "delta_pi", "adjusted_pvalue", "significant"]
    # )


def generate_test_data(num_genes=20):
    np.random.seed(42)
    data = []
    for gene_id in range(1, num_genes):
        num_isoforms = np.random.randint(2, 16)  # Between 2 and 15 isoforms
        for isoform_id in range(1, num_isoforms + 1):
            count_A = np.random.randint(0, 100)
            count_B = np.random.randint(0, 100)
            data.append([f"gene{gene_id}", f"isoform{isoform_id}", count_A, count_B])

    return pd.DataFrame(data, columns=["gene_id", "isoform_id", "count_A", "count_B"])


def run_test():

    logger.info("** Test data:\n")
    df = generate_test_data()
    print(df)

    logger.info("** Test results:\n")
    DE_results = differential_isoform_tests(df)
    print(DE_results)

    return


if __name__ == "__main__":
    run_test()
