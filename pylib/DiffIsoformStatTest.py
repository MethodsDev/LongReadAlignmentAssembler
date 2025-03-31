#!/usr/bin/env python3

import os, sys, re
import logging
import argparse
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact, MonteCarloMethod
from statsmodels.stats.multitest import multipletests
import glob
import statsmodels.api as stats


def differential_isoform_tests(
    df,
    min_reads_per_gene=25,
    min_delta_pi=0.1,
    top_isoforms_each=5,
    test="chi2",
    show_progress_monitor=True,
):

    logger = logging.getLogger(__name__)
    logger.debug("Running differential_isoform_tests()")

    ## Requirements:
    ## df had columns 'gene_id', 'count_A', and 'count_B' with transcripts as separate rows.

    ## test can be 'chi2' or 'fisher'

    ## method initially written by chatgpt based on methods description in:
    ## https://www.nature.com/articles/s41467-020-20343-5
    ## A spatially resolved brain region- and cell type-specific isoform atlas of the postnatal mouse brain
    ## Jogelkar et al., Nature Communications volume 12, Article number: 463 (2021)
    ## modified by bhaas

    assert test in (
        "chi2",
        "fisher",
    ), "Error, not recognizing test type: {}, must be chi2 or fisher".format(test)

    results = []
    grouped = df.groupby("gene_id")

    debug_mode = logger.level == logging.DEBUG

    num_groups = len(grouped)

    # progress monitor.
    group_counter = 0

    for gene_id, group in grouped:

        group_counter += 1
        if show_progress_monitor and group_counter % 1000 == 0:
            frac_done = "{:.2f}% done".format(group_counter / num_groups * 100)
            print(
                f"\r[{group_counter}/{num_groups}] = {frac_done}   ",
                file=sys.stderr,
                end="",
            )

        if debug_mode:
            logger.debug(str(group))

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
            if debug_mode:
                logger.debug("matrix:\n" + str(matrix))
                logger.debug(
                    "skipping due to insuficcient min_delta_pi: {} or {}".format(
                        abs(positive_changes), abs(negative_changes)
                    )
                )
            continue

        # Determine the dominant direction
        if abs(positive_changes) > abs(negative_changes):
            dominant_delta_pi = positive_changes
        else:
            dominant_delta_pi = negative_changes

        if debug_mode:
            logger.debug("testing matrix:\n" + str(matrix))

        if test == "chi2":
            # Perform chi-squared test
            chi2, pvalue, _, _ = chi2_contingency(matrix)
        elif test == "fisher":
            rng = np.random.default_rng()
            method = MonteCarloMethod(rng=rng)
            res = fisher_exact(matrix, method=method)
            pvalue = res.pvalue

        results.append([gene_id, pvalue, dominant_delta_pi])

    # Multiple testing correction
    if results:
        results_df = pd.DataFrame(results, columns=["gene_id", "pvalue", "delta_pi"])
        results_df["test"] = test

        if debug_mode:
            logger.debug("result:\n" + str(results_df))

        return results_df

    logger.debug("-didnt meet requirements to test.")

    return None


def FDR_mult_tests_adjustment(df, signif_threshold=0.001, min_abs_delta_pi=0.1):

    assert "pvalue" in df.columns
    assert "delta_pi" in df.columns

    df["adj_pvalue"] = multipletests(df["pvalue"], method="fdr_bh")[1]

    # Identify significant differential splicing events
    df["significant"] = (df["adj_pvalue"] <= signif_threshold) & (
        df["delta_pi"].abs() >= min_abs_delta_pi
    )

    return df


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

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s : %(levelname)s : %(message)s",
        datefmt="%H:%M:%S",
    )

    logger = logging.getLogger(__name__)

    logger.info("** Test data:\n")
    df = generate_test_data()
    print(df)

    logger.info("** Test results (chi2):\n")
    DE_results = differential_isoform_tests(df, test="chi2")
    DE_results = FDR_mult_tests_adjustment(DE_results)
    print(DE_results)

    logger.info("** Test results (fisher):\n")
    DE_results = differential_isoform_tests(df, test="fisher")
    DE_results = FDR_mult_tests_adjustment(DE_results)
    print(DE_results)

    return


if __name__ == "__main__":
    run_test()
