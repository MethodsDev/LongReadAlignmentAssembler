#!/usr/bin/env python3

import os, sys, logging
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests


def round_to_significant_figures(x, sig_figs=3):
    """Round a number to specified significant figures"""
    if x == 0:
        return 0
    return round(x, -int(np.floor(np.log10(abs(x)))) + (sig_figs - 1))


def differential_isoform_tests(
    df,
    group_by_token="gene_id",
    min_reads_per_gene=25,
    min_delta_pi=0.1,
    reciprocal_delta_pi=False,
    top_isoforms_each=5,
    min_reads_DTU_isoform=25,
    show_progress_monitor=True,
    delta_pi_precision=3,
):

    logger = logging.getLogger(__name__)
    logger.debug("Running differential_isoform_tests()")

    results = []
    grouped = df.groupby(group_by_token)
    debug_mode = logger.getEffectiveLevel() == logging.DEBUG
    num_groups = len(grouped)

    for group_counter, (group_by_id, group) in enumerate(grouped, 1):

        if show_progress_monitor and group_counter % 1000 == 0:
            frac_done = f"{group_counter / num_groups * 100:.2f}% done"
            print(
                f"\r[{group_counter}/{num_groups}] = {frac_done}   ",
                file=sys.stderr,
                end="",
            )

        if debug_mode:
            logger.debug(str(group))

        if len(group) < 2:
            continue

        if group[["count_A", "count_B"]].sum().min() == 0:
            logger.debug(
                f"Either count_A or count_B is zero for {group_by_token} {group_by_id}, skipping."
            )
            continue

        # Store original total counts for reporting
        original_total_counts_A = group["count_A"].sum()
        original_total_counts_B = group["count_B"].sum()

        if original_total_counts_A == 0 or original_total_counts_B == 0:
            logger.debug(
                f"Total counts in condition A or B is zero for {group_by_token} {group_by_id}, skipping."
            )
            continue

        if (
            original_total_counts_A < min_reads_per_gene
            or original_total_counts_B < min_reads_per_gene
        ):
            logger.debug(
                f"{group_by_token} {group_by_id} has insufficient reads (A: {original_total_counts_A} or B: {original_total_counts_B}), skipping."
            )
            continue

        # Calculate pi_A and pi_B on the ORIGINAL group data (before filtering)
        pi_A = group["count_A"] / original_total_counts_A
        pi_B = group["count_B"] / original_total_counts_B
        delta_pi = pi_B - pi_A

        # Add these calculations to the original group
        group = group.copy()
        group["pi_A"] = pi_A
        group["pi_B"] = pi_B
        group["delta_pi"] = delta_pi

        # Now filter to top isoforms for downstream analysis
        top_countA = group.nlargest(top_isoforms_each, "count_A")
        top_countB = group.nlargest(top_isoforms_each, "count_B")
        filtered_group = pd.concat([top_countA, top_countB]).drop_duplicates()

        filtered_group["total"] = filtered_group["count_A"] + filtered_group["count_B"]
        filtered_group = filtered_group.sort_values(by="total", ascending=False).head(
            10
        )

        # Extract delta_pi values for the filtered isoforms only
        filtered_delta_pi = filtered_group["delta_pi"]

        # Now work with delta_pi from filtered isoforms (but calculated on original proportions)
        positive_indices = (
            filtered_delta_pi[filtered_delta_pi > 0]
            .sort_values(ascending=False)
            .index[:2]
        )
        negative_indices = (
            filtered_delta_pi[filtered_delta_pi < 0].sort_values().index[:2]
        )

        positive_sum = filtered_delta_pi.loc[positive_indices].sum()
        negative_sum = filtered_delta_pi.loc[negative_indices].sum()

        pass_delta_pi = False
        if reciprocal_delta_pi:
            pass_delta_pi = (
                abs(positive_sum) > min_delta_pi and abs(negative_sum) > min_delta_pi
            )
        else:
            pass_delta_pi = (
                abs(positive_sum) > min_delta_pi or abs(negative_sum) > min_delta_pi
            )

        if not pass_delta_pi:
            logger.debug(f"{group_by_token} {group_by_id} failed delta_pi threshold.")
            continue

        if abs(positive_sum) > abs(negative_sum):
            dominant_delta_pi = positive_sum
            dominant_indices = positive_indices
            alternate_delta_pi = negative_sum
            alternate_indices = negative_indices
        else:
            dominant_delta_pi = negative_sum
            dominant_indices = negative_indices
            alternate_delta_pi = positive_sum
            alternate_indices = positive_indices

        # Use FILTERED group data for transcript information and counts
        dominant_transcript_ids_str = ",".join(
            filtered_group.loc[dominant_indices, "transcript_id"].tolist()
        )
        dominant_counts_A = filtered_group.loc[dominant_indices, "count_A"].sum()
        dominant_counts_B = filtered_group.loc[dominant_indices, "count_B"].sum()
        dominant_total_reads = dominant_counts_A + dominant_counts_B

        # Get pi values for dominant transcripts
        dominant_pi_A_values = filtered_group.loc[dominant_indices, "pi_A"].tolist()
        dominant_pi_B_values = filtered_group.loc[dominant_indices, "pi_B"].tolist()
        dominant_pi_A_str = ",".join(
            [f"{pi:.{delta_pi_precision}f}" for pi in dominant_pi_A_values]
        )
        dominant_pi_B_str = ",".join(
            [f"{pi:.{delta_pi_precision}f}" for pi in dominant_pi_B_values]
        )

        # Extract gene_ids and splice_hashcodes for dominant transcripts
        dominant_gene_ids = []
        dominant_splice_hashcodes = []
        if "gene_id" in filtered_group.columns:
            dominant_gene_ids = filtered_group.loc[dominant_indices, "gene_id"].tolist()
        if "splice_hashcode" in filtered_group.columns:
            dominant_splice_hashcodes = filtered_group.loc[
                dominant_indices, "splice_hashcode"
            ].tolist()

        dominant_gene_ids_str = (
            ",".join(map(str, dominant_gene_ids)) if dominant_gene_ids else ""
        )
        dominant_splice_hashcodes_str = (
            ",".join(map(str, dominant_splice_hashcodes))
            if dominant_splice_hashcodes
            else ""
        )

        alternate_transcript_ids_str = ",".join(
            filtered_group.loc[alternate_indices, "transcript_id"].tolist()
        )
        alternate_counts_A = filtered_group.loc[alternate_indices, "count_A"].sum()
        alternate_counts_B = filtered_group.loc[alternate_indices, "count_B"].sum()
        alternate_total_reads = alternate_counts_A + alternate_counts_B

        # Get pi values for alternate transcripts
        alternate_pi_A_values = filtered_group.loc[alternate_indices, "pi_A"].tolist()
        alternate_pi_B_values = filtered_group.loc[alternate_indices, "pi_B"].tolist()
        alternate_pi_A_str = ",".join(
            [f"{pi:.{delta_pi_precision}f}" for pi in alternate_pi_A_values]
        )
        alternate_pi_B_str = ",".join(
            [f"{pi:.{delta_pi_precision}f}" for pi in alternate_pi_B_values]
        )

        # Extract gene_ids and splice_hashcodes for alternate transcripts
        alternate_gene_ids = []
        alternate_splice_hashcodes = []
        if "gene_id" in filtered_group.columns:
            alternate_gene_ids = filtered_group.loc[
                alternate_indices, "gene_id"
            ].tolist()
        if "splice_hashcode" in filtered_group.columns:
            alternate_splice_hashcodes = filtered_group.loc[
                alternate_indices, "splice_hashcode"
            ].tolist()

        alternate_gene_ids_str = (
            ",".join(map(str, alternate_gene_ids)) if alternate_gene_ids else ""
        )
        alternate_splice_hashcodes_str = (
            ",".join(map(str, alternate_splice_hashcodes))
            if alternate_splice_hashcodes
            else ""
        )

        # Check minimum read count requirements for DTU isoforms
        if dominant_total_reads < min_reads_DTU_isoform:
            logger.debug(
                f"{group_by_token} {group_by_id} dominant isoforms have insufficient reads ({dominant_total_reads}), skipping."
            )
            continue

        if reciprocal_delta_pi and alternate_total_reads < min_reads_DTU_isoform:
            logger.debug(
                f"{group_by_token} {group_by_id} alternate isoforms have insufficient reads ({alternate_total_reads}), skipping."
            )
            continue

        # Use filtered group for chi-squared test matrix
        matrix = filtered_group[["count_A", "count_B"]].values

        pvalue = None
        status = "OK"
        try:
            chi2, pvalue, _, _ = chi2_contingency(matrix)
        except Exception as e:
            logger.debug(f"Chi2 failed for {group_by_token} {group_by_id}: {e}")
            status = "failed"

        # Round delta_pi values to specified precision
        dominant_delta_pi_rounded = round_to_significant_figures(
            dominant_delta_pi, delta_pi_precision
        )
        alternate_delta_pi_rounded = (
            round_to_significant_figures(alternate_delta_pi, delta_pi_precision)
            if reciprocal_delta_pi
            else None
        )

        # Build results array based on reciprocal_delta_pi setting
        result_row = [
            group_by_id,
            pvalue,
            dominant_delta_pi_rounded,
            dominant_transcript_ids_str,
            dominant_pi_A_str,
            dominant_pi_B_str,
            dominant_gene_ids_str,
            dominant_splice_hashcodes_str,
            original_total_counts_A,  # Keep original totals for reference
            original_total_counts_B,  # Keep original totals for reference
            dominant_counts_A,
            dominant_counts_B,
        ]

        if reciprocal_delta_pi:
            result_row.extend(
                [
                    alternate_delta_pi_rounded,
                    alternate_transcript_ids_str,
                    alternate_pi_A_str,
                    alternate_pi_B_str,
                    alternate_gene_ids_str,
                    alternate_splice_hashcodes_str,
                    alternate_counts_A,
                    alternate_counts_B,
                ]
            )

        result_row.append(status)
        results.append(result_row)

    if results:
        columns = [
            group_by_token,
            "pvalue",
            "delta_pi",
            "dominant_transcript_ids",
            "dominant_pi_A",
            "dominant_pi_B",
            "dominant_gene_ids",
            "dominant_splice_hashcodes",
            "total_counts_A",
            "total_counts_B",
            "dominant_counts_A",
            "dominant_counts_B",
        ]
        if reciprocal_delta_pi:
            columns += [
                "alternate_delta_pi",
                "alternate_transcript_ids",
                "alternate_pi_A",
                "alternate_pi_B",
                "alternate_gene_ids",
                "alternate_splice_hashcodes",
                "alternate_counts_A",
                "alternate_counts_B",
            ]
        columns += ["status"]
        results_df = pd.DataFrame(results, columns=columns)
        return results_df

    logger.debug("-didn't meet requirements to test.")
    return None


def generate_test_data(num_genes=20):
    np.random.seed(42)
    data = []
    for gene_id in range(1, num_genes):
        num_isoforms = np.random.randint(2, 16)
        for isoform_id in range(1, num_isoforms + 1):
            count_A = np.random.randint(0, 100)
            count_B = np.random.randint(0, 100)
            # Add splice_hashcode for testing
            splice_hashcode = f"hash_{gene_id}_{isoform_id}"
            data.append(
                [
                    f"gene{gene_id}",
                    f"isoform{isoform_id}",
                    f"transcript{gene_id}_{isoform_id}",
                    splice_hashcode,
                    count_A,
                    count_B,
                ]
            )
    return pd.DataFrame(
        data,
        columns=[
            "gene_id",
            "isoform_id",
            "transcript_id",
            "splice_hashcode",
            "count_A",
            "count_B",
        ],
    )


def FDR_mult_tests_adjustment(df, signif_threshold=0.001, min_abs_delta_pi=0.1):
    print(df)
    assert "pvalue" in df.columns
    assert "delta_pi" in df.columns
    df["adj_pvalue"] = multipletests(df["pvalue"], method="fdr_bh")[1]
    df["significant"] = (df["adj_pvalue"] <= signif_threshold) & (
        df["delta_pi"].abs() >= min_abs_delta_pi
    )
    return df


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
    DE_results = differential_isoform_tests(
        df, reciprocal_delta_pi=True, top_isoforms_each=1
    )
    DE_results = FDR_mult_tests_adjustment(DE_results)
    print(DE_results)


if __name__ == "__main__":
    run_test()
