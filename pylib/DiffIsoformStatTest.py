#!/usr/bin/env python3

import os, sys, logging
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests


def differential_isoform_tests(
    df,
    min_reads_per_gene=25,
    min_delta_pi=0.1,
    reciprocal_delta_pi=False,
    top_isoforms_each=5,
    min_reads_DTU_isoform=25,
    show_progress_monitor=True,
):

    logger = logging.getLogger(__name__)
    logger.debug("Running differential_isoform_tests()")

    results = []
    grouped = df.groupby("gene_id")
    debug_mode = logger.getEffectiveLevel() == logging.DEBUG
    num_groups = len(grouped)

    for group_counter, (gene_id, group) in enumerate(grouped, 1):

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
                f"Either count_A or count_B is zero for gene {gene_id}, skipping."
            )
            continue

        gene_reads_sum = group[["count_A", "count_B"]].sum().sum()
        if gene_reads_sum < min_reads_per_gene:
            logger.debug(
                f"Gene {gene_id} has insufficient reads ({gene_reads_sum}), skipping."
            )
            continue

        top_countA = group.nlargest(top_isoforms_each, "count_A")
        top_countB = group.nlargest(top_isoforms_each, "count_B")
        group = pd.concat([top_countA, top_countB]).drop_duplicates()

        group["total"] = group["count_A"] + group["count_B"]
        group = group.sort_values(by="total", ascending=False).head(10)

        total_counts_A = group["count_A"].sum()
        total_counts_B = group["count_B"].sum()
        if total_counts_A == 0 or total_counts_B == 0:
            logger.debug(
                f"Total counts in condition A or B is zero for gene {gene_id}, skipping."
            )
            continue

        matrix = group[["count_A", "count_B"]].values

        pi_A = group["count_A"] / total_counts_A
        pi_B = group["count_B"] / total_counts_B
        delta_pi = pi_B - pi_A

        positive_indices = delta_pi[delta_pi > 0].sort_values(ascending=False).index[:2]
        negative_indices = delta_pi[delta_pi < 0].sort_values().index[:2]

        positive_sum = delta_pi.loc[positive_indices].sum()
        negative_sum = delta_pi.loc[negative_indices].sum()

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
            logger.debug(f"Gene {gene_id} failed delta_pi threshold.")
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

        dominant_transcript_ids_str = ",".join(
            group.loc[dominant_indices, "transcript_id"].tolist()
        )
        dominant_counts_A = group.loc[dominant_indices, "count_A"].sum()
        dominant_counts_B = group.loc[dominant_indices, "count_B"].sum()
        dominant_total_reads = dominant_counts_A + dominant_counts_B

        alternate_transcript_ids_str = ",".join(
            group.loc[alternate_indices, "transcript_id"].tolist()
        )
        alternate_counts_A = group.loc[alternate_indices, "count_A"].sum()
        alternate_counts_B = group.loc[alternate_indices, "count_B"].sum()
        alternate_total_reads = alternate_counts_A + alternate_counts_B

        # Check minimum read count requirements for DTU isoforms
        if dominant_total_reads < min_reads_DTU_isoform:
            logger.debug(
                f"Gene {gene_id} dominant isoforms have insufficient reads ({dominant_total_reads}), skipping."
            )
            continue

        if reciprocal_delta_pi and alternate_total_reads < min_reads_DTU_isoform:
            logger.debug(
                f"Gene {gene_id} alternate isoforms have insufficient reads ({alternate_total_reads}), skipping."
            )
            continue

        pvalue = None
        status = "OK"
        try:
            chi2, pvalue, _, _ = chi2_contingency(matrix)
        except Exception as e:
            logger.debug(f"Chi2 failed for gene {gene_id}: {e}")
            status = "failed"

        results.append(
            [
                gene_id,
                pvalue,
                dominant_delta_pi,
                dominant_transcript_ids_str,
                total_counts_A,
                total_counts_B,
                dominant_counts_A,
                dominant_counts_B,
                alternate_delta_pi if reciprocal_delta_pi else None,
                alternate_transcript_ids_str if reciprocal_delta_pi else None,
                alternate_counts_A if reciprocal_delta_pi else None,
                alternate_counts_B if reciprocal_delta_pi else None,
                status,
            ]
        )

    if results:
        columns = [
            "gene_id",
            "pvalue",
            "delta_pi",
            "dominant_transcript_ids",
            "total_counts_A",
            "total_counts_B",
            "dominant_counts_A",
            "dominant_counts_B",
        ]
        if reciprocal_delta_pi:
            columns += [
                "alternate_delta_pi",
                "alternate_transcript_ids",
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
            data.append(
                [
                    f"gene{gene_id}",
                    f"isoform{isoform_id}",
                    f"transcript{gene_id}_{isoform_id}",
                    count_A,
                    count_B,
                ]
            )
    return pd.DataFrame(
        data, columns=["gene_id", "isoform_id", "transcript_id", "count_A", "count_B"]
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
    DE_results = differential_isoform_tests(df, reciprocal_delta_pi=True)
    DE_results = FDR_mult_tests_adjustment(DE_results)
    print(DE_results)


if __name__ == "__main__":
    run_test()
