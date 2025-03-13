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


def main():
    parser = argparse.ArgumentParser(
        description="run single cell cluster isoform DE usage tests via pseudobulk chi-square",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        # formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "--sc_cluster_expr_dir",
        type=str,
        required=True,
        help="directory containing LRAA single cell cluster pseudocount expr matrices",
    )

    parser.add_argument(
        "--top_isoforms_each",
        type=int,
        default=5,
        help="max number of top-expressed isoforms from each cluster to compare",
    )

    parser.add_argument(
        "--min_reads_per_gene",
        type=int,
        default=25,
        help="min reads per gene (default: 25)",
    )

    parser.add_argument(
        "--min_delta_pi", type=float, default=0.1, help="min delta pi (default: 0.1)"
    )

    parser.add_argument(
        "--collapse_splice_identical",
        action="store_true",
        default=False,
        help="collapse isoforms that have identical splice patterns",
    )

    parser.add_argument(
        "--ignore_unspliced",
        action="store_true",
        default=False,
        help="exclude unspliced isoforms",
    )

    parser.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="prefix for output files",
    )

    parser.add_argument(
        "--run_demo_test",
        action="store_true",
        default=False,
        help="run functional test using random data",
    )

    parser.add_argument(
        "--summary_matrices_only",
        action="store_true",
        default=False,
        help="just generate count and TPM matrices and then stop - no iso diff usage stats performed",
    )

    args = parser.parse_args()

    if args.run_demo_test:
        run_test()
        sys.exit(0)

    top_isoforms_each = args.top_isoforms_each
    min_reads_per_gene = args.min_reads_per_gene
    min_delta_pi = args.min_delta_pi
    sc_cluster_expr_dir = args.sc_cluster_expr_dir
    collapse_splice_identical = args.collapse_splice_identical
    ignore_unspliced = args.ignore_unspliced
    output_prefix = args.output_prefix

    cluster_names, counts_big_df, TPM_big_df = build_big_sc_cluster_count_df(
        sc_cluster_expr_dir, collapse_splice_identical, ignore_unspliced
    )

    counts_big_df_fname_prefix = f"{output_prefix}.full.counts.tsv"
    TPM_big_df_fname_prefix = f"{output_prefix}.full.TPM.tsv"

    if ignore_unspliced:
        counts_big_df_fname_prefix += ".ignore_unspliced"
        TPM_big_df_fname_prefix += ".ignore_unspliced"

    if collapse_splice_identical:
        counts_big_df_fname_prefix += ".collapse_splice_identical"
        TPM_big_df_fname_prefix += ".collapse_splice_identical"

    counts_big_df.to_csv(f"{counts_big_df_fname_prefix}.tsv", sep="\t", index=False)
    TPM_big_df.to_csv(f"{TPM_big_df_fname_prefix}.tsv", sep="\t", index=False)

    if args.summary_matrices_only:
        logger.info("-only summary matrices being generated. Stopping now.")
        sys.exit(0)

    ############################################################
    ## pairwise compare clusters for diff isoform usage analysis
    ############################################################

    # all_test_results = None
    for i in range(len(cluster_names)):
        cluster_i = cluster_names[i]
        for j in range(i + 1, len(cluster_names)):
            cluster_j = cluster_names[j]

            logger.info("Testing pair: {} vs {}".format(cluster_i, cluster_j))
            test_df = counts_big_df[
                ["gene_id", "transcript_id", cluster_i, cluster_j]
            ].copy()
            test_df.rename(
                columns={cluster_i: "count_A", cluster_j: "count_B"}, inplace=True
            )

            test_df = test_df.loc[((test_df["count_A"] > 0) | (test_df["count_B"] > 0))]

            print(test_df)

            test_df_results = differential_isoform_tests(
                test_df, min_reads_per_gene, min_delta_pi, top_isoforms_each
            )

            if test_df_results is not None:
                test_df_results["cluster_A"] = cluster_i
                test_df_results["cluster_B"] = cluster_j

                # all_test_results = pd.concat([all_test_results, test_df_results])
                test_df_results.to_csv(
                    f"{output_prefix}.{cluster_i}-vs-{cluster_j}.diff_iso.tsv",
                    sep="\t",
                    index=False,
                )

    # print(all_test_results)

    sys.exit(0)


def differential_isoform_tests(
    df, min_reads_per_gene=25, min_delta_pi=0.1, top_isoforms_each=5
):

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


def build_big_sc_cluster_count_df(
    sc_cluster_expr_dir, collapse_splice_identical=False, ignore_unspliced=False
):
    files = glob.glob(f"{sc_cluster_expr_dir}/*quant.expr*")

    count_data_frames = []
    TPM_data_frames = []

    sample_names = list()

    for file in files:
        sample_name = os.path.basename(file)
        sample_name = re.sub("\\.quant.expr.*$", "", sample_name)
        logger.info("-incorporating {} into count and TPM matrices".format(sample_name))
        df = pd.read_csv(file, sep="\t")  # Adjust separator if necessary

        if ignore_unspliced:
            df = df[pd.notna(df["introns"])]

        if collapse_splice_identical:
            df = collapse_splice_identical_isoforms(df)

        df = df[["gene_id", "transcript_id", "all_reads", "TPM"]]

        df_counts = df[["gene_id", "transcript_id", "all_reads"]].rename(
            columns={"all_reads": sample_name}
        )
        df_TPM = df[["gene_id", "transcript_id", "TPM"]].rename(
            columns={"TPM": sample_name}
        )

        count_data_frames.append(df_counts)
        TPM_data_frames.append(df_TPM)
        sample_names.append(sample_name)

    # Merge all dataframes on gene_id and transcript_id
    merged_counts_df = count_data_frames[0]
    for df in count_data_frames[1:]:
        merged_counts_df = pd.merge(
            merged_counts_df, df, on=["gene_id", "transcript_id"], how="outer"
        )

    merged_TPM_df = TPM_data_frames[0]
    for df in TPM_data_frames[1:]:
        merged_TPM_df = pd.merge(
            merged_TPM_df, df, on=["gene_id", "transcript_id"], how="outer"
        )

    return sample_names, merged_counts_df, merged_TPM_df


def collapse_splice_identical_isoforms(df):

    def get_group_key(row):
        return row["introns"] if pd.notna(row["introns"]) else row["exons"]

    df["group_key"] = df.apply(get_group_key, axis=1)

    def process_group(group_df):
        group_key = group_df.name
        total_reads = group_df["all_reads"].sum()
        total_TP = group_df["TPM"].sum()
        first_row = group_df.iloc[
            0
        ].copy()  # Get a copy of the first row to preserve other columns.
        first_row["transcript_id"] = group_key
        first_row["all_reads"] = total_reads
        return first_row[["gene_id", "transcript_id", "all_reads", "TPM"]]

    result_df = df.groupby("group_key").apply(process_group).reset_index(drop=True)

    return result_df


def run_test():

    logger.info("** Test data:\n")
    df = generate_test_data()
    print(df)

    logger.info("** Test results:\n")
    DE_results = differential_isoform_tests(df)
    print(DE_results)

    return


if __name__ == "__main__":
    main()
