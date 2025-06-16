#!/usr/bin/env python3

import os, sys, re
import logging
import argparse
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
import glob


sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../../pylib"])
)

from DiffIsoformStatTest import *

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
        "--sc_cluster_counts_matrix",
        type=str,
        required=True,
        help="sc cluster counts matrix. First two columns should be gene_id and transcript_id",
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
        "--stat_test",
        default="chi2",
        required=False,
        choices=["chi2", "fisher"],
        help="test to use: chi-squared or Fishers exact test (default: chi2)",
    )

    parser.add_argument(
        "--signif_threshold",
        default=0.001,
        help="significance threshold for stat test to mark as signfiicantly DE",
    )

    parser.add_argument(
        "-d",
        "--debug",
        action="store_true",
        default=False,
        help="debug mode - more verbose",
    )

    args = parser.parse_args()

    top_isoforms_each = args.top_isoforms_each
    min_reads_per_gene = args.min_reads_per_gene
    min_delta_pi = args.min_delta_pi
    sc_cluster_counts_matrix = args.sc_cluster_counts_matrix
    output_prefix = args.output_prefix
    stat_test = args.stat_test
    signif_threshold = args.signif_threshold

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, force=True)

    ## begin

    counts_big_df = pd.read_csv(sc_cluster_counts_matrix, sep="\t")

    column_names = list(counts_big_df.columns)

    assert column_names[0] == "gene_id"
    assert column_names[1] == "transcript_id"

    cluster_names = column_names[2:]

    ############################################################
    ## pairwise compare clusters for diff isoform usage analysis
    ############################################################

    all_test_results = None
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

            # print(test_df)

            test_df_results = differential_isoform_tests(
                test_df, min_reads_per_gene, min_delta_pi, top_isoforms_each, stat_test
            )

            if test_df_results is not None:
                test_df_results["cluster_A"] = cluster_i
                test_df_results["cluster_B"] = cluster_j

                if all_test_results is None:
                    all_test_results = test_df_results
                else:
                    all_test_results = pd.concat([all_test_results, test_df_results])

    if all_test_results is not None:
        # all_test_results = pd.concat([all_test_results, test_df_results])

        # perform mult test pvalue adjustment
        all_test_results = FDR_mult_tests_adjustment(
            all_test_results, min_delta_pi, signif_threshold
        )

        all_test_results.to_csv(
            f"{output_prefix}.diff_iso.{stat_test}.tsv",
            sep="\t",
            index=False,
        )
    else:
        logger.info("No results to report")

    sys.exit(0)


if __name__ == "__main__":
    main()
