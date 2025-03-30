#!/usr/bin/env python3

import os, sys, re
import logging
import argparse
import pandas as pd
import numpy as np


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


top_isoforms_each = 1
min_delta_pi = 0.1
min_reads_per_gene = 5


def main():
    parser = argparse.ArgumentParser(
        description="run single cell cluster isoform DE usage tests via pseudobulk chi-square",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        # formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "--pseudotime_expr_matrix",
        type=str,
        required=True,
        help="pseudotime expression matrix",
    )

    args = parser.parse_args()

    pseudotime_expr_matrix = args.pseudotime_expr_matrix

    pstime_df = pd.read_csv(pseudotime_expr_matrix, sep="\t")
    pstime_df_colnames = list(pstime_df.columns)
    pseudotimes = pstime_df_colnames[1:]
    pstime_df_colnames[0] = "transcript_id"
    pstime_df.columns = pstime_df_colnames

    pstime_df["gene_id"] = pstime_df["transcript_id"].apply(
        lambda x: ":".join(x.split(":")[0:-1])
    )

    # walk through pseudotime to see what switch events are happening sequentially
    for next_i in range(1, len(pseudotimes)):
        i = next_i - 1
        pseudotime_A = pseudotimes[i]
        pseudotime_B = pseudotimes[next_i]

        test_df = pstime_df[
            ["gene_id", "transcript_id", pseudotime_A, pseudotime_B]
        ].copy()

        # shorten string
        pseudotime_A = pseudotime_A[0 : min(5, len(pseudotime_A))]
        pseudotime_B = pseudotime_B[0 : min(5, len(pseudotime_B))]

        logger.info("Testing {} vs. {}".format(pseudotime_A, pseudotime_B))

        test_df.columns = ["gene_id", "transcript_id", "count_A", "count_B"]

        # remove cases wehre rowsums = 0
        test_df = test_df[test_df["count_A"] + test_df["count_B"] > 0]

        test_df_results = differential_isoform_tests(
            test_df, min_reads_per_gene, min_delta_pi, top_isoforms_each, "fisher"
        )

        if test_df_results is not None:
            test_df_results["pseudotime_A"] = pseudotime_A
            test_df_results["pseudotime_B"] = pseudotime_B

            # all_test_results = pd.concat([all_test_results, test_df_results])
            test_df_results.to_csv(
                f"pstime.{pseudotime_A}-vs-{pseudotime_B}.diff_iso.tsv",
                sep="\t",
                index=False,
            )


if __name__ == "__main__":
    main()
