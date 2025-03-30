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
    )

    parser.add_argument(
        "--pseudotime_expr_matrix",
        type=str,
        required=True,
        help="pseudotime expression matrix",
    )

    # logging.basicConfig(level=logging.DEBUG, force=True)

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

    # Iterate over each column as a pivot
    for pivot in range(1, len(pseudotimes)):
        pseudotime_A = pseudotimes[:pivot]  # All columns to the left
        pseudotime_B = pseudotimes[pivot:]  # All columns to the right

        test_df = pstime_df[
            ["gene_id", "transcript_id"] + pseudotime_A + pseudotime_B
        ].copy()

        test_df["count_A"] = test_df[pseudotime_A].sum(axis=1)
        test_df["count_B"] = test_df[pseudotime_B].sum(axis=1)
        test_df = test_df[["gene_id", "transcript_id", "count_A", "count_B"]]

        pseudotime_A_label = "_".join([col[:5] for col in pseudotime_A])
        pseudotime_B_label = "_".join([col[:5] for col in pseudotime_B])

        logger.info("Testing {} vs. {}".format(pseudotime_A_label, pseudotime_B_label))

        # Remove rows where total counts are zero
        test_df = test_df[test_df["count_A"] + test_df["count_B"] > 0]

        test_df_results = differential_isoform_tests(
            test_df, min_reads_per_gene, min_delta_pi, top_isoforms_each, "fisher"
        )

        if test_df_results is not None:
            test_df_results["pseudotime_A"] = pseudotime_A_label
            test_df_results["pseudotime_B"] = pseudotime_B_label

            test_df_results.to_csv(
                f"pstime.{pseudotime_A_label}-vs-{pseudotime_B_label}.diff_iso.tsv",
                sep="\t",
                index=False,
            )


if __name__ == "__main__":
    main()
