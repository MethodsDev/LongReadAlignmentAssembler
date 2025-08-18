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

from DiffIsoformStatTest import (
    differential_isoform_tests,
    FDR_mult_tests_adjustment,
)

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
        "--sc_cluster_fraction_matrix",
        type=str,
        required=False,
        default=None,
        help=(
            "Optional matrix of fraction of cells expressing each feature (0..1). "
            "Format should mirror counts matrix: first two columns gene_id and transcript_id, "
            "remaining columns are cluster names."
        ),
    )

    parser.add_argument(
        "--min_cell_fraction",
        type=float,
        default=0.0,
        help=(
            "Minimum fraction of cells that must express a feature in each cluster being compared. "
            "Used only if --sc_cluster_fraction_matrix is provided. Range: 0..1 (default: 0.0 disables filtering)."
        ),
    )

    parser.add_argument(
        "--splice_hashcode_id_mappings",
        type=str,
        required=True,
        help="LRAA.gene_transcript_splicehashcode.tsv.wAnnotIDs file",
    )

    parser.add_argument(
        "--group_by_feature",
        default="gene_id",
        choices=["gene_id", "splice_hashcode"],
        help="organize groupings according to gene_id or splice_hashcode",
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
        "--min_reads_DTU_isoform",
        type=int,
        default=25,
        help="min reads per DTU isoform",
    )

    parser.add_argument(
        "--min_delta_pi", type=float, default=0.1, help="min delta pi (default: 0.1)"
    )

    parser.add_argument(
        "--reciprocal_delta_pi",
        action="store_true",
        default=False,
        help="require the alt direction isoform to have at least --min_delta_pi and --min_reqds_DTU_isoform too",
    )

    parser.add_argument(
        "--ignore_unspliced",
        action="store_true",
        default=False,
        help="exclude unspliced isoforms (only works when used with splice pattern collapsed isoforms)",
    )

    parser.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="prefix for output files",
    )

    parser.add_argument(
        "--signif_threshold",
        type=float,
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
    sc_cluster_fraction_matrix = args.sc_cluster_fraction_matrix
    min_cell_fraction = args.min_cell_fraction
    output_prefix = args.output_prefix
    signif_threshold = args.signif_threshold
    splice_hashcode_id_mappings_file = args.splice_hashcode_id_mappings
    group_by_feature = args.group_by_feature

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, force=True)

    ## begin

    counts_big_df = pd.read_csv(sc_cluster_counts_matrix, sep="\t")

    column_names = list(counts_big_df.columns)

    assert column_names[0] == "gene_id"
    assert column_names[1] == "transcript_id"

    cluster_names = column_names[2:]

    # If using fraction filtering, load and validate the fraction matrix
    fraction_big_df = None
    if sc_cluster_fraction_matrix:
        if not (0.0 <= float(min_cell_fraction) <= 1.0):
            raise ValueError(
                f"--min_cell_fraction must be within [0,1], got: {min_cell_fraction}"
            )
        logger.info(
            f"-loading fraction-expression matrix: {sc_cluster_fraction_matrix} (min_cell_fraction={min_cell_fraction})"
        )
        fraction_big_df = pd.read_csv(sc_cluster_fraction_matrix, sep="\t")

        frac_cols = list(fraction_big_df.columns)
        assert frac_cols[0] == "gene_id"
        assert frac_cols[1] == "transcript_id"

        # Keep only overlapping clusters between counts and fraction matrices
        fraction_cluster_names = frac_cols[2:]
        overlapping_clusters = sorted(
            list(set(cluster_names).intersection(set(fraction_cluster_names)))
        )
        if not overlapping_clusters:
            logger.warning(
                "No overlapping cluster columns between counts and fraction matrices; disabling fraction filtering."
            )
            fraction_big_df = None
        else:
            if set(cluster_names) - set(overlapping_clusters):
                missing = sorted(list(set(cluster_names) - set(overlapping_clusters)))
                logger.warning(
                    "Fraction matrix missing some clusters present in counts: {}. "
                    "Filtering will only apply to pairs where both clusters are present.".format(
                        ",".join(missing)
                    )
                )

    ################################
    ## Incorporate splice hash codes
    ################################

    splice_hashcode_id_mappings_df = pd.read_csv(
        splice_hashcode_id_mappings_file, sep="\t"
    )

    sp_hash_mapping1 = dict(
        zip(
            splice_hashcode_id_mappings_df["transcript_id"],
            splice_hashcode_id_mappings_df["transcript_splice_hash_code"],
        )
    )
    sp_hash_mapping2 = dict()
    if "new_transcript_id" in splice_hashcode_id_mappings_df.columns:
        sp_hash_mapping2 = dict(
            zip(
                splice_hashcode_id_mappings_df["new_transcript_id"],
                splice_hashcode_id_mappings_df["new_transcript_splice_hash_code"],
            )
        )

    # include the splice code direct mappings as well in case the splice pattern collapsing already happened earlier.
    sp_hash_mapping3 = dict()
    if "transcript_splice_hash_code" in splice_hashcode_id_mappings_df.columns:

        sp_hash_mapping3 = dict(
            zip(
                splice_hashcode_id_mappings_df["transcript_splice_hash_code"],
                splice_hashcode_id_mappings_df["transcript_splice_hash_code"],
            )
        )

    sp_hash_mapping4 = dict()
    if "new_transcript_splice_hash_code" in splice_hashcode_id_mappings_df.columns:
        sp_hash_mapping4 = dict(
            zip(
                splice_hashcode_id_mappings_df["new_transcript_splice_hash_code"],
                splice_hashcode_id_mappings_df["new_transcript_splice_hash_code"],
            )
        )

    sp_hash_mapping = {
        **sp_hash_mapping1,
        **sp_hash_mapping2,
        **sp_hash_mapping3,
        **sp_hash_mapping4,
    }

    counts_big_df["splice_hashcode"] = counts_big_df["transcript_id"].map(
        sp_hash_mapping
    )

    #########
    ## Exclude unspliced if indicated
    #########

    if args.ignore_unspliced:
        logger.info("-pruning unspliced isoforms")
        logger.info("before pruning, have {} rows".format(counts_big_df.shape[0]))
        mask_to_exclude = (
            counts_big_df["splice_hashcode"].astype(str).str.contains(":iso-", na=False)
        )
        counts_big_df = counts_big_df[~mask_to_exclude]
        logger.info("after pruning, have {} rows".format(counts_big_df.shape[0]))

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
                ["gene_id", "transcript_id", "splice_hashcode", cluster_i, cluster_j]
            ].copy()
            test_df.rename(
                columns={cluster_i: "count_A", cluster_j: "count_B"}, inplace=True
            )

            test_df = test_df.loc[((test_df["count_A"] > 0) | (test_df["count_B"] > 0))]

            # Optionally filter by cell-fraction expression per cluster
            if (
                fraction_big_df is not None
                and min_cell_fraction > 0.0
                and cluster_i in fraction_big_df.columns
                and cluster_j in fraction_big_df.columns
            ):
                frac_subset = fraction_big_df[
                    ["gene_id", "transcript_id", cluster_i, cluster_j]
                ].copy()
                frac_subset.rename(
                    columns={cluster_i: "frac_A", cluster_j: "frac_B"}, inplace=True
                )
                before_n = test_df.shape[0]
                test_df = test_df.merge(
                    frac_subset, on=["gene_id", "transcript_id"], how="left"
                )
                # Keep rows with sufficient fraction in both clusters
                test_df = test_df.loc[
                    (test_df["frac_A"].fillna(0.0) >= min_cell_fraction)
                    & (test_df["frac_B"].fillna(0.0) >= min_cell_fraction)
                ].copy()
                after_n = test_df.shape[0]
                # Drop helper columns
                test_df.drop(columns=["frac_A", "frac_B"], inplace=True)
                logger.info(
                    f"Applied fraction filter ({cluster_i} vs {cluster_j}): {before_n} -> {after_n} features"
                )
            elif sc_cluster_fraction_matrix and min_cell_fraction > 0.0:
                logger.warning(
                    f"Skipping fraction filter for pair {cluster_i} vs {cluster_j} (missing columns in fraction matrix)."
                )

            # print(test_df)

            test_df_results = differential_isoform_tests(
                df=test_df,
                group_by_token=group_by_feature,
                min_reads_per_gene=min_reads_per_gene,
                min_delta_pi=min_delta_pi,
                top_isoforms_each=top_isoforms_each,
                reciprocal_delta_pi=args.reciprocal_delta_pi,
                min_reads_DTU_isoform=args.min_reads_DTU_isoform,
            )

            if test_df_results is not None:
                test_df_results["cluster_A"] = cluster_i
                test_df_results["cluster_B"] = cluster_j

                # If fraction data is available for both clusters, add fraction columns
                if (
                    fraction_big_df is not None
                    and cluster_i in fraction_big_df.columns
                    and cluster_j in fraction_big_df.columns
                ):
                    # Build quick lookup dicts for transcript -> fraction per cluster
                    frac_lookup_A = dict(
                        zip(
                            fraction_big_df["transcript_id"],
                            fraction_big_df[cluster_i],
                        )
                    )
                    frac_lookup_B = dict(
                        zip(
                            fraction_big_df["transcript_id"],
                            fraction_big_df[cluster_j],
                        )
                    )

                    def _format_frac_list(transcript_ids_str, lookup):
                        if pd.isna(transcript_ids_str) or transcript_ids_str == "":
                            return ""
                        tids = [
                            t.strip()
                            for t in str(transcript_ids_str).split(",")
                            if t.strip()
                        ]
                        vals = []
                        for t in tids:
                            v = lookup.get(t, np.nan)
                            if pd.isna(v):
                                vals.append("NA")
                            else:
                                try:
                                    vals.append(f"{float(v):.3f}")
                                except Exception:
                                    vals.append("NA")
                        return ",".join(vals)

                    # Dominant isoform fractions
                    test_df_results["dominant_frac_A"] = test_df_results[
                        "dominant_transcript_ids"
                    ].apply(lambda s: _format_frac_list(s, frac_lookup_A))
                    test_df_results["dominant_frac_B"] = test_df_results[
                        "dominant_transcript_ids"
                    ].apply(lambda s: _format_frac_list(s, frac_lookup_B))

                    # Alternate isoform fractions (if present)
                    if "alternate_transcript_ids" in test_df_results.columns:
                        test_df_results["alternate_frac_A"] = test_df_results[
                            "alternate_transcript_ids"
                        ].apply(lambda s: _format_frac_list(s, frac_lookup_A))
                        test_df_results["alternate_frac_B"] = test_df_results[
                            "alternate_transcript_ids"
                        ].apply(lambda s: _format_frac_list(s, frac_lookup_B))

                    # Record the threshold used (same across this pair)
                    test_df_results["min_cell_fraction"] = min_cell_fraction

                if all_test_results is None:
                    all_test_results = test_df_results
                else:
                    all_test_results = pd.concat([all_test_results, test_df_results])

    if all_test_results is not None:
        # all_test_results = pd.concat([all_test_results, test_df_results])

        # perform mult test pvalue adjustment and set significance status
        all_test_results = FDR_mult_tests_adjustment(
            all_test_results, signif_threshold, min_delta_pi
        )

        all_test_results.to_csv(
            f"{output_prefix}.diff_iso.tsv",
            sep="\t",
            index=False,
        )
    else:
        logger.info("No results to report")

    sys.exit(0)


if __name__ == "__main__":
    main()
