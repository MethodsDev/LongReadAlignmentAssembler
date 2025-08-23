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

    # Removed min_cell_fraction filtering; fraction matrix now only used for annotation.


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
        "--save_annotated_isoform_details",
        action="store_true",
        default=False,
        help=(
            "If set, also output an isoform-level annotated table including per-isoform pi values, "
            "delta_pi, and flags indicating whether the isoform/gene was tested, used in the chi2 test, "
            "and membership in dominant/alternate sets. File: <output_prefix>.diff_iso.annotated_isoforms.tsv"
        ),
    )

    parser.add_argument(
        "--limit_clusters",
        type=str,
        default=None,
        help=(
            "Comma-separated list of EXACTLY two cluster names to restrict the comparison to. "
            "Useful for quick testing. Example: --limit_clusters clusterX,clusterY"
        ),
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
    output_prefix = args.output_prefix
    signif_threshold = args.signif_threshold
    splice_hashcode_id_mappings_file = args.splice_hashcode_id_mappings
    group_by_feature = args.group_by_feature
    save_annot = args.save_annotated_isoform_details

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, force=True)

    ## begin

    counts_big_df = pd.read_csv(sc_cluster_counts_matrix, sep="\t")

    column_names = list(counts_big_df.columns)

    assert column_names[0] == "gene_id"
    assert column_names[1] == "transcript_id"

    cluster_names = column_names[2:]

    # Optional restriction to a specific pair of clusters
    if args.limit_clusters:
        restricted = [c.strip() for c in args.limit_clusters.split(",") if c.strip()]
        if len(restricted) != 2:
            raise ValueError("--limit_clusters must specify exactly two distinct cluster names separated by a comma.")
        if restricted[0] == restricted[1]:
            raise ValueError("--limit_clusters requires two DIFFERENT cluster names.")
        missing = [c for c in restricted if c not in cluster_names]
        if missing:
            raise ValueError(f"--limit_clusters: cluster(s) not found in matrix columns: {','.join(missing)}")
        cluster_names = restricted
        logger.info(f"Restricting analysis to cluster pair: {cluster_names[0]} vs {cluster_names[1]}")

    # If using fraction filtering, load and validate the fraction matrix
    fraction_big_df = None
    if sc_cluster_fraction_matrix:
        logger.info(
            f"-loading fraction-expression matrix (annotation only): {sc_cluster_fraction_matrix}"
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
    all_annotated_isoforms = []  # collect annotated per-pair if requested
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

            # Build pairwise fraction dataframe when available
            pair_fraction_df = None
            if (
                fraction_big_df is not None
                and cluster_i in fraction_big_df.columns
                and cluster_j in fraction_big_df.columns
            ):
                frac_subset = fraction_big_df[["gene_id", "transcript_id", cluster_i, cluster_j]].copy()
                # Rename cluster-specific fraction columns to new descriptive names (fraction of cells expressing isoform)
                frac_subset.rename(columns={cluster_i: "cell_detect_frac_A", cluster_j: "cell_detect_frac_B"}, inplace=True)
                # Provide legacy column names for backward compatibility
                frac_subset["frac_A"] = frac_subset["cell_detect_frac_A"]
                frac_subset["frac_B"] = frac_subset["cell_detect_frac_B"]
                test_df = test_df.merge(
                    frac_subset, on=["gene_id", "transcript_id"], how="left"
                )
                # Provide the pairwise fraction dataframe to downstream testing (annotation only)
                pair_fraction_df = frac_subset[["gene_id", "transcript_id", "cell_detect_frac_A", "cell_detect_frac_B", "frac_A", "frac_B"]].copy()
                # fraction reporting and filtering handled in differential_isoform_tests

            # Run DTU tests for this pair
            test_df_results = None
            annotated_df = None
            ditsu_return = differential_isoform_tests(
                df=test_df,
                group_by_token=group_by_feature,
                min_reads_per_gene=min_reads_per_gene,
                min_delta_pi=min_delta_pi,
                top_isoforms_each=top_isoforms_each,
                reciprocal_delta_pi=args.reciprocal_delta_pi,
                min_reads_DTU_isoform=args.min_reads_DTU_isoform,
                fraction_df=pair_fraction_df,
                return_annotated_df=save_annot,
            )
            if save_annot:
                # ditsu_return is a tuple (results_df or None, annotated_df)
                test_df_results, annotated_df = ditsu_return
            else:
                test_df_results = ditsu_return

            if test_df_results is not None:
                test_df_results["cluster_A"] = cluster_i
                test_df_results["cluster_B"] = cluster_j

                if all_test_results is None:
                    all_test_results = test_df_results
                else:
                    all_test_results = pd.concat([all_test_results, test_df_results])

            # Handle annotated isoforms (collect even if no significant test rows produced)
            if save_annot and annotated_df is not None:
                annotated_subset = annotated_df.copy()
                annotated_subset["cluster_A"] = cluster_i
                annotated_subset["cluster_B"] = cluster_j
                # Keep ALL isoforms (even if not tested) for complete transparency.
                # Make grouping_id unique per cluster comparison by appending cluster pair.
                if "grouping_id" in annotated_subset.columns:
                    annotated_subset["grouping_id"] = (
                        annotated_subset["grouping_id"].astype(str).fillna("NA")
                        + f"|{cluster_i}|{cluster_j}"
                    )
                all_annotated_isoforms.append(annotated_subset)

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
        if save_annot and all_annotated_isoforms:
            annot_df = pd.concat(all_annotated_isoforms, ignore_index=True)
            # enforce 3-decimal rounding on all float columns prior to write
            float_cols = [c for c in annot_df.columns if annot_df[c].dtype.kind in ("f", "d")]
            if float_cols:
                annot_df[float_cols] = annot_df[float_cols].round(3)
            annot_outfile = f"{output_prefix}.diff_iso.annotated_isoforms.tsv"
            annot_df.to_csv(annot_outfile, sep="\t", index=False, float_format="%.3f")
            logger.info(f"Wrote annotated isoform details to: {annot_outfile}")
    else:
        logger.info("No results to report")

    sys.exit(0)


if __name__ == "__main__":
    main()
