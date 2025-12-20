#!/usr/bin/env python3

## Method for differential isoform usage testing inspired by:
## Jogelkar et al. (2021) "A spatially resolved brain region- and cell type-specific isoform 
## atlas of the postnatal mouse brain"
## Nature Communications, volume 12, Article number: 463
## https://www.nature.com/articles/s41467-020-20343-5
## Extended with reciprocal delta-pi testing, cell detection fraction filtering, and 
## additional quality control criteria.

import os, sys, re
import logging
import argparse
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
from concurrent.futures import ProcessPoolExecutor, as_completed


sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../../../pylib"])
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


def process_cluster_pair(params):
    """
    Worker function to process a single cluster pair comparison.
    Must be at module level for multiprocessing pickling.
    
    Args:
        params: tuple of (cluster_i, cluster_j, counts_big_df, fraction_big_df, 
                         group_by_feature, min_reads_per_gene, min_delta_pi, 
                         top_isoforms_each, reciprocal_delta_pi, min_reads_DTU_isoform,
                         save_annot, min_cell_fraction)
    
    Returns:
        tuple: (test_df_results, annotated_subset)
    """
    (cluster_i, cluster_j, counts_big_df, fraction_big_df, group_by_feature,
     min_reads_per_gene, min_delta_pi, top_isoforms_each, reciprocal_delta_pi,
     min_reads_DTU_isoform, save_annot, min_cell_fraction) = params
    
    logger.info("Testing pair: {} vs {}".format(cluster_i, cluster_j))
    # Build column list: always include gene_id, transcript_id, splice_hashcode, and cluster columns
    cols_to_select = ["gene_id", "transcript_id", "splice_hashcode", cluster_i, cluster_j]
    # Add gene_symbol if it exists in the dataframe
    if "gene_symbol" in counts_big_df.columns:
        cols_to_select.insert(2, "gene_symbol")  # Insert after gene_id and transcript_id
    
    test_df = counts_big_df[cols_to_select].copy()
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
        reciprocal_delta_pi=reciprocal_delta_pi,
        min_reads_DTU_isoform=min_reads_DTU_isoform,
        fraction_df=pair_fraction_df,
        return_annotated_df=save_annot,
        min_cell_fraction=min_cell_fraction,
    )
    if save_annot:
        # ditsu_return is a tuple (results_df or None, annotated_df)
        test_df_results, annotated_df = ditsu_return
    else:
        test_df_results = ditsu_return

    if test_df_results is not None:
        test_df_results["cluster_A"] = cluster_i
        test_df_results["cluster_B"] = cluster_j

    # Handle annotated isoforms (collect even if no significant test rows produced)
    annotated_subset = None
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
    
    return test_df_results, annotated_subset


def main():
    parser = argparse.ArgumentParser(
        description="run single cell cluster isoform DE usage tests via pseudobulk chi-square",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        # formatter_class=argparse.RawTextHelpFormatter,
    )

    # New option to use resources directory
    parser.add_argument(
        "--resources_dir",
        type=str,
        required=False,
        default=None,
        help=(
            "Path to resources directory (e.g., <sparseM_dir>.resources/) generated by "
            "prep_sparse_matrix_for_diff_iso_usage_analysis.py. When specified, automatically "
            "locates required input files. Overrides individual --sc_cluster_counts_matrix, "
            "--sc_cluster_fraction_matrix file specifications."
        ),
    )

    parser.add_argument(
        "--sc_cluster_counts_matrix",
        type=str,
        required=False,
        default=None,
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
        help="gene_transcript_splicehashcode file (same as used with --gene_transcript_splicehash in prep script)",
    )

    parser.add_argument(
        "--group_by_feature",
        default="gene_id",
        choices=["gene_id", "splice_hashcode", "gene_symbol"],
        help="organize groupings according to gene_id, splice_hashcode, or gene_symbol (extracted from gene_id by splitting on '^')",
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
        "--min_cell_fraction",
        type=float,
        default=0.05,
        help=(
            "Minimum cell-detection fraction required for the dominant isoform set in condition A. "
            "If --reciprocal_delta_pi is set, also require alternate isoform set in condition B to meet this."
        ),
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
        "--CPU",
        type=int,
        default=1,
        help="number of processes to use for parallel cluster comparisons (default: 1)",
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
    min_cell_fraction = args.min_cell_fraction
    output_prefix = args.output_prefix
    signif_threshold = args.signif_threshold
    group_by_feature = args.group_by_feature
    save_annot = args.save_annotated_isoform_details

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, force=True)

    ## Handle resources_dir option to auto-locate input files
    if args.resources_dir:
        resources_dir = os.path.abspath(args.resources_dir)
        if not os.path.exists(resources_dir):
            logger.error(f"Resources directory not found: {resources_dir}")
            sys.exit(1)
        
        logger.info(f"Using resources directory: {resources_dir}")
        
        # Determine the basename from resources_dir (e.g., "foo.resources" -> "foo")
        # The prep script creates: <sparseM_dir>.resources/
        # So resources_dir basename should end with ".resources"
        resources_basename = os.path.basename(resources_dir)
        if resources_basename.endswith(".resources"):
            sparseM_basename = resources_basename[:-len(".resources")]
        else:
            # If it doesn't end with .resources, just use the basename as-is
            sparseM_basename = resources_basename
        
        # Build expected filenames based on naming convention from prep script
        expected_counts_file = os.path.join(resources_dir, f"{sparseM_basename}.cluster_pseudobulk.matrix.for_diff_iso_usage")
        expected_fractions_file = os.path.join(resources_dir, f"{sparseM_basename}.cell_fractions_expressed.matrix.for_diff_iso_usage")
        
        # Use explicit file if provided, otherwise use expected file
        sc_cluster_counts_matrix = args.sc_cluster_counts_matrix or expected_counts_file
        if not os.path.exists(sc_cluster_counts_matrix):
            logger.error(f"Expected counts matrix not found: {sc_cluster_counts_matrix}")
            sys.exit(1)
        logger.info(f"Using counts matrix: {sc_cluster_counts_matrix}")
        
        # Fractions matrix is optional
        sc_cluster_fraction_matrix = args.sc_cluster_fraction_matrix or expected_fractions_file
        if not os.path.exists(sc_cluster_fraction_matrix):
            logger.warning(f"Expected fractions matrix not found: {sc_cluster_fraction_matrix}")
            sc_cluster_fraction_matrix = None
        else:
            logger.info(f"Using fractions matrix: {sc_cluster_fraction_matrix}")
        
        splice_hashcode_id_mappings_file = args.splice_hashcode_id_mappings
    else:
        # Traditional mode: require explicit file paths
        sc_cluster_counts_matrix = args.sc_cluster_counts_matrix
        sc_cluster_fraction_matrix = args.sc_cluster_fraction_matrix
        splice_hashcode_id_mappings_file = args.splice_hashcode_id_mappings
        
        if not sc_cluster_counts_matrix:
            logger.error("--sc_cluster_counts_matrix is required when --resources_dir is not provided")
            sys.exit(1)

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

    # Extract gene_symbol if needed for grouping
    if group_by_feature == "gene_symbol":
        logger.info("-extracting gene_symbol from gene_id by splitting on '^'")
        counts_big_df["gene_symbol"] = counts_big_df["gene_id"].str.split("^").str[0]
        logger.info(f"Extracted {counts_big_df['gene_symbol'].nunique()} unique gene symbols")

    #########
    ## Exclude unspliced if indicated
    #########

    if args.ignore_unspliced:
        logger.info("-pruning unspliced isoforms")
        logger.info("before pruning counts matrix, have {} rows".format(counts_big_df.shape[0]))
        mask_to_exclude = (
            counts_big_df["splice_hashcode"].astype(str).str.contains(":iso-", na=False)
        )
        counts_big_df = counts_big_df[~mask_to_exclude]
        logger.info("after pruning counts matrix, have {} rows".format(counts_big_df.shape[0]))
        
        # Also prune unspliced from fraction matrix if present
        if fraction_big_df is not None:
            logger.info("before pruning fraction matrix, have {} rows".format(fraction_big_df.shape[0]))
            # Add splice_hashcode to fraction_big_df for filtering
            fraction_big_df["splice_hashcode"] = fraction_big_df["transcript_id"].map(sp_hash_mapping)
            mask_to_exclude_frac = (
                fraction_big_df["splice_hashcode"].astype(str).str.contains(":iso-", na=False)
            )
            fraction_big_df = fraction_big_df[~mask_to_exclude_frac]
            logger.info("after pruning fraction matrix, have {} rows".format(fraction_big_df.shape[0]))

    ############################################################
    ## pairwise compare clusters for diff isoform usage analysis
    ############################################################

    # Generate list of cluster pairs to process
    cluster_pairs = []
    for i in range(len(cluster_names)):
        for j in range(i + 1, len(cluster_names)):
            cluster_pairs.append((cluster_names[i], cluster_names[j]))
    
    logger.info(f"Processing {len(cluster_pairs)} cluster pairs using {args.CPU} process(es)")
    
    all_test_results = None
    all_annotated_isoforms = []  # collect annotated per-pair if requested
    
    # Process cluster pairs with multiprocessing for true parallelism
    if args.CPU > 1:
        # Prepare parameters for each cluster pair
        params_list = [
            (cluster_i, cluster_j, counts_big_df, fraction_big_df, group_by_feature,
             min_reads_per_gene, min_delta_pi, top_isoforms_each, args.reciprocal_delta_pi,
             args.min_reads_DTU_isoform, save_annot, min_cell_fraction)
            for cluster_i, cluster_j in cluster_pairs
        ]
        
        with ProcessPoolExecutor(max_workers=args.CPU) as executor:
            # Submit all tasks
            future_to_pair = {
                executor.submit(process_cluster_pair, params): (params[0], params[1])
                for params in params_list
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_pair):
                cluster_i, cluster_j = future_to_pair[future]
                try:
                    test_df_results, annotated_subset = future.result()
                    
                    if test_df_results is not None:
                        if all_test_results is None:
                            all_test_results = test_df_results
                        else:
                            all_test_results = pd.concat([all_test_results, test_df_results])
                    
                    if annotated_subset is not None:
                        all_annotated_isoforms.append(annotated_subset)
                        
                except Exception as e:
                    logger.error(f"Error processing pair {cluster_i} vs {cluster_j}: {e}")
                    raise
    else:
        # Single-threaded mode (original behavior)
        for cluster_i, cluster_j in cluster_pairs:
            params = (cluster_i, cluster_j, counts_big_df, fraction_big_df, group_by_feature,
                     min_reads_per_gene, min_delta_pi, top_isoforms_each, args.reciprocal_delta_pi,
                     args.min_reads_DTU_isoform, save_annot, min_cell_fraction)
            test_df_results, annotated_subset = process_cluster_pair(params)
            
            if test_df_results is not None:
                if all_test_results is None:
                    all_test_results = test_df_results
                else:
                    all_test_results = pd.concat([all_test_results, test_df_results])
            
            if annotated_subset is not None:
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
