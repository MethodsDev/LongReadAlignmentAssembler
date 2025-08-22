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
    # Optional fraction filtering with pairwise dataframe having 'frac_A' and 'frac_B'
    fraction_df=None,
    min_cell_fraction=0.0,
    # Return the annotated per-isoform dataframe in addition to summary results
    return_annotated_df=False,
):

    logger = logging.getLogger(__name__)
    logger.debug("Running differential_isoform_tests()")

    results = []

    # Prepare annotated dataframe if requested
    annotated_df = None
    if return_annotated_df:
        annotated_df = df.copy()
        # Initialize annotation columns
        for col, dtype, default in [
            ("pi_A", float, np.nan),
            ("pi_B", float, np.nan),
            ("delta_pi", float, np.nan),
            ("total_counts_A_gene", float, np.nan),
            ("total_counts_B_gene", float, np.nan),
            ("gene_tested", bool, False),
            ("evaluated_isoform", bool, False),  # part of filtered_group used in chi2
            ("dominant_set", bool, False),
            ("alternate_set", bool, False),
            ("candidate_top_isoform", bool, False),  # top isoforms considered even if test aborted
            ("skip_reason", object, ""),  # reason gene not fully tested / no result row
        ]:
            if col not in annotated_df.columns:
                annotated_df[col] = default
        # Optionally add fraction columns (copy over if present externally later)
        if fraction_df is not None:
            for frac_col in ["frac_A", "frac_B"]:
                if frac_col in fraction_df.columns and frac_col not in annotated_df.columns:
                    # Map fractions by transcript_id (assumes transcript_id uniqueness)
                    frac_map = dict(zip(fraction_df["transcript_id"], fraction_df[frac_col]))
                    annotated_df[frac_col] = annotated_df["transcript_id"].map(frac_map)

    # Prepare fraction lookups if provided
    # Determine if fraction data is available for reporting vs. for filtering
    have_fraction_data = False
    use_fraction_filter = False
    frac_lookup_A = None
    frac_lookup_B = None
    if fraction_df is not None:
        required_cols = {"gene_id", "transcript_id", "frac_A", "frac_B"}
        missing = required_cols - set(fraction_df.columns)
        if missing:
            logger.warning(
                f"fraction_df is missing required columns: {','.join(sorted(missing))}. Fraction reporting/filtering disabled."
            )
        else:
            have_fraction_data = True
            # Build per-transcript lookups for fast access
            frac_lookup_A = dict(
                zip(fraction_df["transcript_id"], fraction_df["frac_A"])
            )
            frac_lookup_B = dict(
                zip(fraction_df["transcript_id"], fraction_df["frac_B"])
            )
            # Enable filtering only if threshold > 0
            use_fraction_filter = (
                isinstance(min_cell_fraction, (int, float))
                and float(min_cell_fraction) > 0.0
            )
            # Always apply 'both' logic for simplicity and biological justification
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
        # Store original total counts for reporting (can be zero leading to skip)
        original_total_counts_A = group["count_A"].sum()
        original_total_counts_B = group["count_B"].sum()

        # Calculate pi only if totals > 0 to avoid division by zero
        if original_total_counts_A > 0 and original_total_counts_B > 0:
            pi_A = group["count_A"] / original_total_counts_A
            pi_B = group["count_B"] / original_total_counts_B
            delta_pi = pi_B - pi_A
        else:
            pi_A = pd.Series([np.nan]*len(group), index=group.index)
            pi_B = pd.Series([np.nan]*len(group), index=group.index)
            delta_pi = pd.Series([np.nan]*len(group), index=group.index)

        # Add these calculations to the original group
        group = group.copy()
        group["pi_A"] = pi_A
        group["pi_B"] = pi_B
        group["delta_pi"] = delta_pi

        # If annotating, store per-isoform metrics now
        if return_annotated_df:
            annotated_df.loc[group.index, "pi_A"] = pi_A
            annotated_df.loc[group.index, "pi_B"] = pi_B
            annotated_df.loc[group.index, "delta_pi"] = delta_pi
            annotated_df.loc[group.index, "total_counts_A_gene"] = original_total_counts_A
            annotated_df.loc[group.index, "total_counts_B_gene"] = original_total_counts_B

        # Early annotation for single-isoform or zero-count cases
        if len(group) < 2:
            if return_annotated_df:
                annotated_df.loc[group.index, "skip_reason"] = "single_isoform"
            continue

        if original_total_counts_A == 0 or original_total_counts_B == 0:
            logger.debug(
                f"Total counts in condition A or B is zero for {group_by_token} {group_by_id}, skipping."
            )
            if return_annotated_df:
                annotated_df.loc[group.index, "skip_reason"] = "zero_total_counts"
            continue

        # Select top isoforms by counts from the full group
        top_countA = group.nlargest(top_isoforms_each, "count_A")
        top_countB = group.nlargest(top_isoforms_each, "count_B")
        filtered_group = pd.concat([top_countA, top_countB]).drop_duplicates()

        filtered_group["total"] = filtered_group["count_A"] + filtered_group["count_B"]
        filtered_group = filtered_group.sort_values(by="total", ascending=False).head(10)

        # Mark candidate top isoforms regardless of downstream filtering outcome
        if return_annotated_df:
            annotated_df.loc[filtered_group.index, "candidate_top_isoform"] = True

        # Now enforce min_reads_per_gene after pi computation so skipped genes retain pi annotation
        if (
            original_total_counts_A < min_reads_per_gene
            or original_total_counts_B < min_reads_per_gene
        ):
            logger.debug(
                f"{group_by_token} {group_by_id} has insufficient reads (A: {original_total_counts_A} or B: {original_total_counts_B}), skipping."
            )
            if return_annotated_df:
                annotated_df.loc[group.index, "skip_reason"] = "insufficient_gene_reads"
            continue

        # Presence-based fraction gating on the selected top isoforms
        if use_fraction_filter:
            thr = float(min_cell_fraction)

            def _ge(v):
                try:
                    return (not pd.isna(v)) and float(v) >= thr
                except Exception:
                    return False

            has_A = any(
                _ge(frac_lookup_A.get(tid, np.nan))
                for tid in top_countA["transcript_id"].tolist()
            )
            has_B = any(
                _ge(frac_lookup_B.get(tid, np.nan))
                for tid in top_countB["transcript_id"].tolist()
            )

            # Require presence in both clusters' top sets
            pass_fraction = has_A and has_B

            if not pass_fraction:
                logger.debug(
                    f"{group_by_token} {group_by_id} failed presence-based fraction gating (A:{has_A}, B:{has_B})."
                )
                if return_annotated_df:
                    annotated_df.loc[group.index, "skip_reason"] = "fraction_gating_fail"
                continue

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
            if return_annotated_df:
                annotated_df.loc[group.index, "skip_reason"] = "delta_pi_fail"
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
            if return_annotated_df:
                annotated_df.loc[group.index, "skip_reason"] = "dominant_isoform_low_reads"
            continue

        if reciprocal_delta_pi and alternate_total_reads < min_reads_DTU_isoform:
            logger.debug(
                f"{group_by_token} {group_by_id} alternate isoforms have insufficient reads ({alternate_total_reads}), skipping."
            )
            if return_annotated_df:
                annotated_df.loc[group.index, "skip_reason"] = "alternate_isoform_low_reads"
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

        # Prepare fraction reporting if available
        dominant_frac_A_str = None
        dominant_frac_B_str = None
        alternate_frac_A_str = None
        alternate_frac_B_str = None
        if have_fraction_data:
            def _format_frac_list(transcript_ids_str, lookup):
                if transcript_ids_str is None or transcript_ids_str == "":
                    return ""
                tids = [t.strip() for t in str(transcript_ids_str).split(",") if t.strip()]
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

            dominant_frac_A_str = _format_frac_list(dominant_transcript_ids_str, frac_lookup_A)
            dominant_frac_B_str = _format_frac_list(dominant_transcript_ids_str, frac_lookup_B)
            if reciprocal_delta_pi:
                alternate_frac_A_str = _format_frac_list(alternate_transcript_ids_str, frac_lookup_A)
                alternate_frac_B_str = _format_frac_list(alternate_transcript_ids_str, frac_lookup_B)

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

        # Append fraction reporting for dominant isoforms if available
        if have_fraction_data:
            result_row.extend([
                dominant_frac_A_str,
                dominant_frac_B_str,
            ])

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
            if have_fraction_data:
                result_row.extend([
                    alternate_frac_A_str,
                    alternate_frac_B_str,
                ])

        # Add status and optional min_cell_fraction for transparency
        result_row.append(status)
        if have_fraction_data:
            result_row.append(min_cell_fraction)
        results.append(result_row)

        # Update annotation flags for this gene if requested
        if return_annotated_df:
            annotated_df.loc[group.index, "gene_tested"] = True
            # Indices in filtered_group correspond to filtered_group's own index (original df indices retained during concat/drop_duplicates)
            filtered_indices = filtered_group.index
            annotated_df.loc[filtered_indices, "evaluated_isoform"] = True
            annotated_df.loc[dominant_indices, "dominant_set"] = True
            if reciprocal_delta_pi:
                annotated_df.loc[alternate_indices, "alternate_set"] = True
            # Clear skip_reason (if any) since test succeeded
            annotated_df.loc[group.index, "skip_reason"] = ""

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
        if have_fraction_data:
            columns += [
                "dominant_frac_A",
                "dominant_frac_B",
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
            if have_fraction_data:
                columns += [
                    "alternate_frac_A",
                    "alternate_frac_B",
                ]
        columns += ["status"]
        if have_fraction_data:
            columns += ["min_cell_fraction"]
        results_df = pd.DataFrame(results, columns=columns)
        if return_annotated_df:
            return results_df, annotated_df
        return results_df

    logger.debug("-didn't meet requirements to test.")
    if return_annotated_df and annotated_df is not None:
        return None, annotated_df
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
    DE_results, annotated = differential_isoform_tests(
        df, reciprocal_delta_pi=True, top_isoforms_each=1, return_annotated_df=True
    )
    DE_results = FDR_mult_tests_adjustment(DE_results)
    print(DE_results)
    print("\nAnnotated isoform dataframe (subset):")
    print(annotated.head())


if __name__ == "__main__":
    run_test()
