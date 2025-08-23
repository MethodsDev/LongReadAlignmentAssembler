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
    # Optional per-isoform fraction reporting dataframe. Historically columns named 'frac_A','frac_B'
    # representing FRACTION OF CELLS with detected expression for that isoform in each condition.
    # New preferred column names: 'cell_detect_frac_A','cell_detect_frac_B'. Either naming is accepted.
    fraction_df=None,
    min_cell_fraction=0.0,
    # Return the annotated per-isoform dataframe in addition to summary results
    return_annotated_df=False,
    # Fixed decimal places for numeric outputs in returned tables
    output_decimal_places=3,
    # When True add clearer, more descriptive duplicate columns (leaving legacy names for compatibility)
    add_descriptive_fraction_columns=True,
):
    """Perform differential isoform usage (DTU) testing between two conditions.

    Two fraction concepts:
      (a) Isoform expression fractions within a gene (pi_A / pi_B).
      (b) Cell-detection fractions across cells (cell_detect_frac_A / B).
    Legacy names (frac_A/frac_B, delta_pi) are retained as aliases.
    """

    logger = logging.getLogger(__name__)
    logger.debug("Running differential_isoform_tests()")

    results = []

    # Annotated DF preparation
    annotated_df = None
    if return_annotated_df:
        annotated_df = df.copy()
        add_cols = [
            ("pi_A", float, np.nan),
            ("pi_B", float, np.nan),
            ("delta_pi", float, np.nan),
            ("isoform_expr_frac_A", float, np.nan),
            ("isoform_expr_frac_B", float, np.nan),
            ("delta_isoform_expr_frac", float, np.nan),
            ("total_counts_A_gene", float, np.nan),
            ("total_counts_B_gene", float, np.nan),
            ("gene_tested", bool, False),
            ("evaluated_isoform", bool, False),
            ("dominant_set", bool, False),
            ("alternate_set", bool, False),
            ("candidate_top_isoform", bool, False),
            ("skip_reason", object, "."),
            ("grouping_id", object, None),
        ]
        for col, _, default in add_cols:
            if col not in annotated_df.columns:
                annotated_df[col] = default
        # Map cell-detection fractions if provided
        if fraction_df is not None:
            legacy_to_new = {"frac_A": "cell_detect_frac_A", "frac_B": "cell_detect_frac_B"}
            for src in ["frac_A", "frac_B", "cell_detect_frac_A", "cell_detect_frac_B"]:
                if src in fraction_df.columns:
                    dest = legacy_to_new.get(src, src)
                    if dest not in annotated_df.columns:
                        m = dict(zip(fraction_df["transcript_id"], fraction_df[src]))
                        annotated_df[dest] = annotated_df["transcript_id"].map(m)

    # Fraction lookups
    have_fraction_data = False
    frac_lookup_A = frac_lookup_B = None
    if fraction_df is not None:
        present = set(fraction_df.columns)
        if {"frac_A", "frac_B"}.issubset(present) and not {"cell_detect_frac_A", "cell_detect_frac_B"}.issubset(present):
            fraction_df = fraction_df.copy()
            fraction_df["cell_detect_frac_A"] = fraction_df["frac_A"]
            fraction_df["cell_detect_frac_B"] = fraction_df["frac_B"]
        required = {"gene_id", "transcript_id", "cell_detect_frac_A", "cell_detect_frac_B"}
        missing = required - set(fraction_df.columns)
        if missing:
            logger.warning("fraction_df missing columns: %s -- disabling fraction annotation", ",".join(sorted(missing)))
        else:
            have_fraction_data = True
            frac_lookup_A = dict(zip(fraction_df["transcript_id"], fraction_df["cell_detect_frac_A"]))
            frac_lookup_B = dict(zip(fraction_df["transcript_id"], fraction_df["cell_detect_frac_B"]))

    grouped = df.groupby(group_by_token)
    num_groups = len(grouped)
    debug_mode = logger.getEffectiveLevel() == logging.DEBUG

    for group_counter, (group_by_id, group) in enumerate(grouped, start=1):

        if show_progress_monitor and group_counter % 1000 == 0:
            pct = 100 * group_counter / num_groups
            print(f"\r[{group_counter}/{num_groups}] = {pct:.2f}% done   ", file=sys.stderr, end="")

        if debug_mode:
            logger.debug("Processing %s=%s", group_by_token, group_by_id)

        original_total_counts_A = group["count_A"].sum()
        original_total_counts_B = group["count_B"].sum()

        if original_total_counts_A > 0 and original_total_counts_B > 0:
            pi_A = group["count_A"] / original_total_counts_A
            pi_B = group["count_B"] / original_total_counts_B
            delta_pi = pi_B - pi_A
        else:
            pi_A = pd.Series([np.nan] * len(group), index=group.index)
            pi_B = pd.Series([np.nan] * len(group), index=group.index)
            delta_pi = pd.Series([np.nan] * len(group), index=group.index)

        group = group.copy()
        group["pi_A"] = pi_A
        group["pi_B"] = pi_B
        group["delta_pi"] = delta_pi

        if return_annotated_df:
            annotated_df.loc[group.index, ["pi_A", "isoform_expr_frac_A"]] = np.column_stack([pi_A, pi_A])
            annotated_df.loc[group.index, ["pi_B", "isoform_expr_frac_B"]] = np.column_stack([pi_B, pi_B])
            annotated_df.loc[group.index, ["delta_pi", "delta_isoform_expr_frac"]] = np.column_stack([delta_pi, delta_pi])
            annotated_df.loc[group.index, "total_counts_A_gene"] = original_total_counts_A
            annotated_df.loc[group.index, "total_counts_B_gene"] = original_total_counts_B
            annotated_df.loc[group.index, "grouping_id"] = group_by_id

        if len(group) < 2:
            if return_annotated_df:
                annotated_df.loc[group.index, "skip_reason"] = "single_isoform"
            continue

        if original_total_counts_A == 0 or original_total_counts_B == 0:
            if return_annotated_df:
                annotated_df.loc[group.index, "skip_reason"] = "zero_total_counts"
            continue

        top_countA = group.nlargest(top_isoforms_each, "count_A")
        top_countB = group.nlargest(top_isoforms_each, "count_B")
        filtered_group = pd.concat([top_countA, top_countB]).drop_duplicates()

        filtered_group["total"] = filtered_group["count_A"] + filtered_group["count_B"]
        filtered_group = filtered_group.sort_values("total", ascending=False).head(10)

        if return_annotated_df:
            annotated_df.loc[filtered_group.index, "candidate_top_isoform"] = True

        if original_total_counts_A < min_reads_per_gene or original_total_counts_B < min_reads_per_gene:
            if return_annotated_df:
                annotated_df.loc[group.index, "skip_reason"] = "insufficient_gene_reads"
            continue

        # Fraction-based gating disabled; we only annotate fractions now.

        filtered_delta_pi = filtered_group["delta_pi"]

        positive_indices = filtered_delta_pi[filtered_delta_pi > 0].sort_values(ascending=False).index[:2]
        negative_indices = filtered_delta_pi[filtered_delta_pi < 0].sort_values().index[:2]

        positive_sum = filtered_delta_pi.loc[positive_indices].sum()
        negative_sum = filtered_delta_pi.loc[negative_indices].sum()

        if reciprocal_delta_pi:
            pass_delta_pi = (abs(positive_sum) > min_delta_pi) and (abs(negative_sum) > min_delta_pi)
        else:
            pass_delta_pi = (abs(positive_sum) > min_delta_pi) or (abs(negative_sum) > min_delta_pi)

        if not pass_delta_pi:
            if return_annotated_df:
                annotated_df.loc[group.index, "skip_reason"] = "delta_pi_fail"
            continue

        if abs(positive_sum) > abs(negative_sum):
            dominant_delta_pi, dominant_indices = positive_sum, positive_indices
            alternate_delta_pi, alternate_indices = negative_sum, negative_indices
        else:
            dominant_delta_pi, dominant_indices = negative_sum, negative_indices
            alternate_delta_pi, alternate_indices = positive_sum, positive_indices

        dominant_transcript_ids_str = ",".join(filtered_group.loc[dominant_indices, "transcript_id"].tolist())
        dominant_counts_A = filtered_group.loc[dominant_indices, "count_A"].sum()
        dominant_counts_B = filtered_group.loc[dominant_indices, "count_B"].sum()
        dominant_total_reads = dominant_counts_A + dominant_counts_B

        dominant_pi_A_values = filtered_group.loc[dominant_indices, "pi_A"].tolist()
        dominant_pi_B_values = filtered_group.loc[dominant_indices, "pi_B"].tolist()
        dominant_pi_A_str = ",".join([f"{pi:.{delta_pi_precision}f}" for pi in dominant_pi_A_values])
        dominant_pi_B_str = ",".join([f"{pi:.{delta_pi_precision}f}" for pi in dominant_pi_B_values])

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

        alternate_transcript_ids_str = ",".join(filtered_group.loc[alternate_indices, "transcript_id"].tolist())
        alternate_counts_A = filtered_group.loc[alternate_indices, "count_A"].sum()
        alternate_counts_B = filtered_group.loc[alternate_indices, "count_B"].sum()
        alternate_total_reads = alternate_counts_A + alternate_counts_B

        alternate_pi_A_values = filtered_group.loc[alternate_indices, "pi_A"].tolist()
        alternate_pi_B_values = filtered_group.loc[alternate_indices, "pi_B"].tolist()
        alternate_pi_A_str = ",".join([f"{pi:.{delta_pi_precision}f}" for pi in alternate_pi_A_values])
        alternate_pi_B_str = ",".join([f"{pi:.{delta_pi_precision}f}" for pi in alternate_pi_B_values])

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

        if dominant_total_reads < min_reads_DTU_isoform:
            if return_annotated_df:
                annotated_df.loc[group.index, "skip_reason"] = "dominant_isoform_low_reads"
            continue
        if reciprocal_delta_pi and alternate_total_reads < min_reads_DTU_isoform:
            if return_annotated_df:
                annotated_df.loc[group.index, "skip_reason"] = "alternate_isoform_low_reads"
            continue

        # Use filtered group for chi-squared test matrix
        matrix = filtered_group[["count_A", "count_B"]].values

        pvalue = None
        status = "OK"
        try:
            chi2, pvalue, _, _ = chi2_contingency(matrix)
        except Exception:
            status = "failed"

        dominant_delta_pi_rounded = round(dominant_delta_pi, output_decimal_places)
        alternate_delta_pi_rounded = round(alternate_delta_pi, output_decimal_places) if reciprocal_delta_pi else None

        dominant_cell_detect_frac_A_str = dominant_cell_detect_frac_B_str = None
        alternate_cell_detect_frac_A_str = alternate_cell_detect_frac_B_str = None
        if have_fraction_data:
            def _fmt(ids_str, lookup):
                if not ids_str:
                    return ""
                vals = []
                for t in [t.strip() for t in ids_str.split(",") if t.strip()]:
                    v = lookup.get(t, np.nan)
                    if pd.isna(v):
                        vals.append("NA")
                    else:
                        vals.append(f"{float(v):.3f}")
                return ",".join(vals)

            dominant_cell_detect_frac_A_str = _fmt(dominant_transcript_ids_str, frac_lookup_A)
            dominant_cell_detect_frac_B_str = _fmt(dominant_transcript_ids_str, frac_lookup_B)
            if reciprocal_delta_pi:
                alternate_cell_detect_frac_A_str = _fmt(alternate_transcript_ids_str, frac_lookup_A)
                alternate_cell_detect_frac_B_str = _fmt(alternate_transcript_ids_str, frac_lookup_B)

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
            original_total_counts_A,
            original_total_counts_B,
            dominant_counts_A,
            dominant_counts_B,
        ]

        # Append fraction reporting for dominant isoforms if available
        if have_fraction_data:
            result_row.extend([
                dominant_cell_detect_frac_A_str,
                dominant_cell_detect_frac_B_str,
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
                    alternate_cell_detect_frac_A_str,
                    alternate_cell_detect_frac_B_str,
                ])

        # Add status and optional min_cell_fraction for transparency
        result_row.append(status)
        if have_fraction_data:
            result_row.append(min_cell_fraction)
        results.append(result_row)

        # Update annotation flags for this gene if requested
        if return_annotated_df:
            annotated_df.loc[group.index, "gene_tested"] = True
            annotated_df.loc[filtered_group.index, "evaluated_isoform"] = True
            annotated_df.loc[dominant_indices, "dominant_set"] = True
            if reciprocal_delta_pi:
                annotated_df.loc[alternate_indices, "alternate_set"] = True
            annotated_df.loc[group.index, "skip_reason"] = "."

    if not results:
        if return_annotated_df and annotated_df is not None:
            return None, annotated_df
        return None

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
        columns += ["dominant_cell_detect_frac_A", "dominant_cell_detect_frac_B"]
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
            columns += ["alternate_cell_detect_frac_A", "alternate_cell_detect_frac_B"]
    columns += ["status"]
    if have_fraction_data:
        columns += ["min_cell_fraction"]
    results_df = pd.DataFrame(results, columns=columns)

    # Add legacy aliases
    if have_fraction_data:
        if "dominant_cell_detect_frac_A" in results_df.columns:
            results_df["dominant_frac_A"] = results_df["dominant_cell_detect_frac_A"]
            results_df["dominant_frac_B"] = results_df["dominant_cell_detect_frac_B"]
        if reciprocal_delta_pi and "alternate_cell_detect_frac_A" in results_df.columns:
            results_df["alternate_frac_A"] = results_df["alternate_cell_detect_frac_A"]
            results_df["alternate_frac_B"] = results_df["alternate_cell_detect_frac_B"]
    # Expression fraction aliases
    results_df["isoform_expr_frac_delta"] = results_df["delta_pi"]
    results_df["dominant_isoform_expr_frac_A"] = results_df["dominant_pi_A"]
    results_df["dominant_isoform_expr_frac_B"] = results_df["dominant_pi_B"]
    if reciprocal_delta_pi:
        results_df["alternate_isoform_expr_frac_A"] = results_df["alternate_pi_A"]
        results_df["alternate_isoform_expr_frac_B"] = results_df["alternate_pi_B"]

    # Rounding
    float_cols = [c for c in results_df.columns if results_df[c].dtype.kind in ("f", "d")]
    if float_cols:
        results_df[float_cols] = results_df[float_cols].round(output_decimal_places)
    if return_annotated_df and annotated_df is not None:
        float_cols_ann = [c for c in annotated_df.columns if annotated_df[c].dtype.kind in ("f", "d")]
        if float_cols_ann:
            annotated_df[float_cols_ann] = np.round(annotated_df[float_cols_ann], output_decimal_places)
        return results_df, annotated_df
    return results_df


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
    # Standardize numeric precision (3 decimals)
    for col in ["pvalue", "adj_pvalue", "delta_pi"]:
        if col in df.columns:
            df[col] = df[col].round(3)
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
