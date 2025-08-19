# Long Read Alignment Assembler (LRAA)

Isoform Discovery and/or Quantification from Long Read RNA-Seq

Visit the [LRAA wiki](https://github.com/MethodsDev/LongReadAlignmentAssembler/wiki) for user documentation


## Optional single-cell fraction-of-cells support for DTU (sc utilities)

Some single-cell workflows benefit from restricting differential transcript usage (DTU) tests to transcripts that are sufficiently expressed across cells. The utility `util/sc/sc_pseudobulk_test_isoform_DiffUsage.py` supports an optional fraction-expression matrix to annotate results and, if desired, gate which tests proceed.

Inputs
- Counts matrix (`--sc_cluster_counts_matrix`): TSV with columns
	- `gene_id`, `transcript_id`, then one column per cluster (integer counts)
- Fraction matrix (`--sc_cluster_fraction_matrix`, optional): TSV with columns
	- `gene_id`, `transcript_id`, then one column per cluster with the fraction of cells expressing the transcript (0..1)

CLI flags (subset)
- `--min_cell_fraction <float>` (default 0.0): When > 0 and a fraction matrix is provided, enforces a presence gate before testing (see below). When 0.0, the fractions are only reported.

Behavior
1) For each cluster pair (A vs B), the driver builds a pairwise fraction dataframe with columns `frac_A` and `frac_B` from the provided fraction matrix.
2) Inside `pylib/DiffIsoformStatTest.py:differential_isoform_tests`, we:
	 - Check total read count thresholds first.
	 - Select top isoforms independently by counts for A (`top_countA`) and B (`top_countB`).
	 - Presence gate (applied only if `min_cell_fraction` > 0 and fraction data is present): require that
		 - at least one transcript in `top_countA` has `frac_A >= min_cell_fraction`, and
		 - at least one transcript in `top_countB` has `frac_B >= min_cell_fraction`.
		 If this gate fails, the group is skipped before significance testing.
	 - Compute proportions and run the chi-square DTU test on selected isoforms as before.

Output columns (when a fraction matrix is provided)
- `dominant_frac_A`, `dominant_frac_B`: comma-separated fractions (3 decimals) corresponding to `dominant_transcript_ids` order in clusters A and B.
- If reciprocal mode is used:
	- `alternate_frac_A`, `alternate_frac_B`: same for `alternate_transcript_ids`.
- `min_cell_fraction`: the threshold applied for that run (0.0 if annotation-only).

Notes
- Missing fraction values are reported as `NA` and do not pass the presence gate.
- The presence gate is intentionally strict (both clusters) for biological interpretability.
- If the fraction matrix lacks a cluster, that pair falls back to annotation-only for that comparison.
