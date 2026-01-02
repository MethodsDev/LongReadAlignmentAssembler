# Differential Isoform Usage Testing Methods

## Overview

Differential isoform usage (DIU) analysis identifies genes where the relative proportions of isoforms change between conditions or cell populations, even when overall gene expression remains constant. LRAA provides a pseudobulk chi-squared testing framework specifically designed for single-cell long-read RNA-seq data to detect reciprocal isoform switching events between cell clusters.

Our approach was inspired by methods described in Jogelkar et al. (2021), "A spatially resolved brain region- and cell type-specific isoform atlas of the postnatal mouse brain" (Nature Communications, volume 12, Article number: 463; https://www.nature.com/articles/s41467-020-20343-5). We incorporated extensions including reciprocal delta-π testing to identify symmetric isoform switches, cell detection fraction filtering to ensure isoforms are broadly expressed within clusters, and additional quality control criteria for read depth and effect size thresholds.

## Pseudobulk Aggregation

Read counts for each isoform are aggregated across all cells within each cluster to generate pseudobulk count matrices. This aggregation strategy increases statistical power by pooling sparse single-cell observations while preserving biological variability at the cluster level. Pairwise comparisons are performed between all cluster combinations to comprehensively assess DIU across the cellular landscape.

## Isoform Expression Fractions and Delta-π

For each gene, isoform expression fractions (π) are calculated within each cluster as the proportion of the gene's total read counts attributed to each isoform:

$$\pi_{\text{isoform}} = \frac{\text{count}_{\text{isoform}}}{\sum \text{count}_{\text{all isoforms in gene}}}$$

The change in isoform usage between clusters is quantified using delta-π (Δπ), representing the difference in isoform expression fractions:

$$\Delta \pi = \pi_{\text{cluster B}} - \pi_{\text{cluster A}}$$

Positive Δπ values indicate increased isoform usage in cluster B relative to cluster A, while negative values indicate decreased usage.

## Reciprocal Isoform Switching

To identify biologically meaningful isoform switching events, LRAA focuses on genes exhibiting reciprocal changes in isoform usage—where one or more isoforms increase in proportion while alternative isoforms decrease. Isoforms are partitioned into two sets based on the direction and magnitude of Δπ:

- **Dominant set**: The top isoforms (by total read count across clusters) with the largest absolute Δπ in one direction (either positive or negative). These represent the isoforms whose usage changes most substantially between clusters.

- **Alternate set**: Isoforms with substantial Δπ in the opposite direction from the dominant set. These represent the reciprocal isoforms whose usage compensates for the dominant set changes.

The reciprocal testing mode requires both sets to exhibit biologically meaningful changes, ensuring that detected DIU events reflect genuine isoform switching rather than spurious differences in low-abundance transcripts.

## Filtering Criteria

Genes are tested for DIU only if they meet stringent quality and biological relevance criteria:

### Read Depth Requirements
- **Minimum gene-level reads**: At least 25 total reads per gene in each cluster being compared. This ensures sufficient statistical power for chi-squared testing.
- **Minimum isoform-level reads**: Each isoform set (dominant and alternate) must be supported by at least 25 reads across the two clusters. This prevents testing of lowly-expressed isoforms prone to technical noise.

### Delta-π Thresholds
- **Minimum delta-π magnitude**: Both dominant and alternate isoform sets must exhibit |Δπ| ≥ 0.1 (10 percentage points). This threshold ensures that detected switches represent substantial shifts in isoform proportions rather than minor fluctuations.
- **Reciprocal requirement**: When enabled, the alternate set must meet the same delta-π threshold as the dominant set, enforcing symmetric isoform switching.

### Cell Detection Fraction
Cell detection fraction quantifies the proportion of cells within a cluster where an isoform is detected (has non-zero read count). This metric distinguishes broadly-expressed isoforms from those present in only a small subset of cells.

For DIU testing, the dominant isoform set must be detected in at least a minimum fraction of cells (e.g., $$5\%$$) in the cluster where it is enriched (higher $$\pi$$ value). When reciprocal testing is enabled, the alternate isoform set must also meet this criterion in its enriched cluster. This filter ensures that detected DIU events are driven by consistent isoform expression patterns across the cell population rather than stochastic detection in rare cells. If the pseudobulk resources lack a cell detection fraction matrix, this filter is skipped automatically, and the detection-fraction columns in the output remain empty.

### Transcript Selection
To focus analysis on major isoform switches and reduce computational burden, only the most highly-expressed isoforms from each cluster are considered. The `--top_isoforms_each` parameter controls how many top-ranked isoforms (by read count) per cluster enter the comparison. Setting this to 1 restricts testing to the single most abundant isoform from each cluster, identifying dominant isoform switches while excluding minor variants.

### Splice Pattern Filtering
Unspliced isoforms can be excluded from analysis (`--ignore_unspliced` flag) to focus exclusively on spliced transcript variants. This is particularly useful when working with splice pattern-collapsed isoform sets, where the primary interest is alternative splicing rather than intron retention or unprocessed transcripts.

## Statistical Testing

DIU significance is assessed using chi-squared contingency tests. For each gene meeting the filtering criteria, a 2×N contingency table is constructed, where rows represent the two clusters being compared and columns represent the N isoforms included in the test (typically from the filtered top isoforms and their reciprocal counterparts). Cell values are the pseudobulk read counts for each isoform in each cluster.

The chi-squared test evaluates whether the distribution of reads across isoforms differs significantly between clusters:

H<sub>0</sub>: Isoform proportions are identical between clusters  
H<sub>A</sub>: Isoform proportions differ between clusters

The test statistic follows a chi-squared distribution with (rows-1) × (columns-1) degrees of freedom. P-values are computed directly from this distribution.

## Multiple Testing Correction

Because DIU testing is performed for potentially thousands of genes across multiple cluster pairs, raw p-values are adjusted for multiple testing using the Benjamini-Hochberg false discovery rate (FDR) procedure. This controls the expected proportion of false discoveries among rejected null hypotheses.

Genes are considered significantly differentially used if they meet both statistical and effect size criteria:
- FDR-adjusted $$p$$-value < $$0.001$$ (default threshold)
- $$|\Delta \pi| \geq 0.1$$ for the dominant isoform set

## Grouping and Annotation

Testing can be organized at different feature levels depending on the biological question:

- **Gene ID**: Tests isoform usage differences within individual genes defined by their LRAA-assigned gene IDs.
- **Gene Symbol**: Tests are grouped by gene symbol, aggregating all isoforms sharing the same symbol. This is useful when analyzing data with gene symbol annotations and enables comparison across studies.
- **Splice hashcode**: Groups isoforms by their splice junction pattern (splice hashcode), enabling analysis of splice variant usage independent of other transcript features.

When splice pattern-collapsed isoform sets are used, transcript IDs are mapped to their representative splice hashcodes via a mapping file (`--splice_hashcode_id_mappings`). This ensures that isoforms differing only in transcript start/end positions but sharing splice structure are analyzed together.

## Output Files

The analysis produces detailed results tables:

### Main DIU Results (`<prefix>.diff_iso.tsv`)
Contains one row per significantly differentially-used gene per cluster pair, including:
- Group identifier (gene ID, gene symbol, or splice hashcode)
- Raw and FDR-adjusted p-values
- Delta-π values for dominant and alternate isoform sets
- Isoform IDs comprising each set
- Isoform expression fractions (π) in each cluster
- Read counts for each isoform set in each cluster
- Cell detection fractions for each isoform set (when available)

### Annotated Isoform Details (`<prefix>.diff_iso.annotated_isoforms.tsv`)
Optional detailed table containing per-isoform information for all tested genes:
- Per-isoform π values and Δπ
- Boolean flags indicating whether each isoform was tested, included in the chi-squared test, and membership in dominant/alternate sets
- Reasons for exclusion (if not tested)
- Cell detection fractions

This comprehensive output enables detailed investigation of isoform switching patterns and validation of detected DIU events.

## Computational Performance

DIU testing supports parallel processing across cluster pairs via the `--CPU` parameter. Because pairwise cluster comparisons are independent, they can be distributed across multiple CPU cores for substantial speedup in datasets with many clusters. The analysis uses Python's multiprocessing framework with process-based parallelism to avoid global interpreter lock limitations.

## Implementation

The DIU testing framework is implemented in:
- `util/sc/diff_iso_usage/sc_pseudobulk_test_isoform_DiffUsage.py`: Main command-line interface and orchestration
- `pylib/DiffIsoformStatTest.py`: Core statistical testing functions and filtering logic

Data preparation utilities for single-cell sparse matrices are provided in `util/sc/` to generate the pseudobulk count and cell fraction matrices required as input.

## Example Command

```bash
sc_pseudobulk_test_isoform_DiffUsage.py \
  --resources_dir LP00714.LRAA.sc^splice_pattern-sparseM.resources \
  --top_isoforms_each 1 \
  --reciprocal_delta_pi \
  --ignore_unspliced \
  --min_cell_fraction 0.05 \
  --splice_hashcode_id_mappings ../data/LP00714.gene_transcript_splicehashcode.withGeneSymbols.tsv \
  --output_prefix LP00714.LRAA.sc^splice_pattern.top_only_w_recip_delta_pi.minCF0.05 \
  --CPU 11 \
  --group_by_feature "gene_symbol"
```

This command:
- Uses splice pattern-collapsed isoforms with gene symbol grouping
- Restricts to the top-expressed isoform from each cluster
- Requires reciprocal Δπ ≥ 0.1 in both directions
- Excludes unspliced isoforms
- Requires ≥5% cell detection fraction for tested isoform sets
- Runs pairwise comparisons across 11 parallel processes

## Manuscript Methods Text for DIU

Differential isoform usage (DIU) was assessed by aggregating isoform-level read counts across all cells within each cluster to form pseudobulk matrices, followed by pairwise comparisons across every cluster pair; unspliced transcripts were excluded from all analyses. For splicing-focused analyses, isoforms sharing the same ordered intron chain were collapsed via splice-pattern hashcodes, summing their counts per cluster before DIU testing so that splice-junction usage changes are measured independently of transcript termini. For termini-focused analyses, isoforms were tested individually, and significant results were filtered to retain pairs that share splice patterns but differ in transcription start and/or end sites, isolating alternative start/stop usage.

Within each analysis mode, isoform (or splice-pattern group) usage fractions for a given gene were computed as

$$
\pi_i = \frac{c_i}{\sum_j c_j},
$$

where $$c_i$$ is the pseudobulk read count for isoform (or group) $$i$$. Usage differences between clusters A and B were summarized as

$$
\Delta \pi_i = \pi_{i,B} - \pi_{i,A}.
$$

Isoforms were partitioned into dominant and alternate sets based on the magnitude and sign of $$\Delta \pi_i$$, enabling detection of reciprocal isoform switches in which one set increases while the other decreases.

This framework adapts the pseudobulk strategy of Jogelkar et al. (2021) with several extensions: (1) reciprocal $$\Delta \pi$$ testing that requires both dominant and alternate sets to exceed effect-size thresholds (default $$|\Delta \pi| \ge 0.1$$); (2) optional filtering on cell detection fraction so that, when a cell-fraction matrix is provided, each enriched isoform set must be detected in at least a specified fraction of cells (default 5%) in its enriched cluster; and (3) additional quality-control filters on minimum gene-level and isoform-set read depths (default ≥25 reads per cluster and per set). Genes passing these filters were evaluated by chi-squared contingency tests on 2×$$N$$ pseudobulk count tables, where rows correspond to clusters and columns to the tested isoforms/sets. P-values were adjusted with the Benjamini–Hochberg procedure, and DIU calls required both the FDR criterion (default FDR < $$0.001$$) and the effect-size thresholds described above.
