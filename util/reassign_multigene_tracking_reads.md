# reassign_multigene_tracking_reads.py Walkthrough

## Purpose

`reassign_multigene_tracking_reads.py` is a post-processing step for merged LRAA
quantification outputs. It corrects reads that were assigned to transcripts from
multiple genes, usually because secondary alignments created cross-gene
candidates.

It takes:

```text
--quant_expr       pre-correction quant.expr
--tracking         pre-correction quant.tracking or .gz
--output_expr      corrected quant.expr
--output_tracking  corrected quant.tracking or .gz
```

The main invariant is that only connected cross-gene read components get
re-EM'd. Everything else keeps its existing tracking fraction.

## Input Helpers

`open_maybe_gzip()` handles plain or gzipped tracking files.

`iter_non_comment_lines()` skips LRAA version comments beginning with `#`.

`read_expr_rows()` reads `quant.expr` with `csv.DictReader`.

`read_tracking_fieldnames()` grabs the first non-comment tracking header.

## Main Setup

`main()` reads expression rows and tracking field names, then builds:

```python
expr_by_tx[(gene_id, transcript_id)] = expr_row
gene_to_transcripts[gene_id] = [(gene_id, transcript_id), ...]
```

It then creates a temporary working directory with intermediate files for sorted
tracking, cross-read gene pairs, component candidates, and EM fractions.

## Step 1: Sort Tracking By Read

The tracking file is externally sorted by `read_name`. This lets the rest of the
script stream groups of rows per read instead of loading the whole tracking file
into memory.

`sort_tracking_body()` shells out to GNU `sort`.

`iter_tracking_groups()` then yields:

```python
read_name, [tracking_row, tracking_row, ...]
```

from that sorted file.

## Step 2: Deduplicate Tracking Rows

`dedupe_tracking_group()` collapses duplicate rows with the same:

```python
(read_name, gene_id, transcript_id)
```

If duplicates have `read_weight`, it keeps the maximum weight. This matters
because merged contig or shard outputs can create repeated candidates for the
same read/transcript.

## Step 3: Find Cross-Gene Reads

`scan_read_sorted_tracking()` is the first major pass over sorted tracking.

For each read, it gets the set of genes represented in that read's candidate
rows:

```python
genes = sorted({row["gene_id"] for row in deduped_rows})
```

If `len(genes) > 1`, the read is affected. The script writes pairs like:

```text
readA    gene1
readA    gene2
```

to `cross_read_genes.tsv`.

It also uses `DisjointSet` to union genes touched by the same cross-gene read.
That creates connected gene components.

Example:

```text
read1 -> geneA, geneB
read2 -> geneB, geneC
```

means `geneA`, `geneB`, and `geneC` are one component and must be EM'd together.

If no affected reads are found, the script copies `quant.expr` unchanged and
writes deduplicated tracking.

## Step 4: Assign Reads To Components

`write_read_components()` converts the cross-read gene file into:

```text
read_name    component_id
```

and returns:

```python
component_genes[component_id] = {gene_id, ...}
```

Then `main()` inverts that into:

```python
gene_to_component[gene_id] = component_id
```

## Step 5: Split Fixed Reads From EM Reads

`write_component_candidate_rows()` does another pass over sorted tracking.

For reads not touching any affected component, it accumulates fixed counts
directly from existing `frac_assigned`:

```python
fixed_counts[tx_key] += frac
```

For reads in an affected component, it writes candidate rows:

```text
component_id    read_name    gene_id    transcript_id    read_weight
```

A key guard checks whether a read somehow spans multiple disconnected affected
components. That should not normally happen if the disjoint-set logic captured
the graph correctly.

The candidate rows are then sorted by component and read.

## Step 6: Run EM Per Component

`expr_init_counts` is built from original `all_reads`.

`run_component_ems_to_fraction_file()` streams sorted candidate components. For
each component, it builds:

```python
component_reads = [read_name, ...]
component_transcripts = all transcripts from all genes in the component
read_candidates[read_name][(gene_id, transcript_id)] = read_weight
```

Then it calls `run_component_em()`.

Inside `run_component_em()`:

1. Initializes transcript counts from original `quant.expr` `all_reads`.
2. Adds alpha regularization for transcripts participating in ambiguous reads.
3. Runs the E-step: assigns each read fractionally by
   `read_weight * current_abundance`.
4. Runs the M-step: updates transcript counts from expected read support plus
   alpha plus `min_expr`.
5. Stops when L1 count movement is below `tol`.
6. Recomputes final reported counts from read fractions only, excluding alpha
   and `min_expr`.

The output from component EM is written to `read_fractions.unsorted.tsv` as:

```text
read_name    gene_id    transcript_id    frac_assigned    read_weight
```

Then that fractions file is sorted by read, gene, and transcript.

## Step 7: Recompute quant.expr

The script recomputes expression fields from:

```python
updated_counts = fixed_counts + EM reassigned counts
updated_uniq_reads = fixed unique reads + EM unique-ish reads
```

It calculates:

```python
gene_total = sum transcript counts for gene
isoform_fraction = all_reads / gene_total
unique_gene_read_fraction = uniq_reads / gene_total
TPM = all_reads / total_reported_read_count * 1e6
RPM_total_reads = all_reads * rpm_scale
```

`rpm_scale` is derived from the median original
`RPM_total_reads / all_reads`. That preserves whatever total-read scaling LRAA
had already applied.

Finally it writes corrected `quant.expr`.

## Step 8: Rewrite Tracking

`write_final_tracking_from_read_sorted()` streams the read-sorted original
tracking and the read-sorted corrected fraction file together.

For each deduplicated original tracking row:

1. If `(gene_id, transcript_id)` has a corrected fraction for that read, it
   replaces `frac_assigned`.
2. If `read_weight` exists, it also writes the formatted weight from the EM
   candidate file.
3. Otherwise, it writes the original row unchanged.

Then `main()` writes the final corrected tracking file.

## Mental Model

The full script can be summarized as:

```text
1. Read quant.expr metadata.
2. Sort tracking by read.
3. Detect reads that touch multiple genes.
4. Build connected gene components from those reads.
5. Keep all unaffected read fractions fixed.
6. For affected components only:
   run per-read weighted EM.
7. Recompute expr totals from fixed + corrected read fractions.
8. Rewrite tracking with corrected fractions only where needed.
```

The reason this exists separately from normal LRAA EM is that normal
quantification happens locally inside LRAA's contig/strand processing. This
script operates after merge time, when secondary alignments can reveal
cross-gene multimapping relationships that need to be reconciled globally.
