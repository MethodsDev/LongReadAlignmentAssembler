# BAM Correction Tags Documentation

## Overview

The `correct_bam_alignments.py` utility produces four output BAM files with different purposes:

1. **`{prefix}.corrected.bam`** - Full dataset with corrections applied
2. **`{prefix}.faulty_only.bam`** - Only original alignments that were corrected (with XC tags using `~` delimiters)
3. **`{prefix}.corrected_only.bam`** - Only the corrected versions (with original metadata tags)
4. **`{prefix}.annotated_originals.bam`** - Full dataset in original form with correction metadata tags

## Tag Schemes

### Corrected Alignments (in `corrected.bam` and `corrected_only.bam`)

These BAMs contain the **corrected** alignment positions and CIGARs, with tags preserving the original values:

- **XC:Z** - Correction type: `"left"`, `"right"`, or `"both"` 
  - Indicates which end of the alignment was corrected
- **OC:Z** - Original CIGAR string
  - The CIGAR before correction was applied
- **OA:i** - Original alignment position (`reference_start`)
  - The genomic coordinate before correction

**Example:**
```
Position: 2047
CIGAR: 28S49M1150N117M295N22M1D79M...195N6M
XC:Z:right
OC:Z:28S49M1150N117M295N22M1D79M...6S
OA:i:2048
```
This read was corrected on the right end. The current alignment (position 2047, CIGAR ending in `195N6M`) is the corrected version. The original was at position 2048 with CIGAR ending in `6S` (soft clip).

### Faulty Alignments (in `faulty_only.bam`)

This BAM contains only the **original** alignments that were successfully corrected. These are the "before" versions. Each read has an XC tag with `~` delimiters to differentiate from corrected versions:

- **XC:Z** - Correction type with `~` delimiters: `"~left~"`, `"~right~"`, or `"~both~"`
  - The `~` characters distinguish faulty originals from corrected versions
  - Indicates which end was correctable

**Example:**
```
Position: 2048
CIGAR: 28S49M1150N117M295N22M1D79M...6S
XC:Z:~right~
```
This is the original faulty alignment (position 2048, soft clip ending). The `~right~` tag indicates the right end was correctable and was successfully corrected in the other output BAMs.

### Annotated Originals (in `annotated_originals.bam`)

This BAM contains **all reads** from the input BAM in their original alignment form. Reads that were successfully corrected have XC tags with `~` delimiters (matching `faulty_only.bam`):

- **XC:Z** - Correction type with `~` delimiters: `"~left~"`, `"~right~"`, or `"~both~"`
  - The `~` characters mark original faulty alignments
  - Only present on reads that were successfully corrected

**Example:**
```
Position: 2048
CIGAR: 28S49M1150N117M295N22M1D79M...6S
XC:Z:~right~
```
This is the **original** alignment (position 2048, CIGAR ending in `6S`). The `~right~` tag indicates it was corrected on the right end in the corrected output BAMs.

Reads that were not corrected appear in this BAM without any XC tags.

## Use Cases

### `corrected.bam` and `corrected_only.bam`
Use these when you want to work with the corrected alignments for downstream analysis (quantification, visualization, etc.). The OC/OA tags let you trace back to the original alignment if needed.

### `faulty_only.bam`
Use this to examine just the original faulty alignments that were successfully corrected. The `~delimiters~` in the XC tag distinguish these from corrected versions. Good for:
- Analyzing what types of alignment errors were present
- Comparing before/after correction
- Understanding the distribution of left vs. right vs. both-end corrections

### `annotated_originals.bam`
Use this when you need to:
- Identify which reads were corrected within the context of the full dataset
- Filter corrected vs. uncorrected reads (by presence of XC tag)
- Preserve original alignments for compatibility while marking which ones were corrected
- The XC tags (with `~` delimiters) match the `faulty_only.bam` format, making it easy to identify original faulty alignments

## Validation

All output BAMs:
- Pass `samtools quickcheck`
- Are coordinate-sorted and indexed
- Preserve all indels (insertions and deletions) correctly
- Have validated CIGAR strings (query length matches sequence length)

## Technical Details

Corrections are made by:
1. Building a splice graph from the GTF annotation
2. Identifying soft-clipped read ends that can be extended across known introns
3. Only correcting **one end** (left OR right) per alignment, never both simultaneously
4. Preserving the original CIGAR operations for the uncorrected portion of the alignment

Tags use standard SAM format types:
- `Z` (string) for XC (with or without `~` delimiters), OC
- `i` (signed integer) for OA

## Tag Summary Table

| BAM File | Contains | XC Tag Format | Other Tags |
|----------|----------|---------------|------------|
| `corrected.bam` | All reads with corrections applied | `"left"`, `"right"`, `"both"` | OC, OA (on corrected reads) |
| `corrected_only.bam` | Only corrected reads | `"left"`, `"right"`, `"both"` | OC, OA |
| `faulty_only.bam` | Only original faulty reads | `"~left~"`, `"~right~"`, `"~both~"` | None |
| `annotated_originals.bam` | All reads in original form | `"~left~"`, `"~right~"`, `"~both~"` (on corrected reads only) | None |
