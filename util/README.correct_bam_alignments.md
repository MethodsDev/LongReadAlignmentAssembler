# BAM Alignment Correction Utility

## Overview

`correct_bam_alignments.py` is a standalone utility that leverages LRAA's splice graph-guided alignment correction to fix soft-clipped alignments in BAM files. It uses reference transcript annotations and genome sequence to identify and correct alignment errors where soft-clipped sequences can be realigned across known splice junctions.

## How It Works

1. **Builds Splice Graph**: Constructs a splice graph from normalized read alignments and reference GTF annotations
2. **Identifies Correction Candidates**: Finds alignments with soft-clipping (5-20bp by default) that overlap known introns
3. **Corrects Alignments**: Extends alignments across splice junctions where the soft-clipped sequence exactly matches the genome
4. **Outputs Corrected BAMs**: Generates three BAM files with corrections applied and tracking tags

## Correction Method

The tool corrects alignments using the same method as LRAA's main pipeline:

- **Left-side correction**: When left soft-clipping overlaps an intron's right boundary, attempts to realign the soft-clipped sequence to the left of the intron
- **Right-side correction**: When right soft-clipping overlaps an intron's left boundary, attempts to realign the soft-clipped sequence to the right of the intron
- **Exact matching**: Only corrects when the soft-clipped sequence exactly matches the genomic sequence (case-insensitive)

## Usage

### Basic Usage

```bash
python util/correct_bam_alignments.py \
  --bam input.bam \
  --gtf reference.gtf \
  --genome genome.fa \
  --output_prefix corrected
```

### Skip Normalization (for small datasets)

```bash
python util/correct_bam_alignments.py \
  --bam input.bam \
  --gtf reference.gtf \
  --genome genome.fa \
  --output_prefix corrected \
  --skip_normalization
```

### Process Specific Contigs with Parallelization

```bash
python util/correct_bam_alignments.py \
  --bam input.bam \
  --gtf reference.gtf \
  --genome genome.fa \
  --output_prefix corrected \
  --contig chr1,chr2,chr3 \
  --CPU 8
```

## Arguments

### Required

- `--bam`: Input BAM file (must be coordinate-sorted; will be indexed if needed)
- `--gtf`: Reference GTF file with transcript annotations
- `--genome`: Reference genome FASTA file (will be indexed if needed)
- `--output_prefix`: Prefix for all output files

### Optional

- `--normalize_max_cov_level`: Coverage normalization level for splice graph construction (default: 1000)
- `--skip_normalization`: Skip BAM normalization (useful for small/targeted datasets)
- `--CPU`: Number of parallel processes (default: 4)
- `--contig`: Comma-separated list of contigs to process (default: all contigs in BAM)
- `--config_update`: JSON file with LRAA config updates (advanced users)
- `--debug`: Enable verbose debug logging
- `--keep_temp`: Keep temporary per-contig files for inspection

## Outputs

### BAM Files

1. **`{prefix}.corrected.bam`**: Full BAM file with corrections applied
   - Contains all alignments from input BAM
   - Corrected alignments have updated CIGAR strings and positions
   - Tagged with correction metadata (see below)

2. **`{prefix}.faulty_only.bam`**: Original faulty alignments
   - Contains only the original (uncorrected) versions of alignments that were corrected
   - Useful for comparing before/after

3. **`{prefix}.corrected_only.bam`**: Corrected alignments only
   - Contains only the corrected versions of alignments
   - Useful for analyzing corrections

All BAM files are coordinate-sorted and indexed (`.bai` files created automatically).

### Correction Tags

Corrected alignments in the full BAM are tagged with:

- `XC:Z:{type}`: Correction type
  - `left`: Left soft-clip was corrected
  - `right`: Right soft-clip was corrected
  - `both`: Both sides were corrected
- `OC:Z:{cigar}`: Original CIGAR string before correction
- `OA:i:{pos}`: Original alignment start position (0-based)

### Statistics Report

**`{prefix}.correction_stats.tsv`**: Tab-separated statistics per contig/strand

Columns:
- `contig`: Chromosome/contig name
- `strand`: `+` or `-`
- `total_reads`: Total primary alignments processed
- `corrected_reads`: Number of alignments corrected
- `pct_corrected`: Percentage corrected
- `left_corrections`: Number with left-side correction
- `right_corrections`: Number with right-side correction
- `both_corrections`: Number with both sides corrected

The file includes a summary row with totals across all contigs.

## Alignment Filtering

The tool processes **only primary alignments**:
- Skips unmapped reads
- Skips secondary alignments (SAM flag 0x100)
- Skips supplementary alignments (SAM flag 0x800)

## Performance Considerations

### Normalization

By default, the tool normalizes the BAM to max coverage of 1000 reads per position for splice graph construction. This:
- Improves splice graph quality by reducing noise
- Speeds up processing for high-coverage regions
- Can be disabled with `--skip_normalization` for small datasets

### Parallelization

- Uses multiprocessing to process contigs in parallel
- Set `--CPU` to number of available cores
- Each contig/strand is processed independently
- Results are merged at the end

### Memory Usage

- Processes one contig/strand at a time (low memory footprint per worker)
- Temporary files written to disk during processing
- Final merge uses `samtools merge` for efficiency

## Configuration

The tool uses LRAA's default configuration values for:
- Soft-clip correction thresholds (`min_softclip_realign_test=5`, `max_softclip_realign_test=20`)
- Splice site consensus sequences
- Intron length limits
- And other splice graph construction parameters

Advanced users can override these with `--config_update` pointing to a JSON file:

```json
{
  "min_softclip_realign_test": 3,
  "max_softclip_realign_test": 25,
  "min_alt_splice_freq": 0.01
}
```

## Examples

### Example 1: Basic correction with default settings

```bash
python util/correct_bam_alignments.py \
  --bam sample.bam \
  --gtf gencode.v38.annotation.gtf \
  --genome hg38.fa \
  --output_prefix sample_corrected
```

Output:
- `sample_corrected.corrected.bam`
- `sample_corrected.faulty_only.bam`
- `sample_corrected.corrected_only.bam`
- `sample_corrected.correction_stats.tsv`

### Example 2: Targeted correction (chr1 only, no normalization)

```bash
python util/correct_bam_alignments.py \
  --bam chr1_reads.bam \
  --gtf gencode.chr1.gtf \
  --genome chr1.fa \
  --output_prefix chr1_corrected \
  --contig chr1 \
  --skip_normalization
```

### Example 3: High-throughput with parallelization

```bash
python util/correct_bam_alignments.py \
  --bam large_sample.bam \
  --gtf annotation.gtf \
  --genome genome.fa \
  --output_prefix corrected \
  --CPU 16 \
  --normalize_max_cov_level 500
```

## Troubleshooting

### "No splice graph for contig"
- Check that your GTF contains annotations for that contig
- Verify contig names match between BAM, GTF, and genome FASTA
- Try running with `--debug` for more details

### High memory usage
- Reduce `--CPU` to use fewer parallel processes
- Enable normalization if not already active
- Process specific contigs one at a time with `--contig`

### Slow performance
- Increase `--CPU` for more parallelization
- Consider reducing `--normalize_max_cov_level` to speed up splice graph construction
- Use `--contig` to process only specific regions of interest

## Technical Details

### Dependencies

The tool requires the same dependencies as LRAA:
- Python 3.7+
- pysam
- networkx
- intervaltree

### Workflow

1. Parse command-line arguments and initialize configuration
2. Index BAM and genome FASTA if needed
3. Normalize BAM by strand (optional)
4. For each contig/strand in parallel:
   - Parse GTF transcripts
   - Build splice graph from normalized BAM + GTF
   - Process alignments from original BAM
   - Identify correction candidates
   - Apply corrections using splice graph
   - Write temporary BAM files
5. Merge temporary BAMs into final outputs
6. Generate statistics report
7. Clean up temporary files (unless `--keep_temp`)

### Code Integration

The utility integrates with LRAA's core modules:
- `Splice_graph`: Builds splice graph from alignments and GTF
- `Pretty_alignment`: Handles alignment correction
- `Pretty_alignment.to_corrected_pysam_alignment()`: Converts corrected segments back to SAM format
- `LRAA_Globals`: Configuration management

## Citation

If you use this utility, please cite LRAA:

*(Citation details to be added)*
