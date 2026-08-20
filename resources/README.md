# Bundled resources

## `human_rDNA_cassette.fa`

The human ribosomal DNA complete repeating unit, GenBank
[U13369.1](https://www.ncbi.nlm.nih.gov/nuccore/U13369.1), 42,999 bp. Contains the
45S precursor rRNA transcription unit (5' external transcribed spacer, 18S, ITS1,
5.8S, ITS2, 28S, 3' external transcribed spacer) plus the intergenic spacer.

Used as the default `--rdna_mask_fasta` reference: at startup LRAA aligns this
sequence against `--genome` (whole genome, or one chunk's mini-genome under
`--chunk`) with `minimap2`, and reads whose alignment overlaps a hit are excluded
from splice-graph construction, coverage normalization, and quantification (see
`--no_rdna_mask` to disable, `--rdna_mask_fasta` to substitute a different
organism's cassette).

Rationale: the reference genome's rDNA loci (GRCh38 acrocentric short arms, plus
several unplaced/decoy scaffolds) carry a handful of collapsed, partial repeat-unit
copies. A locus's true in-vivo copy number is far higher and mostly absent from the
assembly, so the reads that do align there are extreme, low-mapping-quality
multi-mapping piles with almost no isoform content of their own -- Ensembl/GENCODE
annotates only the small mature-rRNA fragments (e.g. `5.8S`) as genes at these
positions, leaving the much larger 18S/28S/spacer footprint unannotated and
therefore unmaskable by GTF biotype alone. Sequence homology to the repeat unit,
not genome annotation, is what actually identifies these reads.
