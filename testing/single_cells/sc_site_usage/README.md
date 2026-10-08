# Site-first differential TSS / PolyA usage across single-cell clusters

An execution test, and a walk-through, of LRAA's site-usage pipeline
(`util/sc/site_usage/`): does a gene start (TSS) or end (PolyA) its transcripts at
different sites in different cell clusters, and when it does, is that a change of
terminus alone or a switch that comes with alternative splicing?

The inputs are seven genes from the PBMC PacBio Kinnex single-cell data of the LRAA
paper (all ~9,900 clustered cells, 14 Seurat clusters), chosen for clear switches in the
full, genome-wide analysis:

| gene    | site kind | what switches                                    | full analysis (best switch)         |
|---------|-----------|--------------------------------------------------|-------------------------------------|
| SELENOH | TSS       | two TSSs 79 nt apart on one first exon           | naive CD4 T (7) vs non-classical monocyte (11), delta 0.70 |
| EMP3    | TSS       | tandem TSSs, 183 nt apart                        | cluster 0 vs 13, delta 0.68          |
| CIAO2A  | TSS       | tandem TSSs, 105 nt apart                        | cluster 2 vs 1, delta 0.75           |
| AIF1    | TSS       | alternative first exon                           | naive CD4 T (7) vs classical monocyte (1), delta 0.78 |
| POLR2K  | PolyA     | tandem 3' UTR                                     | cluster 3 vs 11, delta 0.32          |
| CMPK1   | PolyA     | tandem 3' UTR                                     | cluster 10 vs 1, delta 0.44          |
| ELOVL5  | PolyA     | intronic PolyA (a different last exon)           | naive B (6) vs classical monocyte (1), delta 0.45 |

## Running it

```
make test          # the pipeline, then check_results.py (about 3 min with CPU=4)
make test CPU=8
make figures       # read-track figures for the showcase events (read_tracks.*.pdf)
make clean
```

Requirements: python3 with pysam, numpy, scipy, pandas; R with DEXSeq, stageR,
data.table, argparse, BiocParallel and Matrix (`make figures` also needs tidyverse and
cowplot). This test is part of `testing/single_cells/Makefile`'s `make test`, to be run
in the single-cell (sc) image, which carries DEXSeq and stageR; the lean core image does
not.

`check_results.py` checks what is robust at this scale: every gene above is tested and
stable; the switches for SELENOH, EMP3, CIAO2A, AIF1 (TSS) and ELOVL5 (PolyA) are found
between the clusters named above, gain the expected site, have the full analysis' delta
(within 0.05), and are classified as above (tandem / different terminal exon); and
read-track data is built for the showcase events. With seven genes rather than
thousands, two things legitimately differ from the genome-wide run and are not checked:
p-values (DEXSeq shares information across genes to estimate dispersion) and the
expression-based switch class ("reciprocal" etc.; it is computed from read ends per
million site-ending reads in each cluster, and these seven genes are the whole library).
So POLR2K and CMPK1 are stable genes here but no single cluster pair reaches
significance, and some switches are "concordant" here that are "reciprocal" genome-wide.

## Inputs (`data/`)

Cut from the full data by `data/prepare_test_data.py` (provenance only; it needs the
full data set). Everything is restricted to the seven genes' spans +/- 1 kb.

| file | what |
|------|------|
| `reads.bam` | aligned reads (minimap2 splice:hq) with `CB` cell barcodes; primary alignments, base qualities and unused tags dropped |
| `cell_clusters.tsv` | cell barcode -> Seurat cluster |
| `cluster_cell_types.tsv` | cluster -> majority cell type (figure labels only) |
| `sites.TSS.bed`, `sites.PolyA.bed` | LRAA's integrated TSS / PolyA site beds (cluster-guided sites, plus initial-run sites farther than 50 nt from every cluster-guided site; `util/integrate_TSS_PolyA_sites.py`) |
| `models.gtf` | LRAA isoform models, transcript ids prefixed with the gene symbol (`SELENOH^t:chr11:+:comp-569:iso-10`) |
| `gene_trans_map.tsv` | LRAA transcript id -> symbol-prefixed gene / transcript ids |
| `cluster_quant.tar.gz` | per-cluster `quant.expr` (for unique full-splice-match read counts per isoform) |
| `tracking.tsv.gz` | per-cluster read-to-isoform assignments (`quant.tracking`), for the read-track figures |

## Walk-through

Each step is a Makefile target; the scripts live in `util/sc/site_usage/` unless noted.

### 1. Site table (`prep_site_table.py`)

Reads the TSS and PolyA site beds and assigns each site to a gene: through the gene
symbols of the isoforms carrying it (cluster-guided sites), else the symbol whose
transcript span covers it on its strand. A site assigned to more than one symbol
(read-through or cis-fusion models) or to none is kept, so that its reads are not
counted toward a neighbouring site, but is left out of the test (`competing = False`).
LRAA's sites are used as they are: LRAA already absorbs read ends within 50 nt
(`max_dist_between_alt_{TSS,polyA}_sites`) into one site, so no further merging is done.

Output: `test.sites.tsv` (one row per site: id, kind, position, gene, carrying
isoforms, whether tested) and `test.gene_spans.tsv`.

### 2. Per-cell read support (`util/sc/site_read_support_to_sparse_matrix.py`)

Counts, per cell, the reads whose 5' end supports each TSS site and whose 3' end
supports each PolyA site, with LRAA's own rules for site support (the counting reuses
LRAA's code, not a reimplementation):

- **Which reads.** Those LRAA's quantification keeps (`Util_funcs.quant_discard_reason`):
  primary alignments, not duplicates or QC failures, no implausibly long intron
  (> 200 kb), percent identity >= 97 (`--HiFi`), mapping quality >= 0 (the production
  setting). Reads without a cell barcode are skipped.
- **Which strand.** The strand of the transcript the read came from: its alignment
  orientation, flipped when minimap2's `ts` tag marks the read antisense to its transcript.
- **Where its ends are.** As LRAA places them (`Pretty_alignment`): a soft-clipped polyA
  tail at the 3' end (>= 7 bases, mostly A) and a few untemplated G's at the 5' end (up
  to 3, or a run of >= 3 next to the alignment; the reverse transcriptase's mark of the
  cap) are stripped first. The TSS end is then the alignment's 5' end, the PolyA end its
  3' end.
- **When an end counts.** Only with no residual soft clip at that end
  (`max_soft_clip_at_TSS` / `max_soft_clip_at_PolyA` = 0: an end with unexplained
  clipped bases is not trusted as a terminus), and only within 25 nt
  (`int(max_dist_between_alt_*_sites / 2)`) of a site of that kind on the read's strand.
  Sites are more than 50 nt apart, so an end can match at most one site; the nearest is
  taken. One count per end: a read can support one TSS and one PolyA site.
- **Counts are reads**, not LRAA's coverage-normalization weights.

Output: `support/test.{TSS,PolyA}-sparseM/` (sites x cells, LRAA's sparse-matrix
layout), per-site tables, sites x clusters pseudobulk counts, and a summary of reads
used / discarded by reason and of ends at a site, rejected for soft clipping, or at no
site.

### 3. Count files for the test (`site_counts_from_lraa_support.py`)

Keeps the clustered cells, orders sites as in the site table (sites without reads are
zero rows), and writes `test.site_counts.mtx.gz`, `test.barcodes.tsv.gz`,
`test.sites.tsv.gz`, `test.summary.tsv`.

### 4. Differential site usage (`site_usage_dexseq.R`)

DEXSeq, with a gene's sites standing in for exon bins and genes for groups, run
separately for TSS and PolyA sites.

**Attribution.** Testing differential site usage between single-cell populations with
DEXSeq on pseudo-bulk replicates is not new here. It is the approach of Sierra (Patrick
et al. 2020, *Genome Biology* 21:167, doi:10.1186/s13059-020-02071-7), which aggregates
each cell population's cells into pseudo-bulk profiles used as DEXSeq replicates to call
differential polyA-peak usage, and of SCAPE (Zhou et al. 2022, *Nucleic Acids Research*
50:e66, doi:10.1093/nar/gkac167), which shuffles cells into six pseudo-replicates for
DEXSeq. We follow that approach. What differs is the input (LRAA's TSS and PolyA sites
and per-cell counts from long-read ends, rather than peaks called from short-read 3'-tag
coverage), the inclusion of TSS sites, and two additions (stageR site confirmation,
stability over several random dealings of cells). DEXSeq: Anders, Reyes & Huber 2012,
*Genome Research* 22:2008, doi:10.1101/gr.133744.111. stageR: Van den Berge et al. 2017,
*Genome Biology* 18:151, doi:10.1186/s13059-017-1277-0.

- **Sites tested:** a site must take >= 10% of its gene's read ends in at least one
  cluster with >= 20 gene read ends there, and have >= 10 read ends overall; a gene
  needs >= 2 such sites (and cross-gene / unassigned sites are excluded).
- **Pseudo-replicates:** the data are one sample, so each cluster's cells are dealt at
  random into 3 pseudo-replicates and read ends summed per site (clusters need >= 30
  cells).
- **Test:** `~ sample + exon + cluster:exon` against `~ sample + exon`, an omnibus test of
  whether site usage differs across all clusters; `perGeneQValue` for the gene FDR;
  stageR (`method = "dtu"`) to confirm which sites carry the change.
- **Seeds:** dealing cells is random, so the test is repeated under 5 seeds; a gene is
  **stable** when its q < 0.05 under at least 4.
- **Pairwise contrasts:** for stable genes, each pair of clusters is tested on that
  pair's pseudo-replicates (first seed), reusing the dispersions fitted across all
  clusters; BH over all pairwise site tests. Output includes each site's usage (share of
  the gene's read ends at its tested sites) in each cluster.
- **Size factors:** DEXSeq's default needs a site with no zero in any pseudo-replicate;
  with few genes there may be none (as here for PolyA), and DESeq2's "poscounts"
  geometric means are used instead (logged). Genome-wide the default is always used.

Caveat: pseudo-replicates measure cell-to-cell sampling noise within one sample, not
variation between donors, so p-values are optimistic in absolute terms; seed stability
guards against calls that depend on one dealing of cells. Effect size, reciprocity and
full-length read support carry as much weight as the FDR.

Output: `test.dexseq.{TSS,PolyA}.genes.tsv` (q per seed, seeds significant, stable),
`.sites.tsv`, `.cluster_usage.tsv.gz`, `.pairwise.tsv.gz`.

### 5. Switch events (`annotate_site_usage_events.py`)

An **event** is a stable gene and a pair of clusters with a site at pairwise padj < 0.05
and |delta usage| >= 0.2, with >= 20 of the gene's read ends (at tested sites) in each
cluster. It is oriented so its largest significant change is a gain: the **gained site**
rises from cluster A to cluster B, the **lost site** is the one that falls most. Each
event is annotated with:

- `event_type`, from the isoform models carrying the two sites: TSS `tandem_TSS` (both on
  one first exon) or `alt_first_exon`; PolyA `tandem_3UTR`, `intronic_PolyA` or
  `alt_last_exon`;
- unique full-splice-match (FSM) reads of the best-supported isoform carrying each site,
  summed over the cluster quantifications;
- `switch_class`, from each site's read ends per million site-ending reads in each
  cluster (+1), at 1.5-fold: reciprocal (gained site up, lost site down), concordant,
  one site changes, neither;
- flags: `monoexonic` (a site carried only by single-exon models), `downstream_TSS_no_FSM`
  (a gained TSS downstream of the lost one with no isoform starting there holding 5 FSM
  reads, which is what 5'-truncated reads look like), `close_sites` (< 30 nt apart),
  `A_rich_downstream` (PolyA, if a genome is given);
- **high confidence** = reciprocal, >= 5 unique FSM reads at both sites, no flags.

Output: `test.dexseq.{TSS,PolyA}.events.tsv`.

### 6. Terminal usage or alternative splicing (`classify_site_pairs_by_splicing.py`)

Classifies each event's site pair from the reads at the two sites (read ends assigned to
the nearest site within 25 nt; pooled over all cells; up to 3,000 reads per site) and
their introns. Of the two sites, the **inner** one is nearer the gene body (proximal
PolyA, downstream TSS); reads from the **outer** site pass the inner site's position on
their way into the gene, so they show whether the two sites share a terminal exon.

- **alternative splicing: terminal exon** -- at least half of the outer site's
  informative reads (reaching past the inner site, or carrying an intron between the
  sites) splice out a stretch between the two sites, or fewer than half of those spanning
  the inner site's terminal intron carry it: an alternative first or last exon, an
  intronic PolyA, a retained intron.
- **alternative splicing: internal** -- the same terminal exon, but some intron further
  in is carried by shares of the two sites' reads differing by >= 0.25.
- **alternative terminal usage** -- the same terminal exon and the same splicing: tandem
  TSSs, tandem 3' UTR PolyA sites.
- **unspliced site** -- fewer than 10 spliced reads at the inner site (a monoexonic
  model, reads inside an intron): read ends alone cannot tell splicing from pre-mRNA.
- **unresolved** -- too few outer reads reach the inner site's terminal intron (long
  3' UTRs whose distal-site reads start inside the UTR).

Intron shares are always taken over the reads *spanning* the intron or position in
question: reads at a distal PolyA site are longer and lose more of their 5' ends, so raw
intron frequencies would fall at the distal site and look like a splicing difference.

Output: `test.dexseq.site_pairs.splicing.tsv`.

### 7. Showcase events and read tracks

`select_showcase_events.py` picks, per site kind and splicing group, the high-confidence
events between clusters of enough cells and reads (here >= 50 cells and >= 20 gene read
ends, scaled down from 200 / 50 genome-wide) that are **dominant switches**, judged on
the read ends: the gained site is the gene's most-used site of its kind in the cluster
gaining it, and the lost site the most-used in the other cluster. The best such event
per gene, ranked by |delta|. Here SELENOH, CIAO2A, EMP3 and AIF1 (TSS) qualify; AIF1's
PolyA switch is real, but another PolyA site leads in one of its clusters.
Dominance is not judged on isoform quantifications: isoforms that share their introns
and differ only at a terminus fit the same reads, so the quantification spreads reads
among them whatever their ends (EMP3's upstream-TSS isoform stays the top isoform in T
cells, where only ~7% of read starts are at its TSS). For each site the isoform drawn is
the one carrying it, with >= 3 unique FSM reads, that has the most reads in the cluster
favouring the site; both are written to `events.tsv` (`gained_tx`, `lost_tx`).

`build_site_event_read_tracks.py` then gathers, for each, what a read-track figure
draws:

- the isoform shown for each site: of those carrying it with >= 5 unique FSM reads, the
  one with the most reads assigned in the cluster that favours the site (cluster B for
  the gained site, cluster A for the lost one), so the pair drawn carries the switch.
  Ranking by unique FSM reads alone picks short fragment models, since full-length
  reads are shared among near-identical full-length models. The manifest records each
  isoform's reads in both clusters and the gained isoform's share of the pair
  (`gained_pair_frac_A/B`), showing whether the pair itself switches;
- 30 reads per cluster from the event's two clusters, split between the two isoforms in
  proportion to their reads there: unique FSM reads first, then, if there are too few,
  "compatible" reads -- reads whose 5' (TSS) or 3' (PolyA) end is within 25 bp of the
  isoform's site and whose alignment fits the model (introns a consecutive run of the
  model's, no block reaching into a model intron); drawn lighter in the figures;
- each cluster's read-end density: the 5' (TSS) or 3' (PolyA) ends of all of the
  cluster's reads in the region, on the gene's strand.

Output: `read_tracks/manifest.tsv` and per event `<tag>.reads.tsv`, `.ends.tsv`,
`.totals.tsv`. `make figures` draws them (`figures.R`, using
`util/sc/notebook_templates/site_usage_funcs.R`): read-end density per cluster on top,
the two isoforms (gained site blue, lost site orange), the sampled reads drawn as aligned
below, and a zoom on the varying terminus (or one per site when they are far apart).
