# Differential TSS / PolyA Site Usage Across Single-Cell Clusters

## Overview

This is a **site-first** analysis of alternative transcription start sites (TSS) and
alternative polyadenylation (PolyA) in single-cell long-read data. For each gene it asks
whether the share of its reads that start (TSS) or end (PolyA) at each of its sites
differs between cell clusters, and, for each switch found, whether it is a change of
terminus alone (the same transcript structure starting or ending elsewhere) or a switch
that comes with alternative splicing (a different first or last exon, an intronic PolyA,
a different internal splice).

It complements the isoform-level differential usage test
([differential_isoform_usage_methods.md](differential_isoform_usage_methods.md)). That test
compares isoform fractions from LRAA's EM quantification; for isoforms that differ only at
a terminus, the EM apportions reads compatible with both, and a significant shift can rest
on that apportioning rather than on reads that reach either terminus. Here the evidence is
the reads' own 5' and 3' ends at LRAA's sites, and isoforms are brought in only afterwards,
to describe and illustrate each switch.

Code: `util/sc/site_usage/`, `util/sc/site_read_support_to_sparse_matrix.py`,
`util/integrate_TSS_PolyA_sites.py`, and the plotting helpers in
`util/sc/notebook_templates/site_usage_funcs.R`. An execution test with a walk-through on
seven PBMC genes is in `testing/single_cells/sc_site_usage/` (part of the single-cell
test suite; needs the sc image for DEXSeq / stageR). The reference application is the
PBMC Kinnex analysis of the LRAA paper (LRAA-paper repo,
`PBMCs_kinnex/LRAA_PBMCs_eval/__Site_Usage_Analysis/`: Makefile + notebook
`sc_site_usage_analysis.Rmd`). Numbers quoted below are from that run.

## Attribution

Testing differential site usage between single-cell populations with **DEXSeq on
pseudo-bulk replicates** is not new here. It is the approach of **Sierra** (Patrick et al.
2020, *Genome Biology* 21:167, doi:10.1186/s13059-020-02071-7), which aggregates each cell
population's cells into pseudo-bulk profiles used as DEXSeq replicates to call
differential polyA-peak usage, and of **SCAPE** (Zhou et al. 2022, *Nucleic Acids
Research* 50:e66, doi:10.1093/nar/gkac167), which shuffles cells into six pseudo-replicates
for DEXSeq. We follow that approach. What differs is the input (LRAA's TSS and PolyA sites
and per-cell counts from long-read ends, instead of peaks called from short-read 3'-tag
coverage), the inclusion of TSS sites, stageR site confirmation, stability over several
random dealings of cells, and everything around the test (site support rules, splicing
classification, read-level illustration). DEXSeq: Anders, Reyes & Huber 2012, *Genome
Research* 22:2008, doi:10.1101/gr.133744.111. stageR: Van den Berge et al. 2017,
*Genome Biology* 18:151, doi:10.1186/s13059-017-1277-0.

## Pipeline at a glance

```
 LRAA site beds (cluster-guided + initial run)
        │  0. integrate_TSS_PolyA_sites.py        one site collection
        ▼
 1. prep_site_table.py                 sites -> genes, testable or not
        ▼
 2. site_read_support_to_sparse_matrix.py   per-cell read support (LRAA's rules)
    site_counts_from_lraa_support.py        -> sites x cells counts
        ▼
 3. site_usage_dexseq.R                 DEXSeq on pseudo-replicates, 5 seeds,
        │                               stable genes, pairwise cluster contrasts
        ▼
 4. annotate_site_usage_events.py       switch events + structure, FSM, switch class
        ▼
 5. classify_site_pairs_by_splicing.py  alternative terminal usage vs alt splicing
        ▼
 6. select_showcase_events.py           showcase events
    build_site_event_read_tracks.py     read-track data
        ▼
    notebook (site_usage_funcs.R)       funnel, class tables, figures (presentation only)
```

All computation runs in the Makefile (`make`); the notebook only reads results and plots
(`make report` knits it). TSS and PolyA are handled in parallel, the same way.

---

## Step 0. One site collection (`util/integrate_TSS_PolyA_sites.py`)

The single-cell workflow reports sites twice: from the initial whole-sample run (the
"basic" catalog) and from the cluster-guided run (per-cluster runs, merged). Every
cluster-guided site is kept; a basic site is added only if it lies **more than 50 nt**
(`max_dist_between_alt_{TSS,polyA}_sites`) from every cluster-guided site of the same kind
and strand. 50 nt is the distance over which LRAA's own site definition absorbs read ends
into one site, so no two sites of one run are closer; a basic site within that distance
is the same site called twice. (Before LRAA f283c326 the cutoff was half of that, 25 nt,
which let 26-50 nt near-duplicates through.)

Output: `<prefix>.integrated.TSS.bed`, `<prefix>.integrated.PolyA.bed` (LRAA's site-bed
columns plus `source` = cluster_guided | basic). PBMC: 24,848 TSS and 57,593 PolyA sites.

## Step 1. Site table (`prep_site_table.py`)

**Inputs:** the two integrated beds; LRAA's final gtf with gene symbols in the transcript
ids (`SELENOH^t:chr11:+:comp-569:iso-10`); the gene / transcript id map
(`gene_transcript_splicehashcode.withGeneSymbols.tsv`).

**What it does, per site:**

- **Gene assignment:** the gene symbols of the isoforms carrying the site (cluster-guided
  sites); for basic sites, whose transcript ids refer to the initial catalog, the gene
  whose transcript span covers the site on its strand.
- **Testable or not (`competing`):** a site assigned to exactly one gene is testable. A
  site assigned to several genes (read-through / cis-fusion models: `cross_gene`) or none
  (`no_gene`) is kept in the table, so its reads are counted toward it and cannot spill
  onto a neighbouring site, but is not tested.
- **No merging** (`--merge_dist_* 0`, the default): LRAA's sites are used as they are.
- **Read window:** 25 nt (used by the step-5 classifier; the counting in step 2 applies
  LRAA's own 25 nt rule).

**Example, SELENOH TSS:** 57,741,377 (290 reads; carried by 5 models, some labelled
SELENOH and some TMX2, the neighbouring gene LRAA placed in the same locus) is
`cross_gene` and not tested; 57,741,491 (3,708) and 57,741,570 (3,565), each carried by
SELENOH models only, are tested.

**Outputs:** `<prefix>.sites.tsv` (site_id, kind, chrom, strand, pos, span, window,
source, support, pas, internal_priming, gene_symbol, gene_key, competing, exclusion,
transcript_ids) and `<prefix>.gene_spans.tsv`. PBMC: 22,569 testable TSS sites (2,279
not), 55,706 testable PolyA sites (1,887 not).

## Step 2. Per-cell read support (`util/sc/site_read_support_to_sparse_matrix.py`, then `site_counts_from_lraa_support.py`)

**Inputs:** the aligned reads (BAM with `CB` cell-barcode tags and minimap2 `ts` tags),
the integrated site beds, the cell clusters; `--HiFi` for HiFi data (the PBMC run).

Counts, per cell, the reads whose 5' end supports each TSS site and whose 3' end supports
each PolyA site, **with LRAA's own rules for site support** (by calling LRAA's code, not
a reimplementation):

1. **Which reads:** those LRAA's quantification keeps (`Util_funcs.quant_discard_reason`):
   primary alignments, not duplicates or QC failures, no intron > 200 kb, percent identity
   >= 97 (`--HiFi`; else LRAA's default), mapping quality >= `min_mapping_quality` (0).
   Reads without a cell barcode are skipped.
2. **Which strand:** the transcript's: the alignment orientation, flipped when minimap2's
   `ts:A:-` marks the read antisense to its transcript.
3. **Where its ends are:** as LRAA places them (`Pretty_alignment`): a soft-clipped polyA
   tail at the 3' end (>= 7 bases, mostly A) and untemplated G's at the 5' end (up to 3,
   or a run of >= 3 next to the alignment: reverse transcriptase's mark of the cap) are
   stripped first. The TSS end is the alignment's 5' end, the PolyA end its 3' end.
4. **When an end counts:** only with **no residual soft clip** at that end
   (`max_soft_clip_at_TSS` / `max_soft_clip_at_PolyA` = 0) and **within 25 nt**
   (`int(max_dist_between_alt_*_sites / 2)`) of a site of that kind on the read's strand
   (sites are > 50 nt apart, so at most one matches; the nearest is taken). One count per
   end: a read can support one TSS site and one PolyA site.
5. **Counts are reads**, not LRAA's coverage-normalization weights (`--weighted` sums those
   instead).

PBMC: 84.7M reads seen, 80.5M used (3.1M supplementary, 1.0M below 97% identity, a few with
long introns). TSS: 26.9M ends at a site, 10.1M rejected for soft clip, 43.4M at no site.
PolyA: 18.0M at a site, 31.9M rejected for soft clip, 30.5M at no site. Ends at no site
are mostly 5'-truncated reads (TSS) and internally primed ends LRAA deliberately has no
site for (PolyA). (Most of the PolyA soft-clip rejections are reads keeping a 1-2 base
non-genomic 3' clip after tail stripping; allowing it was tested on chr19 and not adopted:
more false than true sites gained, and unique-FSM counts shift.)

`site_counts_from_lraa_support.py` then keeps the clustered cells, orders sites as in the
site table (sites without reads are zero rows) and writes the step-3 inputs.

**Outputs:** `lraa_site_support/<prefix>.{TSS,PolyA}-sparseM/` (sites x cells, LRAA's
sparse-matrix layout), per-site and per-cluster tables, a summary of reads and ends; and
for step 3:

```
<prefix>.site_counts.mtx.gz   MatrixMarket, sites (rows) x cells (columns), integer
    %%MatrixMarket matrix coordinate integer general
    82441 9936 15337825
    302 1 18                    <- site 302, cell 1: 18 read ends
<prefix>.sites.tsv.gz         site id per row          (PolyA:chr1:184928:+)
<prefix>.barcodes.tsv.gz      barcode, cluster per column  (AAACAACGACAGTCTA  Cluster_3)
```

**Example, SELENOH:** read starts in naive CD4 T cells (cluster 7) / non-classical
monocytes (cluster 11): 57,741,491: 96 / 270; 57,741,570: 407 / 34. Over the two tested
sites the upstream TSS has 19% of the gene's read starts in cluster 7 and 89% in
cluster 11.

## Step 3. Differential site usage (`site_usage_dexseq.R`)

Run separately per site kind (`--kind TSS|PolyA`).

**1. Sites tested.** Testable sites of the kind, in clusters of >= 30 cells
(`--min_cluster_cells`). A site is kept if, in at least one cluster where the gene has
>= 20 read ends at its testable sites (`--min_gene_cluster_reads`), it takes >= 10% of them
(`--min_site_usage`), and it has >= 10 read ends overall (`--min_site_reads`); a gene needs
>= 2 kept sites. The share in this filter is over all the gene's testable sites; afterwards
usage is recomputed over the kept sites only, so a gene's tested sites' usages sum to 1 in
each cluster. Rationale: sites that are a tiny share everywhere carry little information
and their near-zero counts are where spurious significance comes from; "in at least one
cluster" keeps sites important in only one cell type; the 20-read floor keeps a 10% share
from resting on a handful of reads. PBMC: TSS 22,569 testable sites in 5,030 multi-site
genes -> 10,520 sites in 3,945 genes; PolyA 55,706 in 9,417 -> 23,367 in 6,361.

**2. Usage per cluster:** each kept site's share of its gene's read ends, per cluster
(`cluster_usage.tsv.gz`); the "usage" of all later tables and figures.

**3. Pseudo-replicates:** the data are one sample, so each cluster's cells are dealt at
random into 3 groups (`--pseudoreps`) and each group's read ends summed per site: 14
clusters x 3 = 42 "samples". SELENOH, first seed:

| site | C1.r1 | C1.r2 | C1.r3 | C7.r1 | C7.r2 | C7.r3 | C11.r1 | C11.r2 | C11.r3 |
|---|---|---|---|---|---|---|---|---|---|
| TSS 57,741,491 | 308 | 343 | 330 | 28 | 33 | 35 | 94 | 76 | 100 |
| TSS 57,741,570 | 67 | 57 | 55 | 165 | 127 | 115 | 8 | 10 | 16 |

**4. DEXSeq** (sites in the role of exon bins, genes as groups):

```r
DEXSeqDataSet(pb, sampleData, design = ~ sample + exon + cluster:exon, featureID = site, groupID = gene)
estimateSizeFactors -> estimateDispersions -> testForDEU -> DEXSeqResults(independentFiltering = FALSE)
```

Each site is modelled as "this site vs. the rest of its gene" in every pseudo-replicate
(negative binomial GLM). `sample` absorbs each pseudo-replicate's amount of the gene, so
only relative site use is tested, not expression; `cluster:exon` lets the site's share
depend on the cluster. The likelihood-ratio test of `~ sample + exon + cluster:exon`
against `~ sample + exon` asks whether the site's share differs among any of the clusters
(13 df with 14 clusters). Dispersions are estimated per site and shrunk towards a trend
fitted across all sites (as in DESeq2), which is why many genes are needed.
`perGeneQValue` gives the gene-level FDR ("some site of this gene changes"); **stageR**
(`stageRTx`, `method = "dtu"`) screens genes on that q and confirms which sites carry the
change. Size factors: DEXSeq's default needs a site with no zero in any pseudo-replicate;
with few genes there may be none, and DESeq2's poscounts geometric means are used instead
(logged; never needed genome-wide).

**5. Seeds:** the dealing is random, so all of this runs under 5 seeds (`--n_seeds`); a
gene is **stable** when q < 0.05 under >= 4 (`--min_stable_seeds`). PBMC: 1,315 of 3,945
TSS genes and 619 of 6,361 PolyA genes. SELENOH: 5 of 5.

**6. Pairwise contrasts** (stable genes, first seed): for each pair of clusters (91 pairs),
DEXSeq on that pair's 6 pseudo-replicates, each site taking its dispersion from the
all-cluster fit (6 samples are too few to estimate it); p-values BH-adjusted over all
pairwise site tests; each joined to the site's usage in the two clusters and the change.

**Caveat:** pseudo-replicates from one donor capture cell-to-cell sampling noise, not
variation between donors, so p-values are optimistic in absolute terms; seed stability
guards against calls that depend on one dealing; effect size, reciprocity and full-length
read support (step 4) carry as much weight as the FDR.

**Outputs** (`<prefix>.dexseq.<KIND>.*`): `genes.tsv` (q per seed, seeds significant,
stable), `sites.tsv` (stageR padj per seed, seeds confirmed), `seed<N>.sites.tsv.gz`,
`cluster_usage.tsv.gz`, `pairwise.tsv.gz` (site x cluster pair: p, padj, usage in each,
delta).

## Step 4. Switch events (`annotate_site_usage_events.py`)

**Inputs:** step-3 outputs, the site table, LRAA's gtf, the per-cluster `quant.expr`
files (unique FSM reads per isoform); optionally a genome (A-richness past PolyA sites)
and the isoform-level DTU table (cross-reference).

**Event:** a stable gene and a cluster pair with >= 20 gene read ends (at tested sites) in
each cluster (`--min_gene_reads`; with fewer, shares jump to 0 or 1 and the pairwise test,
on all-cluster dispersions, can still call them) and a site at pairwise padj < 0.05 with
|delta usage| >= 0.2 (`--fdr`, `--min_delta`). Oriented so its largest significant change is
a gain: the **gained site** rises from cluster A to B, the **lost site** falls most.

**Annotations:**

- `event_type` from the exon structures of the isoforms carrying the two sites: TSS
  `tandem_TSS` (both on one first exon) / `alt_first_exon`; PolyA `tandem_3UTR` /
  `intronic_PolyA` / `alt_last_exon`. (Model-based; step 5 decides from the reads.)
- FSM support: per site, the number of isoforms carrying it and the largest unique-FSM
  read count among them (summed over the cluster quantifications).
- `switch_class`: each site's read ends per million site-ending reads in each cluster
  (+1), 1.5-fold: **reciprocal** (gained site up, lost site down), concordant (both up /
  both down), one site changes, neither. A share can flip with only one site moving; a
  reciprocal change in expression is the more interesting switch.
- flags: `monoexonic` (a site carried only by single-exon models), `downstream_TSS_no_FSM`
  (gained TSS downstream of the lost one with no isoform starting there holding 5 unique
  FSM reads: what 5'-truncated reads look like), `close_sites` (< 30 nt),
  `A_rich_downstream` (PolyA, >= 12 A in the 20 genomic bases past the gained site).
- **high_confidence** = reciprocal, >= 5 unique FSM reads at both sites (`--min_FSM`), no
  flags.
- isoform-level DTU on the same gene x cluster pair, if given.

**Example, SELENOH, cluster 7 -> 11:** gained 57,741,491 (19% -> 89%, 96 -> 270 reads),
lost 57,741,570 (81% -> 11%, 407 -> 34), |delta| 0.70, padj ~1e-60, `tandem_TSS`, 79 nt
apart, unique FSM 3,394 / 3,287, gained site log2FC +1.88, lost -3.19: reciprocal, high
confidence. (Not called by the isoform-level DTU for this pair; to be investigated.)

**Output:** `<prefix>.dexseq.<KIND>.events.tsv`, one row per event. PBMC: 581 TSS genes
and 206 PolyA genes with a high-confidence event.

## Step 5. Alternative terminal usage or alternative splicing (`classify_site_pairs_by_splicing.py`)

Classifies each event's site pair (gene + gained + lost site; once, pooled over all cells)
**from the reads at the two sites**: up to 3,000 reads per site whose 5' end (TSS) or 3'
end (PolyA) lies within 25 nt of it, with their aligned spans and introns (CIGAR `N`).

Of the two sites the **inner** one is nearer the gene body (downstream TSS; proximal
PolyA); reads from the **outer** site pass the inner site's position on their way into the
gene, so they show whether the two sites share a terminal exon. Two tests on the outer
site's reads:

**(a) Is there splicing between the two sites?** Among outer reads that reach past `I` (or
carry an intron between the sites), the share with an intron between `O` and `I` (or over
`I`).

**(b) Do they join the gene where the inner site's reads do?** `[t]` is the inner site's
most common intron next to the varying end (its first intron for TSS, last for PolyA).
Among outer reads spanning `[t]`, the share carrying it.

Legend for the diagrams: `===` exon / aligned read, `---` intron / spliced gap, `O` outer
site, `I` inner site; TSS orientation (transcription left to right).

```
(a)             O        I
                |        |
outer read:     =========================--------=======          no intron between O and I
outer read:     ======----------------------=====--------=======  intron covers I: yes

(b)             O        I
                |        |
inner reads:             ===============[--t--]=======
outer read:     ========================[--t--]=======             carries [t]: yes
outer read:     ========================-----------=======        spans [t], other intron: no
```

**Classes:**

- **alt_terminal_usage** (same terminal exon, same splicing; tandem TSSs / tandem 3' UTR):
  (a) < 50% and (b) >= 50%, and no intron further in differs by >= 0.25.
- **alt_splicing:terminal_exon**: (a) >= 50% or (b) < 50%: an alternative first or last
  exon, an intronic PolyA, a retained intron.
- **alt_splicing:internal**: same terminal exon, but an intron further in is carried by
  shares of the two sites' reads differing by >= 0.25.
- **unspliced_site**: < 10 spliced reads at the inner site (monoexonic models, reads inside
  an intron): read ends alone cannot tell splicing from pre-mRNA.
- **unresolved**: < 10 outer reads reach the inner site's terminal intron (long 3' UTRs whose
  distal-site reads start inside the UTR).

Intron shares always use as denominator only the reads **spanning** the intron or
position: distal-PolyA reads are longer and lose more of their 5' ends, so raw intron
frequencies would fall at the distal site and fake a splicing difference.

```
alt_terminal_usage (tandem TSS): SELENOH
                O        I
inner reads:             ===============[--t--]=====------=====
outer reads:    ========================[--t--]=====------=====
  (a) 0.2%   (b) 98%   further in: 0.014

alt_splicing:terminal_exon (alternative first exon): AIF1
                O                          I
inner reads:                               ======[--t--]=====
outer reads:    =======---------------------------===[--t--]=====
  (a) 100%   -> O has its own first exon

alt_splicing:internal
                O        I
inner reads:             ===============[--t--]====-----------------=====
outer reads:    ========================[--t--]====------====-------=====
  (a) no   (b) yes   exon 3 skipped by I-reads, included by O-reads

unspliced_site
inner reads:             ==================          (no introns)

PolyA (mirrored: the varying end is on the right)
tandem 3' UTR: POLR2K, CMPK1
                                         I            O
inner reads:   =====-----[--t--]=============
outer reads:   =====-----[--t--]==========================

intronic PolyA / alternative last exon: ELOVL5
                         I                               O
inner reads:   ====[--t--]=====
outer reads:   ====----------------------------------------=========
  (a) 100%   (b) 0%

unresolved (long 3' UTR)
                         I                                          O
inner reads:   ====[--t--]=====
outer reads:                       ==================================   start inside the UTR
```

PBMC (site pairs): TSS 727 alternative terminal usage, 1,394 terminal exon, 27 internal,
489 unspliced site, 49 unresolved; PolyA 61 / 683 / 13 / 283 / 82. Among high-confidence
genes (best event), about a third of TSS switches and a tenth of PolyA switches are
alternative terminal usage; the rest come with alternative splicing, almost all a
different terminal exon.

**Output:** `<prefix>.dexseq.site_pairs.splicing.tsv` (class and the read-level numbers
behind it: reads, spliced fractions, inner adjacent intron, outer share carrying it, outer
share spliced between sites, splicing divergence and its intron).

## Step 6. Showcase events and read tracks

**6a. Selection (`select_showcase_events.py`)** -- the same rule as the notebook's showcase
tables: high-confidence events; both clusters >= 200 cells (`--min_cluster_cells`) with
>= 50 gene read ends each (`--min_gene_reads`); split by site kind and splicing group
(alternative terminal usage; alternative splicing = terminal exon or internal); the best
event per gene (largest |delta|, ties by gene); ranked by |delta|; the top 6
(`--n_terminal_usage`) and 15 (`--n_alt_splicing`) per kind. PBMC: 42 events. Output:
`events.tsv` (tag, e.g. `SELENOH.TSS.terminal_usage`; gene; kind; sites; clusters).

**6b. Read-track data (`build_site_event_read_tracks.py`)**, per event:

1. the isoform drawn for each site: of those carrying it, the one with the most unique
   FSM reads (SELENOH: iso-10, 3,394; iso-17, 3,287);
2. sampled reads: from the event's two clusters, reads assigned uniquely and as full-splice
   matches to either isoform (per-cluster `quant.tracking`), 30 per cluster sampled from the
   two isoforms together, so each isoform's share of the cluster's reads is kept; their
   aligned blocks from the BAM;
3. read-end density: each cluster's 5' (TSS) or 3' (PolyA) ends of **all** its reads in
   the region, on the gene's strand;
4. totals: unique FSM reads per isoform and cluster.

The tracking file is read once for all events. Outputs: `manifest.tsv` and per event
`<tag>.reads.tsv`, `<tag>.ends.tsv`, `<tag>.totals.tsv`.
`extract_isoform_read_tracks.py` does the same for chosen isoforms and clusters of one gene.

**6c. Figures (`site_usage_funcs.R`)**, two per showcase event:

- the site-usage figure (`plot_site_event`): per-cell site-usage UMAPs (smoothed over the
  cells' SNN neighbourhood), each site's share across clusters, each site's read ends per
  million;
- the read-track figure (`plot_site_event_read_tracks` / `plot_isoform_read_tracks`):

```
 ┌──────────── whole gene ─────────────┐ ┌──── TSS (or PolyA) zoom ────┐
 │ read-start density, cluster A  ▁▂█▁ │ │  ▁ ▁█▁      ▂▅█▂            │
 │ read-start density, cluster B  █▃▁▁ │ │  ▃██▃▁       ▁▂▁            │
 │ isoform at gained site  ■■──■■──■■  │ │ ■■■■■■■■■■■                 │
 │ isoform at lost site     ■──■■──■■  │ │       ■■■■■■■               │
 │ cluster A label (reads: 113 / 396)  │ │                             │
 │   sampled reads, as aligned ■──■──■ │ │  (same reads, zoomed)       │
 │ cluster B label (reads: 349 / 41)   │ │                             │
 │   sampled reads             ■──■──■ │ │                             │
 └─────────────────────────────────────┘ └─────────────────────────────┘
   dashed lines mark the two sites; blue = gained-site isoform, orange = lost
```

  Sites more than 400 nt apart get one zoom each; the whole-gene density bins scale with
  the gene's length, the zooms use 4 bp bins. Existing figure files are kept on re-runs
  unless `options(site_usage.overwrite_figures = TRUE)` (`OVERWRITE_PDFS` in the notebook).

The figures let a switch be checked at three levels, none depending on the EM's
apportioning of shared reads: where the reads start or end (density), that the two
isoforms are real and full-length supported (models, FSM counts), and that individual
reads switch between the cell types (sampled reads).

## Running it

```
# LRAA-paper: PBMCs_kinnex/LRAA_PBMCs_eval/__Site_Usage_Analysis/
make            # steps 0-6 (inputs named at the top of the Makefile)
make report     # also knit sc_site_usage_analysis.Rmd

# LRAA: the execution test on seven genes (~3-5 min)
cd testing/single_cells/sc_site_usage && make test [CPU=4] && make figures
```

Runtime on the full PBMC data (16-core workstation): site support ~14 min (14 cores);
DEXSeq ~25 min per site kind; events and splicing classification a few minutes; showcase
read tracks ~3 min; knitting ~5 min.

## Settings

| step | option | default | meaning |
|---|---|---|---|
| 0 | `--TSS_window` / `--PolyA_window` | 50 (LRAA config) | a basic site within this of a cluster-guided site is dropped |
| 1 | `--merge_dist_TSS` / `--merge_dist_PolyA` | 0 | extra merging of sites (off) |
| 1 | `--window_TSS` / `--window_PolyA` | 25 / 25 (Makefile) | read window for step 5 |
| 2 | `--HiFi` | off | 97% identity floor |
| 2 | `--weighted` | off | sum XW weights instead of reads |
| 3 | `--pseudoreps` | 3 | pseudo-replicates per cluster |
| 3 | `--n_seeds` / `--min_stable_seeds` | 5 / 4 | stability |
| 3 | `--min_cluster_cells` | 30 | clusters tested |
| 3 | `--min_site_usage` / `--min_gene_cluster_reads` / `--min_site_reads` | 0.1 / 20 / 10 | site filter |
| 3 | `--fdr` | 0.05 | gene q, stageR, pairwise padj |
| 4 | `--min_delta` / `--min_gene_reads` / `--min_FSM` | 0.2 / 20 / 5 | event, high confidence |
| 5 | `--min_spliced` / `--min_spanning` / `--min_adjacent_share` / `--min_divergence` / `--max_reads_per_site` | 10 / 10 / 0.5 / 0.25 / 3000 | splicing classes |
| 6 | `--min_cluster_cells` / `--min_gene_reads` / `--n_terminal_usage` / `--n_alt_splicing` / `--max_reads` | 200 / 50 / 6 / 15 / 30 | showcase selection, reads drawn |

LRAA's site parameters (`max_dist_between_alt_*_sites` 50, `min_alignments_define_*_site`
5, `max_soft_clip_at_*` 0) were evaluated on chr22 / chr19 for human (narrower PolyA
window, higher TSS minimum, a 3' soft-clip allowance) and deliberately left at their
defaults for now.

## Caveats and known limitations

- **One donor:** pseudo-replicates measure sampling noise within one sample; p-values are
  optimistic; read as "reproducible across random dealings", not population-level.
- **Sites LRAA does not call are invisible:** about half of read ends fall at no site. On
  the PolyA side these are mostly internally primed ends, which LRAA vetoes by design
  (vetoed peaks are 6-16% atlas-supported and ~5% PAS-bearing); a few real sites next to
  genomic A-runs (e.g. TXNIP's main site) are lost with them. Not rescued on a PAS alone,
  by choice (precision over sensitivity).
- **PolyA support excludes ~40% of reads** (1-2 base non-genomic 3' clips; see step 2).
- **Component-wide site pruning in LRAA:** `min_{TSS,PolyA}_iso_fraction` is computed over a
  splice-graph component, which can span neighbouring genes, so a gene next to a highly
  expressed neighbour can lose its main site (e.g. DDAH2's TSS next to CLIC1). Narrow (about
  170 TSS / 110 PolyA genes in PBMC, mostly immune-receptor segments); set aside.
- **5' truncation** can mimic downstream TSSs; the `downstream_TSS_no_FSM` flag and the FSM
  requirement guard against it.
- **Small inputs** (the seven-gene test): DEXSeq's dispersion trend and the per-cluster
  library totals behind `switch_class` differ from a genome-wide run; the test checks only
  what is robust at that scale (see its README).

## References

- Patrick R, et al. (2020) Sierra: discovery of differential transcript usage from
  polyA-captured single-cell RNA-seq data. *Genome Biology* 21:167.
  doi:10.1186/s13059-020-02071-7
- Zhou R, et al. (2022) SCAPE: a mixture model revealing single-cell polyadenylation
  diversity and cellular dynamics during cell differentiation and reprogramming. *Nucleic
  Acids Research* 50:e66. doi:10.1093/nar/gkac167
- Anders S, Reyes A, Huber W (2012) Detecting differential usage of exons from RNA-seq data.
  *Genome Research* 22:2008-2017. doi:10.1101/gr.133744.111
- Van den Berge K, Soneson C, Robinson MD, Clement L (2017) stageR: a general stage-wise
  method for controlling the gene-level false discovery rate in differential expression and
  differential transcript usage. *Genome Biology* 18:151. doi:10.1186/s13059-017-1277-0
- Squair JW, et al. (2021) Confronting false discoveries in single-cell differential
  expression. *Nature Communications* 12:5692. doi:10.1038/s41467-021-25960-2
