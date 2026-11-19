# Gene counts that include intronic reads (single-nucleus)

Status: design plus a working prototype on branch `sn-intronic-gene-counts` (not on devel).
Second revision, 2026-11-18, which applies the user's answers to the first round's open
questions (section 6). It is opt-in via `--config_update '{"gene_body_counts": true}'`.
With the key off, LRAA's outputs are unchanged. With it on, `quant.expr` and
`quant.tracking` stay identical to a key-off run (tested, chrM included since the
`--oversimplify` tie-break became deterministic; section 5).

Third revision, 2026-11-19: the key now works in discovery mode too (de novo and
ref-guided), on the final quant of the assembled models, and the single-cell workflow
has an opt-in single-nucleus mode, `sc_nucseq`, that clusters on gene_full (section 4,
"Single-nucleus mode in the WDLs").

## The problem

In quant-only mode a read is assigned to an isoform only if it is structurally compatible
with it: at least `fraction_read_align_overlap` (0.75) of the read must lie on the isoform's
exons. A pre-mRNA read lying in an intron matches no isoform. It gets no tracking row, so it
is missing from the isoform matrix and also from the gene matrix that the `util/sc` builders
derive from the tracking file. Single-nucleus libraries are mostly pre-mRNA. In the ALS
LP00714 nucseq test window below, 87% of reads are unspliced, and the standard gene matrix
holds 4,307 reads of the 11,550 that lie in annotated genes on their own strand.

Cell Ranger ("include introns"), STARsolo `GeneFull` and alevin-fry's USA mode all count
these reads toward their gene. This design adds the same thing to LRAA.

## What "chr-simplify" turned out to be

The repo has no `chr_simplify` or `simplify` method. The nearest match is **`--oversimplify`**
(`LRAA_Globals.config["oversimplify_enabled"/"oversimplify_contigs"]`, default `chrM,M` in
the single-cell workflow). Its quant-only implementation, `LRAA:_run_oversimplify_best_overlap`:

- works per contig;
- assigns each read to the single transcript with the most exonic bases of overlap, breaking
  ties with `_oversimplify_tie_break` (sorted tied ids indexed by an md5 of the read name, so
  the pick does not depend on run, chunking or read order; it was an unseeded `random.choice`
  until 2026-11-18);
- does no EM;
- skips any read with zero exonic overlap.

The new mode generalizes that idea in three ways. It works on genes (exon union plus gene
span) rather than transcripts. It scores intronic overlap as well as exonic. It runs
genome-wide and *adds to* the normal isoform quantification instead of replacing it. Its
assignment logic does not change; the path now also writes gene-body rows (section 4).

## Design in brief

During the streaming quant pass, every retained read already gets an isoform verdict. When
`gene_body_counts` is on, each read is also written to a **sidecar**
`<prefix>.quant.gene_body.tracking.gz`:

- **Isoform-assigned reads** keep their isoform assignment. The EM fractions are summed per
  gene, so these rows reproduce what the standard gene matrix counts (verified below).
- **Reads the isoform pass leaves unassigned** (no graph path, spacer path, or a path that
  matches no isoform) go to the most appropriate same-strand gene by the rules below.
- **Every row** carries a splice class: S (spliced), U (unspliced) or A (ambiguous).
- **Matrices.** A new builder, `util/sc/gene_body_tracking_to_sparse_matrix.py`, makes
  `gene_full` = S + U + A plus one matrix per class. These sit beside the existing gene
  matrix; they do not replace it.

## 1. Read classes and strand

**Strand.** Strand comes from LRAA's existing read-strand logic: `Util_funcs.transcribed_strand`
uses the `ts` tag, corroborated by splice motifs, and falls back to the aligned strand. Each
read is processed only in the pass for its own contig-strand, so a gene on the other strand
is never a candidate. An unspliced read has no `ts` evidence and takes its aligned strand.
**Inputs are expected to be strand-specific** (oriented 10x / Kinnex reads). In an
unstranded library, half the unspliced reads would land on the wrong strand and count as
antisense or intergenic. There is no special handling for that, by decision
(user, 2026-11-18). In the ALS window every spliced read carries `ts:A:+`, so that library
is oriented.

A log-only check flags it. A spliced read's transcribed strand is fixed by its junctions,
so in an oriented library its aligned strand matches the unit's strand almost always,
while unstranded cDNA gives about half the other way. `GeneBodyRecorder` counts spliced
reads whose aligned strand disagrees with the unit's strand and, at `finish()`, logs a
WARNING for the contig-strand when they are more than 20% of at least 200 spliced reads
(`STRAND_CHECK_*` in `GeneBodyAssigner`). A reverse-oriented library trips it too, and
should: its unspliced reads all land on the wrong strand. Nothing written changes.

| class | definition | handled by |
|---|---|---|
| exonic-compatible | isoform-compatible (unchanged: graph path plus 0.75 overlap rule, EM) | existing quant; sidecar `source=isoform` |
| intronic / pre-mRNA / partial | retained read, no isoform verdict, inside a same-strand gene span | `GeneBodyIndex`; sidecar `source=gene_body:<why>` |
| intergenic / antisense | no same-strand candidate gene | counted (`reads_no_gene`), not written |
| filtered | rejected by `quant_discard_reason` (MAPQ, per-id, rDNA, ...) | excluded, as for isoforms |
| cross-component ambiguous | path spans isoforms of two EM components (existing rule) | excluded, as rescue excludes it |

`<why>` is one of:

- `no_path`, `spacer_path`, `unassignable`: the streaming pass's unassigned exits;
- `no_quant_unit`: the unit made no isoform quant at all (section 4);
- `oversimplify_no_exon_overlap`: an `--oversimplify` contig read with no exonic overlap.

Transcriptome rescue, when on, gets first claim on these reads. A rescued read counts as
`isoform`.

**Splice class** is computed against the gene the read is credited to. The gene's intronic
sequence is its span minus the union of its exons.

- **U**: the read covers at least `gene_body_min_feature_bp` (20) bases of that intronic
  sequence. This holds whether the read is spliced or not, so partially spliced molecules
  with a retained intron count as U.
- **S**: not U, and the read has at least one junction.
- **A**: not U and unspliced, so it lies entirely on exons. It fits mature mRNA and pre-mRNA
  equally well. A is kept as its own category (user, 2026-11-18), and gene_full = S + U + A.

An isoform-assigned read can still be U: the 0.75 rule lets up to 25% of a read run into an
intron.

## 2. Choosing "the most appropriate gene"

Rules for a read with no isoform assignment, applied in order
(`pylib/GeneBodyAssigner.py`). Only genes on the read's own strand are considered.

0. **Candidate genes.**
   - A gene is a candidate if its span (first exon start to last exon end, over all
     isoforms) holds at least `gene_body_min_read_frac_in_span` (0.5) of the read's aligned
     bases.
   - There is no extra allowance for unannotated 3' extensions (user, 2026-11-18). The 0.5
     fraction already tolerates an overhang of up to half the read.
   - **Readthrough genes are excluded** (`gene_body_exclude_readthrough_genes`, default on).
     A gene counts as readthrough when every one of its transcripts carries the GTF tag
     `readthrough_transcript` (GENCODE), or its gene line carries `readthrough_gene`. It is
     never a gene-body candidate. Its isoform-assigned reads still count, and still get a
     class.
   - This needed a parser change: the GTF parser (`Transcript.GTF_contig_to_transcripts`) now
     keeps every `tag` value in `Transcript.get_gtf_tags()`. GENCODE repeats the `tag` key on
     a line, and the meta dict keeps only the last value. The meta dict itself is unchanged,
     so GTF re-export is unaffected (tested).
1. **Junctions.** Prefer the gene with the most read junctions that exactly match one of its
   annotated introns. A spliced read is placed by its splicing.
2. **Exonic bases.** Prefer the gene with the most exonic bases (overlap with its exon
   union). Exonic overlap under 20 bp counts as zero, so a read end grazing an exon does not
   outrank the gene whose intron holds the read. This handles nested genes (a read on a
   nested gene's exon beats the host's intron) and tandem genes (the gene holding more of the
   read wins). It follows the precedent of Cell Ranger and STARsolo
   `GeneFull_ExonOverIntron` / `GeneFull_Ex50pAS`: exonic beats intronic.
3. **Containment.** Prefer a gene whose span contains the whole read.
4. **Isoform-read weighting** (`gene_body_tie_weighting = isoform_reads`). Weight each tied
   gene by its isoform-assigned reads in this contig-strand unit, counted over the full
   library.
   - **Dominance.** If one gene holds at least `gene_body_dominance_frac` (0.9) of the total
     weight, it takes the whole read.
   - **Otherwise** the read is split in proportion to the weights; genes with no weight get
     nothing.
   - The weights exist only once the stream has ended. Tied reads are therefore spilled to a
     temp file and resolved in `GeneBodyRecorder.finish()`, so memory stays bounded by the
     number of genes.
5. **3' end.** Applies only when no tied gene has any isoform reads. Prefer the gene whose
   span holds the read's 3' end (strand-aware). It comes after rule 4 (user, 2026-11-18):
   with internal priming, the 3' end is a priming site rather than a transcript end.
6. **Innermost gene.** Prefer the smallest span.
7. **Still tied.** `gene_body_ambiguous_policy`: `drop` (default) writes rows with frac 0,
   so the read is recorded but counted nowhere; `split` gives each gene 1/n. The outcome is
   deterministic.

`GeneBodyRecorder.counts` (logged per contig-strand) records which rule settled each tie:
`tie_dominance`, `tie_weight_split`, `tie_three_prime`, `tie_innermost`, `tie_ambiguous`.

**Readthrough history (ALS chr8 test).** In round 1, before any readthrough handling and with
"innermost" ranked ahead of expression, the GENCODE readthrough gene MRPS28-TPD52
(ENSG00000276418, 162 kb) took 553 intronic reads from TPD52 (196 kb). Its span is smaller
and covers TPD52's 3' region. Weighting by isoform reads fixed that (ENSG00000276418 down to
19 reads). Excluding readthrough genes, now the default, takes it to 0 and removes every tie
in the window (section 5).

### Comparison with other tools (as documented by each tool)

| | exonic vs intronic | multi-gene read | antisense | S/U/A |
|---|---|---|---|---|
| Cell Ranger >= 7 (include introns) | exonic if >= 50% on exons; otherwise intronic if inside a gene | not counted | not counted (3' sense chemistries) | no |
| STARsolo GeneFull | any overlap with the gene's full span | not counted | no | Velocyto feature |
| STARsolo GeneFull_ExonOverIntron / _Ex50pAS | exonic overlap prioritized; Ex50pAS mimics Cell Ranger's 50% rule | not counted | per variant | Velocyto feature |
| alevin-fry USA (splici/spliceu) | S targets = transcripts, U = gene bodies | `cr-like`: discard; `cr-like-em` / `parsimony-em`: EM across genes | no | S, U, A layers |
| velocyto | intron overlap makes a read U | discarded | no | spliced / unspliced / ambiguous |
| **LRAA gene_body (this)** | junction match, then exonic bases, then containment | isoform-read weight: dominant gene (>= 0.9) takes it, else proportional split; with no weights, 3' end, innermost, then drop/split | never (strand-split passes) | S, U, A per row |

Like Cell Ranger's references, LRAA leaves readthrough genes out by default, but through a
config key rather than by editing the GTF. Long reads make S/U/A cleaner than with short
reads: a junction is observed directly in the read, not inferred from which k-mers a fragment
hits.

## 3. Outputs

**Sidecar** `<prefix>.quant.gene_body.tracking.gz` starts with LRAA's version and command
comment lines, then this header:

```
gene_id  read_name  frac_assigned  splice_class  source  n_genes  aligned_bp  exonic_bp  intronic_bp
```

- `gene_id` is the output gene id (`symbol^ENSG...`), the same as in the tracking file.
- `read_name` is `barcode^umi^name`.
- There is one row per (read, gene). Fractions sum to 1 per read. The exception is a read
  tied under `drop`, whose rows have frac 0.
- `source` is `isoform` or `gene_body:<why>` (section 1).
- Rows for tied reads come after the rest of their unit, so the file is not in read order.

**Matrices.** `util/sc/gene_body_tracking_to_sparse_matrix.py --gene_body_tracking <sidecar>
--output_prefix P` writes `P^gene_full-sparseM` and `P^gene_{spliced,unspliced,ambiguous}-sparseM`.

- The layout is 10x (features x barcodes), the same as the existing builder.
- All four share one feature list and one barcode list, so S + U + A = gene_full entry for
  entry.
- Reads without both cell tags are skipped.
- `--exclude_gene_body` gives the isoform-only reads. Their sum equals the existing
  `^gene-sparseM` cell by cell (tested).
- `gene_full` sits **beside** the existing gene matrix, which is unchanged.

**Config keys.** All are set through `--config_update`; no new CLI flags.

| key | default | meaning |
|---|---|---|
| `gene_body_counts` | false | write the sidecar; needs `--stream_reads` (and `--gtf` under `--quant_only`), else LRAA exits with an error at startup. In discovery mode the genes are the models the run reports |
| `gene_body_min_read_frac_in_span` | 0.5 | candidate-gene threshold (rule 0) |
| `gene_body_exclude_readthrough_genes` | true | readthrough-tagged genes are never candidates (rule 0) |
| `gene_body_min_feature_bp` | 20 | exonic evidence (rule 2) and U threshold |
| `gene_body_tie_weighting` | isoform_reads | rule 4 on (`isoform_reads`) or off (`none`) |
| `gene_body_dominance_frac` | 0.9 | rule 4 dominance cut, must be in (0.5, 1]; 1.0 means always split |
| `gene_body_ambiguous_policy` | drop | rule 7 |

## 4. Where it lives, and what it interacts with

**Code:**

- `pylib/GeneBodyAssigner.py` (new):
  - `GeneBodyIndex`: gene models per contig-strand, readthrough flags, overlap query,
    `assign` (rules 0-3) and `resolve` (rules 4-7);
  - `GeneBodyRecorder`: sidecar rows, spill and `finish`;
  - `segments_from_cigar` and `introns_from_segments`.
- `pylib/Transcript.py`: the parser keeps GTF `tag` values in `get_gtf_tags()`.
- `pylib/StreamingQuant.py:stream_assign` takes a new `gene_body_recorder` argument.
  - `to_gene_body` is called at the three unassigned exits after rescue, at the assigned
    exit with the rows just written, and for rescued low-per-id reads.
  - Supporting changes: `assign_path` gains a `sink` argument, and `offer_to_rescue` returns
    its rows.
  - Without a recorder, these paths behave exactly as before.
- `LRAA:run_quant_only` takes a new `gene_body_fh` argument.
  - It builds the index from `reference_transcripts` (the annotation as supplied), so a gene
    whose models all failed to map still owns its reads.
  - It calls `finish()` after the stream.
  - **Early exits.** On its two early exits with an annotation present ("no splice graph
    created", "no input transcripts mapped into the splice graph"), it now runs
    `_gene_body_geometry_only_pass`. That pass applies the same retention filter under the
    final-quant thresholds, so those units' reads still get gene-body rows
    (`gene_body:no_quant_unit`). Such a unit has no isoform reads, so rule 4 has no weights
    there. This closes the round-1 gap where those units wrote a header-only sidecar.
- `LRAA:_run_oversimplify_best_overlap` (chrM by default) takes `gene_body_fh`.
  - Each best-overlap read gets an `isoform` row for its chosen model's gene, so MT- genes
    reach gene_full and `percent.mt`.
  - A read with no exonic overlap goes to `GeneBodyIndex`
    (`gene_body:oversimplify_no_exon_overlap`).
  - The read segments come from the Pretty_alignment the path already reads.
  - The ref-guided-discovery oversimplify branch passes it too: there the reference models
    are carried forward and quantified the same way.
- `LRAA:_run_oversimplify_aggregate` (de novo `--oversimplify`, one pseudo-gene per
  contig-strand) takes `gene_body_fh`: each read gets an `isoform` row for the aggregate
  gene, so a de novo chrM reaches gene_full and `percent.mt`.
- **Discovery mode.** `run_transcript_assembly` takes `gene_body_fh` and hands it to its
  final `run_quant_only` call, with a new `gene_body_index_transcripts` argument set to
  the models it is about to report. The index is built from those, not from
  `reference_transcripts` (only the reference subset there, nothing in de novo), so the
  sidecar's gene ids are the ones `quant.tracking` carries. Only the final quant records;
  the draft quant does not stream.
- `LRAA:_contig_job_runner` opens the per-unit sidecar whenever the key is on.
- `LRAA:_run_jobs_schedule_and_merge`:
  - concatenates the per-unit sidecars;
  - removes them during cleanup;
  - refuses to reuse a finished shard that lacks a sidecar when the key is on (tested).
- `pylib/ChunkedRun.py:merge_gene_body_tracking` merges the units' sidecars (all or none;
  some-but-not-all is an error). In discovery it renames gene ids through the same id plan
  the tracking merge applies (`merged_id_plan` / `mapped_id`), so an id two chunks minted
  for different models is renamed identically in both files. `LRAA:_run_chunked_mode`
  publishes the result and exits with an error if the key is on and no sidecar came back.
  Chunked mode is the default.
- `LRAA` main validates the settings right after `--config_update` is applied.

**Interactions:**

- **EM.** None. Gene-body reads never enter `mp_counter`, the assignment table or theta, so
  `quant.expr` and the tracking file are identical with the key on or off.
- **Coverage normalization.** The sidecar comes from the second, streaming pass over the
  full `--bam`, so the counts and rule-4 weights are full-library numbers.
- **Mismapping filter.** It runs on the merged discovery outputs and rewrites only the gtf
  and `quant.expr`. Like `quant.tracking`, the sidecar keeps the rows of a model it drops,
  so gene_full and the standard gene matrix stay consistent with each other.
- **Transcriptome rescue.** Rescue runs first. A rescued read is `isoform`.
- **Batch path (`--no_stream_reads`).** Refused at startup.
- **Discovery mode's final quant** (`run_transcript_assembly` -> `run_quant_only`,
  `report_quants=False`). Wired, since that pass already streams through
  `StreamingQuant.stream_assign`: no second pass is needed. De novo models carry no GTF
  tags, so readthrough exclusion applies only to reference genes in ref-guided mode.

### Single-nucleus mode in the WDLs (`sc_nucseq`)

`LRAA-singlecell.wdl` has `Boolean sc_nucseq = false`, never set from the data. When true:

1. **Gene-body counts in both whole-library passes.** `LRAA_runner.wdl`,
   `LRAA_chunk_scatter.wdl` and `LRAA.wdl` take `Boolean gene_body_counts`. The runner
   task merges `"gene_body_counts": true` into its own config file (the chunk leaf into its
   worker config), so the caller's `config_update_json` passes through untouched. The init
   pass gets it from `sc_nucseq` in every mode (discovery via the hook above, quant-only
   directly); the cluster-guided final quant gets it too; the per-cluster discovery runs
   do not. `init_config_update_json` stays unset, so `reuse_init_as_final_quant` (quant-only
   cluster-guided) still applies, and it hands the init pass's gene-body sidecar and
   matrices down with the rest.
2. **Merging the sidecar.** The runner emits `File? LRAA_gene_body_tracking`. `LRAA.wdl`
   concatenates the per-chromosome sidecars (`mergeGeneBodyTracking`,
   `util/misc/concat_gene_body_tracking.py`; it fails when a shard's is missing), the
   by_chunk merge stages and merges the units' sidecars, and `scattering = off` publishes
   its own. `mergedGeneBodyTracking` is the result.
3. **One matrix build, either way the standard matrices are built.**
   `sc_build_gene_body_sparse` (in `LRAA-build_sparse_matrices_from_tracking.wdl`) runs the
   builder over the merged sidecar and `util/sc/gene_body_cell_stats.py` over its output.
   It does not care whether the standard matrices came from shards
   (`sc_sparse_from_shards`) or the library-wide build. Outputs: gene_full and S/U/A
   tarballs, and a per-barcode table (`gene_full`, `spliced`, `unspliced`, `ambiguous`,
   `n_genes_full`, `spliced_frac`, `unspliced_frac`).
4. **Cell calling and clustering on gene_full.** `FilterGoodCells` takes
   `gene_full_sparse_tar_gz`: emptyDrops then runs on gene_full, and the standard gene,
   isoform and splice-pattern matrices are filtered to its cells
   (`filter_good_cells.R --gene_matrix_dir`). Seurat clusters the (filtered) gene_full
   with `min_features_nucseq` (1000) and `--nucseq`: variable genes are ranked within 10
   quantile bins of log10 gene span by the z-score of vst's standardized variance, as
   mean.var.plot z-scores dispersion within bins of mean. Spans come from the GTF the genes
   were named from (`initial_annot_gtf` in quant-only, the init GTF in discovery) through
   `util/sc/gtf_gene_spans.py`, which uses LRAA's `symbol^id` feature naming. The chosen
   genes and their bins are written out. With `sc_nucseq` off, the R scripts take their old
   path.
5. **Final quant and symbols.** `LRAA-cell_cluster_guided.wdl` (`sc_nucseq` input) builds
   the same matrices and table from the final quant. `Incorporate_gene_symbols.wdl` takes
   the four gene-body matrices as optional inputs and names them like the gene matrix.

## 5. Prototype results (small data only, workstation, every LRAA run under a 20 GB memory cap)

**Unit tests (new files).**

- `pylib/test_gene_body_assigner.py` (27 tests):
  - S / U / A classes;
  - antisense genes excluded;
  - the span threshold;
  - nested, readthrough, tandem and grazing-exon geometries;
  - junction-match precedence;
  - the dominance cut (at, above and below 0.9, and via the config key) and the
    proportional split;
  - weighting outranks the 3' end; the 3' end decides only without weights, on both
    strands; innermost; drop and split; determinism;
  - readthrough exclusion (on, off, mixed-tag genes, gene-level tag);
  - the parser keeps repeated GTF tags and leaves meta unchanged;
  - CIGAR parsing, recorder rows, spill-then-`finish`, and the chunked sidecar merge.
- `pylib/test_gene_body_lraa_paths.py` (2 tests): `--oversimplify` rows (isoform plus
  intronic), and the geometry-only pass for units with no isoform quant.

**Existing tests that still pass.** Oversimplify, streaming (resolver, rescue,
cross-component, served fraction, table fractions, XW single-cell), config resolution,
quant.expr schema, chunked tracking merge, read-assignment summaries, GTF re-export / region
restriction / tabix, reference provenance and gene-id output. In total 134 + 30, all passing.

**PBMC Kinnex chr19 fixture** (cells, `--no_chunk`, round-2 code):

- `quant.expr` (655 rows) and `quant.tracking` (7,329 rows) are identical with the key on
  and off.
- Sidecar `isoform` rows match the tracking file's per-(read, gene) fractions exactly (4,878
  reads).
- gene_full / standard = 1.07.

**ALS LP00714 nucseq, chr8:79.0-80.5 Mb** (14,175 reads; GENCODE v39 genes wholly inside
the window):

- `quant.expr` (90 rows) and `quant.tracking` (15,343 rows) are identical with the key on
  and off, and the `isoform` rows match tracking exactly.

*Rule 4 (old order) before the move.* With the round-1 order, 3,052 reads were still tied
after rules 0-3. The 3'-end rule narrowed **none** of them: every tie was between genes
whose spans all held the read's 3' end. Moving it after the weighting therefore changes
nothing in this window. After the move, `tie_three_prime` = 0.

*Before and after:*

- round 1: no readthrough exclusion, proportional split, 3' rule before weighting;
- no-RT-excl: round-2 code with `gene_body_exclude_readthrough_genes: false`;
- default: round-2 defaults;
- chunked: round-2 defaults with `--approx_MB_per_cut 1` (145 chunks).

| gene | standard gene matrix | round 1 | no-RT-excl | default | chunked |
|---|---|---|---|---|---|
| TPD52 | 1,536 | 4,079 | 4,079 | 4,080 | 4,080 |
| MRPS28 | 538 | 3,352 | 3,354 | 3,360 | 3,360 |
| STMN2 | 1,473 | 2,794 | 2,794 | 2,794 | 2,794 |
| HEY1 | 541 | 553 | 553 | 553 | 553 |
| MRPS28-TPD52 readthrough (ENSG00000276418) | 0 | 19 | 17 | 0 | 0 |
| all genes | 4,307 | 11,560 | 11,560 | 11,550 | 11,550 |

What changed between the columns:

- **no-RT-excl.** All 3,052 ties were MRPS28 or TPD52 against the readthrough gene, and all
  were settled by dominance (`tie_dominance` = 3,052). The readthrough gene got some EM
  leakage weight in round 1's proportional split, but never 10% of it.
- **default.** No ties at all. Compared with no-RT-excl, the readthrough gene's 17
  gene-body reads are gone: its spliced readthrough molecules move to MRPS28 (+6) and TPD52
  (+1), and 10 reads that only it spanned become `reads_no_gene`.
- **Splice classes (default):** gene_full = 11,550 = S 881 + U 7,041 + A 3,628.
- **Chunked parity.** The chunked sidecar has the same rows as the unchunked one.
- **Resume parity.** A key-off run followed by a key-on rerun under the same prefix
  recomputes the shards and gives the same rows as a fresh key-on run.
- **Builder parity.** `--exclude_gene_body` equals the existing gene matrix cell by cell
  (3,609 entries; round 1, unaffected by round 2).

**chrM via `--oversimplify chrM`** (the same window plus 2% of the sample's chrM reads,
11,436 reads):

- 13 MT- genes now appear in gene_full, for example MT-ND3 1,463, MT-ATP6 1,410 and
  MT-CO2 1,362, all `isoform` rows, class A. The same rows come back from the chunked run.
- Sidecar `isoform` rows match tracking for all 13,407 isoform-assigned reads.
- **Deterministic tie-break (2026-11-18, user-approved).** Before it, on/off
  `quant.tracking` differed in exactly 1 read: a chrM read tied between two MT models
  (`is_unique=0`) that the unseeded `random.choice` broke differently from run to run, and
  three key-off runs differed from each other. With `_oversimplify_tie_break`, three key-off
  runs give identical `quant.expr` (103 rows) and `quant.tracking` (24,443 rows), key-on
  equals key-off on both, and the chunked run equals the unchunked one (`quant.expr` up to
  row order; tracking compared without the per-process `mp_id` column).
  Test: `pylib/test_oversimplify_deterministic_ties.py`.

**Cost.** In round 1, every read through `assign()` against all GENCODE v39 chr8 genes took
24-39 us. Round 2: ALS window wall time 96.5 s on versus 98.0 s off (`--no_chunk`, 4 CPUs).

**Discovery mode and `sc_nucseq` (2026-11-19, chr19 2 Mb PBMC Kinnex fixture,
`testing/single_cells/sc_full_pipe`, whole cells, so few intronic reads).**

- LRAA CLI, chunked at 0.2 Mb, key on vs off: de novo and ref-guided gtf, `quant.expr` and
  `quant.tracking` byte-identical. The sidecar's `isoform` rows equal the tracking file per
  (read, gene): 342,139 pairs de novo, 344,807 ref-guided, 0 mismatches, with the chunk id
  plan applied in both. gene_full / tracking = 1.06 (S 87%, U 4%, A 9% de novo).
  Strandedness check: 0 antisense spliced reads (oriented library), no warning.
- Workflow, local images built from the branch, `sc_nucseq` off vs the branch's base
  commit: quant-only basic on the shard path (17 outputs) and on the library-wide path with
  emptyDrops (22), cluster-guided quant-only (32) and cluster-guided de novo (53) all
  identical, except a chunk plan's `bam_identity` (a path fingerprint) and the
  cluster pseudobulk matrices' row order (`util/build_LRAA_expr_matrices.py` iterates a set;
  same rows).
- `sc_nucseq` on (min_features_nucseq 1 for the fixture): init gtf, tracking and standard
  matrices unchanged; gene_full 90 genes x 12,655 barcodes (standard 82 x 12,427),
  S + U + A = gene_full entry for entry; per-barcode table for every barcode; emptyDrops on
  gene_full kept 3,003 cells (2,861 on the gene matrix, `lower_threshold` 20); clusters
  from gene_full. Cluster-guided quant-only reused the init pass's gene-body outputs as the
  final ones (no second quant); de novo built them from the final quant (177 genes). In de
  novo, 4 of 165 genes had no span: the mismapping filter had dropped their models from the
  init gtf, while the tracking file keeps their reads.

### Known limitations

- Classification and the junction rule read the raw CIGAR, while the isoform verdict uses
  LRAA's soft-clip-corrected, terminal-intron-pruned alignment. A terminal intron that LRAA
  prunes still counts as a junction for rule 1 and for S.
- Rule 4 weights come from the unit's own isoform reads. A unit routed through the
  geometry-only pass has none, so its ties fall to the 3' end, innermost and drop rules.
- The matrix builder holds rows as integer (gene, barcode) triples and folds them into
  sparse matrices every 5 M rows, so memory follows the nonzeros (2026-11-19; it was a
  string-keyed dict). Not yet run at genome scale.

## 6. Decisions on the round-1 questions (user, 2026-11-18)

1. A (unspliced, entirely exonic) stays a separate category; gene_full = S + U + A. Done.
2. Ties are split by isoform-read weight, with a dominance cut: a gene with at least 90% of
   the weight takes the whole read (`gene_body_dominance_frac`). Done.
3. No downstream 3' allowance; candidates are span-only with the 0.5 fraction. Unchanged.
4. Unstranded libraries get no special handling; inputs are expected to be strand-specific.
   Documented (section 1 and the module docstring).
5. chrM: the `--oversimplify` path writes gene-body rows. Done.
6. gene_full sits beside the gene matrix. In single-nucleus mode clustering uses gene_full.
   Done: `sc_nucseq` in the single-cell workflow (section 4), with the discovery-mode hook.
7. Readthrough-tagged genes are excluded as candidates, on by default
   (`gene_body_exclude_readthrough_genes`), with GTF tags now parsed. Done.
8. The 3'-end rule moved after the weighting. It decided 0 reads in the ALS window, before
   and after. Done.

## 7. Still open

- A Terra or fleet-VM run of `sc_nucseq` on a full single-nucleus sample (only the 2 Mb
  test fixture has been run through the workflow).
- Sizing the matrix-build task (`memoryGBbuildSparseMatrices`) on a genome-scale sidecar.
- Optionally, a per-unit gene-body summary TSV (the counts are logged per contig-strand now).
- Running on the full ALS sample on a fleet VM or Terra (never on this workstation), and
  comparing with Cell Ranger include-introns counts if they exist for LP00714.
