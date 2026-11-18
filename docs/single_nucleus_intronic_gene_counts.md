# Gene counts that include intronic reads (single-nucleus)

Status: design plus a working prototype (branch `sn-intronic-gene-counts`, not on devel).
Opt-in via `--config_update '{"gene_body_counts": true}'`. With the key off, LRAA's outputs
are byte-for-byte unchanged.

## The problem

In quant-only mode a read is assigned to an isoform only if it is structurally compatible
with it: at least `fraction_read_align_overlap` (0.75) of the read must lie on the isoform's
exons. A pre-mRNA read lying in an intron matches no isoform. It gets no tracking row, so it
is missing from the isoform matrix and also from the gene matrix that the `util/sc` builders
derive from the tracking file. Single-nucleus libraries are mostly pre-mRNA. In the ALS
LP00714 nucseq test window below, 87% of reads are unspliced, and the standard gene matrix
holds 4,307 reads of the 11,560 that lie in annotated genes on their own strand.

Cell Ranger ("include introns"), STARsolo `GeneFull` and alevin-fry's USA mode all count
these reads toward their gene. This design adds the same thing to LRAA.

## What "chr-simplify" turned out to be

The repo has no `chr_simplify` or `simplify` method. The nearest match is **`--oversimplify`**
(`LRAA_Globals.config["oversimplify_enabled"/"oversimplify_contigs"]`, default `chrM,M` in
the single-cell workflow). Its quant-only implementation, `LRAA:_run_oversimplify_best_overlap`:

- works per contig;
- assigns each read to the single transcript with the most exonic bases of overlap, breaking
  ties at random;
- does no EM;
- skips any read with zero exonic overlap, so a deep-intronic read is dropped there too.

The new mode generalizes that idea in three ways. It works on genes (exon union plus gene
span) rather than transcripts. It scores intronic overlap as well as exonic. It runs
genome-wide and *adds to* the normal isoform quantification instead of replacing it.
`--oversimplify` itself is untouched.

## Design in brief

During the streaming quant pass, every retained read already gets an isoform verdict. When
`gene_body_counts` is on, each read is also written to a **sidecar**
`<prefix>.quant.gene_body.tracking.gz`:

- **Isoform-assigned reads** keep their isoform assignment. The EM fractions are summed per
  gene, so these rows reproduce what the standard gene matrix counts (verified below).
- **Reads the isoform pass leaves unassigned** (no graph path, spacer path, or a path that
  matches no isoform) go to the most appropriate same-strand gene by the rules below.
- **Every row** carries a splice class: S (spliced), U (unspliced) or A (ambiguous).

`quant.expr` and `quant.tracking` are not touched. A new builder,
`util/sc/gene_body_tracking_to_sparse_matrix.py`, turns the sidecar into gene matrices: total
(`gene_full`) and S / U / A.

## 1. Read classes

Strand comes from LRAA's existing read-strand logic. `Util_funcs.transcribed_strand` uses the
`ts` tag, corroborated by splice motifs, and falls back to the aligned strand. Each read is
processed only in the pass for its own contig-strand, so a gene on the other strand is never
a candidate.

| class | definition | handled by |
|---|---|---|
| exonic-compatible | isoform-compatible (unchanged: graph path plus 0.75 overlap rule, EM) | existing quant; sidecar `source=isoform` |
| intronic / pre-mRNA / partial | retained read, no isoform verdict, inside a same-strand gene span | `GeneBodyIndex.assign`; sidecar `source=gene_body:<why>` |
| intergenic / antisense | no same-strand gene span holds >= 50% of the read | counted (`reads_no_gene`), not written |
| filtered | rejected by `quant_discard_reason` (MAPQ, per-id, rDNA, ...) | excluded, as for isoforms |
| cross-component ambiguous | path spans isoforms of two EM components (existing rule) | excluded, as rescue excludes it |

`<why>` is `no_path`, `spacer_path` or `unassignable`. Transcriptome rescue, when on, gets
first claim on these reads. A rescued read counts as `isoform`.

**Splice class** is computed against the gene the read is credited to. The gene's intronic
sequence is its span minus the union of its exons.

- **U**: the read covers at least `gene_body_min_feature_bp` (20) bases of that intronic
  sequence. This holds whether the read is spliced or not, so partially spliced molecules
  with a retained intron count as U.
- **S**: not U, and the read has at least one junction.
- **A**: not U and unspliced, so it lies entirely on exons. It fits mature mRNA and pre-mRNA
  equally well, which is alevin-fry's "ambiguous".

An isoform-assigned read can still be U: the 0.75 rule lets up to 25% of a read run into an
intron. That is reported honestly.

## 2. Choosing "the most appropriate gene"

Rules for a read with no isoform assignment, applied in order
(`pylib/GeneBodyAssigner.py`). Only genes on the read's own strand are considered.

0. **Candidate genes.** A gene is a candidate if its span (first exon start to last exon
   end, over all isoforms) holds at least `gene_body_min_read_frac_in_span` (0.5) of the
   read's aligned bases. This allows unannotated 3' UTR extensions up to half the read and
   rejects reads that are mostly intergenic.
1. **Junctions.** Prefer the gene with the most read junctions that exactly match one of its
   annotated introns. A spliced read is placed by its splicing.
2. **Exonic bases.** Prefer the gene with the most exonic bases (overlap with its exon
   union). Exonic overlap under 20 bp counts as zero, so a read end grazing an exon does not
   outrank the gene whose intron holds the read. This rule handles nested genes (a read on a
   nested gene's exon beats the host's intron) and tandem readthrough (the gene holding more
   of the read wins). It is the precedent STARsolo `GeneFull_ExonOverIntron` /
   `GeneFull_Ex50pAS` and Cell Ranger set: exonic beats intronic.
3. **Containment.** Prefer a gene whose span contains the whole read.
4. **3' end.** Prefer a gene whose span holds the read's 3' end (strand-aware). This matters
   for 3' libraries.
5. **Expression-weighted split.** If genes are still tied, split the read among them in
   proportion to each gene's isoform-assigned reads in this contig-strand unit, counted over
   the full library. A gene with no isoform reads gets nothing. The isoform counts exist
   only at the end of the stream, so tied reads go to a temp file and are resolved in
   `GeneBodyRecorder.finish()`. Memory stays bounded by the number of genes.
   Config `gene_body_tie_weighting` takes `isoform_reads` (default) or `none`.
6. **Innermost gene.** If no tied gene has isoform reads, prefer the gene with the smallest
   span: a gene nested in a host's intron, or a constituent gene inside an annotated
   readthrough gene.
7. **Still tied.** `gene_body_ambiguous_policy`: `drop` (default) writes rows with frac 0,
   so the read is recorded but counted nowhere; `split` gives each gene 1/n. The outcome is
   deterministic. There is no `random.choice`, unlike `--oversimplify`.

**Why rule 5 comes before rule 6 (found on real data).** On the first ALS run, with
innermost-span as rule 5, the GENCODE readthrough gene MRPS28-TPD52 (ENSG00000276418,
protein_coding, 162 kb) took 553 intronic reads from TPD52. Its span is *smaller* than
TPD52's (196 kb) and covers TPD52's 3' region, so "innermost" picked the wrong gene. That
readthrough gene has 0 isoform-assigned reads, against 1,536 for TPD52. With the
expression-weighted rule, ENSG00000276418 drops to 19 reads, mostly genuine spliced
MRPS28-TPD52 molecules that rule 1 places by their junctions, and TPD52 rises from 3,525 to
4,079. Cell Ranger's references avoid this by removing readthrough genes from the GTF.
LRAA's GTF parser keeps no `tag` attribute, and the rule here does not depend on the
organism.

### Comparison with other tools (as documented by each tool)

| | exonic vs intronic | multi-gene read | antisense | S/U/A |
|---|---|---|---|---|
| Cell Ranger >= 7 (include introns) | exonic if >= 50% on exons; otherwise intronic if inside a gene | not counted | not counted (3' sense chemistries) | no |
| STARsolo GeneFull | any overlap with the gene's full span | not counted | no | Velocyto feature |
| STARsolo GeneFull_ExonOverIntron / _Ex50pAS | exonic overlap prioritized; Ex50pAS mimics Cell Ranger's 50% rule | not counted | per variant | Velocyto feature |
| alevin-fry USA (splici/spliceu) | S targets = transcripts, U = gene bodies | `cr-like`: discard; `cr-like-em` / `parsimony-em`: EM across genes | no | S, U, A layers |
| velocyto | intron overlap makes a read U | discarded | no | spliced / unspliced / ambiguous |
| **LRAA gene_body (this)** | junction match, then exonic bases, then containment, then 3' end | split by the genes' isoform-assigned reads; else innermost; else drop/split | never (strand-split passes) | S, U, A per row |

Long reads make S/U/A cleaner than with short reads. A junction is observed directly in the
read, not inferred from which k-mers a fragment hits.

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
- `source` is `isoform` or `gene_body:{no_path,spacer_path,unassignable}`.
- Rows for tied reads come after the rest of their unit, so the file is not in read order.

**Matrices.** `util/sc/gene_body_tracking_to_sparse_matrix.py --gene_body_tracking <sidecar>
--output_prefix P` writes:

- `P^gene_full-sparseM`
- `P^gene_spliced-sparseM`
- `P^gene_unspliced-sparseM`
- `P^gene_ambiguous-sparseM`

The layout is 10x (features x barcodes), the same as the existing builder, which it reuses
for writing. The four share one feature list and one barcode list, so S + U + A equals
gene_full entry for entry. Reads without both cell tags are skipped, as in the existing
builder. `--exclude_gene_body` gives the isoform-only reads split by class. The prototype
accumulates in a dict (bounded by non-zero entries, not by reads). At genome scale it should
reuse `singlecell_tracking_to_sparse_matrix.LevelAccumulator`'s chunked path.

Nothing existing changes. The tracking schema, `_append_without_first_line`'s field-count
padding, ChunkedRun's header-equality checks, the direct reader's column positions and the
quant.expr TPM denominators all assume the current tracking schema. Adding sentinel rows or
columns there would have broken them, or inflated isoform TPMs.

**Config keys.** All are set through `--config_update`; no new CLI flags.

| key | default | meaning |
|---|---|---|
| `gene_body_counts` | false | write the sidecar (quant-only with `--gtf` and `--stream_reads`; otherwise LRAA exits with an error at startup) |
| `gene_body_min_read_frac_in_span` | 0.5 | candidate-gene threshold (rule 0) |
| `gene_body_min_feature_bp` | 20 | exonic evidence (rule 2) and U threshold |
| `gene_body_tie_weighting` | isoform_reads | rule 5 on/off |
| `gene_body_ambiguous_policy` | drop | rule 7 |

## 4. Where it lives, and what it interacts with

**Code:**

- `pylib/GeneBodyAssigner.py` (new): `GeneBodyIndex` (gene models per contig-strand, overlap
  query, `assign`, `resolve`), `GeneBodyRecorder` (sidecar rows, spill and `finish`),
  `segments_from_cigar`.
- `pylib/StreamingQuant.py:stream_assign` takes a new `gene_body_recorder` argument.
  `to_gene_body` is called at the three unassigned exits after rescue, at the assigned exit
  with the rows just written, and for rescued low-per-id reads. `assign_path` gains a `sink`
  argument and `offer_to_rescue` returns its rows. Without a recorder these paths are
  unchanged.
- `LRAA:run_quant_only` takes a new `gene_body_fh` argument. It builds the index from
  `reference_transcripts` (the annotation as supplied), not from the graph-filtered
  `input_transcripts`, so a gene whose models all failed to map still owns its reads. It
  calls `finish()` after the stream.
- `LRAA:_contig_job_runner` opens the per-unit sidecar when quant-only and the key is on. The
  file has a header only on `--oversimplify` contigs.
- `LRAA:_run_jobs_schedule_and_merge`:
  - concatenates the per-unit sidecars;
  - removes them from the cleanup list;
  - refuses to reuse a shard (resume) that lacks the sidecar when the key is on. This bug was
    found while testing: a rerun silently reused shards written without the sidecar.
- `pylib/ChunkedRun.py:merge_gene_body_tracking` merges the units' sidecars (all or none;
  some-but-not-all is an error). `LRAA:_run_chunked_mode` publishes the result and exits with
  an error if the key is on and no sidecar came back. Chunked mode is the default.
- `LRAA` main checks the settings right after `--config_update` is applied.

**Interactions:**

- **EM.** None. Gene-body reads never enter `mp_counter`, the assignment table or theta.
  Isoform fractions, TPM, `fraction_read_align_overlap` and `quant.expr` are identical with
  the key on or off (tested).
- **Coverage normalization.** The first (EM) pass reads the normalized BAM. The sidecar is
  written by the second, streaming pass over the full `--bam`. Gene-body counts and the
  rule-5 weights are therefore full-library numbers, in the same units as the tracking file.
- **Mismapping filter.** `_maybe_filter_mismappings` runs only in discovery mode. The single-
  cell final quant is quant-only against an already-filtered merged GTF, so there is no
  interaction.
- **Transcriptome rescue.** Rescue runs first. A rescued read is `isoform`.
- **`--oversimplify` contigs (chrM).** The prototype writes nothing for them (header only),
  so `gene_full` has no MT- genes. Open item: emit `source=isoform` rows from
  `_run_oversimplify_best_overlap`. It already has the read segments and the chosen
  transcript.
- **Batch path (`--no_stream_reads`).** Not supported; LRAA refuses the combination at
  startup. Supporting it would need a second hook after `q.report_quant_results`.
- **Discovery mode's final quant** (`run_transcript_assembly` -> `run_quant_only`,
  `report_quants=False`): not wired. The single-cell workflow's final quant is quant-only,
  so this is not needed there.
- **Single-cell WDL** (`WDL/LRAA-cell_cluster_guided.wdl`). The final quant is
  `LRAA_final_quant` (quant-only). Passing
  `final_quant_config_update_json = {"gene_body_counts": true}` reaches the workers: chunk
  workers receive `--config_update` keys (`_collect_explicit_config_overrides`). The rest is
  design only, not prototyped:
  - `WDL/subwdls/LRAA_runner.wdl` needs an optional output for the sidecar;
  - `LRAA.wdl` needs to gather it across `by_chromosome` scatters (concatenation; reads do
    not span chromosomes);
  - `LRAA_chunk_scatter.wdl` / `util/misc/merge_chunk_outputs.py` need the same, since
    `merge_gene_body_tracking` is file-keyed and callable from there;
  - the cluster-guided workflow needs to output the sidecar and run the new matrix builder
    next to `LRAA-build_sparse_matrices_from_tracking.wdl`;
  - the shard path (`sc_sparse_from_shards`, `tracking_to_shard_sparse.py`) needs a
    per-shard sidecar or a post-merge build.
- **Cell clustering.** Whether downstream Seurat clustering (`LRAA-gene_sparseM_to_seurat_clusters.wdl`)
  should use `gene_full` for nuclei is a user choice. It would change cluster calls.
- **Size.** The sidecar adds roughly one row per retained read, about 120 bytes before gzip,
  compared with about 2.2 rows per assigned read in the tracking file.

## 5. Prototype results (small data only, on the workstation, memory-capped)

The unit tests are in `pylib/test_gene_body_assigner.py` (23 tests). They cover:

- U / S / A classes, including a partially spliced read and exon-edge wobble under 20 bp;
- antisense genes excluded;
- the span threshold;
- nested genes (exon versus host intron, the nested gene's own intron, host intron outside);
- junction-match precedence;
- annotated readthrough genes;
- tandem genes and a grazing exon;
- the 3' end rule on both strands;
- the expression-weighted split beating innermost (the TPD52 geometry);
- drop and split ties, and determinism;
- CIGAR parsing;
- recorder rows and per-gene EM sums;
- spill-then-`finish` resolution;
- the chunked sidecar merge (all, none, partial).

The existing streaming, rescue, config-resolution, chunked tracking-merge and summary tests
also pass (90 + 55).

**PBMC Kinnex chr19 fixture** (`testing/single_cells/data`, cells, not nuclei; `--no_chunk`):

- `quant.expr` (655 rows) and `quant.tracking` (7,329 rows, ignoring the per-process
  `mp_id`) are identical with the key on and off.
- Sidecar `isoform` rows match the tracking file's per-(read, gene) fractions exactly (4,878
  reads, 0 mismatches).
- 342 more reads get a gene, so `gene_full` / standard = 1.07. Of the gene-body reads, 214
  are U, 121 are S (novel or partial junctions) and 7 are A.

**ALS LP00714 nucseq, chr8:79.0-80.5 Mb** (14,175 reads, 13% spliced, all spliced reads
`ts:A:+`, so the library is oriented; GENCODE v39 genes wholly inside the window):

- `quant.expr` and `quant.tracking` are again identical on and off, and the `isoform` rows
  match the tracking file exactly.
- Gene totals rise from 4,307 to 11,560 (x2.68): S 879, U 7,053, A 3,628. On the minus
  strand, 3,052 reads were tied and resolved by rule 5.
- 2,294 reads had no same-strand gene. Many lie in genes cut by the window edge, so they are
  not in the test GTF.

| gene | standard gene matrix | gene_full | S | U | A |
|---|---|---|---|---|---|
| TPD52 | 1,536 | 4,079 | 121 | 2,567 | 1,391 |
| MRPS28 | 538 | 3,352 | 469 | 2,604 | 279 |
| STMN2 | 1,473 | 2,794 | 248 | 1,315 | 1,231 |
| HEY1 | 541 | 553 | 22 | 0 | 531 |
| ENSG00000249328 (lncRNA) | 11 | 312 | 8 | 297 | 7 |
| MRPS28-TPD52 readthrough | 0 | 19 | 3 | 16 | 0 |

**Chunked parity.** With chunked mode at `--approx_MB_per_cut 1` (145 chunks, 144 cuts, 0
reads severed), the sidecar has the same 14,633 rows as the unchunked run, ignoring order.

## 6. Open questions for the user

1. **Should the S / U / A split be kept as built?** Long-read U includes partially spliced
   molecules, and A counts toward neither layer in velocity tools by default. Should A count
   as spliced (velocyto's convention for exonic reads), stay separate, or be split?
2. **Rule-5 weighting** splits a tied intronic read by the genes' isoform-assigned
   (mature-RNA) reads. An alternative is "winner takes all", or dropping such reads as Cell
   Ranger and STARsolo do. Splitting gives fractional UMIs; is that acceptable?
3. **Candidate threshold 0.5 and the 20 bp feature threshold.** A read that runs well past an
   annotated 3' end is called intergenic. Should a downstream-extension allowance be added?
4. **Orientation.** For unspliced reads the strand is the aligned strand. That is right for
   oriented Kinnex / 10x reads (all spliced reads in the ALS window carry `ts:A:+`), but it
   would split an unstranded library across strands. Should LRAA refuse, or warn, when
   `gene_body_counts` meets a library with many `ts:A:-` spliced reads?
5. **chrM.** Write `isoform`-source rows from the `--oversimplify` path so `percent.mt` works
   on `gene_full`? (Probably yes.)
6. **Which matrix drives clustering** for nucleus data in the cluster-guided workflow:
   `gene_full` or the standard gene matrix? And should the new matrices sit beside the
   existing ones or replace the gene level when the key is on?
7. **Readthrough genes.** Rule 5 handled MRPS28-TPD52. Should LRAA also offer to ignore
   `readthrough_transcript`-tagged genes (Cell Ranger's reference practice)? That would need
   the GTF parser to keep tags.

## 7. Next steps if adopted

- chrM rows from the oversimplify path;
- WDL plumbing (section 4);
- a genome-scale matrix builder (`LevelAccumulator` reuse);
- optionally a per-unit gene-body summary TSV (the counts are logged per contig-strand now);
- run on the full ALS sample on a fleet VM or Terra, never on this workstation, and compare
  with Cell Ranger include-introns counts if they exist for LP00714.
