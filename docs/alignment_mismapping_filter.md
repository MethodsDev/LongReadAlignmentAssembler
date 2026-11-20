# Whole-genome alignment-mismapping filter (v0.40.0)

Added 2026-09-24 on branch `devel`. ON BY DEFAULT. This document describes what
the filter removes, why it is safe, where it runs, and how to tune or disable it.
Measured results below are marked MEASURED; the corpus is the ONT cDNA Sequin
spike-in depth series (`chrIS`) and the PacBio/simulated benchmark corpora
(SIRVs, mouse, arabidopsis).

## What it removes

Some assembled isoforms are not transcripts: they are **alignment/strand
mismapping artifacts of a much-higher-expressed transcript**. Two forms dominate
in de novo / reference-free output:

- **Wrong-strand near-mirrors.** On unstranded/ambiguous data, a minority of a
  gene's reads get a cryptic canonical splice on the OPPOSITE strand a few bp
  from the real junction (the motif is genuinely GT-AG on the other strand at a
  shifted position, so it is not a motif violation). At depth, enough such reads
  accumulate for LRAA to assemble opposite-strand models that shadow the real
  gene. MEASURED (ONT Sequin, depth 2M): gffcompare labels 86 antisense
  exon-overlap (`x`) + 25 antisense intron-match (`s`) such models; they grow
  ~6x with depth (19 at 150k → 111 at 2M) while StringTie3/IsoQuant produce
  essentially none. v0.39.0's `ts` tag fix reduced this but does not remove it
  at depth.

- **Run-ons and chimeras.** 2-exon models with a single spurious long intron
  (MEASURED: up to ~106 kb, spanning genes), and same-strand mismappings, whose
  cDNA is near-identical to a real, higher-expressed transcript elsewhere.

All of these carry a small minority of the source transcript's reads (MEASURED
Sequin 2M: median ~0.01% of the opposing sense gene). Removing them and
renormalizing changes ~0.07% of the library (MEASURED) — negligible, which is
why the quant is repaired by row removal rather than requantification.

## The two detectors (unioned)

A model is removed if EITHER detector flags it.

### 1. MIRROR (coordinate; needs no cDNA)

A multi-exon model `M` is flagged when there is an opposite-strand model `N` with:

- exonic base overlap `>= mismap_min_base_overlap` (default 0.5), and
- EVERY internal splice site of `M` within `mismap_junction_tolerance` bp
  (default 20) of an exon boundary of `N`, and
- `expr(M) < mismap_max_expr_fraction * expr(N)` (default 1%), taking the
  highest-expressed qualifying `N`.

Catches the wrong-strand `s` near-mirrors. MEASURED: 0 real transcripts removed
on SIRVs (all 12, incl. E2), mouse, arabidopsis, because a genuine antisense
isoform has its OWN splice structure (its sites do not sit on the sense gene's
boundaries), and on PacBio HiFi there are no cryptic-strand reads to build such
a model in the first place.

### 2. SEQUENCE (minimap2 cDNA all-vs-all)

Spliced cDNA is extracted for every model (from the genome fasta) and aligned
all-vs-all with minimap2 (`-x map-ont`; identity from the PAF `de` tag). A model
`M` is flagged when its best hit to a **different-gene** model `N` is:

- `>= mismap_min_seq_identity` % identity (default 99) over
  `>= mismap_min_seq_coverage` of `M`'s length (default 0.85, merged over HSPs),
  and
- `expr(M) < mismap_max_expr_fraction * expr(N)`.

Catches the `x` run-ons, chimeras, and same-strand mismappings.

Two design points that make this both effective and safe:

- **minimap2 (nucleotide), not diamond.** The artifacts are the reverse
  complement of a coding transcript; a translated/protein search (diamond,
  blastx) has no matching ORF and would MISS them.
- **The cross-gene requirement is load-bearing.** A genuine MINOR isoform is
  also near-identical and lower-expression than its gene's dominant isoform, but
  that match is SAME-gene and excluded. Only a mismapping matches a DIFFERENT
  highly-expressed gene. MEASURED: paralog exposure is confined to the
  low-expression tail by the `<1%` gate — on mouse the net effect is F1-positive
  at every identity threshold, and raising `mismap_min_seq_identity` toward 100
  cuts it further (paralogs have diverged; a mismapping is ~100% identical).

## Why no requant

The `<1%` gate guarantees what is removed is tiny (MEASURED: 0.0736% of reads on
Sequin 2M). The filter therefore rewrites the quant.expr directly:

- drop the removed rows,
- renormalize surviving `TPM` to sum to 1e6,
- recompute `isoform_fraction` and `unique_gene_read_fraction` ONLY for genes
  that lost a member (unchanged genes keep their exact original values),
- leave `RPM_total_reads` unchanged (its denominator is total sequenced reads).

Dropping (rather than reassigning) the removed reads is the correct semantics:
they are mismapping artifacts owed to no surviving model. A requant that
reassigned them could misattribute artifact reads to a neighbor.

## Where it runs

It is a WHOLE-GENOME, POST-MERGE stage. It requires the merged genome-wide model
set (the mirror partner or the sequence match can be anywhere), so it cannot run
per chunk.

### Direct `LRAA` (no WDL)

`LRAA` runs it in-process at the top level, after the per-contig models are
merged into `<prefix>.gtf` + `<prefix>.quant.expr`, before BED and splice
collapse. It rewrites both files in place and writes
`<prefix>.gtf.mismapping_filter.log`. A chunk worker (`LRAA --no_chunk` with the
worker env set) NEVER runs it on its slice; the chunked top level runs it once.
Disable with `--no_filter_mismappings`. No-op under `--quant_only`.

### WDL (`WDL/LRAA.wdl`)

Every per-shard `LRAA_runner` passes `--no_filter_mismappings`, and the workflow
runs the filter ONCE as a dedicated task, `alignment_mismapping_filter`, on the
merged genome-wide gtf + quant (whichever of the three scattering arms produced
them), before `splice_pattern_collapse`. Its filtered gtf/quant become the
`mergedGTF`/`mergedQuantExpr` outputs; the removal log is `mismappingFilterLog`.
Toggle with the workflow input `filter_mismappings` (default true).

```
 merge (by_chunk | by_chromosome | off)
        └─ merged whole-genome gtf + quant.expr
              └─ alignment_mismapping_filter  (gtf + quant + genome
                    → filtered gtf + filtered quant.expr + log)
                    └─ splice_pattern_collapse → outputs
```

### Oversimplify contigs are exempt (v0.42.0)

Models on `--oversimplify` contigs (e.g. chrM) are carried forward from the
reference whatever the reads did, so the filter never removes them. Before this
exemption the filter's coverage dependence made sparsely and densely covered
cluster-guided inputs lose different chrM models, and the final merge of the
cluster-guided single-cell pipeline refused the disagreeing record sets. The
standalone util and the WDL task take `--oversimplify` / `oversimplify` for this.

## Configuration (`pylib/LRAA_Globals.py`)

| key | default | meaning |
|---|---|---|
| `filter_mismappings` | `True` | master switch (`--no_filter_mismappings` sets False) |
| `mismap_min_seq_identity` | `99.0` | % identity of the cDNA-vs-cDNA match |
| `mismap_min_seq_coverage` | `0.85` | fraction of the query cDNA the match must cover |
| `mismap_max_expr_fraction` | `0.01` | remove only if expr < this fraction of the match's expr |
| `mismap_junction_tolerance` | `20` | bp tolerance, mirror splice site vs opp-strand boundary |
| `mismap_min_base_overlap` | `0.5` | mirror: min exonic base-overlap fraction |

The standalone util exposes the same as `--min_seq_identity`,
`--min_seq_coverage`, `--max_expr_fraction`, `--junction_tolerance`,
`--min_base_overlap`.

## Components

- `pylib/AlignmentMismappingFilter.py` — the detectors and the gtf/quant rewrite.
- `util/filter_LRAA_isoforms_by_mismapping.py` — standalone CLI (the WDL task
  runs this): `--gtf --quant_expr --genome --output_gtf --output_quant_expr`.
- `LRAA` — `--no_filter_mismappings`; `_maybe_filter_mismappings` called
  post-merge in `_run_jobs_schedule_and_merge` (top level only) and
  `_run_chunked_mode`.
- `WDL/LRAA.wdl` — `filter_mismappings` input, `alignment_mismapping_filter`
  task; `WDL/subwdls/LRAA_runner.wdl` — `no_filter_mismappings` (default true).

## Validation (MEASURED)

Standalone on Sequin 2M (`chrIS`), defaults: 93 of 329 models removed = 76 of the
111 gffcompare-labeled antisense artifacts (`x`+`s`) plus 17 other artifacts, and
**0 of the 163 true-positive models**; filtered quant.expr sums to exactly 1e6.
Harmlessness on PacBio/sim: SIRVs (all 12, incl E2), mouse, arabidopsis remove 0
real transcripts with the mirror detector; the sequence detector's paralog
exposure is bounded by the `<1%` gate and net-F1-positive on mouse.

## Design history: rules tested and rejected (2026-09-23)

Before v0.40.0, three simpler "antisense filter" rules were scored offline by
re-scoring the v0.39.0 de novo benchmark GTFs (SIRVs, mouse and arabidopsis
simulations, cell lines). None improved identification F1:

- **Relative abundance alone** (remove an antisense model below 1% of the
  opposite-strand gene). Net harmful: it removed about 40 real chains per
  artifact. Real antisense isoforms (SIRV E2, simulations) are often below 1%,
  while the artifacts on those corpora sat above it, so no threshold on this
  axis helps. This is why expression is only ever a gate on top of a structural
  or sequence match, never a rule on its own.
- **Exon-footprint mirror without junction concordance** (high exon overlap with
  an opposite-strand model at low read fraction). Safer, since real antisense
  has its own splice structure, but still about 20 real for 1 artifact on the
  simulations, because real overlapping antisense gene pairs share exon
  footprint. Adding the requirement that EVERY splice site sit within 20 bp of
  the partner's exon boundaries is what made the MIRROR detector harmless.
- **Exact wrong-strand intron-chain duplicate.** Zero cases on any corpus: the
  artifact is not the same intron relabelled but a cryptic canonical splice a
  few bp away on the other strand (median shift 4 bp, max 15 bp on Sequin `s`
  models; 53 of 56 sites within 10 bp). Hence the 20 bp tolerance.

The Sequin ONT spike-ins (no true antisense) then supplied a labelled artifact
set. Read-level checks showed the `x` run-on reads are the resident gene's own
sense-oriented, uniquely mapped (MAPQ 60), spliced reads that minimap2 spliced on
the opposite strand at a novel cryptic canonical site: not contamination, not
genomic DNA, not ligation chimeras (about 0.5% supplementary alignments). The
read-level cryptic opposite-strand splice rate is not specific to cDNA (it was
higher in direct RNA at matched depth), so the artifact is a general
splice-aligner plus depth effect. The mirror rule caught 25 of 25 `s` models
but only 9 of 86 `x` models, because run-ons use new splice sites rather than
the sense gene's boundaries; a blastn check showed both classes are
near-identical copies (`s` median 99.0% identity at full coverage, `x` 100% at
0.85) of an isoform with far higher expression (median ratio 0.0014 and
0.0004), whereas real isoforms with similar matches are far MORE expressed than
their match (median 203x). That became the SEQUENCE detector. Two variants were
rejected along the way: restricting it to antisense matches (on mouse it flagged
8, 7 of them real; mismapping artifacts there are mostly same-strand), and
identity 95-98% (more paralog exposure on a real genome; 99% kept the mouse
F1 gain while losing fewer paralogs).

## Real-data benchmark (v0.40.0, 2026-09-24)

LRAA v0.40.0 de novo on whole GRCh38, filtered vs the same run with
`--no_filter_mismappings`, scored by the paper benchmark:

| sample | platform | models removed | effect |
|---|---|---|---|
| BT474 | PacBio HiFi | 35 of 87,707 (0.04%) | intron-chain TP -2, FP -3; F1 unchanged. Both TP losses are same-strand near-identical (>=99.6%) copies below 1% relative expression, the sequence detector's paralog cost |
| SG-NEx MCF7 | ONT cDNA | 74 of 31,081 (0.24%) | TP 0 lost, FP -51; novel FDR 0.274 -> 0.271 |
| SG-NEx A549 | ONT cDNA | 154 (0.43%), de novo and ref-guided | about -127 FP in each regime, known TPR unchanged |
| SIRVs x12 (incl. E2), MORF PacBio + ONT | both | 0 in all 28 runs (de novo and ref-guided) | identical scores |

Nearly all ONT removals came from the SEQUENCE detector (MCF7: 72 sequence, 1
mirror, 1 both). The filter is self-gating: on clean HiFi data there are few
artifacts for it to find, so no platform switch is needed.
