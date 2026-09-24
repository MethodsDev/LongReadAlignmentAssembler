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
