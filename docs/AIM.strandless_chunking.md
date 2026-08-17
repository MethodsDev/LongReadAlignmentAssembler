# Research aim: chunk BEFORE strand separation

Proposed 2026-08-16. Feasibility below is read out of the code and verified by
execution where marked; the performance figures are projections and say so.

## The proposal

Today the pipeline strand-separates the WHOLE bam first, then cuts each
contig-strand into chunks:

```
  strand split (serial, whole bam)  ->  cuts per contig-strand  ->  per chunk: extract | normalize | quant  ->  merge
```

Instead, cut the raw bam and do strand separation INSIDE each chunk, in
parallel with every other chunk:

```
  cuts per contig (strandless)  ->  per chunk: extract | STRAND SPLIT | normalize | quant  ->  merge
```

## Why it matters more than anything else on the list

Strand split is the largest serial phase in the pipeline and it was not on the
previous study's defect list at all. Measured on this box, A0/A1, chr1:

| phase | whole-genome bam | chr1-only bam (projected) |
|---|---|---|
| strand split | **1150.8 s** | ~120 s |
| cut selection | 19-24 s | ~10 s |
| extraction (tabix) | 129 s | 129 s |
| chunked makespan @8 | 283 s | 283 s |

It is 64% of the A0 arm's wall clock and 68% of A1's. It is also the phase that
`--contig chr1` cannot restrict, because `separate_bam_by_strand.py` has no
contig option -- so a single-contig run strand-splits all 18.1 M records to use
1.9 M of them.

Parallelising it composes with fixing serial extraction (defect 7):

| serial prep | p | Amdahl ceiling |
|---|---|---|
| today, chr1 inputs: 120 + 10 + 129 | 0.912 | 11.4x |
| + strandless chunking (split moves into the parallel phase) | 0.940 | 16.7x |
| + extraction parallelised 8-wide as well | 0.988 | 84x |

`[INFERENCE]` on every row but the first; they assume per-chunk work divides
cleanly at 8-wide, which is what the study has to measure.

## Feasibility: three things checked, all favourable

**1. Strand assignment is per-read local.** `split_bam_by_strand` decides a
read's strand from that read alone: splice dinucleotides read off the genome
sequence (`infer_spliced_orient`), falling back to annotation overlap
(`infer_transcribed_orient_via_annotation_mapping`). There is no cross-read
aggregation and no library-level statistic; the counters are reporting only.

So a read's strand is IDENTICAL computed genome-wide or inside a chunk. This is
the property the whole idea rests on, and it holds.

**2. A strandless region already parses and already keeps both strands.**
Verified by execution:

```
  parse_region("chr1+:1-1000")  -> Region(chrom='chr1', strand='+', ...)
  parse_region("chr1:1-1000")   -> Region(chrom='chr1', strand='',  ...)
  _strand_matches(forward, "")  -> True
  _strand_matches(reverse, "")  -> True
```

`_strand_matches` opens with `if not strand: return True`, so the extractor
emits both strands for a strandless region rather than needing a new mode.

**3. Island blocking already unions the strands when strand is empty.**
`find_islands` skips its strand filter on a falsy strand, so a strandless cut
selection blocks on the union of both strands' gene loci -- which is exactly the
constraint a shared cut point has to satisfy.

## What actually has to be built

Less than the above suggests. The primitives exist; the wiring does not.

1. Cut selection over a raw bam with an empty strand. Depth and spanning-
   alignment counts become strand-agnostic, which is the right objective for a
   cut that must serve both strands. `select_contig_cut_points.py` does no
   strand filtering of its own today -- it assumes a pre-split bam -- so this
   needs checking rather than assuming.
2. A per-chunk strand-split step between extract and normalize, producing two
   normalize/quant work units per chunk instead of one.
3. Merge accounting: severed reads are currently counted per contig-strand, and
   would become per chunk across both strands.

## Risks, in the order I would test them

**Fewer admissible cut positions.** Union islands are larger and more numerous
than either strand's alone, so the gaps between them are fewer. The selector
already degrades gracefully and reports it: at 10 Mb chr1 gives 24 cuts and 0
severed reads per strand, but at 2.5 Mb it reports 13 unplaced targets, 19
severed alignments and 4 constraint compromises PER STRAND. Union will hit that
wall sooner. Measure the admissible-position count under union blocking before
building anything -- it is a cheap selector-only run and it can kill the idea.

**Strand flipping must precede any strand filter.** `separate_bam_by_strand`
REWRITES `read.is_reverse` when the inferred strand disagrees with the aligner
(`num_records_strand_flipped`). The extractor's `_strand_matches` reads the raw
flag. So a chunk must be extracted strandlessly and split afterwards; extracting
with `chr1+:...` from a raw bam would silently assign every flipped read to the
wrong strand. This is the one ordering constraint that will produce a plausible
wrong answer if violated.

MEASURED, 2026-08-16. "Plausible wrong answer" is now a number. Extracting a
strand-suffixed region from the raw bam on purpose
(`StrandlessParity.py ordering-cost`) puts **56 of 245 emitted records (23%) on
the wrong strand** at `minigenome+:1014001-1502500`, and **21 of 25,757** at
`chr1+:9919601-20000000`. Every one of those reads is then quantified against the
opposite strand's transcripts, and nothing downstream revisits the flag, so the
run completes and its tables look ordinary. The 300x spread between the two
figures is the substrate, not the geometry: the minigenome's bam is not oriented
with transcription (47.5% of records flipped by the split) while HG002 Kinnex
largely is (1.2%). Neither is zero, and the larger one is what justifies guarding
this ordering in four places rather than documenting it once.

**Per-chunk genome sequence.** `split_bam_by_strand` caches one contig's
sequence via `prev_chrom`. Inside a chunk that becomes the mini-contig, which is
smaller and better -- but it is a behaviour change worth confirming rather than
assuming.

CONFIRMED, and it is the same run that confirms the invariant: the in-chunk split
reads its dinucleotides out of the mini-contig, and over 1,064,804 chr1 records
it reached the same verdict as the split reading whole chr1, including for all
12,471 records whose flag it rewrote. The mini-contig sequence is the same
substring at shifted coordinates, and the measurement says so rather than the
argument.

## How it gets validated

The parity harness already answers this exactly: `--arm both` runs the chunked
arm and a whole-contig control over the same substrate, and the strandless
variant has to reproduce the current chunked arm's `quant.expr` and
`quant.tracking`, not merely come close. Anything else is a different pipeline
wearing the same name.

The tiny corpus gate applies first and costs 23 s: `555 / 541 / 2266 /
severed 5`.

### Result, 2026-08-16

Both gates pass. `pylib/StrandlessParity.py`, driven from
`testing/strandless_parity/`; the reasoning and the figures are in
`docs/chunked_quant_parity_evaluation.md`.

- **Strand-assignment invariant, real chr1, no quant.** 1,064,804 of 1,064,804
  records agree between the whole-contig split and a split run inside the chunk
  holding the read. Under `--infer_read_orient` -- the only regime in which a
  read's strand can disagree with its alignment flag, hence the only one in which
  this could have failed -- the split rewrote 12,471 flags genome-wide and the
  same 12,471 in-chunk, with zero disagreement. The aligner-flag figure is the
  weaker of the two and should not be quoted alone: there a record's strand is a
  property of the record and cannot disagree with itself.
- **Tiny corpus, 0.5 Mb cuts.** Strand-first reproduced all four pinned values
  (`555 / 541 / 2266 / severed 5`) and strandless gave `555 / 541 / 2266` with
  **`quant.expr` byte-identical**. `quant.tracking` is 541 rows in both arms, in
  the same order, differing only in `mp_id`, a per-process MultiPath counter.
- **The severed figure does differ, legitimately: 6 strandless against 5.** The
  extra read is `m160629_191727_..._s1_p0/69114/ccs`, reverse-oriented, spanning
  928,251-1,020,755, severed by the shared union cut at 1,014,000 that the minus
  strand's own cut at 1,046,300 misses. It reaches neither arm's `quant.tracking`,
  which is why 5-vs-6 leaves the tables byte-identical. Chunk inputs differ by
  exactly that one name, 2,261 against 2,260.
- **One defect found and fixed.** The first strandless run emitted 1,110 expr rows
  against 555: both per-orientation quant units were handed the chunk's
  both-strand mini GTF, so every transcript was quantified twice, each copy
  against the other strand's reads. Stage 3b now splits the mini GTF as well as
  the mini bam. The gate keeps a regression test and refuses to key on
  `transcript_id` when it is not unique rather than silently taking one row of
  each pair.
