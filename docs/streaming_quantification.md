# Two-Pass Streaming Quantification

## Overview

Quantification's memory cost is dominated by holding reads, not by the annotation. The
single-pass path materializes each shard's alignments, collapses them into multipaths, runs EM,
and emits tracking from the retained read associations: 1,266 MB of alignment objects for 173k
reads on a chr20 HiFi shard, plus ~90 bytes per read name for tracking
(`pylib/StreamingQuant.py:5-9`). That is affordable at millions of reads and impossible at
billions.

Streaming splits the work in two. **Pass 1** quantifies normally against the coverage-normalized
BAM, whose read count is bounded by the depth target rather than by library size, and hands over
the abundances it settled on. **Pass 2** streams the full BAM and, for each read, recomputes its
splice-graph path, looks up the answer pass 1 already computed for that path, writes the tracking
row, and forgets the read. Resident memory in pass 2 is the lookup table plus per-isoform
accumulators — thousands of entries, independent of how many reads stream past.

This is the consumer coverage normalization was built toward: the thinned BAM carries the
structural evidence and the abundance estimate, and the full library is then attributed against
that estimate without ever being held. See `docs/coverage_normalization.md` for the thinning
itself.

`--stream_reads` is **on by default** since v0.25.0 (`LRAA_Globals.config["stream_reads"]`,
`pylib/LRAA_Globals.py:446`). `--no_stream_reads` restores the single-pass path.

## What makes the second pass possible

Two properties measured on chr20 (`pylib/StreamingQuant.py:17-25`):

- **A read's path through the splice graph is a function of the read and the graph alone**,
  identical across runs (86,655/86,655 reads). Pass 2 can therefore recompute each read's path
  without anything carried over from pass 1 except the graph.
- **99.97% of full-BAM reads land on a path pass 1 already resolved** — provided the table also
  records the paths that matched *nothing*. Caching only successes drops that to 89%, because
  ~28% of paths have no compatible isoform and would otherwise be re-tested for every read
  carrying them.

## The assignment table

`StreamingQuant.AssignmentTable` maps a canonical splice-graph path to the per-isoform split a
completed quantification settled on (`pylib/StreamingQuant.py:48-192`). A lookup returning
`None` means the path was never seen and the caller must resolve it; a lookup returning an empty
list means it was seen and matched no isoform, which is a real answer and is cached as one.

Fractions are recomputed from the converged theta rather than reused from the ones EM returned.
The two now agree, because EM runs a final E-step at the theta it returns, so its fractions and
counts describe the abundances it hands back (`docs/methods.md`, "Quantification via EM").
Recomputing is therefore equivalent for a path pass 1 saw, and necessary for one it did not — a
resolver has no returned fraction to reuse and must derive it from theta. Keeping one expression
for both is what makes seen and resolved paths comparable at all.

**The denominator is the read-sharing component, not the gene.** That is the unit EM runs on and
therefore the unit theta is normalized within. A path compatible with isoforms of two genes in
one component gets ONE split across all of them; per-gene groups would each sum to 1 and count
the read twice (`pylib/StreamingQuant.py:91-94`).

With `run_EM` off there is no theta and the split is equal across a gene's compatible isoforms,
matching `Quantify._estimate_isoform_read_support` in that mode. The table refuses to infer which
regime it is in: `em_was_run` is passed explicitly, and theta-without-EM or EM-without-theta both
raise, because either mismatch would silently emit a table that is internally consistent and
built from the wrong abundances (`pylib/StreamingQuant.py:108-124`).

## Paths the first pass never saw

Reads surviving coverage normalization do not cover every path in the full BAM, so pass 2 meets
paths absent from the table. `StreamingQuant.make_path_resolver` answers them as pass 1 would
have: the same anchoring, the same compatibility cascade, the same 3'-end weighting, against the
same splice graph, the same transcript set and the same theta
(`pylib/StreamingQuant.py:364-445`).

Every step is a call into `Quantify` rather than a reimplementation. The cascade in particular is
ten ordered tests, and order decides which isoform a read lands on, so a second copy would not
fail loudly — it would quietly assign some reads differently.

Each resolution is cached, so the read at hand and every later read on that path are assigned
from one resolution. **No read is dropped for being unfamiliar**; completeness of the tracking
file is why the pass exists.

One exception, and it is accounted rather than hidden: a path resolved in-stream can turn out to
be compatible with isoforms from more than one independently-normalized component — thinning
removed the read that would have joined those components in pass 1, so pass 1 never saw reason to
quantify them together. `rows_for_multipath` refuses to invent a split and raises
`CrossComponentAmbiguousPath`; `assign_path` counts the read under
`totals.reads_cross_component_ambiguous` and writes no tracking row, exactly as for a path
matching no isoform. This is collateral of the same normalization that already drops severed
reads at chunk boundaries (`pylib/StreamingQuant.py:648-655`).

## Outputs and what they mean

Pass 2 writes both quantification outputs, and they are consistent with each other because both
come from the same streamed sums:

- **`LRAA.quant.tracking.gz`** — one row per (read, compatible isoform) over the **full** BAM, not
  over the reads normalization retained. This is strictly more complete than the single-pass
  default's tracking file.
- **`LRAA.quant.expr`** — counts are sums of per-read fractional assignments over the whole BAM
  (`StreamingQuant.write_expr`, `pylib/StreamingQuant.py:847-932`). TPM is normalized over the
  transcripts passed in, i.e. this contig/strand's set, matching the default path's per-shard
  behaviour, which the merge renormalizes genome-wide. `RPM_total_reads` requires
  `num_total_reads` and raises rather than emitting a zero-filled column.

**These are not the counts the single-pass path reports.** Streaming performs one expectation
step at the abundances pass 1 estimated, where the default path would go on to re-estimate them
from the full read set. How close the two land is a measured quantity, not an identity: on ONT
chr20, gene totals are preserved, while 22 of 1785 expressed transcripts differ by more than 5%
at the isoform level, holding 8% of read mass (`pylib/StreamingQuant.py:856-862`).

`transcripts` may be a superset of what the totals cover; every lookup defaults to zero, so an
isoform that never remapped onto the final-quant splice graph still gets a reported zero-count
row instead of vanishing from the file.

## The served fraction

How much of a unit pass 1 answered directly is **reported, never gated on**
(`StreamingQuant.served_read_fraction`, `pylib/StreamingQuant.py:935-957`):

```
served = (reads_streamed - reads_on_stream_resolved_paths) / reads_streamed
```

`None` when the unit streamed nothing, since a rate over zero reads is undefined and reporting
0.0 would read as "the table answered none of it".

A low served fraction costs time, not correctness — every read is assigned either way — and says
pass 1 precomputed little of *this* unit's work. The log line also reports reads per resolve
(near 1.0 means the resolution cache never paid, a different complaint) and the unseen fraction
over the *accounted* denominator, `reads_assigned + reads_unassignable`, which excludes reads with
no graph path and reads on a spacer path (`pylib/StreamingQuant.py:1016-1052`).

There is deliberately no threshold. An earlier version aborted above 25% unseen — after the
stream had already written every tracking row, so the abort could only destroy a finished, correct
result.

## Flags and required combinations

| flag | default | meaning |
|---|---|---|
| `--stream_reads` / `--no_stream_reads` | on | two-pass streaming, or the pre-v0.25.0 single-pass in-memory path |
| `--stream_reads_rescue_unassigned` | tracks transcriptome rescue (so: on) | run transcriptome rescue inside the streaming pass against a resident `mappy` index |
| `--stream_reads_rescue_unassigned_to_targets` | off | additionally offer rescue the reads that mapped to a graph path matching no target |

**Pass 1 must read a thinner BAM than pass 2**, and this is enforced, not assumed. `--stream_reads`
exits unless coverage normalization will run (`--normalize_max_cov_level > 0` without `--no_norm`)
or a distinct `--bam_for_sg` is supplied — the composition the chunked pipeline uses, normalizing
per chunk and then running with `--bam_for_sg <normalized> --no_norm`
(`LRAA:1874-1896`). If both passes resolved to one file, pass 1 would materialize every read pass
2 streams, which is the cost the mode exists to avoid.

Under stock defaults the invariant holds by normalization: `normalize_max_cov_level` is 1000
(`pylib/LRAA_Globals.py:278`) and `--no_norm` is off, so pass 1 reads the normalized BAM
(`bam_file_for_pass1 = bam_file_for_sg`, `LRAA:5702-5703`).

**Transcriptome rescue.** The streaming loop maps every record genomically only. Left to itself it
would rescue nothing, while pass 1 rescues only the reads it saw, so a rescue-enabled run would
report a rescue summary covering a fraction of the library.
`--stream_reads_rescue_unassigned` closes that by rescuing candidates inside the stream against a
resident `mappy` index over this contig-strand's models; without it, the combination is refused
(`LRAA:1897-1918`). A rescued read is not accounted specially — its verdict is a splice-graph path
that goes through the same table lookup, resolution and tracking write as a genomic one — except
that rescue resolutions are kept out of the served/unseen accounting, whose denominator excludes
rescued reads and would otherwise report a ratio above 1
(`pylib/StreamingQuant.py:685-699`).

`--stream_reads_rescue_unassigned_to_targets` is off by default because it is the one candidate
category where the two paths cannot target the same reads: the batch path derives it from its own
first pass, which reads the normalized BAM, while the stream reads the full one. Measured on ONT
chr20: batch 3,442 against streaming 11,196 with nothing batch-only, i.e. 43% more reads than the
batch path targets. An extension, not a reproduction (`LRAA:926-941`).

**Read retention is untouched.** The streaming loop's discard test branches rather than dropping
outright, and narrowly: a read rejected *only* for low percent identity is offered to rescue,
every other discard reason still discards, and the read is counted filtered either way.
`Util_funcs.quant_discard_reason` is the single retention policy shared by quantification,
coverage normalization and the chunked pipeline, so changing what it effectively retains here
would move depth measurement and `XW` sampling weights with it
(`pylib/StreamingQuant.py:657-665`).

## Where streaming runs

- **Quantification-only** (`--quant_only`): `run_quant_only`'s own call reports here, writing
  `quant.expr` from the streamed totals.
- **Discovery**: the final-quant tail streams too, but reports afterward over the full settled
  transcript set rather than the graph-filtered `input_transcripts`, so an isoform that failed to
  remap onto the final-quant splice graph still gets a zero-count row (`LRAA:5857-5893`).

Pass 2 runs whenever streaming is on, independent of whether the call reports quantities: it is
the accounting pass that writes the tracking file, not a reporting step.

## Interaction with XW read weights

By default pass 1's theta is estimated from *retained* read counts. `EM.py:87` feeds
`mp.get_read_weight()` into EM unconditionally, but that equals the read count unless
`--use_XW_read_weights_for_quant` populated the weight registry, so thinning's depth-dependent
acceptance is not divided back out of the abundance estimate. Final counts are still on the scale
of the full library, because they are streamed sums over every read; what the thinning affects is
the *split* within a read-sharing component, which comes from pass 1's theta.

`--use_XW_read_weights_for_quant` corrects that split, and streaming is what makes the flag safe
for per-read consumers: the merged tracking file covers the full BAM rather than only retained
reads, so the incompleteness marker is not emitted and single-cell matrices and `--tag_bam` are
permitted (`LRAA:1975-2056`, `LRAA:4114-4126`). See `docs/coverage_normalization.md`, "XW weights
during quantification".
