# Chunked DISCOVERY: why the guard was relaxed, and what replaced it

Proposed 2026-08-16 on branch `denovo-chunking`, branched from `devel`. Revised
2026-08-17 after review: the first revision made zero severed a HARD constraint in
discovery, and that was REJECTED. This is a proposal for review, not a merge.
Everything below is either MEASURED, with the measurement named, or marked
`[INFERENCE]`.

## What changed

`LRAA --chunk` used to refuse anything but `--quant_only --gtf`:

```python
    if not args.quant_only:
        sys.exit(
            "Error: --chunk requires --quant_only. Isoform discovery cannot be "
            "split across chunk boundaries."
        )
```

It now serves three configurations, and cut selection is IDENTICAL in all three:

| invocation | mode | cut rule |
|---|---|---|
| `--chunk --quant_only --gtf A` | quantification | severing is PRICED: minimise it, drop/count/name what is severed |
| `--chunk --gtf A` | ref-guided discovery | same rule |
| `--chunk` | de novo discovery | same rule |

The price is weighted: a severed monoexonic alignment costs 1, a severed multi-exon
one costs `--chunk_severed_multiexon_weight` (default 10). The search around each
target widens progressively up to the wiggle window and stops at the first position
that severs nothing; on reaching the maximum it takes the best position available.
The ANNOTATION remains a hard block and is the only thing that can decline a cut.

`--gtf` is now required only by `--quant_only`. Quant-only behaviour is unchanged;
see "Quant mode is unchanged" below, which is a measurement, not a claim.

## The rule this proposal first shipped, and why it was rejected

The first revision promoted severing from a cost to a REFUSAL in discovery:
`--require_zero_severed` struck every severing position from the candidate set, and
a target whose window held none was DECLINED. The reasoning was that a severed read
costs quantification one read but can cost discovery a whole locus, which is true
and is still the reason severing is priced at all.

The USER rejected it, and the reason is about where the data is going rather than
about this corpus:

> I don't think zero-severed should ever be a hard constraint. It's a preferred
> situation but it will be impossible to achieve as data sets get deeper. At some
> point, every base will be covered by a read. We will find reasonable cut points
> though and collateral severed reads is an acceptable cost.

That is decisive, and the failure mode it names is the worst kind. A rule that
forbids severing does not degrade into "fewer cuts"; on a contig where every base is
covered it declines EVERY target, so chunking silently switches itself off on
precisely the inputs it exists for, while the run reports success. Nor can it be
tuned back on by widening the window: the problem is not a window that is too small,
it is a genome with no clean position anywhere in it.

Measured on the shipped corpus, at a geometry tight enough to force the choice --
`testing/single_contig`, de novo, strandless, 0.2 Mb spacing, 0.02 Mb window, 2,266
retained primary alignments:

| contract | cuts placed / 16 | declined | chunks | alignments severed |
|---|---|---|---|---|
| REJECTED, `--require_zero_severed` | 11 | 5, for severing | 12 | 0 |
| this revision, severing priced | **16** | 0 | **17** | 44 (0 monoexonic, 44 multi-exon) |

The hard rule gave up 5 of 16 cuts -- and 5 chunks of parallelism -- to avoid
severing 44 alignments out of 2,266, or 1.94%. Reproduced by running the selector
from `483aaa4` (the last commit carrying the hard rule) beside the current one on the
identical bam.

## The theory, and what it does and does not license

A chunk boundary at position `b` can damage discovery only by cutting through a
locus. In de novo mode a locus IS a connected component of the splice graph, whose
nodes come from alignment blocks and whose edges come from `N` CIGAR ops. If no
retained alignment spans `b`, then no node and no edge crosses `b`: the graph is
already disconnected there, and since isoforms are reconstructed per connected
component, the partition cannot change what is reconstructed. If an alignment DOES
span `b`, the graph is cut through a locus and both halves are reconstructed
separately, as truncated or spurious models.

What that licenses is the SHAPE of the cost -- severing is scored on an alignment's
whole REFERENCE span, introns included, and a spliced alignment costs more than a
monoexonic one because the `N` op is an edge. What it does NOT license is a veto,
for the reason above. The theory says which positions are better; it cannot say that
the best available position is unacceptable, because something has to be chosen.

### Why the cost is the REFERENCE span, and what a zero costs buys

Worth spelling out, because the obvious cheap cost function misses a whole damage
class. `Chr21DenovoParity` revised their own diagnostic after finding that
`comp-2240:iso-10` was damaged while sitting 87 kb from the nearest cut: the cut at
40,005,600 fell inside its COMPONENT's span (39,953,924-40,137,010), and cutting a
component repartitions models that never come near the boundary. Their conclusion,
correctly, is that a gate on "no model within X bp of a cut differs" misses this, and
the right granularity is the component.

The spanning cost is at that granularity by construction. It counts an alignment as
severed when its full REFERENCE span contains the position — `aln.reference_start` to
`aln.reference_end`, introns included, not its aligned blocks. So a position scoring
ZERO has every retained alignment wholly left or wholly right of it. Nodes are blocks
and edges are `N` ops, and both lie inside an alignment's reference span, so nothing
on the left shares a node or an edge with anything on the right, and connectivity —
which is transitive through shared nodes and edges — cannot cross `b`. No component
can straddle a position no alignment spans.

That is what a zero-cost cut buys, and it is why the search hunts for one rather than
settling near the target. What it is NOT is a promise the selector can keep on
demand: whether a zero exists is a fact about the data.

Measured independently by `Chr21DenovoParity` over the 535 components of the
unchunked chr21 de novo run, using LRAA's `gene_id` as the component id, so nothing
is inferred:

| partition | alignments severed | components straddling a cut |
|---|---|---|
| 4 default cuts | 0 | **0 of 535** |
| 22 stress cuts | 940 | **3**, at 26,000,200 / 40,005,600 / 41,992,400 |

Those three cuts sever 861 / 30 / 47 alignments. The two stress cuts that sever
exactly ONE alignment each straddle NO component — so zero spanning alignments is
SUFFICIENT for zero straddling components without being tight: a single severed
alignment did not happen to connect two loci at either cut.

Read correctly, that table is an argument about GEOMETRY, not about a gate. At the
shipped 10 Mb spacing the cuts sever nothing and no component straddles one, so the
partition is provably harmless there — and the selector reaches that outcome by
searching, not by refusing. At the stress geometry, whose window was starved 40-fold
below what the data need, three cuts sever and three components straddle. The fix for
that is a window wide enough to find the gaps that exist, which the 1 Mb default
already is: at 2 Mb spacing with the full 1 Mb window, chr21 places all 22 cuts at
ZERO severed (measured here, section "Progressive expansion").

The reference-span detail is load-bearing and is pinned by
`test_a_severed_spliced_read_is_counted_by_its_reference_span`, whose fixture has no
aligned base anywhere in the window — only an intron. Fail-verified by mutating the
cost to count every alignment at weight 1: the spliced read's cost falls from 10 to 1
and that test fails, along with six others. Restored, 73 pass.

## The evidence

From `FINDINGS.chr21_denovo_parity.md` (LRAA_PAPER_Analyses `05af6d8`), measured by
a sibling slice on chr21, HG002 PacBio Kinnex, 199,153 retained primary alignments.
It is not re-derived here.

**At the shipped 10 Mb spacing the cuts sever ZERO alignments**, and chunked de novo
discovery differs from the unchunked run by ONE model out of 1462, none gained. That
one model sits 2,724,898 bp from the nearest cut, so it is not a boundary artifact.
The monoexonic bucket — the one that would carry severed-read damage — is identical,
120 against 120, every span matching.

**Forced to a 2 Mb spacing with a 0.02 Mb window the selector has to sever 940
alignments, and the damage is exactly the predicted shape.** 17 of the 20 lost
models SPAN a cut. At `chr21(-):41,986,832-42,010,348` a gene with eight isoforms
whose first exon crosses the cut at 41,992,400 loses ALL EIGHT, and the chunked arm
emits two SPURIOUS MONOEXONIC models at 634 bp and 797 bp from the cut whose left
ends coincide with the originals'. The same shape recurs at 26,000,200 and
40,005,600.

Those two figures are from different arms and different diagnostics, and conflating
them would overstate the case, so: the 20 lost models are the stress arm with rescue
ON, keyed by model span. With rescue OFF the stress arm loses 18, and under the
component-level diagnostic all 18 are attributed to a severed component rather than
17 of 18 by span. The severing damage and the rescue damage are separable, and both
counts are of severing once rescue is taken out of the picture.

Independently, on chr1 at the 10 Mb default, all 24 cut targets are placeable with
zero spanning retained alignments, de novo and ref-guided alike, at a median search
radius of 25 kb against a 500 kb maximum. And on chr21 with GENCODE
(`Chr21RefGuidedParity`), annotation-constrained selection places 4 cuts at 0
severed, at coordinates that differ from the annotation-free arm's — so the two
discovery flavours reach zero by different routes, and both reach it.

Read together: the guard was well founded in KIND, and wrong in DEGREE. The
mechanism it feared is real, is local, and is caused by severing — which the
selector measures per candidate position and can therefore SEARCH AWAY FROM. The
first revision of this note ended that sentence with "and can therefore refuse",
which is the step the user rejected.

## The design: search outward, price what is left, never refuse

### Progressive expansion, with the wiggle window as a maximum

The search starts at a small radius around the target and widens through
`DEFAULT_EXPANSION_RADII` — 5 kb, 25 kb, 100 kb, 250 kb, 500 kb — clipped to the
wiggle half-width, stopping as soon as an annotation-compliant position severs
NOTHING. On reaching the maximum it takes the cheapest position available and the
alignments it severs are collateral: counted, named, and reported per cut.

Terminating early is not an approximation, and this matters because it is what makes
the whole scheme free. If any compliant position at radius R severs nothing then the
whole-window minimum is also zero, and among zero-cost positions the tie-break is
distance from the target, so every candidate that could win already lies inside R.
Where windows overlap and the anti-sliver floor binds, a target the joint solve had
to compromise is expanded another rung and the solve repeated, which carries the
property over to the joint problem. Verified on chr21, strandless de novo, by running
the laddered search and a single-rung search at the maximum radius (which IS the flat
scan this module used to do) and comparing: **identical cut positions at both 10 Mb
and 2 Mb spacing.**

### Each rung fetches only the new ANNULUS

This is the part that is easy to get wrong in the direction of a regression. The
ladder's radii sum to more than the window: re-scanning `target ± radius` from
scratch at every rung costs `2*(5+25+100+250+500)` kb = 1.76 Mb per target against
the flat 1 Mb of scanning the whole window once. So in the regime where nothing
terminates early — which is the deep regime this design exists to serve — naive
expansion is 1.76x WORSE than never expanding at all. Each rung therefore reads only
the two intervals the previous rung did not cover, and carries the running best
forward.

Measured on chr21 (46.7 Mb, 199,153 retained primary alignments), strandless de novo,
1 Mb maximum window, counting bases requested from `pysam.fetch`:

| spacing | annulus | naive re-scan | flat, no expansion |
|---|---|---|---|
| 10 Mb, 4 cuts | **1.03 Mb** in 12 intervals | 1.79 Mb in 8 | 4.00 Mb in 4 |
| 2 Mb, 22 cuts | **2.55 Mb** in 50 intervals | 4.17 Mb in 36 | 22.00 Mb in 22 |

All three choose the same positions. Annulus fetching is 1.74x cheaper than naive
re-scanning at 10 Mb and 1.64x at 2 Mb, and 3.9x / 8.6x cheaper than not expanding.

One honest qualification, because the two numbers above are easy to conflate. On THIS
corpus naive re-scanning is still cheaper than flat (1.79 against 4.00 Mb), because
early termination is nearly universal here — 3 of the 4 chr21 targets stop at the
first 5 kb rung. The 1.76x-worse-than-flat figure is the WORST CASE, reached when
nothing terminates early. Both are true; the annulus form is the one that is never
worse than flat, which is why it is what is built.

### The top of the ladder is load-bearing on real data

The terminal search radii on chr21 at 10 Mb spacing are **5 kb, 5 kb, 5 kb and
500 kb**. Three targets find a clean cut almost on top of themselves; the fourth has
to travel to the maximum, and lands at 40,382,800 — the +382.8 kb offset
`FINDINGS.index_span_probe.md` reports independently. At 2 Mb spacing, 15 of 22
targets stop at 5 kb, four at 25 kb, one at 100 kb and two at 500 kb.

That is a second, independent argument for the wiggle window staying ABSOLUTE: one
chr21 target in four needs a search radius of more than 250 kb, and a window derived
as 10% of a 2 Mb spacing would offer a 100 kb half-width and could not reach it.

### The severing cost is weighted by exon structure

A severed monoexonic alignment costs 1; a severed multi-exon one costs
`--chunk_severed_multiexon_weight`, default 10. The asymmetry is the same mechanism
as the reference-span rule: an `N` op is an edge of the splice graph, so severing a
spliced read removes structure while severing a monoexonic read removes depth.

Recorded honestly, because a later reader running a sweep will otherwise conclude the
feature is broken: **on real data the weight is nearly inert.** Of the alignments
severed at 2 Mb spacing, monoexonic against spliced — chr21 20 kb window: 14 / 926;
chr21 200 kb: 0 / 743; chr1 20 kb: 5 / 2593; chr1 200 kb: 2 / 85. 98-100% spliced at
every geometry measured, because spanning probability scales with genomic span and a
monoexonic read is a kb or two long against tens of kb for a spliced one. It shows up
in the shipped corpus too: at 0.5 Mb spacing on `minigenome` every one of the 5
severed alignments is multi-exon, so `best_spanning_ignoring_annotation` goes from 1
to 10 and from 2 to 20 when the weight is applied.

It is not inert everywhere, and where it acts is the interesting part. At 2 Mb spacing
with a 0.02 Mb window on chr21, K=10 against K=1 moves exactly ONE of the 22 cuts —
and it is **41,992,400 → 41,990,700**, the cut this whole design is motivated by, the
one that replaced an eight-isoform gene with two spurious monoexonic models:

| | position | severed | monoexonic | multi-exon | weighted |
|---|---|---|---|---|---|
| K=1 | 41,992,400 | 47 | 12 | 35 | 47 |
| K=10 | 41,990,700 | 49 | 25 | 24 | 265 |

K=10 gives up 13 more monoexonic severings to save 11 spliced ones. Whether that
particular trade rescues that particular gene is not measured here — it needs a
discovery arm, and the geometry is a stress arm nobody should run — but the direction
is the one the cost function is meant to express, and it fires at the worst cut in
the corpus rather than at an irrelevant one.

**What NOT to do with the inertness.** It is not a licence to count multi-exon
alignments only and let monoexonic severing be free. At depth the mix is exactly what
is in motion, and a free class of severing is an unbounded cost waiting for the data
that produces it.

### The annotation is still a hard block

Unchanged, and the only reason a target can now be DECLINED. A read can be dropped
and accounted for; a locus cannot — `genes_contained` emits a gene whole or not at
all, so a locus straddling a boundary is contained by neither neighbour and both omit
it. A target whose whole window is annotation-blocked is declined, the two chunks it
would have separated stay joined, and the run continues. That is accepted behaviour,
not a defect to be fixed by relaxing the block.

On `testing/single_contig` at 0.2 Mb spacing with a 0.002 Mb window that fires hard:
14 of 16 targets per orientation are declined, every one of them for the annotation,
and the two that place sever 2 alignments between them. `devel` leaves exactly the
same targets unplaced at that geometry, with the same reason and the same wording
intent — only the message text differs.

## Silent degradation was the real risk, so it is reported

Two things a run must not do quietly: produce fewer chunks than its geometry implies,
and sever reads. The second is now the EXPECTED outcome rather than a refused one, so
the placement report is the only place anyone sees what the boundaries cost. Every
run — quant-only as well as discovery — prints it before the expensive phase and
stores it in `timing.json` and `outputs.json` under `cut_placement`. Actual output
from the de novo smoke run below:

```
CUT PLACEMENT (discovery mode): severing is a COST, never a veto -- a severed read is dropped,
counted and named, and a multi-exon one costs more than a monoexonic one. Only the annotation
can decline a target.

  minigenome+    6 requested, 6 placed, 0 declined for annotation, 0 otherwise unplaced,
                 0 tail-merged -> 7 chunk(s); 1 alignment(s) severed (0 monoexonic, 1 multi-exon)
      cut at 520300 (target 500000, +20300) severs 1 alignment(s): 0 monoexonic, 1 multi-exon;
      searched 50000 bp
  minigenome-    6 requested, 6 placed, 0 declined for annotation, 0 otherwise unplaced,
                 0 tail-merged -> 7 chunk(s); 0 alignment(s) severed (0 monoexonic, 0 multi-exon)

  TOTAL 12 cut(s) requested, 12 placed, 0 declined for annotation, 0 otherwise unplaced,
  0 tail-merged -> 14 chunk(s); 1 alignment(s) severed (0 monoexonic, 1 multi-exon)
  1 alignment(s) were severed and are dropped, counted and named
  (<cuts>.dropped_reads.tsv, <cuts>.severed_reads.bam). 1 of them carried junctions,
  which is the part that can cost a model.
```

Only cuts that cost something get a line, so the report stays readable on a
chromosome; a clean cut is visible in the per-contig-strand totals. Each declined
target carries how many grid positions its window held and the radius the search
reached, so a reader can tell "the annotation blocked everything" from "the window was
narrow" without rerunning anything.

## Where the cost is computed and reported

1. **`util/misc/select_contig_cut_points.py`.** One objective, `spanning_counts`,
   weighted by `--severed_multiexon_weight`, minimised jointly by the existing DP.
   Nothing is struck from the candidate set. The per-cut monoexonic/multi-exon split
   is counted in the same pass that gathers the dropped read names, so it costs no
   extra I/O and cannot disagree with the names.
2. **`ChunkedRun.verify_severed_accounting`.** Unchanged in kind: the reads selection
   NAMED must be the reads extraction DROPPED. The first revision additionally raised
   whenever discovery dropped anything; that raise is gone, because a nonzero drop is
   now the expected outcome. The identity check remains, because the parity
   comparison's pruned baseline depends on it.
3. **`ChunkedRun.cut_placement_report`.** One function, both modes, printed and
   stored. There is no mode-conditional cut rule left to test, so the unit test
   asserts the opposite of what it used to: that `--require_zero_severed` appears on
   NEITHER command line and that both modes pass the same weight.

## What discovery chunking does NOT guarantee

A zero-cost cut is provably harmless (section "Why the cost is the REFERENCE span"),
and at the shipped geometry every cut is zero-cost on both measured chromosomes. What
is NOT guaranteed is that a zero-cost cut exists — that is a fact about the data, and
where none does the selector takes the cheapest position and says so. On a deep enough
library there will be no clean position anywhere, every cut will sever, and the
reported counts are the only honest account of what that cost. This is the trade the
user accepted explicitly, and it replaces a rule that would have declined every cut
instead.

Separately, and independent of severing:

`pylib/IsoformReadRescue.py` is scoped to a whole contig-strand and is position-blind:
it maps read sequences against a FASTA of ALL models on that contig-strand with
`minimap2 -a --secondary=yes -N 50`. In the unchunked chr21 run, a read whose only
primary alignment is at `chr21:6,500,980` was assigned to a paralogous model at
`chr21:43,109,302` — 36.6 Mb away. Chunking confines rescue to a chunk, that
assignment cannot happen, and the model it supported drops from 2 reads to 1 and is
not reported. That is the single default-spacing difference, and 17 of 230,864
assignment rows in that run cross a default chunk boundary (8 reads, 12 transcripts).

**Measured since, by `Chr21DenovoParity`: with transcript-space rescue ablated in
BOTH arms, the 5-chunk de novo arm reproduces the whole-contig arm EXACTLY** — 1461
against 1461, zero differences in intron chains, in exact exon structures and in
monoexonic spans. So the entire default-spacing divergence was rescue scope, and at
the shipped 10 Mb spacing on chr21 there is no residual partition effect at all once
it is removed. Note also that the config key governs the de novo path despite its CLI
help saying quant-only (`pylib/LRAA.py` reads it at 282, 1233, 1416, 1470), which is
its own small documentation defect.

Two consequences, and neither is implemented here.

1. Exact reproducibility is REACHABLE, not merely close: disable rescue inside chunks,
   or better, scope it to the whole contig by running one rescue pass over the MERGED
   model set after stage 6. The second keeps the feature and removes the scope
   dependence, so it is the one I would build; it needs a design decision about where
   a post-merge pass belongs, which is why it is not in this proposal.
2. The ablation is a ready-made acceptance test for whoever does it: rescue on, expect
   exactly one divergent model on chr21 at 10 Mb spacing; rescue off, expect zero.

What remains genuinely out of reach of any CUT rule is the choice itself — whether
position-blind rescue, which can rest a model's support on a locus 36.6 Mb from the
read's own primary alignment, is desirable. That is a question for the authors, and
this slice does not answer it.

## Quant mode is unchanged

Verified by running `LRAA --chunk --quant_only` on `testing/single_contig` from this
branch and from `devel` (`fdd6da7`, the main worktree), same corpus, same geometry,
at two geometries:

| geometry | severed | merged `quant.expr` | merged `quant.tracking` | cut positions | dropped read names | drop detail | severed bam |
|---|---|---|---|---|---|---|---|
| 0.5 Mb / 0.1 Mb wiggle | 5 | md5 identical | md5 identical | identical, both orientations | identical | identical | 1 = 1, 3 = 3 |
| 0.2 Mb / 0.002 Mb wiggle | 2 | md5 identical | md5 identical | identical, both orientations | identical | identical | 2 = 2, 0 = 0 |

The second row is the one that exercises the annotation-decline path and the first the
cost-minimising selection with real severing; both produce the same cuts and the same
quantification as `devel`.

On a real chromosome as well as the synthetic one. Selector-only, chr21 strandless de
novo, 1 Mb window: **4 cuts at 10 Mb spacing and 22 at 2 Mb, positions identical to
`devel`'s**, including the 40,382,800 offset the published findings report.

A key-by-key comparison of the flattened cut manifests against `devel` reports **zero
REMOVED keys** and two kinds of changed value:

* `cuts[].best_spanning_ignoring_annotation`, from 1 to 10 and from 2 to 20. This is
  the same alignments priced with the multi-exon weight — a reported diagnostic, not a
  coordinate — and the exact 10x confirms every severed alignment on this corpus is
  spliced.
* `unplaced_targets[].reason`, reworded to name the annotation as the decline reason.

Everything else is ADDED: `params.severed_multiexon_weight`,
`counts.targets_declined_annotation`, `counts.alignments_dropped_monoexonic`,
`counts.alignments_dropped_multiexon`, `counts.severed_weighted_cost_at_cuts`, and per
cut `severed_monoexonic`, `severed_multiexon`, `severed_weighted_cost`,
`search_radius`; per unplaced target `declined_annotation`, `best_spanning_in_window`,
`search_radius`. `severed_read_accounting.zero_severed_required`, which the first
revision added to `timing.json`, is GONE.

**Where design and measurement pull in different directions, stated rather than
resolved silently.** The weighted cost is one objective serving both modes, so in
principle it can reorder two DIRTY positions in a quant run and move a quant cut. It
cannot move a cut whenever some position in the window severs nothing, because
weighted zero and unweighted zero are the same set — and that is every target measured
at shipped geometry on chr1 and chr21. It did not move one anywhere on this corpus, at
either geometry, in either orientation. The one place it does move a cut is chr21 at a
deliberately starved 0.02 Mb window (above), which is a stress geometry, not a quant
default. `--severed_multiexon_weight 1` reproduces `devel`'s ranking exactly if that
is ever wanted. The alternative — a weight that applies to discovery only — was
rejected: it would be two objectives that have to agree on cut coordinates, which is
two chances to disagree, and quant's own manifest already prices a severed read at
what it costs.

## Smoke evidence for the new mode

`testing/single_contig`, `minigenome`, 3,421,379 bp, 2,507 alignments, `--HiFi`,
`--cpu_budget 4`, strand-first, de novo (`--chunk` with no `--gtf`), 0.5 Mb spacing /
0.1 Mb window:

**Exit 0. 12 cuts requested, 12 placed, 6 per contig-strand, 0 declined, 14 chunks.**
Six cuts on each of `minigenome+` and `minigenome-`, so the partition is genuinely
multi-cut per contig-strand rather than the vacuous one-chunk-per-orientation. Wrote
`dn.LRAA.ref-free.gtf`, `.quant.expr`, `.quant.tracking.gz`; merged GTF has 49
transcripts and 49 distinct transcript ids.

One cut, `minigenome+` at 520,300, severs one MULTI-EXON alignment after searching to
its 50 kb maximum radius. That is the new contract exercised end to end: under the
rejected rule that target would have been declined and the run would have produced 13
chunks with a silently doubled one. The severed read is in `dropped_reads.tsv` and in
`severed_reads.bam`, and the placement report names it above.

Terminal search radii across the 12 cuts: 5 kb (4 cuts), 25 kb (4), 50 kb = the
maximum (4). Median 25 kb, which is the figure the chr1 measurement predicted.

The id namespacing is load-bearing rather than decorative: on the first run, stripping
the per-unit prefix collapses those 49 transcript ids to 46, i.e. three real
collisions on a 3.4 Mb contig. LRAA names a model after its contig and a per-run
component index, and every chunk's mini contig carries the same contig name, so every
chunk emits a `comp-1`. Unpatched concatenation fuses unrelated models; that is what
produced 37 spurious chromosome-crossing "models" in the chr21 work before it was
diagnosed.

## How much does the window width actually buy?

Measured while the hard rule was still in place, which makes the table more useful
than it looks: the "declined" column is what the REJECTED contract would have thrown
away, and the last column is what this revision severs instead. Swept on `minigenome`
at 0.2 Mb spacing, selector only, strandless:

| wiggle | placed / requested under the REJECTED rule | it declined | alignments THIS revision severs |
|---|---|---|---|
| 0.002 Mb | 7 / 16 | 9 | 378 |
| 0.01 Mb | 8 / 16 | 8 | 331 |
| 0.02 Mb | 11 / 16 | 5 | 44 |
| 0.06 Mb | 13 / 16 | 3 | 24 |
| 0.08 Mb | 14 / 16 | 2 | 16 |
| 0.14 Mb | 15 / 16 | 1 | 1 |
| 0.2 Mb | 15 / 16 | 1 | 1 |

Read as a statement about the WINDOW, which is what it is, this is the strongest
argument in the note for widening rather than refusing. Every row that declines a cut
is a row where a wider window would have found a clean one: going from 0.002 Mb to
0.2 Mb takes the severing cost from 378 alignments to 1 and the declines from 9 to 1,
on the same reads, with no change to the rule. The damage is a property of the search
radius, not of the contract.

It also disposes of the argument the first revision made from this table — that the
refusal is safe because its cost concentrates where the alternative is worst. True as
far as it goes, and irrelevant to the case the user named: at depth the last column
never reaches zero at any width, so the "declined" column never reaches zero either,
and a rule whose cost is concentrated where the alternative is worst becomes a rule
that declines everything.

### The requirement is ABSOLUTE, not proportional to spacing

The table above holds spacing FIXED, so wiggle in bases and the wiggle-to-spacing
ratio move together and cannot be told apart in it. An earlier draft of this note
read it as the ratio being operative. **That reading was wrong**, and it matters,
because it is the premise of a change that would reintroduce the exact damage this
mode exists to prevent.

Deconfounded by varying spacing at FIXED absolute wiggle, same corpus:

| spacing | wiggle | ratio | targets | declined | decline fraction |
|---|---|---|---|---|---|
| 0.2 Mb | 0.04 Mb | 0.200 | 16 | 5 | 0.312 |
| 0.4 Mb | 0.04 Mb | 0.100 | 8 | 2 | 0.250 |
| 0.8 Mb | 0.04 Mb | 0.050 | 3 | 1 | 0.333 |
| 0.2 Mb | 0.1 Mb | 0.500 | 16 | 2 | 0.125 |
| 0.4 Mb | 0.1 Mb | 0.250 | 8 | 1 | 0.125 |
| 0.8 Mb | 0.1 Mb | 0.125 | 3 | 1 | 0.333 |

Hold absolute wiggle at 0.04 Mb and the decline fraction stays near 0.3 while the
ratio falls 4x. Hold it at 0.1 Mb and the fraction stays near 0.125 over the same 4x
ratio range. Compare across the two blocks instead and the fraction halves when
absolute wiggle rises — 0.312 to 0.125 at spacing 0.2, 0.250 to 0.125 at spacing 0.4.
The ratio predicts nothing: ratio 0.200 gives 0.312 while the LARGER ratio 0.250 gives
0.125, and ratio 0.125 gives a higher decline fraction than the smaller ratio 0.100.
(The 0.8 Mb rows carry only 3 targets; the 16- and 8-target rows are the load-bearing
ones.)

So `SpanIndexProbe` is right and this corpus agrees with it: the distance you may have
to travel to reach a read-free gap is a property of the sequence and the library, in
BASES. Measured at chromosome scale, max nearest-zero distance 431.4 kb on chr1 and
382.8 kb on chr21, which is why ~0.7-0.8 Mb of total wiggle suffices at 2 Mb spacing
and why the shipped 1 Mb default already covers it. Spacing decides how many windows
there are, not how likely each one is to succeed.

Two points about the FORM of that evidence, because they decide what may be built on
it.

First, the inversions are what settle it, not the flat rows. Decline fraction being
flat under a 4x change in ratio is consistent with absolute; ratio 0.200 declining MORE
than the larger ratio 0.250, and ratio 0.125 more than the smaller 0.100, is
INCONSISTENT WITH RATIO. A monotone relationship cannot survive an inversion. Reporting
it as "inconsistent with ratio" rather than "consistent with absolute" is the stronger
and the honest form.

Second, and this is the part a tidy story gets wrong: ratio is not operative on this
corpus EITHER. The tempting reconciliation — that `minigenome`'s gaps scale with its
own length so ratio governs synthetic corpora while absolute distance governs real
chromosomes — reaches the right conclusion by a false mechanism, and a false mechanism
that agrees with the answer is the most durable kind. It is false here: absolute wiggle
predicts the decline fraction ON MINIGENOME, across a 4x range of ratios, in both
blocks of the table above. Ratio was never operative anywhere; it was a confound in a
sweep that held spacing fixed. Anyone who writes down "ratio governs small or synthetic
contigs" has licensed the proportional-wiggle rule the next section exists to prevent.

The mechanism is unremarkable once stated: a window of width W offers W/depth_window
candidate positions, and whether any of them is read-free depends on the local
distribution of read-free gaps — a property of the sequence and the library, measured
in bases. Spacing changes how many windows a contig has. It does not change what is
inside any one of them.

### The misreading to prevent

**Do not make the default wiggle proportional to spacing.** A "keep the ratio at 0.1"
rule looks harmless and is calibrated by the 10 Mb / 1 Mb default, but at 2 Mb spacing
it yields a 0.2 Mb window — below the 382.8 kb that chr21 actually needs. That is not
a hypothetical, and it has now been measured directly at 2 Mb spacing on HG002 PacBio
Kinnex, de novo, counting the retained primary alignments the cuts would sever:

| contig | 20 kb window | 200 kb window | 1 Mb window |
|---|---|---|---|
| chr21 | 940 | **743** | **0** |
| chr1 | 2598 | 87 | **0** |

The 200 kb column IS the proportional reading at 2 Mb spacing, and on chr21 it still
severs 743 alignments where the absolute 1 Mb default severs NONE. The formula does
not give less margin; it gives a window that fails. And chr1 and chr21 disagree
8.5-fold at identical parameters — 87 against 743 — because chr1 offers 6,258 read-free
gap runs to chr21's 1,295. The requirement tracks the GENOME, not the geometry.

The progressive search adds a third, independent confirmation: at 10 Mb spacing on
chr21 the four terminal search radii are 5 kb, 5 kb, 5 kb and **500 kb**. One target
in four needs more than half the default window, and a 200 kb window would give it a
100 kb half-width.

**No proportional term was added, not even as a floor.** An earlier draft of this
revision implemented `max(1 Mb, 0.10 * spacing)`, on the reading that a floor binds
only at coarse spacing and so is harmless. That was withdrawn: the 10% figure was
about what to pass when TESTING at 2 Mb spacing, not a rule to encode, and a
production default of 1 Mb absolute already exists and is correct. Coupling the two
parameters at all — even in the direction that only ever enlarges the window — invents
a relationship the data does not support and gives the next reader a formula to
generalise. `approx_MB_per_cut_wiggle_window` is an absolute default that a caller may
override, and `pylib/test_select_contig_cut_points.py` asserts it stays 1 Mb across a
500x range of spacings so that anyone who couples them fails a test.

A caller testing a finer spacing should pass a smaller window and expect severing.
That run still exercises cut placement as it exists — the selector still searches for
the best position it can reach — and the reads it severs are reported. What it is not
is a parity arm.

Under an absolute rule the smoke section's decline demonstration is a starved window
in absolute terms — 0.002 Mb on a corpus whose gaps need tens of kb — which is the
same failure the chr21 stress case had at 0.02 Mb against a genome needing hundreds of
kb. Both starve the search. (Those two geometries also happen to share a
wiggle-to-spacing ratio of 0.01; on the evidence above that is a coincidence of the
arithmetic, not the mechanism, and it should not be quoted as one.)

### What NOT to conclude from this corpus

That `minigenome` still severs at 0.2 Mb wiggle is not evidence against the
chromosome-scale result. `minigenome` is a synthetic 3.4 Mb contig carrying 2,507
alignments packed far more densely per Mb than a real chromosome, so its read-free gaps
are scarcer in ABSOLUTE terms and a window of any given width in bases is likelier to
be crowded; and 16 targets at 0.2 Mb spacing is 50x finer than the shipped 10 Mb
default. At the smoke section's own sane geometry — 0.5 Mb spacing, 0.1 Mb window — it
places 12 of 12 and severs 5 alignments out of 2,507.

## Why only the strandless arm needs the raw bam indexed

Stated because the next person to add a chunked mode will face the same question, and
because "I tested it and it was fine" does not survive a new mode. The raw bam is
fetched by coordinate only when there is no stage-1 split output to read, and
`--strandless_chunks` is precisely the flag that skips stage 1. Those are not two
conditions that happen to agree today -- they are one condition named twice, so
`ensure_bam_index`'s `if args.strandless_chunks:` guard is co-extensive with the need
by construction rather than by coincidence.

Discovery does not change that: de novo strand-first reads the stage-1 strand-split
bams for cut selection and for extraction exactly as quant strand-first does.
Verified on a deliberately unindexed copy of the corpus bam, both arms: strand-first
de novo exits 0, places 12 of 12 cuts, and leaves the directory holding `reads.bam`
and no `.bai`; strandless de novo over the same unindexed copy exits 0 and logs
`missing index for .../reads.bam, building it`. One arm proves the guard fires when
needed, the other that it correctly does not fire when it is not.

## Deliberately not done

* **No WDL wiring.** `885a39f` has just landed WDL chunk parameters and a sibling is
  validating them. Adding an unreviewed mode to that surface now would tangle two
  reviews.
* **No `--arm baseline` for discovery.** The control arm is a whole-contig
  quantification against a supplied annotation, which is not the comparison a
  discovery run wants. `run()` refuses the combination rather than producing two
  artifacts that look comparable and are not.
* **No change to the strandless path**, beyond it inheriting the same discovery
  switch as strand-first.
* **The chr21 parity science was not re-run.** It is done, committed and cited.
* **No change to rescue scoping**, although the measurement now shows that is what
  stands between this mode and exact reproducibility on chr21. Scoping
  `IsoformReadRescue` to the merged model set is a design decision about where a
  post-merge pass belongs, and folding it into a proposal about cut placement would
  put two independent decisions behind one review. The acceptance test for it already
  exists: rescue on, one divergent model; rescue off, zero.
* **The CLI help for
  `--no_rescue_unassigned_reads_via_transcriptome_alignment` was not corrected**,
  though it says quant-only while `pylib/LRAA.py` reads the key on the de novo path
  too (282, 1233, 1416, 1470). Reported rather than fixed: it is not this slice's
  file, and a sibling found it.
* **No discovery arm at the geometry where K moves a cut.** K=10 moves chr21's
  41,992,400 cut to 41,990,700 at a 0.02 Mb window, and whether that recovers the
  eight-isoform gene is the interesting question. Answering it needs a full de novo
  discovery arm at a geometry nobody should run in production, which is a science
  slice rather than a change to this one. The selector-level trade is measured and
  reported above; the model-level consequence is not.
* **No expansion-ladder tuning.** 5/25/100/250/500 kb was chosen to make the common
  case pay one 5 kb rung and the rare case still reach the maximum in five steps, and
  it is measurably adequate — median terminal radius 5 kb on chr21, 25 kb on
  `minigenome`. It is not optimised, and it does not need to be: the ladder cannot
  change which position is chosen, only how much is read to find it.
* **No `--severed_multiexon_weight` sweep for a best value.** 10 is the value the user
  asked for. It is a policy statement about what a junction is worth relative to a
  read of depth, not a parameter with an optimum, and on today's data it is nearly
  inert either way (see the mix measurement above).
