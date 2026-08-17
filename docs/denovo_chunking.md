# Chunked DISCOVERY: why the guard was relaxed, and what replaced it

Proposed 2026-08-16 on branch `denovo-chunking`, branched from `devel`. This is a
proposal for review, not a merge. Everything below is either MEASURED, with the
measurement named, or marked `[INFERENCE]`.

## What changed

`LRAA --chunk` used to refuse anything but `--quant_only --gtf`:

```python
    if not args.quant_only:
        sys.exit(
            "Error: --chunk requires --quant_only. Isoform discovery cannot be "
            "split across chunk boundaries."
        )
```

It now serves three configurations:

| invocation | mode | cut rule |
|---|---|---|
| `--chunk --quant_only --gtf A` | quantification | severing is PRICED: minimise it, drop/count/name what is severed |
| `--chunk --gtf A` | ref-guided discovery | severing is FORBIDDEN: a target with no clean position is DECLINED |
| `--chunk` | de novo discovery | same as above |

`--gtf` is now required only by `--quant_only`. Quant-only behaviour is byte-for-byte
unchanged; see "Quant mode is unchanged" below, which is a measurement, not a claim.

## The theory the change rests on

A chunk boundary at position `b` can damage discovery only by cutting through a
locus. In de novo mode a locus IS a connected component of the splice graph, whose
nodes come from alignment blocks and whose edges come from `N` CIGAR ops. If no
retained alignment spans `b`, then no node and no edge crosses `b`: the graph is
already disconnected there, and since isoforms are reconstructed per connected
component, the partition cannot change what is reconstructed. If an alignment DOES
span `b`, the graph is cut through a locus and both halves are reconstructed
separately, as truncated or spurious models.

So severing is not a cost to be minimised in discovery. It is the whole failure
mode, and the correct treatment of it is refusal.

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

Independently, on chr1 at the 10 Mb default, all 24 cut targets are placeable with
zero spanning retained alignments, de novo and ref-guided alike, at a median search
radius of 25 kb against a 500 kb maximum. And on chr21 with GENCODE
(`Chr21RefGuidedParity`), annotation-constrained selection places 4 cuts at 0
severed, at coordinates that differ from the annotation-free arm's — so the two
discovery flavours reach zero by different routes, and both reach it.

Read together: the guard was well founded in KIND, and wrong in DEGREE. The
mechanism it feared is real, is local, and is caused by severing — which the
selector already measures per candidate position and can therefore refuse.

## The design decision: skip and widen, not fail the run

When no position in a target's wiggle window severs zero alignments, discovery
**declines that cut**. The two chunks it would have separated stay joined as one
larger chunk. The run continues.

The alternative was to fail the run. Rejected, for three reasons:

1. **Asymmetric consequences.** A larger chunk is SLOWER and uses more memory. A bad
   cut is WRONG — it deletes real models and invents fake ones. Trading correctness
   for speed at a cut is never the right trade; trading speed for correctness at a
   cut always is.
2. **Chunk spans are already uneven, by design.** The tail-merge rule produces a
   14.4 Mb final chunk on chr20 at the shipped defaults, and the wiggle window means
   a realised span is nominal ± up to the window width. Nothing downstream assumes
   equal chunks, so a wider one is a performance fact, not a structural surprise.
3. **Failing would make the mode unusable on the substrates that need it most.** A
   dense contig is exactly where chunking pays, and exactly where some windows are
   fully covered. A mode that refuses to run there is a mode nobody runs.

A third option — widen the window until a clean position is found — was also
rejected. The window is a promise about chunk geometry, and the selector's docstring
already states that it is never widened; a mode that silently widened it would make
chunk coordinates depend on read depth in a way no caller could predict. Declining
is the same outcome (a bigger chunk) reached by a rule that is stated in advance.

## Silent degradation was the real risk, so it is reported

A run that produced fewer chunks than its geometry implies, for good reasons, is
indistinguishable six months later from a performance regression. So every run —
quant-only as well as discovery — now prints a placement report before the expensive
phase, and stores it in `timing.json` and `outputs.json` under `cut_placement`:

```
CUT PLACEMENT (discovery mode): zero severed is REQUIRED, and a target with no clean position is DECLINED

  minigenome+    16 requested, 9 placed, 7 declined for severing, 0 otherwise unplaced, 1 tail-merged -> 10 chunk(s)
      DECLINED target 600000: DECLINED under the zero-severed requirement: none of the 21 compliant
      position(s) in the window severs zero retained primary alignments, the cheapest severs 2. The cut
      is skipped and the two chunks it would have separated stay joined as one larger chunk
      ...
  minigenome-    16 requested, 7 placed, 9 declined for severing, 0 otherwise unplaced, 1 tail-merged -> 8 chunk(s)
      ...
  TOTAL 32 cut(s) requested, 16 placed, 16 declined for severing, 0 otherwise unplaced, 2 tail-merged -> 18 chunk(s)
  16 cut(s) were DECLINED rather than placed badly. The chunks they would have separated stay joined, so
  this run is slower than its geometry suggests and that is the intended trade: a larger chunk is slower,
  a cut through a locus is wrong.
```

Each declined target carries the number of annotation-compliant positions its window
held and what the cheapest of them would have severed, so a reader can tell "the
window was crowded" from "the window was empty" without rerunning anything.

## Where the rule is enforced

Three places, deliberately, because each catches a different way of getting it wrong.

1. **`util/misc/select_contig_cut_points.py --require_zero_severed`.** Positions with
   nonzero spanning cost are STRUCK from the candidate set, not priced within it. The
   joint solver minimises unplaced targets first, so a severing position left in the
   set would be chosen whenever the window held nothing better — pricing it at
   infinity would not do.
2. **`ChunkedRun.verify_severed_accounting(..., discovery=True)`.** Enforced against
   what EXTRACTION actually dropped, before the expensive phase. Selection promised
   zero; if extraction drops anything, the two tools disagree about which alignments
   are retained, and that is a bug rather than a tight geometry.
3. **The pipeline sets the flag in exactly one place**, `stage_select_cuts`, keyed on
   `args.discovery`, and a unit test asserts the flag appears on the selector command
   line for discovery and does not for quant-only.

## What discovery chunking still does NOT guarantee

Zero severed is necessary, not sufficient, and the measurement says so directly.

`pylib/IsoformReadRescue.py` is scoped to a whole contig-strand and is position-blind:
it maps read sequences against a FASTA of ALL models on that contig-strand with
`minimap2 -a --secondary=yes -N 50`. In the unchunked chr21 run, a read whose only
primary alignment is at `chr21:6,500,980` was assigned to a paralogous model at
`chr21:43,109,302` — 36.6 Mb away. Chunking confines rescue to a chunk, that
assignment cannot happen, and the model it supported drops from 2 reads to 1 and is
not reported. That is the single default-spacing difference, and 17 of 230,864
assignment rows in that run cross a default chunk boundary (8 reads, 12 transcripts).

This is stated here rather than fixed here. It is not a boundary artifact, no cut
rule can address it, and whether position-blind rescue is desirable is a separate
question for the authors.

## Quant mode is unchanged

Verified by running `LRAA --chunk --quant_only` on `testing/single_contig` from this
branch and from its branch point (`f6019a1`, a detached worktree), at two geometries:

| geometry | severed | merged `quant.expr` | merged `quant.tracking` | cut positions | dropped read names | drop detail | severed bam |
|---|---|---|---|---|---|---|---|
| 0.5 Mb / 0.2 Mb wiggle | 0 | byte-identical | md5 identical | identical | identical | identical | 0 = 0 |
| 0.2 Mb / 0.002 Mb wiggle | 2 | byte-identical | md5 identical | identical | identical | identical | 2 = 2 |

The second row is the one that matters: it exercises the cost-minimising selection
and the drop-count-and-name path, and both produce exactly what they produced before.
A key-by-key comparison of the cut manifests, flattened, reports zero changed
values and zero removed keys in both orientations. Five keys are ADDED, and all
five state the contract rather than change it: `params.require_zero_severed` and
`counts.targets_declined_zero_severed` in the manifest, `declined_zero_severed`
and `best_spanning_in_window` on each unplaced target, and
`severed_read_accounting.zero_severed_required` in `timing.json`. All are false or
null on a quant-only run.

## Smoke evidence for the new mode

`testing/single_contig`, `minigenome`, 3,421,379 bp, 2,507 alignments, `--HiFi`,
`--cpu_budget 2`, strand-first:

* **0.5 Mb / 0.2 Mb wiggle.** 12 cuts requested, 12 placed, 0 declined, 0 severed, 14
  chunks (7 per contig-strand). Exit 0. Wrote `dn.LRAA.ref-free.gtf`,
  `.quant.expr`, `.quant.tracking.gz`. Merged GTF: 49 transcripts, 49 distinct
  transcript ids.
* **0.2 Mb / 0.002 Mb wiggle.** 32 requested, 16 placed, **16 declined for severing**,
  0 severed, 18 chunks. Exit 0. Under the soft rule the same geometry severs reads;
  under the hard rule it declines instead and reports every decline.

The id namespacing is load-bearing rather than decorative: on the first run, stripping
the per-unit prefix collapses those 49 transcript ids to 46, i.e. three real
collisions on a 3.4 Mb contig. LRAA names a model after its contig and a per-run
component index, and every chunk's mini contig carries the same contig name, so every
chunk emits a `comp-1`. Unpatched concatenation fuses unrelated models; that is what
produced 37 spurious chromosome-crossing "models" in the chr21 work before it was
diagnosed.

## How often does the refusal actually bind?

A constraint that fires constantly is a footgun, so this was measured rather than
assumed. `SpanIndexProbe` measured on real chromosomes that ~0.7-0.8 Mb of TOTAL
wiggle suffices to give every 2 Mb target a zero-cost position, and the shipped
default is already 1 Mb -- i.e. at sane geometries the constraint should not bind.
Note that the 940-severed stress result came from a 20 kb wiggle OVERRIDE, not from
2 Mb spacing: a wiggle-to-spacing ratio of 0.01.

Swept on `minigenome` at 0.2 Mb spacing, selector only, strandless. The right axis
is the wiggle-to-spacing RATIO, because that is what decides how many grid positions
a window offers:

| ratio | wiggle | placed / requested | declined | alignments the SOFT rule would sever |
|---|---|---|---|---|
| 0.01 | 0.002 Mb | 7 / 16 | 9 | 378 |
| 0.05 | 0.01 Mb | 8 / 16 | 8 | 331 |
| 0.10 | 0.02 Mb | 11 / 16 | 5 | 44 |
| 0.30 | 0.06 Mb | 13 / 16 | 3 | 24 |
| 0.40 | 0.08 Mb | 14 / 16 | 2 | 16 |
| 0.70 | 0.14 Mb | 15 / 16 | 1 | 1 |
| 1.00 | 0.2 Mb | 15 / 16 | 1 | 1 |

Two things to read out of it, and one thing not to.

Declines fall monotonically with the ratio and TRACK what the soft rule would have
severed: at ratio 0.01 the refusal costs 9 cuts and buys back 378 severed
alignments; at ratio 1.0 it costs 1 cut and buys back 1. The constraint is
self-limiting -- it binds hardest exactly where the soft rule does most damage, and
becomes nearly free where the soft rule is nearly harmless.

The decline demonstration in the smoke section is at ratio 0.01, which is the SAME
ratio as the measured chr21 damage case (2 Mb spacing, 0.02 Mb wiggle). It is the
scaled analogue of a geometry known to fragment loci, not an arbitrary stress.

What NOT to read out of it: that the constraint still binds once at ratio 1.0 on this
corpus is not evidence against the chromosome-scale finding. `minigenome` is a
synthetic 3.4 Mb test contig carrying 2,507 alignments packed far more densely per Mb
than a real chromosome, and its 16 targets at 0.2 Mb spacing are 50x finer than the
shipped 10 Mb default. At the smoke section's own sane geometry -- 0.5 Mb spacing with
0.2 Mb wiggle, ratio 0.4 -- it placed 12 of 12 with nothing declined. `[INFERENCE]`
At the shipped 10 Mb / 1 Mb defaults on a real chromosome the refusal should be rare
to absent, which is what chr21 (0 severed at 4 cuts) and chr1 (24 of 24 targets
placeable at 0 severed) both show; that has not been rerun here.

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
