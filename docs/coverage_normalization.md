# Coverage Normalization for Splice-Graph Construction

## Overview

Splice-graph construction scales with read depth, and a handful of loci in any sample carry
depth far beyond what graph construction needs. LRAA thins those loci before building the
graph, leaving the rest of the genome untouched. Reported counts stay on the scale of the original
BAM: with `--no_stream_reads` abundance estimation reads that BAM directly, and under the default
two-pass streaming path the thinned BAM carries the abundance estimate while the original BAM is
streamed to assign every read (`docs/streaming_quantification.md`).

Thinning is performed by `util/normalize_bam_by_strand.py`, invoked from
`_normalize_bam_for_splice_graph` in `LRAA`. Reads are separated by strand, normalized
independently per strand, and merged, because the splice graph is itself built per contig and
strand.

The requirement that shapes the design is not how many reads survive but **which relative
quantities survive**. The splice graph decides what to keep by comparing features against each
other — a junction against the strongest junction sharing its exon island, an unspliced
interval against the intron it bridges, a TSS against the dominant TSS in its component. Those
comparisons are only meaningful if the features being compared were sampled at the same rate,
or if the rates are recorded so they can be divided back out. LRAA does the latter.

## Procedure

Two sequential passes over each strand-specific BAM, over exactly the records the consumer
will read. That set is not redefined here: the normalizer calls
`Util_funcs.quant_discard_reason`, the single retention policy quantification itself consumes,
so the two cannot drift. Unmapped, secondary, supplementary, duplicate and QC-failed records,
improper pairs, and alignments below `--min_mapping_quality` or `--min_per_id` are all
excluded; over-long introns are dropped earlier, at the strand split. Measuring depth over
records the consumer discards would raise measured coverage above anything downstream sees,
lower the acceptance probability, and inflate every surviving read's weight at that locus.

**Pass 1 — measure.** Read depth is accumulated per 100 bp window (`--depth_window`) from
aligned blocks, and junction support is counted exactly. Both come from one CIGAR-only scan.
Both structures are held for all contigs until pass 2 finishes. At 100 bp windows the depth
arrays come to roughly 250 MB for a human genome (20 MB for its largest chromosome), against
2 GB for a per-base array of the same width; junction counts add on the order of 150 MB at the
~0.65 M junctions per strand that chr20 extrapolates to. Strands are processed one after the
other, so those are peak figures rather than a sum.

Windows are offset from the contig's first aligned base, not from coordinate zero, so the
partition travels with the reads. Translate a locus and every read lands in the window it was
in before, giving the same depth estimate and the same acceptance probability. See the origin
section below for why that matters.

That default makes the boundaries a function of which records are in the input, which is right
for a whole input and wrong for a chunk of one: a chunk begins at some other read, so the same
absolute locus bins differently than it did whole. `--window_origin` replaces the implicit
anchor with the caller's grid. Its value is the absolute 0-based reference coordinate that
position 0 of the input maps to — `0` for a whole-contig BAM, the rebase offset for a chunk
extracted by `util/misc/extract_contig_region_inputs.py` (its manifest's `offset`, which is
`region.lend - 1`). Boundaries then sit at absolute multiples of `--depth_window` whatever the
input contains, so a chunk measures the same windows a whole-contig run measures. Two conditions
come with it: cuts must fall on multiples of `--depth_window`, or a window draws its bases from
two chunks; and a whole-contig control run must be given `--window_origin 0`, because on its
default it keeps the first-aligned-base anchor and no chunk grid can match it.

**Pass 2 — sample.** Each read is kept with probability `p` and, if kept, records `1/p`:

| condition | `p` | rationale |
|---|---|---|
| local depth ≤ target | 1 | coverage is already at or below what the graph needs |
| read carries a junction with support < target | 1 | preserve scarce evidence exactly |
| otherwise | `target / local_depth` | thin toward the target |

`local_depth` is the **median** window depth across the read's aligned blocks. The median keeps
a read that merely clips a narrow peak from being judged by that peak.

The uniform draw is `blake2b(seed, read_name)`, so a read's fate follows its name rather than
its coordinate or its position in the file. Runs are reproducible, and shuffling the input
cannot change the outcome — position must not decide, since positional dependence is the defect
this design removes.

The first rule is what makes this a normalization rather than a downsampling: wherever coverage
already sits below the target, every read passes through untouched.

The second rule buys exactness where it matters most. Without it, a junction with 156 reads at a
locus retained at 11% would be judged on about a dozen, and the frequency test that decides
whether to keep it becomes a coin flip; retaining those reads turns its support from an estimate
into an exact count.

**This rule trades the depth target for evidence, and the target is the side that gives.** No
single sub-target junction can be why depth is over target, but nothing stops many of them from
adding up, and alignment noise produces exactly that: a spray of distinct, weakly supported
junctions, each individually exempt. A locus whose excess depth is spread across many minor
junctions rather than concentrated in a few dominant ones will therefore thin less than one
whose isoform structure is simple. Peak retained depth across fifteen high-expression loci ran
to 1.78x the target for this reason, so the target is a target and not a bound.

The exemption is deliberate: over-retention costs graph-construction time, while the counts it
protects decide what the graph contains at all. If a locus is ever slow enough to matter,
lowering `--normalize_max_cov_level` also lowers the support at which a junction stops being
exempt, tightening both halves together.

## The XW tag

Each retained read carries **`XW:f`**, the reciprocal of its acceptance probability: how many
reads it stands for. A read kept unconditionally carries `XW:f:1.0`.

Consumers recover unbiased counts by **summing weights instead of counting reads**.
`Splice_graph._populate_exon_coverage_and_extract_introns` does this at each of its accumulation
sites: intron support, splice-site support, per-base exon coverage, and the TSS and PolyA
position counters. `Pretty_alignment.get_normalization_weight()` exposes the value, captured at
construction because `lighten()` discards the pysam record long before support is tallied.

**A read with no tag weighs 1.** Unnormalized BAMs, and alignment pickles written before the tag
existed, are therefore read exactly as they were.

**Weights compose across passes.** Thinning an already-thinned BAM multiplies: a record kept at
`p1` and then at `p2` stands for `1/(p1*p2)`, and the normalizer reads whatever weight its input
carried and multiplies rather than replacing it (`util/normalize_bam_by_strand.py:592-601`). This
is what keeps the tag's meaning — reads of the ORIGINAL library this record stands for —
independent of how many passes produced it, and it matters whenever `--bam` was thinned upstream
while normalization is left at its default. Replacing the weight instead dropped the earlier
factor outright: on the bundled corpus, chaining a loose pass into a tight one drove a record's
weight from 10.0 down to 3.87 and collapsed the weighted total to 615 against 4,970 input records.

**The three BAMs have roles, and each is checked against its own.** Weighting is not a mode with a
setting; it is a property of the data, and what makes that safe is that each input's role fixes
whether a weight can be there (`LRAA:145-222,2393-2428`):

| input | role | requirement |
|---|---|---|
| `--bam` | the full library — reported counts are scaled by it | must **not** carry `XW` |
| `--bam_for_sg` | the splice-graph evidence, taken as given | must carry `XW` |
| `--bam_for_priors` | the pass-1 abundance evidence, taken as given | must carry `XW` |

`--bam` is what `num_total_reads` counts and what the streaming pass sums assignments over, so a
thinned BAM there would report per-retained-read quantities indistinguishable from full-library
ones and emit a tracking file covering part of the library while looking complete. `--bam_for_sg`
and `--bam_for_priors` are the inputs that exist for thinned data, and an untagged one would be
read as though its survivors were the whole library. All three are refused at startup, because
none of these is detectable downstream — that is what makes each a validation error rather than a
result to be interpreted.

An unweighted priors BAM is the least visible of the three: theta would be estimated
with every read weighing 1, giving a different prior for pass-2 apportionment with no
crash and no anomaly in the output. Hence the startup refusal, which names the flag and
its role (`LRAA:2387-2395`).

One record decides each: normalization tags **every** record it retains, including those kept
whole at $p=1$ (`XW:f:1.0`), so a BAM it produced has the tag on its first aligned record and a
BAM it never touched has it on none. An empty BAM passes any of these checks — it yields an empty
graph, and the chunked pipeline produces empty per-orientation BAMs routinely.

**`--bam_for_sg` is never re-normalized.** Supplying it *is* the decision about what the splice
graph reads, so `--normalize_max_cov_level` does not apply to it. Re-thinning it would compose
acceptance rates for no gain and produce only artifacts the caller could have produced directly.

Weighting is not optional bookkeeping. Any scheme whose acceptance rate varies along the genome
distorts relative support, and every such scheme must either record its rates or restrict
itself to a single rate per unit of comparison.

**The estimator has a name.** Dividing each observation by the probability it was sampled with
is inverse-probability weighting, the Horvitz–Thompson estimator from survey sampling. A kept
read contributes `1/p` and a discarded one contributes 0, so its expected contribution is
`p × 1/p = 1` — exactly what it contributed before thinning. This holds for any `p`, and for any
mixture of different `p` within one sum, which is what lets scarce-junction reads at `p = 1` be
added to thinned reads at `p ≪ 1` and still give a total on the scale of the original BAM.

**Precision is set by the retained count, not the true one.** Unbiasedness is a statement about
the average over draws; a single run's error is not small just because the expectation is right.
For a group of `N` true reads sampled at rate `p`, the weighted estimate has relative standard
deviation

```
sqrt((1 - p) / (N p))  ≈  1 / sqrt(K),    K = N p = reads actually retained
```

At the worked locus below — depth 15,113, target 1000, so `p ≈ 0.066`:

| true reads in the group | retained | relative sd of the weighted estimate |
|---|---|---|
| 45 | 3 | 56% |
| 156 | 10 | 30% |
| 1,500 | 100 | 10% |
| 15,000 | 990 | 3% |

The 24-seed measurement below is consistent with this. Its second row is that locus's 156-read
junction: the measured coefficient of variation on the frequency built from it is 25% (0.0026
against 0.0106) where the rule predicts 30%, and with 24 seeds the sd estimate itself carries
roughly 15% relative uncertainty. One point of agreement, not a validation of the formula across
regimes.

Two consequences worth holding onto. **The scarce-junction rule is variance control, not bias
correction** — the weighted scheme without it is already unbiased, and what the rule removes is
the `1/sqrt(K)` term, by making `K = N`. And **a group that ends up with a handful of retained
reads cannot be rescued by weighting**, however correct the weight is: at ten retained reads the
estimate is 30% noise whatever the target was.

## XW weights during quantification

Quantification consumes the weight through one expression: `EM.py:87` takes each multipath's
support as `mp.get_read_weight()` rather than a read count, and the E-step apportions that
quantity directly. The same call backs the count rollups in `Quantify.py:1824,1846,2212`.

**There is no setting.** `MultiPath.get_read_weight()` sums a per-read registry that
`_populate_read_multi_paths` fills from whichever BAM the pass is reading
(`pylib/LRAA.py:1227-1238`), and an untagged read weighs 1 (`Pretty_alignment.py:61-65`). So
honouring the tag is a no-op on a BAM nobody thinned, and there is nothing for a caller to turn
on. `--use_XW_read_weights_for_quant` and its negation are gone; they were a mode standing in for
the role invariants above, and the invariants are what actually decide whether a weight is there.

The one deliberate exception is a pass, not a flag: `weight_reads=False` tells a single
`_populate_read_multi_paths` call to ignore whatever weights it was handed. Nothing uses it today
— discovery's draft quantification used to, and no longer does. What replaced it is a cleaner
split, because the two families of gate want genuinely different quantities:

- **Proportional / abundance → weighted.** `min_reads_novel_isoform` (via
  `Transcript.get_read_counts_assigned()`), `min_isoform_fraction`, the TPM gates, and
  `min_frac_gene_unique_reads`. These are estimates about the library, so they read weighted
  support and are on the library's scale.
- **Counts of observations → literal reads.** `min_unique_reads_novel_isoform`, the FSM counts,
  and `min_monoexonic_supporting_cells`. "Two unique reads" is a confidence statement about having
  seen a structure twice; a read retained at `p = 1/15` stands for fifteen but was still seen
  once, and re-weighting one observation cannot make it two.

`get_isoform_unique_assigned_read_count` returns both — the literal tally and the weighted support
— and hands each consumer the one it needs. That also fixes a ratio that used to mix them:
`min_frac_gene_unique_reads` divided a literal numerator by a weighted gene total, deflating it by
roughly the acceptance rate at exactly the deep loci thinning touches.

**The draft quantification reads whatever the final quant's estimating pass will read**, by the
same expression (`bam_file_for_sg if stream_reads else bam_file_for_quant`), weighted, at the same
`max_EM_iterations_quant_only` budget. What that shares is the **estimator configuration** — same
BAM, same `XW` treatment, same iteration budget — not the resulting numbers: the draft runs over
the full candidate set while the final pass re-estimates theta over the survivors only, and
dropping a candidate redistributes the reads it held. A filtered isoform's reported abundance will
not equal the value it was filtered on. The point is that both estimates sit on the same scale,
which is what makes a threshold mean the same thing in both places.

Measured on a locus of 3,002 reads thinned to 1,033: the draft EM runs on 1,033 reads,
`min_reads_novel_isoform` is handed **3095.0** — the weighted estimate of the library, 3.1% high,
which is the $1/\sqrt{1033}$ sampling error inverse-probability weighting carries — while the
unique gate is handed **1033**, the literal reads it actually saw.

Two consequences worth stating plainly. Discovery's peak memory drops, because the draft no longer
materializes the full BAM; under streaming the full BAM is now read exactly once, in pass 2. And
novel isoforms get **harder** to call at deep loci: the unique-read requirement now counts retained
reads, so a structure witnessed a handful of times in the library may retain too few to clear it.
That is the intended reading of a confidence threshold, but it is a real sensitivity change.

Everything downstream of that pass weights. What that changes is the isoform **split** within a
read-sharing component; it does not change the scale of the reported counts, because those are
sums over `--bam`, which the role invariant guarantees is the full library. Measured on the
two-isoform fixture in `pylib/test_stream_reads_xw_single_cell_matrix.py`, against a truth of
48.59 reads on the scarce isoform: weighted 48.64, unweighted 58.02 — the unweighted split ran
~19% high on the isoform distinguished by a scarce junction.

**No incompleteness marker is emitted.** LRAA used to stamp a `# WARNING:` line on a tracking file
whose rows covered only retained reads. That could only happen when the quantified BAM was itself
thinned, which the `--bam` role invariant now rejects, so the rows always enumerate the full
library. `util/lraa_merge_header.py` and `util/sc/singlecell_tracking_to_sparse_matrix.py` still
recognize and refuse the marker, deliberately: files written by earlier versions carry it
legitimately and are just as incomplete as they ever were.

Distinct from the other weight in EM: `weight_reads_by_3prime_agreement` (`EM.py:98`) is a
per-(read, transcript) 3'-end agreement weight and is unrelated to `XW`.

## Parameters

| flag | default | meaning |
|---|---|---|
| `--normalize_max_cov_level` | 1000 | target read depth; `0` disables normalization entirely |
| `--depth_window` | 100 | resolution in bases at which depth is measured |
| `--random_seed` | 42 | seed for the per-read draw |
| `--min_per_id` | 0 | percent-identity floor; must match the consumer's `min_per_id` (97 in HiFi mode, else 80). 0 disables. |
| `--min_mapping_quality` | 0 | MAPQ floor; must match whichever value the consumer will actually enforce, which under the quant stage is `--min_mapping_quality_for_final_quant` rather than `--min_mapping_quality`. Bites hardest when set, since multimapping reads carry MAPQ 0 at exactly the paralogous loci where thinning decisions matter. 0 disables. |
| `--max_intron_length` | `config` | alignments containing a longer intron are dropped at the strand split; 0 or negative disables |
| `--window_origin` | unset | absolute coordinate of the input's position 0; pins the window grid to the absolute grid. Unset: anchored on the contig's first aligned base |
| `--input_is_single_strand` | off | the input is already orientation-pure: normalize it directly, skipping the strand split and the merge |

## Caching and method changes

Both the driver and the utility checkpoint their work, and a checkpoint is trusted on sight, so
every name has to carry whatever determines the contents.

`Util_funcs.splice_graph_norm_cache_stem` names the cached BAM and its work directory. The
driver returns as soon as it sees this stem's checkpoint, so it never runs the utility and
never consults the utility's own finer-grained token — anything missing from this name is
invisible:

```
<source>.norm_<target>.maxintron_<cap>.<method>.pid<min_per_id>.mapq<min_mapq>.w<window>.s<seed>.o<origin>.scope<scope>.<identity>.bam
    sample.quant.norm_1000.maxintron_200000.cov5.pid97.mapq0.w100.s42.o0.scopenone.a56fdafac29c.bam
work_<source>.norm_<target>.maxintron_<cap>.<method>.pid<...>.mapq<...>.w<...>.s<...>.o<...>.scope<...>.<identity>/
```

Every component is required rather than defaulted, because a caller that omits one keys a BAM
by a name that does not describe it, and the omission surfaces as a silently reused cache
rather than an error.

`pid` and `mapq` are the percent-identity and mapping-quality floors, and `maxintron_` the
intron cap. Each defines the evidence universe, so two runs differing only in one of them
produce different BAMs from the same source and target. The identity floor is why a HiFi run
must not share a cache with a default one: 97 against 80 admits a materially different read
population. `mapq` is the value the consumer will actually enforce, which under the quant
stage is `--min_mapping_quality_for_final_quant`, not `--min_mapping_quality`.

`w` and `s` are `--depth_window` and `--random_seed`. The driver holds both fixed
(`NORM_DEPTH_WINDOW`, `NORM_RANDOM_SEED` in `LRAA`), which is exactly why omitting them from
the name would go unnoticed: the seed salts the per-read acceptance draw and the window
changes measured depth, so either would silently invalidate every existing cache the day a
default changed.

`o` is `--window_origin`, and it is rendered `onone` when unset rather than `o0`. Unset
anchors the window grid on each contig's first aligned base, which is a different placement
from the absolute grid at coordinate 0, so collapsing the two would let them share one cached
BAM. The driver passes `0` explicitly (`NORM_WINDOW_ORIGIN`) because it normalizes whole
contigs.

`scope` is which part of `<source>` was read: `scopenone` for the whole file, or the sorted,
`+`-joined contig names it was restricted to (a long `--restrict_to_chromosomes` list is bounded
to a count and a digest instead, so the stem cannot exceed a filesystem's path-component limit).
Without it, normalizing one contig directly from a shared whole-genome BAM would collide with a
whole-BAM normalization of that same file — two different outputs, one cache key. `--contig` and
`--restrict_to_chromosomes` are the only callers that pass a restricted scope today; every other
caller normalizes the whole source and reads `scopenone`. A restricted scope is not only a name:
`_normalize_bam_for_splice_graph` pre-filters the utility's INPUT to just those contigs (writing
a `<stem>.contigscope.bam` intermediate, gated behind the same checkpoint as the final output) so
the utility never reads records normalization was never going to keep anyway. This is safe
because it is WHOLE contigs only — `window_bases` in `sift_bam` accumulates depth per contig
already, and the driver's fixed `--window_origin` anchors window boundaries at absolute
coordinates, so which other contigs are present cannot move a kept contig's measured depth or
its per-read accept/reject decisions. `--region`'s sub-contig narrowing is deliberately excluded
from this scope: a read whose alignment spans just outside a narrow region still has to be seen
to measure that region's true depth, a reasoning whole-contig scoping does not need.

**Bump `LRAA_Globals.SPLICE_GRAPH_NORMALIZATION_METHOD` whenever the normalizer changes which
reads it keeps or what it records on them.** No consumer can detect a stale cache for itself: a
BAM from the read-start-binning era carries no `XW` tag, an absent tag legitimately means weight
1, and its distorted counts would be used in silence. The token is the only thing standing
between a superseded cache and a current run.

`<identity>` is `Util_funcs.file_identity_token`: a digest of the resolved path, size, and
modification time. The stem otherwise holds only a basename, so without it a second
`sample.bam` from another directory, or the same path regenerated, would land on the first
one's cache. Contents are not hashed — reading a multi-gigabyte BAM on every startup would cost
more than the step being cached, to cover what the stat pair already catches.

The utility keys its own checkpoints separately, because the two stages depend on different
things. The strand split depends on the input and on `--max_intron_length`, which decides which
records it emits. Sampling, merge, and index are keyed on a digest of the input identity, the
intron cap, target depth, `--depth_window`, `--random_seed`, the method, and the output path —
every flag that changes which reads are kept or where they land. The driver holds the window and
seed fixed, but both are exposed on the command line, and without this a second invocation in
the same directory would silently reuse the first one's sample.

`--window_origin` and `--input_is_single_strand` are in both tokens: the origin decides which
reads survive, and the single-strand flag decides whether the split stage runs at all. Each is
appended to the hashed string only when supplied, so a run passing neither hashes exactly what
it hashed before the options existed and keeps hitting the caches it already has.

The identity and mapping-quality floors are appended to the sampling token only, on the same
when-supplied basis. They change both the measured depth and which records are written, so a
HiFi run must never inherit a default run's sample — but the strand split applies no threshold
of its own, so including them there would invalidate a split that does not depend on them.

## Why read-start binning was replaced

The previous procedure grouped reads into 100 bp bins **by start position** and capped each bin
at 1000 by sampling within the bin. Four consequences followed.

**Its output depended on the coordinate origin.** Bin boundaries came from absolute start
positions, so translating a locus re-partitioned the bins and changed which reads survived. The
control is exact: the same 6,026 reads — identical names, identical CIGAR and sequence, every
position offset by precisely 1,203 bp, and no read starting in the strip between the two origins
— retained 2,804 reads at one origin and 1,980 at the other, and flipped whether a transcript
was emitted. LRAA itself is deterministic; triplicate runs gave byte-identical output. Under the
current scheme the same reads at origins 363 and 1,203 bp apart retain an identical set — 2,837
reads at every origin, symmetric difference zero — because window offsets are measured from the
contig's first aligned base and the acceptance draw is a hash of the read name. Neither depends
on absolute position, so the procedure is translation-invariant by construction rather than
approximately so.

This is why extracted sub-regions could not be trusted to reproduce whole-genome behaviour:
rebasing a window to position 1 is exactly the translation above, so at any deep locus the
extraction changed the answer for reasons unrelated to its read content. `--window_origin` is
what resolves that — given the chunk's rebase offset, a chunk and the whole contig assign the
same absolute locus to the same window, and the reads retained there are identical.

**It did not bound coverage.** Starts were capped per bin while depth accumulates across bins.
Across fifteen high-expression loci, depth after "capping at 1000" ranged to 3,793.

**Its reduction tracked TSS sharpness rather than depth.** Realized reduction varied from 1.0x
to 22.5x at comparable depths, because it depended on how concentrated read starts were.

**It systematically deflated 5' junctions.** A junction near the transcription start site can
only be supported by reads that start there — necessarily the crowded bin, sampled hardest. A
junction further into the gene collects reads starting anywhere upstream, across many sparse
bins that survive whole. Nothing recorded the difference, so the frequency test compared
quantities scaled by different arbitrary factors.

A worked case, `chr7:44,796,793-44,801,286` in HG002: a five-intron chain carried by 147
full-length reads and reported by four other assemblers. Depth at the locus is 15,113, and one
start bin holds 12,327 read starts and survives at 8.1% while every other bin survives whole.
The chain's second intron drew 92% of its support from that bin.

| intron | reads in BAM | after old normalization |
|---|---|---|
| 1 (5') | 320 | 41 |
| **2** | **156** | **25** |
| 3 | 12,980 | 1,672 |
| 4 | 13,500 | 2,198 |
| 5 (3') | 14,597 | 3,326 |

Its frequency against the strongest junction in its exon island is 0.0107 in the BAM and 0.0075
after normalization, against a HiFi `min_alt_splice_freq` of 0.01. The junction was purged, the
chain became unbuildable, and the 147 reads supporting it had nowhere to go. The model was
absent from the earliest candidate artifact — never built, rather than built and filtered.

## Measured behaviour

On that locus, over 24 seeds, measuring the frequency the prune decision actually rests on
(true value 0.0107, threshold 0.01):

| scheme | mean | sd | pruned |
|---|---|---|---|
| read-start binning | 0.0072 | 0.0007 | 100% of seeds |
| depth-targeted, weighted | 0.0106 | 0.0026 | 42% of seeds |
| depth-targeted, weighted, scarce-junction rule | **0.0107** | **0.0003** | **0%** |

The middle row is unbiased but noisy, which is what the scarce-junction rule removes; the
precision rule under "The XW tag" gives the quantity that sets that noise.

Across fifteen high-expression loci: worst relative-support error 2.0%, against 30% previously;
peak depth 1.78x target, against 3.79x; and full retention at every locus already under target.
On chr20 of HG002, 82.8% of reads are retained where the previous scheme retained 80.7%, and
runtime is unchanged — the second pass reads CIGARs only.

## Notes and limitations

**Weighted counts run slightly high** on heavily thinned junctions — about 5% in the worked case
above. The effect is measured; its cause has not been independently established. The standing
attribution is that median window depth underestimates true depth for reads spanning a peak,
making `1/p` overshoot, and it is recorded here as the working hypothesis rather than a settled
mechanism, for two reasons. Underestimating depth *raises* `p`, which *lowers* `1/p`. And
inverse-probability weighting is unbiased for whichever `p` the draw actually used, so a
systematic offset in a weighted sum requires retention to correlate with weight — which a
per-read median could produce, since a read's own measured depth decides both its `p` and whether
it survives, but that link has not been demonstrated. Settling it needs its own experiment:
weighted junction support against unnormalized support, across many junctions and seeds at known
depth, with the per-read median recorded alongside. The 24-seed table above does not answer it —
that measures a junction-frequency ratio, not a weighted count against its true value. Treat the
5% as real and its explanation as open. The direction is conservative either way: it moves ratios
toward the filtering threshold rather than away from it, and it is uncorrected.

**The scarce-evidence exemption is junction-based only.** `_acceptance_probability` inspects
`_read_junctions` and nothing else (`util/normalize_bam_by_strand.py:630-633`), so evidence that
is scarce without being a scarce junction gets no protection: an isoform distinguished from its
neighbours only by TSS position, only by PolyA position, or a monoexonic model has no junction to
claim an exemption with, and its support at a deep locus is estimated from whatever survived at
that locus's own rate — the `1/sqrt(K)` regime above, with `K` in the tens. Extending the
exemption to those features would need a first-pass tally of TSS and PolyA positions, which pass
1 does not keep: it accumulates window depth and junction support only
(`util/normalize_bam_by_strand.py:541-551`).

**Absolute thresholds on support changed meaning.** Support is now on the scale of the original
BAM everywhere. `min_intron_cov_for_filtering`, the one absolute threshold applied to junction
support, previously held high-coverage loci to a far stricter bar than shallow ones because it
compared against sampled counts. It is now consistent across coverage levels.

**Scope of the benefit is narrow but concentrated.** On chr20 of HG002, 45 of 40,323 occupied
start bins exceed the target — 0.11% of bins, but 19.3% of all reads. These cluster into 41
loci, 7 of them cut more than fourfold. Only about 1.3% of recovered reference chains sit at
such a locus, so benchmark-level metrics move little; the correctness gain is at the
high-expression tail.
