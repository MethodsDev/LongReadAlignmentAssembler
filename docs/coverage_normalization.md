# Coverage Normalization for Splice-Graph Construction

## Overview

Splice-graph construction scales with read depth, and a handful of loci in any sample carry
depth far beyond what graph construction needs. LRAA thins those loci before building the
graph, leaving the rest of the genome untouched. Quantification is unaffected: abundance
estimation reads the original BAM, and the thinned BAM is used only as splice-graph evidence.

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

Two sequential passes over each strand-specific BAM.

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

Weighting is not optional bookkeeping. Any scheme whose acceptance rate varies along the genome
distorts relative support, and every such scheme must either record its rates or restrict
itself to a single rate per unit of comparison.

## Parameters

| flag | default | meaning |
|---|---|---|
| `--normalize_max_cov_level` | 1000 | target read depth; `0` disables normalization entirely |
| `--depth_window` | 100 | resolution in bases at which depth is measured |
| `--random_seed` | 42 | seed for the per-read draw |

## Caching and method changes

Both the driver and the utility checkpoint their work, and a checkpoint is trusted on sight, so
every name has to carry whatever determines the contents.

`Util_funcs.splice_graph_norm_cache_stem` names the cached BAM and its work directory:

```
<source>.norm_<target>.<method>.<identity>.bam    sample.quant.norm_1000.cov1.a56fdafac29c.bam
work_<source>.norm_<target>.<method>.<identity>/  the utility runs here
```

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
things. The strand split depends only on the input, so it is keyed on `<identity>` alone and is
reused across settings. Sampling, merge, and index are keyed on a digest of the input identity,
target depth, `--depth_window`, `--random_seed`, the method, and the output path — every flag
that changes which reads are kept or where they land. The driver holds the window and seed
fixed, but both are exposed on the command line, and without this a second invocation in the
same directory would silently reuse the first one's sample.

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
extraction changed the answer for reasons unrelated to its read content.

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

The middle row is unbiased but noisy, which is what the scarce-junction rule removes.

Across fifteen high-expression loci: worst relative-support error 2.0%, against 30% previously;
peak depth 1.78x target, against 3.79x; and full retention at every locus already under target.
On chr20 of HG002, 82.8% of reads are retained where the previous scheme retained 80.7%, and
runtime is unchanged — the second pass reads CIGARs only.

## Notes and limitations

**Weighted counts run slightly high** on heavily thinned junctions — about 5% in the worked
case above. Median window depth somewhat underestimates true depth for reads spanning a peak,
which makes `1/p` overshoot. The effect moves ratios toward the threshold rather than away from
it, so it is conservative, but it is uncorrected.

**Absolute thresholds on support changed meaning.** Support is now on the scale of the original
BAM everywhere. `min_intron_cov_for_filtering`, the one absolute threshold applied to junction
support, previously held high-coverage loci to a far stricter bar than shallow ones because it
compared against sampled counts. It is now consistent across coverage levels.

**Scope of the benefit is narrow but concentrated.** On chr20 of HG002, 45 of 40,323 occupied
start bins exceed the target — 0.11% of bins, but 19.3% of all reads. These cluster into 41
loci, 7 of them cut more than fourfold. Only about 1.3% of recovered reference chains sit at
such a locus, so benchmark-level metrics move little; the correctness gain is at the
high-expression tail.
