#!/usr/bin/env python3

"""Choose the cut points that split a contig-strand into independently
processable chunks.

The pipeline this feeds is: split the BAM by orientation, choose cut points per
(contig, orientation), extract rebased chunks, normalize and run downstream PER
CHUNK in parallel, then merge and translate coordinates back. Normalization moves
inside the chunk because it is ~61 s of serial work that otherwise blocks
everything else, and per-chunk normalization parallelizes it.

Sizing is by SPAN, not by alignment or gene count
-------------------------------------------------
Targets sit at multiples of ``approx_MB_per_cut`` (default 10 Mb) across the
contig, with no cap on chunks per contig: chr20 at 64.4 Mb gets 6, chr1 at 248.9
Mb gets 25. Span-based sizing makes a chunk's coordinates a function of the contig
length alone, so they can be computed without reading the BAM and are stable
across runs and across samples.

Around each target the position is searched over a window whose TOTAL width
``approx_MB_per_cut_wiggle_window`` (default 1 Mb) is a MAXIMUM rather than the
radius actually searched: see "Progressive expansion" below. That default is
ABSOLUTE and is deliberately not derived from the cut spacing -- see the config
comment on the key in ``pylib/LRAA_Globals.py`` for the measurement, which is that
a proportional rule would give 200 kb at 2 Mb spacing and still sever 743
alignments on chr21 where 1 Mb severs none, while chr1 and chr21 disagree 8.5-fold
at identical parameters. A caller may pass a smaller window to test a finer
spacing, and the reads that then get severed are reported rather than refused.

Three conditions on a position, in order of authority
-----------------------------------------------------
1. GRID ALIGNMENT, mechanical. ``(b - grid_origin) % depth_window == 0`` for a
   1-based cut ``b``, where the left chunk is ``[.., b]`` and the right chunk
   begins at ``b + 1``. This is what lets per-chunk normalization reproduce
   whole-contig normalization: if a depth window straddled a boundary it would
   draw aligned bases from two chunks, so its median, the acceptance probability
   derived from it, and the reads kept near that threshold would all differ.
   ``b`` is the LAST base of a window and ``b + 1`` the first base of the next --
   requiring ``b`` to START a window would instead put a window across every cut.
   The pipeline pins ``grid_origin`` to absolute 0, which reduces the condition to
   ``b % depth_window == 0``; the parameter exists so the frame is explicit and
   the relative form is what gets tested.

2. THE ANNOTATION, a hard constraint. When a GTF is supplied a cut may never fall
   inside an annotated locus, and a target with no compliant position in its
   window is REPORTED, never silently skipped. Unlike a read, a locus cannot be
   dropped and accounted for: ``genes_contained`` emits a gene whole or not at
   all, so a locus straddling a boundary is contained by neither neighbour and
   both omit it. Enforced via ``find_islands``/``cut_zones`` from the extractor,
   the same primitives emission itself checks, so a position this module chooses
   cannot be one the extractor then refuses.

3. SPANNING ALIGNMENTS, a weighted cost, and never a veto. Among compliant
   positions, minimise the cost of the retained primary alignments the cut
   severs, since each one is a read that gets dropped. A MONOEXONIC severed
   alignment costs 1 and a MULTI-EXON one costs ``severed_multiexon_weight``
   (default 10): a spliced alignment carries the junction evidence the splice
   graph's edges are built from, a monoexonic one carries none. Ties go to the
   position nearest the target, then to the lower coordinate, so the result is
   deterministic.

Severing is priced, and that is ALL it is
-----------------------------------------
An earlier revision of this module made zero severed a HARD constraint in
discovery: severing positions were struck from the candidate set and a target
whose window held none was declined. That was reviewed and REJECTED, and the
reason is a property of the data rather than a preference.

As libraries get deeper, every base ends up covered by some read. A rule that
forbids severing therefore declines every target on a deep contig, and declining
every target silently disables chunking on exactly the inputs that need it most
while reporting success. Nor can it be tuned out of that regime by widening: it
is not that the window was too small, it is that no clean position exists at any
width.

So a severed read in DISCOVERY is treated as it is in quantification: dropped,
counted and named. What survives from the measurement that motivated the hard
rule is the SHAPE of the cost, not the veto. Measured on chr21, HG002 PacBio
Kinnex (``FINDINGS.chr21_denovo_parity.md``, LRAA_PAPER_Analyses ``05af6d8``): at
the shipped 10 Mb spacing the cuts sever ZERO alignments and de novo discovery
differs from the unchunked run by one model in 1462, that one 2,724,898 bp from
the nearest cut. Forced to 2 Mb spacing with a 0.02 Mb window the selector has to
sever 940 alignments; 17 of the 20 models then lost SPAN a cut, and at
``chr21(-):41,986,832`` a gene of eight isoforms whose first exon crosses the cut
at 41,992,400 loses all eight while the chunked arm emits two SPURIOUS MONOEXONIC
models whose left ends coincide with the originals'. The damage there was to
JUNCTIONS, which is what the multi-exon weight prices, and the geometry that
caused it was a window starved 40-fold below what the data need. Both point at a
better search rather than at a refusal.

Progressive expansion, and why the annulus is not an optimisation
----------------------------------------------------------------
The search starts at a small radius around the target and widens through
``DEFAULT_EXPANSION_RADII`` (5 kb, 25 kb, 100 kb, 250 kb, 500 kb) up to the
wiggle half-width, stopping as soon as an annotation-compliant position severs
NOTHING. At the 10 Mb default on chr1 that terminates at a median radius of 25 kb
against a 500 kb maximum. On reaching the maximum it takes the best position
available, and the alignments that position severs are collateral: counted, named
and reported.

Each rung fetches only the NEW ANNULUS -- two intervals, one either side of what
is already scanned -- and carries the running best forward. That is what keeps
expansion from being a REGRESSION rather than a speed-up. Re-scanning from the
target at every rung costs 2*(5+25+100+250+500) kb = 1.76 Mb against the flat
1 Mb of scanning the whole window once, so naive expansion is 1.76x WORSE than no
expansion at all in the regime where nothing terminates early -- which is the
deep regime this design exists to serve. Annulus fetching makes the worst case
equal to the flat cost and the common case a fraction of it.

Terminating early is not an approximation. If any compliant position at radius R
severs nothing then the whole-window minimum is also zero, and among zero-cost
positions the tie-break is distance from the target, so every candidate that could
win already lies inside radius R. The chosen position is the one a full scan would
choose; only the I/O differs. Where neighbouring windows overlap and the
anti-sliver floor binds, a target the joint solve had to compromise is expanded
another rung and the solve repeated, which carries the same property over to the
joint problem.

The ANNOTATION remains a hard block, and is now the ONLY reason a target can be
declined. A read can be dropped and accounted for; a locus cannot -- see
condition 2. A target whose whole window is annotation-blocked is DECLINED and
reported, and the two chunks it would have separated stay joined. That is
accepted behaviour, not a defect to be fixed by relaxing the block.

Selection is JOINT, not per-target
----------------------------------
Positions must be strictly increasing and no realised chunk may fall below
``minimum_span``, which defaults to half the nominal segment length (floored to a
``depth_window`` multiple, never below one window). Reasoning: with the shipped
defaults adjacent positions are already at least ``segment - wiggle`` = 9 Mb
apart, so the constraint never binds and costs nothing; it becomes active only
when the wiggle window is widened until windows overlap, which is exactly the
regime where two targets could otherwise collapse into a sliver. Half the nominal
span is the natural statement of "no chunk may be less than half the size it was
asked to be", and it is always satisfiable in isolation because placing every
position at its target satisfies it.

The same floor applies to the tail: a trailing target whose residual would be
below ``minimum_span`` is merged into its predecessor rather than left as a
sliver, and the merge is reported. On chr20 that is what turns 6 targets into 5
cuts and 6 chunks, at the cost of a final chunk of 14.4 Mb.

A dynamic program does the joint choice, minimising in lexicographic order:
targets left unplaced, then total severing COST, then total distance from
target, then total coordinate. Ordering "unplaced" first means the selector never
trades a cut away to save a few reads. When windows are disjoint and the span
floor does not bind, the program reduces exactly to the per-target argmin over
``(cost, |offset|, position)`` described above.

Accounting
----------
Every count carries its denominator and nothing is dropped silently. Per cut:
target, chosen position, offset from target, the radius the search stopped at,
the alignments severed split MONOEXONIC against MULTI-EXON and the weighted cost
they add up to, whether the annotation forced a compromise, and how many grid
positions the annotation removed from the window. Per contig-strand: realised
chunk spans, the total dropped against the total retained primary alignments, and
every unplaced target with the reason. The dropped read NAMES are written out,
because the comparison methodology builds a pruned baseline BAM from exactly that
set and a count alone would not let anyone reproduce it.

The split matters more than the total now that severing is expected rather than
exceptional. It is the only place a reader can see what a cut actually cost: two
monoexonic reads and one spliced read are three severed alignments either way, and
the second is the one that can cost a junction.

Measured on chr20, and what it means for wiring the pipeline
------------------------------------------------------------
At the shipped defaults on ``rescue_probe_chr20/mini`` (10 Mb span, 1 Mb wiggle,
depth_window 100, margin 200), both orientations place 5 cuts into 6 chunks and
drop ZERO alignments:

  chr20+  9,985,800 (-14,200) | 20,361,000 (+361,000) | 30 Mb | 40 Mb | 50 Mb
  chr20-  10 Mb | 20 Mb | 30 Mb | 40 Mb | 50 Mb, every one exactly on target

The denominator is 319,698 retained primary alignments = 174,106 (+) + 145,592
(-). That is NOT the 320,240 figure quoted elsewhere: 320,240 counts mapped
primaries BEFORE the intron filter, and 542 of them (340 on +, 202 on -) carry an
intron over ``max_intron_length`` and are discarded at the strand split. Anything
downstream of that split must use 319,698, and a pruned comparison baseline built
against 320,240 would be wrong.

Consequence for the verification arms, stated because it is easy to over-read:
with nothing dropped, pruning removes nothing, so ON THIS DATASET the pruned
baseline is just the ordinary intron-filtered baseline and the "what did dropping
cost" comparison is trivially zero. The chunked-versus-baseline arm can therefore
run at zero tolerance with no pruning step, provided the whole-contig control is
given ``--window_origin 0`` -- its unset default anchors on the first aligned base
per contig and no chunk grid can match that. This is a property of this dataset at
these settings, not of the design: at a 0.02 Mb wiggle the 20 Mb target lands at
20,004,500 and drops 2 named alignments, and at wiggle 0 the 10 Mb and 20 Mb
targets are annotation-blocked outright and come back as unplaced.

That the search earns its keep is measurable rather than assumed. In the 20 Mb
window 7,993 of 10,001 grid positions are annotation-blocked, the worst compliant
position carries 1,861 spanning alignments, and the target itself is both blocked
AND would have severed 46 reads -- yet 1,280 zero-cost compliant positions exist,
the nearest at 20,361,000. Removing the GTF moves that cut to 20,166,800 and the
10 Mb cut to exactly 10,000,000, which isolates the annotation as the sole cause
of both offsets.

The 14.44 Mb final chunk is the anti-sliver floor, not a malfunction: the 60 Mb
target would have left a 4.44 Mb residual against a 5 Mb minimum span, so it is
merged into its predecessor and reported as ``tail_merged``. The tail is the ONLY
place that floor can act at the shipped defaults -- a 1 Mb wiggle leaves adjacent
cuts at least ``span - wiggle`` = 9 Mb apart, comfortably above 5 Mb -- so between
targets the constraint is provably slack and only an enlarged wiggle can engage it.

Strandless selection
--------------------
``--strandless`` chooses ONE set of cuts over the RAW bam, serving both strands,
so that the strand split can move out of the whole-genome serial phase and into
the per-chunk parallel phase. Both conditions above become the union without any
new machinery, and that is worth stating precisely because it is easy to assume
instead of check:

* the annotation. ``find_islands`` skips its strand filter on a falsy strand, so
  the islands are the union of both strands' gene loci and a cut must clear every
  locus on either strand. Measured on chr1 with the shipped margin, the union
  blocks 148.7 Mb -- 40.3% of the contig admissible against 67-69% per strand --
  and at the 10 Mb default all 24 targets still place with a median of 3,726
  admissible positions per wiggle window. At 2.5 Mb it degrades and the selector
  reports it; it is not silently absorbed.
* the alignments. ``spanning_counts`` scores a position through
  ``extractor.retained_for_extraction``, which admits both orientations on a
  falsy strand, so a position is priced by everything it severs on either strand.
  That is the correct objective: the two strand chunks share the boundary, so
  they share its cost.

What ``--strandless`` adds over simply omitting ``--strand`` is the declaration
and its accounting -- strandless region strings, a per-orientation denominator,
and a warning if the bam turns out to hold one orientation, which means a
strand-separated bam was passed and the "strandless" cuts serve one strand.
"""

import argparse
import bisect
import collections
import importlib.util
import json
import os
import sys
from importlib.machinery import SourceFileLoader

import pysam

sys.path.insert(
    0,
    os.path.sep.join(
        [os.path.dirname(os.path.realpath(__file__)), "..", "..", "pylib"]
    ),
)

import LRAA_Globals
import Util_funcs


def _load_extractor():
    """Load the sibling extractor as a module.

    It is a script rather than a package member, so it cannot simply be
    imported. Sharing it matters: the annotation hard constraint and the
    retention filter must be the SAME code the extractor enforces, or this
    module could choose a position the extractor then refuses.
    """

    path = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "extract_contig_region_inputs.py"
    )
    loader = SourceFileLoader("extract_contig_region_inputs", path)
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


extractor = _load_extractor()

ExtractionError = extractor.ExtractionError

MB = 1000000

# The depth-window grid is pinned to absolute 0-based coordinate 0. Chunks are
# rebased, and normalize_bam_by_strand.py is told the rebase offset via
# --window_origin, so chunk-local windows ARE the absolute windows. Anything else
# would make the grid depend on which read happened to come first.
DEFAULT_GRID_ORIGIN = LRAA_Globals.config["chunk_grid_origin"]

DEFAULT_DEPTH_WINDOW = LRAA_Globals.config["chunk_depth_window"]

DEFAULT_MULTIEXON_WEIGHT = LRAA_Globals.config["chunk_severed_multiexon_weight"]

# Search radii in BASES, increasing, each rung fetching only the annulus the
# previous one did not cover. Clipped to the wiggle half-width, which is the
# maximum. Measured at the 10 Mb default on chr1: a median terminal radius of
# 25 kb -- the second rung -- against the 500 kb maximum, because a zero-cost
# position is usually close. The ladder is geometric-ish rather than uniform so
# the common case pays 5 kb and the rare case still reaches 500 kb in five steps;
# a uniform ladder either wastes the common case or takes too many rungs.
DEFAULT_EXPANSION_RADII = (5000, 25000, 100000, 250000, 500000)

# Opening the bam goes through one name so a test can wrap it and observe which
# intervals were fetched. Annulus-only fetching is a property of the FETCHES, not
# of the answer -- a full re-scan at every rung returns the same positions -- so
# the only way to hold the code to it is to count the intervals.
_open_bam = pysam.AlignmentFile


class SelectionError(RuntimeError):
    """Parameters that cannot describe a partitioning."""


CutChoice = collections.namedtuple(
    "CutChoice",
    "index target position offset spanning_dropped compromised "
    "window_lend window_rend grid_positions annotation_blocked "
    "unconstrained_best_spanning "
    # The severed alignments split by structure, the weighted cost they sum to,
    # and how far the progressive search had to travel. Severing is expected
    # rather than exceptional now, so the split is the only place its real cost
    # is visible: a severed spliced read can cost a junction, a monoexonic one
    # cannot.
    "severed_monoexonic severed_multiexon severed_weight search_radius",
)

UnplacedTarget = collections.namedtuple(
    "UnplacedTarget",
    "index target window_lend window_rend grid_positions annotation_blocked reason "
    # The cheapest position in the window by weighted cost, whether the ANNOTATION
    # is what left the window with nothing usable -- the only remaining decline
    # reason -- and the radius the search reached.
    "best_spanning declined_annotation search_radius",
    defaults=(None, False, None),
)

Segment = collections.namedtuple("Segment", "index lend rend span")

Selection = collections.namedtuple(
    "Selection",
    "chrom strand contig_length depth_window grid_origin segment_span wiggle "
    "minimum_span targets cuts segments unplaced tail_merged "
    "total_retained_primary total_dropped dropped_read_names "
    # Declared, and the two orientation counts that substantiate it. Defaulted so
    # a caller that predates strandless selection constructs the same tuple.
    "strandless retained_primary_forward retained_primary_reverse "
    # What a severed multi-exon alignment cost this selection, against 1 for a
    # monoexonic one. Recorded because it decides the cut coordinates: a manifest
    # that does not state it cannot be compared with another run's.
    "multiexon_weight",
    defaults=(False, None, None, DEFAULT_MULTIEXON_WEIGHT),
)


def grid_positions(lend, rend, depth_window, grid_origin=DEFAULT_GRID_ORIGIN):
    """1-based cut positions in ``[lend, rend]`` that end a depth window.

    ``(b - grid_origin) % depth_window == 0``. Returned in increasing order, so
    the list is an arithmetic sequence with step ``depth_window`` and callers may
    index into it arithmetically.
    """

    if depth_window < 1:
        raise SelectionError("depth_window must be >= 1")
    if rend < lend:
        return []
    # smallest b >= lend with (b - grid_origin) divisible by depth_window
    first = lend + (grid_origin - lend) % depth_window
    if first > rend:
        return []
    return list(range(first, rend + 1, depth_window))


def minimum_span_for(segment_span, depth_window):
    """Half the nominal segment, floored to a whole depth window.

    Never below one window, so the floor cannot be degenerate, and always a
    multiple of ``depth_window`` so it cannot fight grid alignment.
    """

    half = (segment_span // 2 // depth_window) * depth_window
    return max(depth_window, half)


def expansion_rungs(half, radii=None):
    """Increasing search radii, ending at ``half``, which is the maximum.

    ``half`` is the wiggle HALF-width. Rungs above it are pointless and a rung
    equal to it ends the ladder, so the result always finishes exactly at ``half``
    and never exceeds it: the window width is a promise, and the last rung is
    where "take the best available" happens. A zero ``half`` gives ``[0]``, the
    target itself.
    """

    if radii is None:
        radii = DEFAULT_EXPANSION_RADII
    rungs = []
    for radius in sorted(radii):
        if radius < 1:
            continue
        clipped = min(radius, half)
        if not rungs or clipped > rungs[-1]:
            rungs.append(clipped)
        if clipped >= half:
            break
    if not rungs or rungs[-1] < half:
        rungs.append(half)
    return rungs


def annulus_intervals(target, inner, outer, window_lend, window_rend):
    """The position ranges radius ``outer`` adds to radius ``inner``.

    ``inner`` is None before anything has been scanned, and the answer is then the
    single interval ``target +/- outer``. Otherwise it is the two intervals either
    side of what is already scanned, which is why widening the search does not
    re-read the middle: ``[target - outer, target - inner - 1]`` and
    ``[target + inner + 1, target + outer]``. Clipped to the window and dropped
    when empty, so a rung the window edge or the contig end already covers costs
    no fetch at all.
    """

    if inner is None:
        spans = [(target - outer, target + outer)]
    else:
        spans = [
            (target - outer, target - inner - 1),
            (target + inner + 1, target + outer),
        ]
    intervals = []
    for lend, rend in spans:
        lend = max(lend, window_lend)
        rend = min(rend, window_rend)
        if lend <= rend:
            intervals.append((lend, rend))
    return intervals


def cut_targets(contig_length, segment_span, minimum_span):
    """Target cut positions, and the trailing targets merged away.

    Targets are the multiples of ``segment_span`` strictly inside the contig. A
    trailing target whose residual tail would be shorter than ``minimum_span`` is
    merged into its predecessor: the same anti-sliver rule the joint constraint
    applies between targets, applied to the end of the contig where there is no
    successor to collide with. Returns ``(targets, merged)``, and ``merged`` is
    reported rather than discarded.
    """

    if segment_span < 1:
        raise SelectionError("segment_span must be >= 1")

    targets = list(range(segment_span, max(contig_length, 1), segment_span))
    merged = []
    while targets and contig_length - targets[-1] < minimum_span:
        merged.append(targets.pop())
    # a first target closer to the contig start than the floor is a sliver too
    while targets and targets[0] < minimum_span:
        merged.append(targets.pop(0))
    merged.sort()
    return targets, merged


def _blocked_by_annotation(positions, zones):
    """Split ``positions`` into (compliant, blocked) against admissible zones.

    ``zones`` are the inclusive ranges of positions no annotated locus forbids,
    from ``extractor.cut_zones``, in increasing order.
    """

    if not zones:
        return [], list(positions)
    starts = [lo for lo, _ in zones]
    compliant = []
    blocked = []
    for position in positions:
        index = bisect.bisect_right(starts, position) - 1
        if index >= 0 and position <= zones[index][1]:
            compliant.append(position)
        else:
            blocked.append(position)
    return compliant, blocked


def is_multiexon(aln):
    """Does this alignment carry a junction?

    An ``N`` CIGAR op is what the splice graph turns into an edge, so one is
    exactly what makes severing this alignment cost more than a read's worth of
    depth. ``retained_for_extraction`` has already discarded alignments whose
    introns are too long to be real, so every ``N`` seen here is one the
    downstream stages would have used.
    """

    return any(op == 3 for op, _ in (aln.cigartuples or ()))


def spanning_counts(
    bam,
    chrom,
    strand,
    positions,
    max_intron_length=None,
    multiexon_weight=1,
):
    """Weighted cost of the retained primary alignments each position severs.

    ``positions`` must be an increasing arithmetic sequence, as returned by
    ``grid_positions``. An alignment spanning 1-based ``[s, e]`` is severed by a
    cut at ``b`` exactly when ``s <= b < e``, so it contributes to the position
    range ``[s, e - 1]``; the contributions are accumulated as a difference array
    over the candidate grid rather than tested position by position, which keeps
    this linear in the alignments overlapping the window rather than quadratic.

    Only the interval ``[first - 1, last + 1)`` is fetched, so a caller may score
    one ANNULUS of a progressively widening search without re-reading the middle:
    an alignment severing a position in this range necessarily overlaps the range,
    so the answer for these positions does not depend on what lies outside it.

    The strand reaches the cost only through ``retained_for_extraction``, which
    reads a falsy strand as "either orientation". So a STRANDLESS selection costs
    a position by every alignment it severs on both strands, which is the correct
    objective for a cut that has to serve both: the two strand chunks share the
    boundary, so they share its price.

    A monoexonic severed alignment costs 1 and a multi-exon one
    ``multiexon_weight``. An earlier revision of this docstring said a weighted
    cost was deliberately NOT built, on the grounds that at the 10 Mb default on
    chr1 no monoexonic alignment spans any sampled cut target and the selector
    already reaches 0 severed, so nothing could be checked against it. The
    measurement stands and the conclusion drawn from it does not: it describes the
    regime where the weight is inert and says nothing about the regime it exists
    for. The weight is built, and it is deliberately built with a MEASURED
    near-inertness on today's corpora, recorded here so nobody later reads a flat
    sweep as a broken feature.

    What the severed alignments actually are, 2 Mb spacing, HG002 PacBio Kinnex,
    monoexonic against spliced:

        chr21, 20 kb window  |  14 mono | 926 spliced
        chr21, 200 kb window |   0      | 743
        chr1,  20 kb window  |   5      | 2593
        chr1,  200 kb window |   2      |   85

    98-100% spliced at every geometry measured, because spanning probability scales
    with genomic span: a monoexonic read covers a kb or two against tens of kb for
    a spliced one, so it almost never reaches a cut at all. K will therefore rarely
    change which position wins on data like this. It is still the right cost --
    severing a junction-bearing read removes an edge the splice graph is built
    from, severing a monoexonic read removes depth -- and the mix is exactly what
    changes as libraries deepen and spacings shrink.

    What NOT to do with that inertness: it is not a licence to count multi-exon
    alignments only and let monoexonic severing be free. At depth the mix is the
    thing in motion, and a free class of severing is an unbounded cost waiting for
    the data that produces it. The default of 1 here is "unweighted", for callers
    scoring a raw count; ``select_cut_points`` passes the configured weight.
    """

    if not positions:
        return []

    first = positions[0]
    last = positions[-1]
    step = positions[1] - positions[0] if len(positions) > 1 else 1
    diff = [0] * (len(positions) + 1)

    # every severing alignment satisfies s <= last and e >= first + 1, hence
    # overlaps 1-based [first, last + 1]; fetch is a superset of that
    for aln in bam.fetch(chrom, max(0, first - 1), last + 1):
        if not extractor.retained_for_extraction(aln, strand, max_intron_length):
            continue
        start = aln.reference_start + 1
        end = aln.reference_end
        lo = max(start, first)
        hi = min(end - 1, last)
        if hi < lo:
            continue
        lo_index = -(-(lo - first) // step)  # ceil
        hi_index = (hi - first) // step
        if hi_index < lo_index:
            continue
        weight = multiexon_weight if is_multiexon(aln) else 1
        diff[lo_index] += weight
        diff[hi_index + 1] -= weight

    counts = []
    running = 0
    for value in diff[:-1]:
        running += value
        counts.append(running)
    return counts


def spanning_alignments(
    bam,
    chrom,
    strand,
    position,
    max_intron_length=None,
    quant_only=False,
    min_per_id=None,
    min_mapping_quality=None,
):
    """The alignments a cut at ``position`` severs.

    ``quant_only=False`` uses the same superset predicate the cost scoring uses,
    so the names reported for a cut match the number reported against it.

    ``quant_only=True`` restricts to alignments quantification would actually use,
    via the shared policy in ``Util_funcs.quant_discard_reason``: mapping quality
    and percent identity as well.  Emission needs that stricter set.  The superset
    is sound for BOUNDING a cut's cost, because over-counting only biases toward
    safer positions -- but a consumer asking which severed reads linked two genes
    is asking about reads quantification sees, and one it would have rejected on
    MAPQ or identity is a false positive there, not a conservative one.

    The thresholds are arguments rather than config reads because the EFFECTIVE
    values are decided by the invoking run, not by the defaults: --HiFi raises
    min_per_id from 80 to 97 inside LRAA, so a selector reading its own config
    would emit reads that run rejects.  Callers must pass what the quant step will
    use.
    """

    for aln in bam.fetch(chrom, max(0, position - 1), position + 1):
        if quant_only:
            if not Util_funcs.retained_for_quantification(
                aln,
                strand,
                max_intron_length=max_intron_length,
                min_per_id=min_per_id,
                min_mapping_quality=min_mapping_quality,
            ):
                continue
        elif not extractor.retained_for_extraction(aln, strand, max_intron_length):
            continue
        start = aln.reference_start + 1
        end = aln.reference_end
        if start <= position < end:
            yield aln


def spanning_read_names(bam, chrom, strand, position, max_intron_length=None):
    """Names of the retained primary alignments a cut at ``position`` severs.

    The set the extractor will drop at this boundary, named so a pruned
    comparison baseline can be built from exactly it.
    """

    return [
        aln.query_name
        for aln in spanning_alignments(bam, chrom, strand, position, max_intron_length)
    ]


def write_severed_alignments_bam(
    header, alignments, path, min_per_id=None, min_mapping_quality=None
):
    """A coordinate-sorted, indexed BAM of the alignments the cuts sever.

    Names alone cannot answer what a severed read was compatible with.  Transcript
    compatibility is a function of exon blocks and junctions, not extent, so two
    reads with identical start and end can match different genes -- and for the
    cross-gene question the junction structure is precisely what decides it.  A
    span would give silent mis-attribution rather than an answer.

    Names also cannot be looked up: a coordinate-indexed BAM cannot fetch by query
    name, so attribution from the manifest alone costs a full pass over the source
    however few names it holds.  Emitting the records while they are in hand is
    what keeps that check local.

    What this is NOT: the set quantification used.  See the CO lines written into
    the output for the scope a consumer has to respect.

    Takes a header rather than an open file because the caller writes ONCE across
    every contig it selected on; a path threaded into per-contig selection would
    leave only the last contig's reads behind.

    Sorted in memory rather than by an external sort: cuts ascend, but reads
    severed by a later cut can start before those severed by an earlier one, and
    the set is small by construction -- a cut severing many reads is one the
    selector rejects.
    """

    by_position = sorted(
        alignments, key=lambda aln: (aln.reference_id, aln.reference_start, aln.query_name)
    )
    # The scope travels with the artifact, because every consumer so far has
    # over-read it.  These records are what EXTRACTION would retain at the stated
    # thresholds, taken from the pre-normalization strand split.  They are NOT the
    # set quantification saw: normalization runs two stages later and samples by
    # local depth, so some of these would have been dropped there.  Intersecting
    # them with a per-chunk normalized bam finds nothing by construction -- chunk
    # extraction discarded them before normalization ran -- so establishing that a
    # severed read dissolved a component needs the whole-library normalized bam,
    # which the chunked arm does not produce.
    stamped = header.to_dict() if hasattr(header, "to_dict") else dict(header)
    notes = list(stamped.get("CO", []))
    notes.append(
        "LRAA severed alignments: extraction-retained, PRE-normalization. "
        "min_per_id={} min_mapping_quality={}".format(min_per_id, min_mapping_quality)
    )
    notes.append(
        "LRAA severed alignments: NOT the quant-visible set. Zero records proves no "
        "component dissolved. Nonzero requires a whole-library normalized bam to "
        "attribute; without one, refuse rather than approximate."
    )
    stamped["CO"] = notes

    with pysam.AlignmentFile(path, "wb", header=stamped) as out:
        for aln in by_position:
            out.write(aln)
    pysam.index(str(path))
    return len(by_position)


def count_retained_primary(bam, chrom, strand, max_intron_length=None):
    """The denominator: retained primary alignments on this contig-strand."""

    total, _, _ = count_retained_primary_by_orientation(
        bam, chrom, strand, max_intron_length
    )
    return total


def count_retained_primary_by_orientation(bam, chrom, strand, max_intron_length=None):
    """``(total, forward, reverse)`` over the same set the denominator counts.

    One pass, because the whole-contig scan is the expensive part and the split
    is what tells a raw bam from a strand-separated one: after
    ``separate_bam_by_strand.py`` every retained record carries the orientation
    it was assigned to, so a zero on either side means the input was already
    split. Reported, never inferred from, since a sparse contig can be one-sided
    for honest reasons.
    """

    total = 0
    forward = 0
    reverse = 0
    for aln in bam.fetch(chrom):
        if not extractor.retained_for_extraction(aln, strand, max_intron_length):
            continue
        total += 1
        if aln.is_forward:
            forward += 1
        else:
            reverse += 1
    return total, forward, reverse


def _solve(candidates, costs, targets, minimum_span, contig_length):
    """Joint choice of one position per target, or fewer.

    Lexicographic minimisation of (unplaced, spanning, distance from target,
    coordinate) over chains that are strictly increasing and leave every realised
    chunk at least ``minimum_span`` long. Returns the list of chosen
    ``(target_index, position, cost)`` in increasing order.

    A dynamic program rather than a per-target argmin because the constraints
    couple neighbours: with overlapping windows the cheapest position for one
    target can crowd out its successor. Each node contributes the additive vector
    ``(-1, cost, |position - target|, position)``, so "fewest unplaced" is
    "most placed" and the key stays additive -- which is what allows the
    predecessor search to be a prefix minimum over positions instead of a scan
    over pairs.
    """

    # nodes flattened, grouped by target, each with the best key of a chain
    # ending there and a backpointer
    best = []  # per target: list of (key, prev_target, prev_index)
    # sorted (position, prefix-min key, target, index) over targets already seen
    seen_positions = []
    seen_keys = []
    seen_refs = []

    for t_index, (positions, position_costs) in enumerate(zip(candidates, costs)):
        row = []
        for p_index, position in enumerate(positions):
            contribution = (
                -1,
                position_costs[p_index],
                abs(position - targets[t_index]),
                position,
            )
            # predecessor must sit at most position - minimum_span
            limit = position - minimum_span
            key = None
            ref = (None, None)
            if limit >= 0:
                cut = bisect.bisect_right(seen_positions, limit) - 1
                if cut >= 0:
                    key = seen_keys[cut]
                    ref = seen_refs[cut]
                elif position >= minimum_span:
                    # start the chain here: the first chunk is [1, position]
                    key = (0, 0, 0, 0)
            # no "start fresh" alternative is considered when a predecessor
            # exists: every predecessor key has placed >= 1, hence a first
            # component <= -1, so it always beats the empty chain's (0, ...)
            if key is None:
                row.append(None)
                continue
            row.append(
                (
                    tuple(a + b for a, b in zip(key, contribution)),
                    ref[0],
                    ref[1],
                )
            )
        best.append(row)

        # fold this target's nodes into the prefix-min structure for later targets
        merged_positions = []
        merged_keys = []
        merged_refs = []
        i = j = 0
        own = [
            (positions[p], row[p][0], (t_index, p))
            for p in range(len(positions))
            if row[p] is not None
        ]
        while i < len(seen_positions) or j < len(own):
            take_own = j < len(own) and (
                i >= len(seen_positions) or own[j][0] <= seen_positions[i]
            )
            if take_own:
                merged_positions.append(own[j][0])
                merged_keys.append(own[j][1])
                merged_refs.append(own[j][2])
                j += 1
            else:
                merged_positions.append(seen_positions[i])
                merged_keys.append(seen_keys[i])
                merged_refs.append(seen_refs[i])
                i += 1
        # running prefix minimum, so a lookup is one bisect
        running_key = None
        running_ref = None
        for k in range(len(merged_positions)):
            if running_key is None or merged_keys[k] < running_key:
                running_key = merged_keys[k]
                running_ref = merged_refs[k]
            merged_keys[k] = running_key
            merged_refs[k] = running_ref
        seen_positions, seen_keys, seen_refs = (
            merged_positions,
            merged_keys,
            merged_refs,
        )

    # finish: the last chunk runs from the final position + 1 to contig_length
    final_key = (0, 0, 0, 0)
    final_ref = (None, None)
    for t_index, row in enumerate(best):
        for p_index, entry in enumerate(row):
            if entry is None:
                continue
            if contig_length - candidates[t_index][p_index] < minimum_span:
                continue
            if entry[0] < final_key:
                final_key = entry[0]
                final_ref = (t_index, p_index)

    chosen = []
    t_index, p_index = final_ref
    while t_index is not None:
        position = candidates[t_index][p_index]
        chosen.append((t_index, position, costs[t_index][p_index]))
        _, prev_t, prev_p = best[t_index][p_index]
        t_index, p_index = prev_t, prev_p
    chosen.reverse()
    return chosen


def select_cut_points(
    bam_filename,
    chrom,
    contig_length,
    strand="",
    gtf=None,
    segment_span=None,
    wiggle=None,
    depth_window=DEFAULT_DEPTH_WINDOW,
    grid_origin=DEFAULT_GRID_ORIGIN,
    margin=None,
    minimum_span=None,
    max_intron_length=None,
    annotation=None,
    gtf_index_cache_dir=None,
    count_denominator=True,
    severed_sink=None,
    # The thresholds the quant step will apply, for the emitted set only.  Not
    # config reads: --HiFi changes min_per_id inside LRAA, so the caller decides.
    min_per_id=None,
    min_mapping_quality=None,
    strandless=False,
    # A severed multi-exon alignment's cost, against 1 for a monoexonic one.
    # None reads the configured default.
    multiexon_weight=None,
    # Increasing search radii in bases; the last rung is always the wiggle
    # half-width. A parameter so a test can pin the ladder.
    expansion_radii=None,
):
    """Choose cut points for one contig-strand. Returns a ``Selection``.

    ``segment_span`` and ``wiggle`` are in BASES here; the megabase values live
    in the config and are converted by the CLI, so this function has one unit.

    ``strandless=True`` declares that the bam holds BOTH orientations and the
    cuts must serve both. It does not change the arithmetic -- a falsy ``strand``
    already unions the annotation islands and costs both orientations -- it
    changes what the run says about itself: the segments are labelled without a
    strand suffix, the denominator is reported per orientation, and a bam that
    turns out to hold only one orientation is reported rather than passed off as
    a strandless selection. Keeping it a declaration rather than a mode is
    deliberate: two code paths that must agree on cut coordinates are two chances
    to disagree.

    ``multiexon_weight`` prices a severed spliced alignment against a monoexonic
    one, default ``chunk_severed_multiexon_weight``. Severing is a cost and never
    a veto: see the module docstring for why a hard zero-severed rule was rejected.

    The search is PROGRESSIVE. ``wiggle`` is the maximum window width, not the
    width searched: rungs from ``expansion_radii`` widen around the target and each
    rung reads only the annulus the previous rung did not cover, stopping as soon
    as an annotation-compliant position severs nothing. Every reported quantity
    except ``grid_positions_examined`` is exact regardless of where the search
    stopped, because it can only stop on a zero cost and zero is the floor: a
    compliant zero inside radius R makes the whole-window minimum zero too, and
    the distance tie-break then puts the winner inside R.
    """

    if strandless and strand:
        raise SelectionError(
            "a strandless selection cannot also filter by orientation, got "
            "strand {!r}. A strandless cut serves BOTH strands; filtering to one "
            "would cost it against half the alignments it will actually "
            "sever.".format(strand)
        )

    if segment_span is None:
        segment_span = int(LRAA_Globals.config["approx_MB_per_cut"] * MB)
    if wiggle is None:
        # The absolute configured default, NOT a function of segment_span. The
        # config comment on the key holds the measurement that says why.
        wiggle = int(LRAA_Globals.config["approx_MB_per_cut_wiggle_window"] * MB)
    if multiexon_weight is None:
        multiexon_weight = DEFAULT_MULTIEXON_WEIGHT
    if margin is None:
        margin = extractor.DEFAULT_MARGIN
    if wiggle < 0:
        raise SelectionError("wiggle window must be >= 0")
    if multiexon_weight < 1:
        raise SelectionError(
            "multiexon weight must be >= 1: a severed spliced alignment cannot "
            "cost less than a monoexonic one, and 0 would make it free"
        )
    if minimum_span is None:
        minimum_span = minimum_span_for(segment_span, depth_window)

    if annotation is None:
        annotation = (
            extractor.load_gtf_for_contig(
                gtf, chrom, strand, cache_dir=gtf_index_cache_dir
            )
            if gtf
            else extractor.Annotation()
        )
    islands = extractor.find_islands(annotation, chrom, strand, margin)
    zones = extractor.cut_zones(islands, 1, contig_length)

    targets, tail_merged = cut_targets(contig_length, segment_span, minimum_span)

    half = wiggle // 2
    windows = []
    for target in targets:
        window_lend = max(1, target - half)
        window_rend = min(contig_length - 1, target + half)
        windows.append((window_lend, window_rend))

    rungs = expansion_rungs(half, expansion_radii)

    # Costed with no BAM access at all, so the two window-wide figures a report
    # quotes stay exact even when the search stops after one rung: how many grid
    # positions the MAXIMUM window holds, and how many of them the annotation
    # forbids. Both are arithmetic over the grid and the admissible zones.
    all_grid_counts = []
    blocked_counts = []
    for window_lend, window_rend in windows:
        on_grid = grid_positions(window_lend, window_rend, depth_window, grid_origin)
        all_grid_counts.append(len(on_grid))
        blocked_counts.append(len(_blocked_by_annotation(on_grid, zones)[1]))

    # Per target: what has been scanned, and the best of it.
    searches = [
        {
            "rung": 0,  # index of the NEXT rung to run
            "radius": None,  # radius already scanned, None before the first rung
            "examined": 0,  # grid positions actually costed
            "scored": {},  # compliant position -> weighted cost
            "best": None,  # cheapest compliant cost seen
            "unconstrained": None,  # cheapest cost seen, annotation ignored
        }
        for _ in targets
    ]

    with _open_bam(bam_filename, "rb") as bam:

        def run_rung(t_index):
            """Cost one more annulus around a target. Returns nothing."""

            state = searches[t_index]
            target = targets[t_index]
            window_lend, window_rend = windows[t_index]
            outer = rungs[state["rung"]]
            state["rung"] += 1
            for lend, rend in annulus_intervals(
                target, state["radius"], outer, window_lend, window_rend
            ):
                on_grid = grid_positions(lend, rend, depth_window, grid_origin)
                if not on_grid:
                    continue
                state["examined"] += len(on_grid)
                # One fetch per annulus interval, covering only that interval.
                # Re-scanning from the target at every rung would cost 1.76x a
                # flat scan of the whole window: the regression this avoids.
                grid_costs = spanning_counts(
                    bam,
                    chrom,
                    strand,
                    on_grid,
                    max_intron_length,
                    multiexon_weight,
                )
                cheapest = min(grid_costs)
                if state["unconstrained"] is None or cheapest < state["unconstrained"]:
                    state["unconstrained"] = cheapest
                compliant = set(_blocked_by_annotation(on_grid, zones)[0])
                for position, cost in zip(on_grid, grid_costs):
                    if position not in compliant:
                        continue
                    state["scored"][position] = cost
                    if state["best"] is None or cost < state["best"]:
                        state["best"] = cost
            state["radius"] = outer

        def exhausted(t_index):
            return searches[t_index]["rung"] >= len(rungs)

        def own_best_key(t_index):
            """The key this target would choose on its own, or None if it cannot."""

            state = searches[t_index]
            if not state["scored"]:
                return None
            target = targets[t_index]
            return min(
                (cost, abs(position - target), position)
                for position, cost in state["scored"].items()
            )

        for t_index in range(len(targets)):
            # Widen until a compliant position severs nothing, or until the
            # maximum. A target whose window is entirely annotation-blocked has no
            # cost at all, so it also searches to the maximum -- which is right:
            # the block may lift further out.
            while not exhausted(t_index) and searches[t_index]["best"] != 0:
                run_rung(t_index)

        while True:
            candidates = []
            costs = []
            for state in searches:
                ordered = sorted(state["scored"].items())
                candidates.append([position for position, _ in ordered])
                costs.append([cost for _, cost in ordered])

            chosen = _solve(candidates, costs, targets, minimum_span, contig_length)
            placed = {t_index: (position, cost) for t_index, position, cost in chosen}

            # A target the JOINT constraints forced off its own optimum is the one
            # case where a wider search can still help: the anti-sliver floor may
            # be satisfiable by a position the truncated search never looked at.
            # Everything else is already the answer a full scan would give.
            compromised = []
            for t_index in range(len(targets)):
                if exhausted(t_index):
                    continue
                if t_index not in placed:
                    compromised.append(t_index)
                    continue
                position, cost = placed[t_index]
                key = (cost, abs(position - targets[t_index]), position)
                if key != own_best_key(t_index):
                    compromised.append(t_index)
            if not compromised:
                break
            for t_index in compromised:
                run_rung(t_index)

        # The severed set per placed position, gathered once: the read NAMES the
        # extractor will drop, and the monoexonic/multi-exon split that says what
        # the cut actually cost. Same predicate as the scoring, so the number
        # reported against a cut is the number of names listed for it.
        dropped_read_names = collections.OrderedDict()  # name -> [cut positions]
        cut_positions = [placed[t_index][0] for t_index in sorted(placed)]
        profiles = {}
        for position in cut_positions:
            monoexonic = multiexon = 0
            for aln in spanning_alignments(
                bam, chrom, strand, position, max_intron_length
            ):
                dropped_read_names.setdefault(aln.query_name, []).append(position)
                if is_multiexon(aln):
                    multiexon += 1
                else:
                    monoexonic += 1
            profiles[position] = (monoexonic, multiexon)

        cuts = []
        unplaced = []
        for t_index, target in enumerate(targets):
            window_lend, window_rend = windows[t_index]
            state = searches[t_index]
            if t_index not in placed:
                best_here = state["best"]
                declined_annotation = not state["scored"]
                if declined_annotation:
                    # Short on purpose: this line is emitted once per declined
                    # target and a contig can decline many, so it carries the FACTS
                    # and the summary carries the rationale once.
                    reason = (
                        "DECLINED for the annotation: all {} grid position(s) in "
                        "the window (searched to radius {} bp, {} bp grid) fall "
                        "inside an annotated locus".format(
                            all_grid_counts[t_index],
                            state["radius"],
                            depth_window,
                        )
                    )
                else:
                    reason = (
                        "{} compliant position(s) exist but none can be used without "
                        "leaving a chunk shorter than the {} bp minimum span".format(
                            len(state["scored"]), minimum_span
                        )
                    )
                unplaced.append(
                    UnplacedTarget(
                        index=t_index,
                        target=target,
                        window_lend=window_lend,
                        window_rend=window_rend,
                        grid_positions=all_grid_counts[t_index],
                        annotation_blocked=blocked_counts[t_index],
                        reason=reason,
                        best_spanning=best_here,
                        declined_annotation=declined_annotation,
                        search_radius=state["radius"],
                    )
                )
                continue
            position, weight = placed[t_index]
            reference = state["unconstrained"]
            monoexonic, multiexon = profiles[position]
            cuts.append(
                CutChoice(
                    index=t_index,
                    target=target,
                    position=position,
                    offset=position - target,
                    # The alignment COUNT, which is what the extractor drops and
                    # what the denominator counts. The weighted figure is beside
                    # it rather than in place of it: they are different questions.
                    spanning_dropped=monoexonic + multiexon,
                    compromised=reference is not None and weight > reference,
                    window_lend=window_lend,
                    window_rend=window_rend,
                    grid_positions=all_grid_counts[t_index],
                    annotation_blocked=blocked_counts[t_index],
                    unconstrained_best_spanning=reference,
                    severed_monoexonic=monoexonic,
                    severed_multiexon=multiexon,
                    severed_weight=weight,
                    search_radius=state["radius"],
                )
            )

        for cut in cuts if severed_sink is not None else ():
            # Emission is stricter: only alignments quantification would use.  A
            # read rejected on MAPQ or identity never reaches a component, so
            # including it would let a consumer attribute a dissolved component to
            # a read no downstream stage ever saw.
            for aln in spanning_alignments(
                bam,
                chrom,
                strand,
                cut.position,
                max_intron_length,
                quant_only=True,
                min_per_id=min_per_id,
                min_mapping_quality=min_mapping_quality,
            ):
                # An alignment spanning several cuts is fetched once per cut, so
                # something has to collapse the repeats.  Not identity: two BAM rows
                # can be byte-identical and still be two retained alignments, so
                # hashing the record -- however completely -- destroys multiplicity
                # and the emitted BAM stops representing every alignment the cuts
                # sever.
                #
                # Ownership instead.  Each alignment is emitted only at the FIRST
                # selected cut inside its span, which is a property of coordinates
                # rather than content: repeats collapse because the later cuts are
                # not the owner, while two identical rows are both yielded at the
                # owning cut and both survive.
                owner = cut_positions[
                    bisect.bisect_left(cut_positions, aln.reference_start + 1)
                ]
                if cut.position == owner:
                    severed_sink.append(aln)

        if count_denominator:
            (
                total_retained,
                retained_forward,
                retained_reverse,
            ) = count_retained_primary_by_orientation(
                bam, chrom, strand, max_intron_length
            )
        else:
            total_retained = retained_forward = retained_reverse = None

    if (
        strandless
        and total_retained
        and not (retained_forward and retained_reverse)
    ):
        # Reported, not refused: the cuts are still valid coordinates. What is
        # not true is the claim the flag makes -- these cuts were costed against
        # one orientation, so they are per-strand cuts wearing a strandless
        # label, and the chunks they produce will have an empty partner after the
        # in-chunk split. Silent here would look exactly like success.
        print(
            "WARNING: --strandless was given but {} holds only one orientation "
            "over {} ({} forward, {} reverse retained primary alignment(s)): this "
            "looks like an already strand-separated bam. Strandless cuts must be "
            "chosen over the RAW bam, with the strand split moved inside the "
            "chunk.".format(
                bam_filename, chrom, retained_forward, retained_reverse
            ),
            file=sys.stderr,
        )

    segments = []
    start = 1
    for cut in cuts:
        segments.append(
            Segment(
                index=len(segments),
                lend=start,
                rend=cut.position,
                span=cut.position - start + 1,
            )
        )
        start = cut.position + 1
    segments.append(
        Segment(
            index=len(segments),
            lend=start,
            rend=contig_length,
            span=contig_length - start + 1,
        )
    )

    return Selection(
        chrom=chrom,
        strand=strand,
        contig_length=contig_length,
        depth_window=depth_window,
        grid_origin=grid_origin,
        segment_span=segment_span,
        wiggle=wiggle,
        minimum_span=minimum_span,
        targets=targets,
        cuts=cuts,
        segments=segments,
        unplaced=unplaced,
        tail_merged=tail_merged,
        total_retained_primary=total_retained,
        total_dropped=sum(cut.spanning_dropped for cut in cuts),
        dropped_read_names=dropped_read_names,
        strandless=strandless,
        retained_primary_forward=retained_forward,
        retained_primary_reverse=retained_reverse,
        multiexon_weight=multiexon_weight,
    )


def selection_to_dict(selection):
    """JSON-shaped manifest. Every count carries its denominator."""

    return {
        "chrom": selection.chrom,
        # null when there is no orientation to report, never "" and never the
        # string "None": these region strings are consumed by the extractor.
        "strand": selection.strand or None,
        "strandless": selection.strandless,
        "contig_length": selection.contig_length,
        "params": {
            "segment_span": selection.segment_span,
            # The MAXIMUM window, which is not what was searched: each cut carries
            # the radius its own search stopped at.
            "wiggle_window": selection.wiggle,
            "depth_window": selection.depth_window,
            "grid_origin": selection.grid_origin,
            "minimum_span": selection.minimum_span,
            "severed_multiexon_weight": selection.multiexon_weight,
        },
        "targets": selection.targets,
        "tail_merged_targets": selection.tail_merged,
        "cuts": [
            {
                "target": cut.target,
                "position": cut.position,
                "offset_from_target": cut.offset,
                "spanning_alignments_dropped": cut.spanning_dropped,
                # What that count is MADE OF, and what it cost the objective.
                # Severing is expected rather than exceptional, so a bare total
                # would hide the only distinction that matters between two cuts
                # dropping the same number of reads.
                "severed_monoexonic": cut.severed_monoexonic,
                "severed_multiexon": cut.severed_multiexon,
                "severed_weighted_cost": cut.severed_weight,
                "search_radius": cut.search_radius,
                "hard_constraint_compromise": cut.compromised,
                "window": [cut.window_lend, cut.window_rend],
                "grid_positions_in_window": cut.grid_positions,
                "grid_positions_blocked_by_annotation": cut.annotation_blocked,
                "best_spanning_ignoring_annotation": cut.unconstrained_best_spanning,
                # what to pass extract_contig_region_inputs.py, and what the
                # chunk starting here must pass the normalizer
                "window_origin_for_next_chunk": cut.position,
            }
            for cut in selection.cuts
        ],
        "segments": [
            {
                "region": "{}{}:{}-{}".format(
                    selection.chrom,
                    selection.strand or "",
                    segment.lend,
                    segment.rend,
                ),
                "lend": segment.lend,
                "rend": segment.rend,
                "span": segment.span,
                "window_origin": segment.lend - 1,
            }
            for segment in selection.segments
        ],
        "unplaced_targets": [
            {
                "target": item.target,
                "window": [item.window_lend, item.window_rend],
                "grid_positions_in_window": item.grid_positions,
                "grid_positions_blocked_by_annotation": item.annotation_blocked,
                # Whether the ANNOTATION is what left this window unusable, which
                # is the only remaining reason a target can be declined, and how
                # dirty the cheapest position on offer was.
                "declined_annotation": item.declined_annotation,
                "best_spanning_in_window": item.best_spanning,
                "search_radius": item.search_radius,
                "reason": item.reason,
            }
            for item in selection.unplaced
        ],
        "counts": {
            "targets": len(selection.targets),
            "cuts_placed": len(selection.cuts),
            "targets_unplaced": len(selection.unplaced),
            "targets_declined_annotation": sum(
                1 for item in selection.unplaced if item.declined_annotation
            ),
            "targets_tail_merged": len(selection.tail_merged),
            "segments": len(selection.segments),
            "retained_primary_alignments": selection.total_retained_primary,
            "retained_primary_forward": selection.retained_primary_forward,
            "retained_primary_reverse": selection.retained_primary_reverse,
            "alignments_dropped_at_cuts": selection.total_dropped,
            "alignments_dropped_monoexonic": sum(
                cut.severed_monoexonic for cut in selection.cuts
            ),
            "alignments_dropped_multiexon": sum(
                cut.severed_multiexon for cut in selection.cuts
            ),
            "severed_weighted_cost_at_cuts": sum(
                cut.severed_weight for cut in selection.cuts
            ),
            "distinct_dropped_read_names": len(selection.dropped_read_names),
        },
    }


def _as_selections(selections):
    return [selections] if isinstance(selections, Selection) else list(selections)


def write_dropped_read_names(selections, path):
    """One read name per line, sorted, deduplicated across every selection.

    Plain names because that is what ``samtools view -N`` and the equivalent
    filters consume: the point is to build a pruned baseline BAM from exactly
    this set, so the file has to be directly usable as the filter input.

    Deduplicated because a read severed by a cut is dropped by BOTH neighbouring
    chunks -- it is contained by neither -- so the same name is reported twice by
    per-chunk accounting and must appear once in a filter set.
    """

    every = set()
    for selection in _as_selections(selections):
        every.update(selection.dropped_read_names)
    with open(path, "wt") as ofh:
        for name in sorted(every):
            print(name, file=ofh)
    return path


def write_dropped_read_detail(selections, path):
    """Attribution for each dropped name: which cut or cuts severed it."""

    with open(path, "wt") as ofh:
        print("read_name\tchrom\tstrand\tcut_positions", file=ofh)
        for selection in _as_selections(selections):
            for name in sorted(selection.dropped_read_names):
                print(
                    "{}\t{}\t{}\t{}".format(
                        name,
                        selection.chrom,
                        selection.strand or ".",
                        ",".join(
                            str(p) for p in selection.dropped_read_names[name]
                        ),
                    ),
                    file=ofh,
                )
    return path


def format_report(selection):
    """Human-readable report. Reports what happened; tunes nothing."""

    lines = []
    label = "{}{}".format(selection.chrom, selection.strand or "")
    if selection.strandless:
        lines.append(
            "# STRANDLESS selection over {}: one set of cuts serving BOTH "
            "strands. Islands are the union of both strands' loci and each "
            "position is costed by the alignments it severs in either "
            "orientation.".format(label)
        )
    lines.append(
        "# {} length {} bp; segment span {} bp; MAXIMUM wiggle window {} bp "
        "(target +/- {}, searched progressively); depth_window {} bp; grid origin "
        "{}; minimum span {} bp; severed multi-exon weight {}".format(
            label,
            selection.contig_length,
            selection.segment_span,
            selection.wiggle,
            selection.wiggle // 2,
            selection.depth_window,
            selection.grid_origin,
            selection.minimum_span,
            selection.multiexon_weight,
        )
    )
    lines.append(
        "# severing is a COST, never a veto: a cut takes the cheapest position its "
        "search reaches and the alignments it severs are dropped, counted and "
        "named. A monoexonic severed alignment costs 1, a multi-exon one {}. Only "
        "the ANNOTATION can decline a target.".format(selection.multiexon_weight)
    )
    declined = [item for item in selection.unplaced if item.declined_annotation]
    lines.append(
        "# {} target(s), {} cut(s) placed, {} unplaced ({} declined for "
        "annotation), {} tail-merged -> {} segment(s)".format(
            len(selection.targets),
            len(selection.cuts),
            len(selection.unplaced),
            len(declined),
            len(selection.tail_merged),
            len(selection.segments),
        )
    )
    if selection.tail_merged:
        lines.append(
            "# tail-merged target(s) {}: the residual would have been shorter than "
            "the {} bp minimum span".format(
                ", ".join(str(t) for t in selection.tail_merged), selection.minimum_span
            )
        )

    lines.append(
        "#target\tposition\toffset\tdropped\tmono\tmultiexon\tweighted\tradius\t"
        "grid_in_window\tannot_blocked\tbest_ignoring_annot\tcompromised"
    )
    for cut in selection.cuts:
        lines.append(
            "{}\t{}\t{:+d}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                cut.target,
                cut.position,
                cut.offset,
                cut.spanning_dropped,
                cut.severed_monoexonic,
                cut.severed_multiexon,
                cut.severed_weight,
                cut.search_radius,
                cut.grid_positions,
                cut.annotation_blocked,
                cut.unconstrained_best_spanning,
                "yes" if cut.compromised else "no",
            )
        )

    for item in selection.unplaced:
        lines.append(
            "# UNPLACED target {} (window {}-{}): {}".format(
                item.target, item.window_lend, item.window_rend, item.reason
            )
        )

    lines.append("#segment\tregion\tspan\twindow_origin")
    for segment in selection.segments:
        lines.append(
            "{}\t{}{}:{}-{}\t{}\t{}".format(
                segment.index,
                selection.chrom,
                selection.strand or "",
                segment.lend,
                segment.rend,
                segment.span,
                segment.lend - 1,
            )
        )

    spans = [segment.span for segment in selection.segments]
    if spans:
        lines.append(
            "# realised spans: min {} bp, max {} bp, mean {:.0f} bp; nominal {} bp "
            "(max deviation {:+.1f}%)".format(
                min(spans),
                max(spans),
                sum(spans) / len(spans),
                selection.segment_span,
                100.0
                * max(
                    (span - selection.segment_span) / selection.segment_span
                    for span in spans
                ),
            )
        )
    denominator = selection.total_retained_primary
    lines.append(
        "# dropped {} alignment(s) at {} cut(s), {} distinct read name(s), out of "
        "{} retained primary alignment(s){}".format(
            selection.total_dropped,
            len(selection.cuts),
            len(selection.dropped_read_names),
            denominator if denominator is not None else "not counted",
            ""
            if not denominator
            else " ({:.4f}%)".format(100.0 * selection.total_dropped / denominator),
        )
    )
    lines.append(
        "# of those, {} monoexonic and {} multi-exon, weighted cost {}".format(
            sum(cut.severed_monoexonic for cut in selection.cuts),
            sum(cut.severed_multiexon for cut in selection.cuts),
            sum(cut.severed_weight for cut in selection.cuts),
        )
    )
    if selection.retained_primary_forward is not None:
        lines.append(
            "# denominator by orientation: {} forward, {} reverse".format(
                selection.retained_primary_forward,
                selection.retained_primary_reverse,
            )
        )
    return "\n".join(lines)


def main(argv=None):

    parser = argparse.ArgumentParser(
        description="choose cut points splitting a contig-strand into chunks: "
        "targets at multiples of --approx_MB_per_cut, each searched outward from "
        "the target up to a MAXIMUM window of --approx_MB_per_cut_wiggle_window, "
        "for a depth-window-aligned position outside every annotated locus whose "
        "severed retained primary alignments cost the least",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bam", type=str, required=True, help="bam aligned reads")
    parser.add_argument(
        "--genome_fa",
        type=str,
        required=False,
        help="genome fasta, for contig lengths. Defaults to the bam header.",
    )
    parser.add_argument("--gtf", type=str, required=False, help="gtf annotation")
    parser.add_argument(
        "--gtf_index_cache_dir",
        type=str,
        default=None,
        help="where to put the tabix index of --gtf when the GTF's own "
        "directory is not writable. The GTF itself is never modified.",
    )
    parser.add_argument(
        "--contig",
        type=str,
        default=None,
        help="restrict to one contig. Default: every contig in the bam header.",
    )
    parser.add_argument(
        "--strand",
        type=str,
        default=None,
        choices=("+", "-"),
        help="orientation to select for. Omit for an already strand-split bam, "
        "in which case every record counts.",
    )
    parser.add_argument(
        "--strandless",
        action="store_true",
        help="the bam holds BOTH orientations and one set of cuts must serve "
        "both: islands are the union of both strands' loci and each position is "
        "costed by what it severs in either orientation. Point --bam at the RAW "
        "bam and move the strand split inside the chunk, AFTER extraction -- "
        "separate_bam_by_strand.py rewrites is_reverse, so any strand filter "
        "applied to a raw bam reads pre-flip flags. Incompatible with --strand.",
    )
    parser.add_argument(
        "--approx_MB_per_cut",
        type=float,
        default=LRAA_Globals.config["approx_MB_per_cut"],
        help="approximate MEGABASES between cut points",
    )
    parser.add_argument(
        "--approx_MB_per_cut_wiggle_window",
        type=float,
        default=LRAA_Globals.config["approx_MB_per_cut_wiggle_window"],
        help="MAXIMUM total width in MEGABASES of the window centred on each "
        "target, not the half-width. The search starts small and widens to this, "
        "stopping early on a position that severs nothing. ABSOLUTE: it is not "
        "derived from --approx_MB_per_cut, because the distance to a read-free "
        "position is a property of the genome and the library rather than of the "
        "chunk size",
    )
    parser.add_argument(
        "--depth_window",
        type=int,
        default=DEFAULT_DEPTH_WINDOW,
        help="resolution in bases at which read depth is measured; cuts land on "
        "this grid so no depth window spans one. Must match "
        "normalize_bam_by_strand.py --depth_window.",
    )
    parser.add_argument(
        "--grid_origin",
        type=int,
        default=DEFAULT_GRID_ORIGIN,
        help="0-based reference coordinate of a depth-window boundary; only its "
        "remainder modulo --depth_window matters. The pipeline pins this to 0.",
    )
    parser.add_argument(
        "--margin",
        type=int,
        default=extractor.DEFAULT_MARGIN,
        help="clearance demanded either side of a cut from annotated loci, in bp. "
        "Must match extract_contig_region_inputs.py --margin.",
    )
    parser.add_argument(
        "--minimum_span",
        type=int,
        default=None,
        help="shortest permitted chunk, in bases. Default: half the segment span, "
        "floored to a depth_window multiple.",
    )
    parser.add_argument(
        "--severed_multiexon_weight",
        type=int,
        default=DEFAULT_MULTIEXON_WEIGHT,
        help="what a severed MULTI-EXON alignment costs, against 1 for a "
        "monoexonic one. Severing is a cost to minimise, never a veto: a spliced "
        "alignment carries the junction evidence the splice graph's edges come "
        "from, a monoexonic one carries depth only. 1 makes the cost a plain "
        "alignment count.",
    )
    parser.add_argument(
        "--max_intron_length",
        type=int,
        default=LRAA_Globals.config["max_intron_length"],
        help="alignments carrying a longer intron are excluded, matching the "
        "strand split. 0 disables.",
    )
    parser.add_argument(
        "--min_per_id",
        type=float,
        default=LRAA_Globals.config["min_per_id"],
        help="percent identity the quant step will require. Affects only which "
        "alignments --severed_reads_bam emits, never cut placement. Must match the "
        "run being selected for: --HiFi raises this to 97 inside LRAA, so leaving "
        "the default of 80 would emit reads that run discards.",
    )
    parser.add_argument(
        "--min_mapping_quality",
        type=int,
        default=int(LRAA_Globals.config["min_mapping_quality"]),
        help="mapping quality the quant step will require. Emission only, as with "
        "--min_per_id.",
    )
    parser.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="writes <prefix>.cuts.json, <prefix>.cuts.txt, "
        "<prefix>.dropped_reads.txt and <prefix>.dropped_reads.tsv",
    )
    parser.add_argument(
        "--severed_reads_bam",
        type=str,
        required=False,
        default=None,
        help="also write the severed primary alignments themselves to this BAM, "
        "coordinate-sorted and indexed. The names in <prefix>.dropped_reads.tsv "
        "cannot be fetched from a coordinate-indexed BAM, and a span cannot answer "
        "what a read was compatible with -- compatibility follows exon blocks and "
        "junctions -- so a consumer asking which severed reads linked two genes "
        "needs the records.",
    )

    args = parser.parse_args(argv)

    segment_span = int(args.approx_MB_per_cut * MB)
    wiggle = int(args.approx_MB_per_cut_wiggle_window * MB)
    if args.strandless and args.strand:
        parser.error(
            "--strandless and --strand are incompatible: a strandless cut serves "
            "both strands, so costing it against one orientation would price it "
            "against half the alignments it severs"
        )
    # Thread argparse's None through rather than coercing it, so the retention
    # predicates receive the value their contract is written against. The report
    # sites that interpolate a strand supply their own placeholder.
    strand = args.strand

    if args.genome_fa:
        with pysam.FastaFile(args.genome_fa) as fasta:
            lengths = dict(zip(fasta.references, fasta.lengths))
    else:
        with pysam.AlignmentFile(args.bam, "rb") as bam:
            lengths = dict(zip(bam.references, bam.lengths))

    # A REFERENCE CONTIG THE BAM DOES NOT CARRY IS SKIPPED, NOT FATAL.
    #
    # The contig list comes from the genome fasta when one is given, but every
    # candidate is priced by FETCHING it from the bam, and pysam raises
    # `ValueError: invalid contig` on a name absent from the bam header. A reference
    # carrying sequences the alignment does not is the normal case rather than the
    # exceptional one -- decoys, EBV, alt scaffolds, anything the aligner was not
    # given or the reads never reached. GRCh38_no_alt.fa has 195 sequences and a
    # 10x SBX PBMC bam aligned against the same assembly has 194: chrEBV, holding
    # zero reads and zero gencode records, killed a whole-genome chunked run at 246 s.
    #
    # Such a contig cannot be partitioned in any case: with no alignments there is no
    # depth to place a cut against and nothing for a chunk to contain. Skipping is the
    # answer rather than refusing, and the count is reported so a reference/bam
    # mismatch big enough to matter is still visible.
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        bam_contigs = set(bam.references)
    absent = [c for c in sorted(lengths) if c not in bam_contigs]
    if absent:
        shown = ", ".join(absent[:5])
        print(
            "NOTE: {} of {} reference sequence(s) are absent from the bam header and "
            "are skipped, having no alignments to partition: {}{}".format(
                len(absent),
                len(lengths),
                shown,
                "" if len(absent) <= 5 else ", ...",
            ),
            file=sys.stderr,
        )
        lengths = {c: n for c, n in lengths.items() if c in bam_contigs}

    if args.contig:
        if args.contig not in lengths:
            raise SelectionError(
                "contig {} is absent{}; known: {}".format(
                    args.contig,
                    " from the bam" if args.contig in bam_contigs or absent else "",
                    ", ".join(sorted(lengths)),
                )
            )
        contigs = [args.contig]
    else:
        contigs = sorted(lengths)

    selections = []
    # One sink for every contig, written once below.  Per-contig writing would
    # leave only the last contig's severed reads on disk.
    severed_sink = [] if args.severed_reads_bam else None
    for chrom in contigs:
        annotation = (
            extractor.load_gtf_for_contig(
                args.gtf, chrom, strand, cache_dir=args.gtf_index_cache_dir
            )
            if args.gtf
            else extractor.Annotation()
        )
        selections.append(
            select_cut_points(
                bam_filename=args.bam,
                chrom=chrom,
                contig_length=lengths[chrom],
                strand=strand,
                segment_span=segment_span,
                wiggle=wiggle,
                depth_window=args.depth_window,
                grid_origin=args.grid_origin,
                margin=args.margin,
                minimum_span=args.minimum_span,
                max_intron_length=args.max_intron_length,
                annotation=annotation,
                severed_sink=severed_sink,
                min_per_id=args.min_per_id,
                min_mapping_quality=args.min_mapping_quality,
                strandless=args.strandless,
                multiexon_weight=args.severed_multiexon_weight,
            )
        )

    if severed_sink is not None:
        with pysam.AlignmentFile(args.bam, "rb") as bam:
            header = bam.header
        written = write_severed_alignments_bam(
            header,
            severed_sink,
            args.severed_reads_bam,
            min_per_id=args.min_per_id,
            min_mapping_quality=args.min_mapping_quality,
        )
        print(
            "-wrote {} severed alignment(s) to {}".format(
                written, args.severed_reads_bam
            )
        )

    report = "\n".join(format_report(selection) for selection in selections)
    print(report)

    with open("{}.cuts.txt".format(args.output_prefix), "wt") as ofh:
        print(report, file=ofh)
    with open("{}.cuts.json".format(args.output_prefix), "wt") as ofh:
        json.dump(
            [selection_to_dict(selection) for selection in selections],
            ofh,
            indent=2,
            sort_keys=True,
        )
        print("", file=ofh)

    write_dropped_read_names(
        selections, "{}.dropped_reads.txt".format(args.output_prefix)
    )
    write_dropped_read_detail(
        selections, "{}.dropped_reads.tsv".format(args.output_prefix)
    )

    unplaced = sum(len(selection.unplaced) for selection in selections)
    declined = sum(
        1
        for selection in selections
        for item in selection.unplaced
        if item.declined_annotation
    )
    if unplaced:
        print(
            "WARNING: {} target(s) could not be placed, {} of them DECLINED "
            "because the annotation left no admissible position in the window; "
            "see the report. These are reported rather than skipped "
            "silently.".format(unplaced, declined),
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
