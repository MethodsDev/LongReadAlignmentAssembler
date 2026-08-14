#!/usr/bin/env python3

"""One CPU budget, DIVIDED across a flat queue of work units. Never multiplied.

What this replaces
------------------
LRAA used to expose three independent parallelism knobs that multiplied rather than
divided: ``--num_parallel_contigs`` (contig/strand workers), ``--concurrency`` on the
chunked-quant orchestrator, and ``--num_threads_per_worker`` (native tool threads,
which also silently doubled as the component-fork count). Nothing checked that their
product fit the host, and the product was the only thing that mattered.

There is now one authoritative number -- the budget -- and every level is derived from
it. The budget is a HARD CAP, not a utilisation target: when the work does not
decompose finely enough to spend it, the honest answer is idle cores, reported as such,
not an allocation contrived to make the number look spent.

The queue is FLAT
-----------------
Two deployment regimes exist and no per-level split serves both:

* large genomes -- human, mouse -- are a few very large contigs, so chunks supply
  essentially all the units and contig-level parallelism supplies almost none;
* small-contig references -- SIRVs, MORFs -- are many contigs each below one chunk
  span, so contig-level parallelism supplies all the units and chunking supplies
  nothing.

A real genome is both at once: chr1 at 248 Mb beside chrM at 16 kb. So the levels are
FLATTENED into a single queue of ``(contig, strand, chunk)`` units and the scheduler
never asks which level a unit came from. One rule then covers both regimes and the
mixed case.

The queue LRAA builds today is one unit per (contig, strand): a whole-genome or
small-contig-reference invocation therefore has many units, and a per-chromosome
sharded invocation has exactly two. Emitting several chunk units per contig-strand is
prepared for -- ``WorkUnit`` carries the region, and ``artifact_token`` keys per-unit
outputs so chunks of one contig-strand cannot overwrite each other -- but it is not
merely a switch to flip. Two integration questions are open and are recorded in the
module notes at the bottom of this file.

Allocation is PHASE-AWARE
-------------------------
Cores that can run assembly work and threads that can only help a native tool step are
not interchangeable, so the two phases are sized by different rules:

* ``allocate_prep`` sizes the serial prep phase -- strand split, index, read count.
  These are samtools invocations that genuinely scale with threads and have nothing
  else competing with them, so the whole budget is spent here and the remainder is
  handed out one core at a time: 16 over 3 steps is 6, 5, 5 and not 5, 5, 5.

* ``allocate`` sizes the main per-unit phase, where the wall time actually is. A unit
  is single-core Python assembly and quant work; it CANNOT consume a second core for
  that work. So the remainder is deliberately NOT distributed here. Twelve units on a
  budget of sixteen means twelve workers and four idle cores, and
  ``budget_shortfall_note`` says so in the log rather than pretending otherwise.

* Component forking is DEMAND-DRIVEN, never reserved: a unit that reaches an eligible
  multipath-graph component draws workers from whatever cores are free at that instant
  and runs in-process when none are. See ``CoreLease``. Reserving would strand cores
  while the queue is full and forfeit the tail once it drains.
"""

import logging
import multiprocessing
import os
import sys
from collections import namedtuple

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import Util_funcs


def default_budget():
    """The budget to use when the user names none.

    ``Util_funcs.available_cpus()`` is affinity-aware, so a cpuset -- every container
    runtime, Slurm ``--cpus-per-task``, taskset -- is respected without the user having
    to restate it on the command line.
    """
    return Util_funcs.available_cpus()


class WorkUnit(namedtuple("WorkUnit", "contig_acc contig_strand chunk_index region cost")):
    """One schedulable piece of LRAA work, and the only element type of the flat queue.

    ``region`` is ``None`` for a unit spanning a whole contig-strand, or an inclusive
    1-based ``(lend, rend)`` for a chunk of one. ``cost`` is a cheap scheduling proxy
    only -- retained read count where available, contig span otherwise -- and is never
    used for anything but launch order.
    """

    __slots__ = ()

    def __new__(cls, contig_acc, contig_strand, chunk_index=0, region=None, cost=0):
        if region is not None:
            lend, rend = int(region[0]), int(region[1])
            if lend >= rend:
                raise ValueError(
                    "WorkUnit region must satisfy lend < rend, got {}-{}".format(
                        lend, rend
                    )
                )
            region = (lend, rend)
        return super().__new__(
            cls,
            str(contig_acc),
            str(contig_strand),
            int(chunk_index),
            region,
            float(cost),
        )

    @property
    def artifact_token(self):
        """Filename stem component that is UNIQUE per unit.

        Every per-unit artifact -- ``quant.expr``, ``quant.tracking``, ``gtf``, the
        ``.ok`` resume checkpoint, the error log, the resource trace -- is keyed on this.

        A whole contig-strand unit keeps the historical ``{contig}.{strand}`` stem, so a
        run resuming from a pre-flat-queue run still matches its own checkpoints. A
        region-restricted unit appends its chunk index and bounds, which is what stops
        two chunks of one contig-strand from overwriting each other's outputs. That also
        fixes a live collision on ``--region``: two different ``--region`` runs sharing an
        ``--output_prefix`` previously resumed from each other's checkpoint, because the
        region was absent from every artifact name.
        """
        stem = "{}.{}".format(self.contig_acc, self.contig_strand)
        if self.region is None:
            return stem
        return "{}.chunk{}_{}-{}".format(
            stem, self.chunk_index, self.region[0], self.region[1]
        )


def order_longest_first(units):
    """Deterministic longest-first launch order for a flat queue.

    Unit cost varies enormously across the queue -- a 10 Mb chunk of chr1 against a
    10 kb SIRV contig -- and span does not bound read count, so a small
    highly-expressed contig can be expensive too. Launching in construction order
    leaves the most expensive unit for the tail of the run; longest-first is the
    standard makespan heuristic and costs nothing given a proxy already to hand.

    This only changes wall time when the queue is longer than the worker pool, so it
    pays off on direct whole-genome and multi-contig invocations, including the
    small-contig SIRV and MORF regime. It is inert on a per-chromosome sharded
    invocation, which has two units.

    Unlike the unseeded shuffle it replaces, the order is deterministic: ties break on
    the unit's own identity, so two runs of the same input launch in the same order.
    """
    return sorted(
        units,
        key=lambda u: (-u.cost, u.contig_acc, u.contig_strand, u.chunk_index),
    )


Allocation = namedtuple(
    "Allocation",
    "budget units unit_workers tool_threads component_workers unallocated_cores",
)


def allocate(budget, num_units, max_unit_workers=None, max_component_workers=None):
    """Size the MAIN per-unit phase from the budget and the flat queue length.

    ``unit_workers`` is how many units run at once, capped by the budget because a unit
    needs a core of its own, and by ``max_unit_workers`` for callers that must serialise
    (``--no_parallelize_contigs``, or ``--debug``/``--include_prelim_transcripts``, whose
    debug artifacts clobber each other across concurrent workers).

    ``tool_threads`` is what each worker may pass to a NATIVE tool step it invokes --
    samtools ``-@``, minimap2 ``-t``. It is one core for the worker's own single-threaded
    work plus its share of whatever the unit count could not claim. It is never zero:
    ``unit_workers <= budget`` by construction, so the floor division cannot round to 0.

    ``component_workers`` is a CEILING on what one unit may ask for when it reaches an
    eligible component, NOT a reservation. Nothing is set aside for it: the actual grant
    is whatever ``CoreLease`` finds free at the moment of the request, which is zero while
    the queue is full. The default ceiling is the spare budget unit parallelism could not
    claim, since that is the most that can ever be free, plus the cores that finished
    units hand back as the queue drains.

    ``unallocated_cores`` is the part of the budget no unit worker can use. It is REPORTED
    rather than allocated away, because a main-phase unit cannot consume a second core for
    its own assembly or quant work -- though a demand-driven component grant can, which is
    exactly the case this leaves room for.
    """
    budget = max(1, int(budget))
    units = max(1, int(num_units))

    unit_workers = min(budget, units)
    if max_unit_workers is not None:
        unit_workers = max(1, min(unit_workers, int(max_unit_workers)))

    tool_threads = budget // unit_workers

    component_workers = budget if max_component_workers is None else int(max_component_workers)
    component_workers = max(0, min(component_workers, budget))
    # A grant of 1 is not parallelism -- it builds the manager and its result transport to
    # run one component serially -- so a ceiling that cannot reach 2 is off.
    if component_workers < 2:
        component_workers = 0

    return Allocation(
        budget=budget,
        units=units,
        unit_workers=unit_workers,
        tool_threads=tool_threads,
        component_workers=component_workers,
        unallocated_cores=budget - unit_workers * tool_threads,
    )


class CoreLease:
    """The budget's live count of free cores, shared across unit-worker processes.

    This is what makes component-level parallelism DEMAND-DRIVEN instead of reserved.
    A static reservation is wrong in both directions: it strands cores while the queue is
    full, and it forfeits the tail when the queue drains. So nothing is reserved. Each
    unit worker holds one permit for its own single-core work, and a unit that reaches an
    eligible multipath-graph component asks for more and is granted whatever is free at
    that instant -- which is zero while the queue is full, and the spare budget once the
    queue has drained. Granted zero means assembling in-process, exactly today's default,
    so this is additive rather than a second code path.

    That the grant arrives precisely when the queue is drained is the point: that is the
    same moment a large component is what bounds the node's latency. Aggregate CPU share
    does not bound makespan -- component reconstruction measured 29.7 CPU-seconds of
    1,862, about 1.6% of run CPU, yet parallelising inside the TAIL unit cuts node
    latency by however long that component takes.

    Three properties this has to hold, all of them easy to get wrong:

    * COMPONENT WORKERS COUNT AGAINST THE BUDGET. They are children of a unit worker, so
      if they did not hold permits the hard cap would be a fiction.
    * NO PREEMPTION, so brief oversubscription is ACCEPTED and reported. Component work
      is not checkpointable, so an outstanding grant cannot be recalled when the ready
      queue refills. The alternative -- granting only when the queue has slack -- needs
      the worker to see the parent's queue depth, which it cannot. So a unit that finds
      no permit free still runs (it would otherwise deadlock behind its own children),
      and says so.
    * GRANTS RELEASE ON ABNORMAL TERMINATION. Every acquisition is paired with a
      ``finally``; a leaked permit would silently shrink the node's capacity for the rest
      of a long run.
    """

    def __init__(self, budget):
        self.budget = max(1, int(budget))
        # Bounded, so a double release raises here rather than inflating the budget.
        self._permits = multiprocessing.BoundedSemaphore(self.budget)

    def take(self, count=1):
        """Acquire up to ``count`` permits without blocking; return how many were taken."""
        taken = 0
        for _ in range(max(0, int(count))):
            if not self._permits.acquire(block=False):
                break
            taken += 1
        return taken

    def give_back(self, count):
        """Return ``count`` permits. Safe to call with 0, and safe from a ``finally``."""
        for _ in range(max(0, int(count))):
            try:
                self._permits.release()
            except ValueError:
                # BoundedSemaphore refuses a release past its ceiling. Reaching here means
                # a bookkeeping bug, not a condition to paper over silently.
                logging.getLogger(__name__).error(
                    "core lease released more permits than it holds; budget accounting is wrong"
                )
                break


def component_worker_grant(lease, ceiling, eligible_components, log=None):
    """Grant component workers from whatever is free NOW; return (granted, release).

    ``eligible_components`` is REQUIRED and is the reason most requests are declined.
    Forking runs DIFFERENT components of one unit concurrently; it does NOT split a single
    component across workers. So a tail unit whose remaining work is one enormous
    component gains nothing from a grant of any width -- the component is indivisible
    either way -- and a grant is worth making only when the unit has a BACKLOG of several
    eligible components. A request that cannot use its grant is declined here rather than
    charged to the budget.

    The ceiling is therefore also capped at the backlog: more workers than components is
    process creation, payload serialisation and shard fsync (measured at roughly 0.18 ms
    CPU and 2.1 ms wall per published unit) charged to the parent whose latency the grant
    was supposed to reduce.

    ``granted`` is 0 or at least 2, and ``release`` MUST be called from a ``finally``:
    a permit leaked on abnormal termination silently shrinks the node's capacity for the
    rest of the run.
    """
    usable = min(int(ceiling or 0), int(eligible_components or 0))
    if lease is None or usable < 2:
        if log is not None and eligible_components == 1:
            log(
                "only one eligible component outstanding, which forking cannot subdivide; "
                "assembling in-process"
            )
        return 0, (lambda: None)

    granted = lease.take(usable)
    if granted < 2:
        lease.give_back(granted)
        if log is not None:
            log(
                "no spare cores free for component workers (the unit queue still holds "
                "the budget); assembling in-process"
            )
        return 0, (lambda: None)

    released = [False]

    def release():
        if not released[0]:
            released[0] = True
            lease.give_back(granted)

    return granted, release


PrepAllocation = namedtuple("PrepAllocation", "budget steps concurrency threads")


def allocate_prep(budget, num_steps):
    """Size the SERIAL PREP phase: strand split, index, read count, chunk extraction.

    Unlike a main-phase unit, a prep step IS a native tool invocation, so a second
    thread is genuinely consumable and the remainder is distributed rather than floored:
    a budget of 16 over 12 steps gives four steps 2 threads and eight steps 1, summing
    to exactly 16. Concurrency is capped at the budget; a caller with more steps than
    that runs them in waves of ``concurrency``.
    """
    budget = max(1, int(budget))
    steps = max(1, int(num_steps))
    concurrency = min(budget, steps)
    base, extra = divmod(budget, concurrency)
    threads = [base + 1] * extra + [base] * (concurrency - extra)
    return PrepAllocation(
        budget=budget, steps=steps, concurrency=concurrency, threads=threads
    )


ALLOCATION_LOG_PREFIX = "CPU allocation:"


def format_allocation(allocation, phase="main"):
    """The single startup line naming units, workers and threads.

    A user who gets unexpected wall time needs to see the division without instrumenting
    anything, so this is one greppable line with every derived number on it.
    """
    return (
        "{} budget={} cores; phase={}; units={}; unit_workers={}; "
        "tool_threads_per_worker={}; component_worker_ceiling={}; "
        "unallocated_cores={}".format(
            ALLOCATION_LOG_PREFIX,
            allocation.budget,
            phase,
            allocation.units,
            allocation.unit_workers,
            allocation.tool_threads,
            allocation.component_workers,
            allocation.unallocated_cores,
        )
    )


SHORTFALL_LOG_PREFIX = "CPU budget note:"


def budget_shortfall_note(allocation):
    """One line stating the part of the budget unit parallelism cannot spend, or None.

    This exists so the flat queue's central limitation is never dressed up as
    utilisation. One chromosome, both strands, budget 16 is TWO units: the other 14
    cores are reachable only by native tool steps inside those two workers, which are a
    small part of wall time. Only more units can put them on assembly work -- and see
    THE HOT-UNIT BOUND below for when more units are not obtainable at any budget.
    """
    if allocation.unit_workers >= allocation.budget:
        return None
    spare = allocation.budget - allocation.unit_workers
    tool_reachable = allocation.unit_workers * (allocation.tool_threads - 1)
    return (
        "{} {} concurrent unit worker(s) cannot spend a budget of {}; {} core(s) sit "
        "beyond unit parallelism -- {} reachable only by native tool steps inside a "
        "worker (samtools, minimap2), {} not reachable at all. Assembly and quant work "
        "is single-core per unit, so only MORE UNITS (finer contig/strand/chunk "
        "partitioning) can put these cores on LRAA work.".format(
            SHORTFALL_LOG_PREFIX,
            allocation.unit_workers,
            allocation.budget,
            spare,
            tool_reachable,
            allocation.unallocated_cores,
        )
    )


# ---------------------------------------------------------------------------------
# THE HOT-UNIT BOUND: the queue length is an INPUT, never a target
# ---------------------------------------------------------------------------------
# `allocate` takes however many units exist and never asks for more. That is deliberate,
# and it is the documented behaviour for the case where a budget-driven refinement WANTS
# to split a unit and cannot.
#
# Chunk boundaries may never fall inside an annotated locus, so a single very long or
# very highly expressed locus defines a unit that is INDIVISIBLE at any budget. On chr20
# the largest such island holds 168 transcripts on plus and 156 on minus, so no uniform
# target below about 168 genes per chunk is exactly achievable there at all.
#
# When refinement cannot reach the budget, REFUSING and accepting fewer units is
# correct. An indivisible hot unit bounds the makespan regardless: with single-core
# units the makespan is bounded below by the longest unit, so cores that the hot unit
# leaves idle could not have shortened the run even if something were scheduled on them.
#
# Two kinds of work CAN legitimately fill cores that drain at the tail, and both reduce
# wall time rather than merely occupying cores. Neither is implemented here; both are
# recorded so the next person does not reach for the wrong tool:
#   * pipeline the EPILOGUE per finished unit rather than after a barrier -- output
#     sorting and compression, per-unit QC, merge preparation, and coordinate
#     back-translation if chunks are ever materialised. Overlapping that with a
#     straggler is free makespan.
#   * PREFETCH per-unit input preparation for queued units. Chunk extraction costs about
#     30 s of serial prep per chromosome and is the Amdahl floor of the measured chunked
#     run; doing it on otherwise-idle cores takes it off the critical path.
#
# Ruled out, so nobody re-litigates them: dynamically re-splitting a RUNNING unit (its
# state is not checkpointable) and speculative duplicate execution (the workload is
# deterministic, so a duplicate is pure waste).
#
# When tail idle is structurally large -- a small chromosome on a sixteen-core node --
# the fix is the DEPLOYMENT SHAPE, packing chr21 and chr22 or a batch of small contigs
# onto one node, not a scheduler trick.


# ---------------------------------------------------------------------------------
# OPEN QUESTIONS for emitting several chunk units per contig-strand
# ---------------------------------------------------------------------------------
# The scheduler and the artifact keying are ready for it. These two are not, and both
# could change the shape of what the unit generator should emit, so neither is a switch:
#
# 1. REGION FILTER SEMANTICS ARE ASYMMETRIC, and the read side carries an OFF-BY-ONE.
#    Both are LATENT today: the cut rule forbids a boundary inside an annotated locus, and
#    on chr20 that left zero of 319,698 alignments severed, so neither filter has anything
#    to act on. They surface on data where a cut MUST sever, which the cutter permits when
#    its wiggle window holds no clean position.
#
#    TRANSCRIPT SIDE IS CONTAINMENT. Transcript.parse_GTF_to_Transcripts filters per
#    FEATURE at pylib/Transcript.py:913 -- `if lend < lend_restrict or rend > rend_restrict:
#    continue` -- so a locus straddling a boundary is NOT double-counted; it is TRUNCATED on
#    both sides, each chunk emitting a mutilated model from the exons on its own side.
#    Containment protects against double counting only while no boundary falls inside a
#    transcript span.
#
#    READ SIDE IS OVERLAP AND DROPS THE REGION'S FIRST BASE.
#    pylib/Bam_alignment_extractor.py:77-79 calls
#        self._pysam_reader.fetch(contig_acc, region_lend, region_rend)
#    where region_lend/region_rend are 1-BASED INCLUSIVE (straight from
#    `--region contig:lend-rend`, carried as restrict_region_lend/rend) while pysam.fetch
#    takes 0-BASED HALF-OPEN [start, stop). So fetch(c, L, R) yields 1-based [L+1, R].
#      - The error is on the LEFT EDGE ONLY. 1-based base L sits at 0-based L-1, below
#        start=L, so it is excluded. The right edge is correct by accident: 1-based R sits
#        at 0-based R-1 < R, so it is included.
#      - EXPOSED BY: an alignment whose entire overlap with the region is that one base,
#        i.e. one ENDING EXACTLY AT region_lend. Demonstrated with a two-record BAM --
#        read A at 1-based 100..200, read B at 201..300; fetch(c, 200, 300) for 1-based
#        [200,300] returns ['B'] where overlap semantics require ['A','B'].
#      - USER-REACHABLE TODAY, not only under chunking: `LRAA --region chr20:L-R` loses any
#        alignment ending exactly at L. Invisible at whole-contig scope, where region_lend
#        is None and the unrestricted fetch(contig_acc) branch is taken.
#      - The fix is `fetch(contig_acc, region_lend - 1, region_rend)`, but it moves quant
#        output on any region-restricted run, so it belongs with boundary accounting and a
#        regression arm rather than here.
#
#    The hazard at a boundary is LOSS, not duplication: overlap does not double-count in
#    practice, because assignment is gated by transcript compatibility and the side holding
#    no compatible transcript drops the read. The converse is real -- at ~10 genes per chunk
#    one read of 82,884 was assigned by NEITHER chunk, since the chunk holding its only
#    compatible transcript could not represent the part of the alignment outside its region.
#    So the invariant a chunked queue needs is READ COMPLETENESS -- no alignment span may
#    cross a boundary -- not coordinate ownership. A 5'-anchor single-owner rule was
#    implemented elsewhere and refuted by measurement: it drops legitimately assigned reads.
#
# 2. MEMORY FOOTPRINT TRADES OFF AGAINST SIMPLICITY, and the win does not extrapolate.
#    Region restriction inside LRAA needs no extraction and no coordinate translation,
#    because restrict_region_lend/rend is already plumbed and absolute -- by far the
#    cheaper route. But the contig SEQUENCE is fetched and upper-cased per contig, so
#    every worker still holds the WHOLE contig FASTA and the memory win is lost.
#    The measured peak-RSS win for chunking came from region-scoped MATERIALISED inputs:
#    0.79x, 2.42 GB SIX WIDE against a 3.08 GB single-process baseline, largest chunk
#    532 MB. Quote it at that concurrency or not at all. It relies on workers not hitting
#    their per-phase memory peaks simultaneously, and homogeneous units launched together
#    will do exactly that; 16-wide is ESTIMATED at 6.5 to 8.5 GB, not measured. The
#    532 MB largest-chunk figure is likewise a property of chr20's expression profile
#    rather than a bound, since an indivisible hot locus sets the per-worker peak.
#    Materialised chunks keep whatever win there is and cost extraction plus translation.
#
# Also inherited by any global merge: TPM is renormalised to the sum of all_reads over
# whatever table is merged, so per-chunk TPM must be recomputed after a global merge.
# num_total_reads is already global regardless of --region and needs nothing.
