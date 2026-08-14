#!/usr/bin/env python3

"""One budget, DIVIDED across a flat queue of work units. Never multiplied.

The bug these guard against is arithmetic, not plumbing: three independent knobs whose
PRODUCT decided how many cores a run took, with nothing checking that the product fit
the host. So every assertion here is against the allocation functions directly. Nothing
in this file spawns a process to check arithmetic.

The two shapes named throughout are the two supported deployment modes:

* CHROMOSOME-SHARD MODE -- the WDL shards per chromosome, so a node gets one contig and
  two strands on about sixteen cores. Two units. Units have to come from chunks here.
* SINGLE-NODE MODE -- a small reference (SIRVs, MORFs) run whole on one node, so every
  contig and both strands are in one invocation. Many units of wildly varying cost.
"""

import re

import pytest

import CpuBudget
import Util_funcs


# --------------------------------------------------------------- the budget itself


def test_default_budget_is_the_affinity_aware_cpu_count(monkeypatch):
    """The default must be what the process may RUN on, not what the machine has.

    Restating the budget on every command line is exactly what a cpuset already knows,
    so the default reads it instead of asking.
    """
    monkeypatch.setattr(Util_funcs.os, "sched_getaffinity", lambda pid: {0, 1, 2, 3},
                        raising=False)
    monkeypatch.setattr(Util_funcs.os, "cpu_count", lambda: 64)
    assert CpuBudget.default_budget() == 4
    assert CpuBudget.default_budget() == Util_funcs.available_cpus()


def test_explicit_budget_is_respected_and_is_a_hard_cap():
    """A named budget is obeyed exactly: workers times threads never exceeds it."""
    for budget in range(1, 33):
        for units in (1, 2, 3, 7, 12, 25, 50, 200):
            alloc = CpuBudget.allocate(budget=budget, num_units=units)
            assert alloc.budget == budget
            assert alloc.unit_workers * alloc.tool_threads <= budget, (budget, units)


# ------------------------------------------------------------------- the two shapes


def test_chromosome_shard_shape_allocates_two_workers_with_eight_tool_threads():
    """CHROMOSOME-SHARD MODE: one contig, two strands, budget 16."""
    alloc = CpuBudget.allocate(budget=16, num_units=2)
    assert alloc.unit_workers == 2
    assert alloc.tool_threads == 8
    assert alloc.unallocated_cores == 0


def test_chromosome_shard_shape_reports_that_fourteen_cores_cannot_do_unit_work():
    """The eight tool threads are NOT sixteen-core utilisation, and must not read as it.

    Native tool threads help samtools and minimap2 steps, which are a small part of wall
    time; the assembly and quant work is single-core per unit. With component forking off
    the run is two busy processes whatever the budget says, so the note has to state that
    only more units can change it.
    """
    alloc = CpuBudget.allocate(budget=16, num_units=2)
    note = CpuBudget.budget_shortfall_note(alloc)
    assert note is not None
    assert "14 core(s) sit beyond unit parallelism" in note
    assert "MORE UNITS" in note


def test_chunked_shape_allocates_twelve_workers_with_one_thread_each():
    """Twelve units on a budget of sixteen: twelve workers, one thread apiece."""
    alloc = CpuBudget.allocate(budget=16, num_units=12)
    assert alloc.unit_workers == 12
    assert alloc.tool_threads == 1


def test_main_phase_does_not_distribute_the_remainder():
    """Four cores idle at twelve units on sixteen, reported rather than allocated away.

    Handing four of the twelve units a second thread would be a lie about the main phase:
    a unit is single-core assembly and quant work and cannot consume one. The budget is a
    hard cap, not a utilisation target.
    """
    alloc = CpuBudget.allocate(budget=16, num_units=12)
    assert alloc.unallocated_cores == 4
    assert alloc.tool_threads == 1
    note = CpuBudget.budget_shortfall_note(alloc)
    assert "0 reachable only by native tool steps" in note
    assert "4 not reachable at all" in note


def test_single_node_shape_uses_every_core_it_has_units_for():
    """SINGLE-NODE MODE: 7 SIRV contigs x 2 strands on 16 cores is 14 busy workers."""
    alloc = CpuBudget.allocate(budget=16, num_units=14)
    assert alloc.unit_workers == 14
    assert alloc.tool_threads == 1
    assert alloc.unallocated_cores == 2


# ------------------------------------------------------------------ the floor at one


@pytest.mark.parametrize("budget,units", [(4, 12), (1, 12), (1, 1), (2, 50), (3, 1000)])
def test_a_budget_below_the_unit_count_still_gives_every_worker_a_thread(budget, units):
    """Never zero threads: floor division cannot round down below one.

    unit_workers is clamped to the budget, so budget // unit_workers is at least 1 by
    construction. A worker allocated 0 native-tool threads would hand `-@ 0` to samtools.
    """
    alloc = CpuBudget.allocate(budget=budget, num_units=units)
    assert alloc.unit_workers == budget
    assert alloc.tool_threads >= 1


def test_budget_and_unit_count_are_floored_at_one_not_rejected():
    """Zero or negative inputs come from arithmetic upstream, not from a user; clamp."""
    alloc = CpuBudget.allocate(budget=0, num_units=0)
    assert alloc.budget == 1 and alloc.unit_workers == 1 and alloc.tool_threads == 1


def test_serialising_hands_the_whole_budget_to_the_single_worker():
    """--no_parallelize_contigs, --debug and --include_prelim_transcripts cap at one
    worker; that worker's native tool steps then get the entire budget."""
    alloc = CpuBudget.allocate(budget=16, num_units=12, max_unit_workers=1)
    assert alloc.unit_workers == 1
    assert alloc.tool_threads == 16
    assert alloc.unallocated_cores == 0


def test_the_shortfall_accounting_is_exact():
    """Tool-reachable plus unreachable must equal the budget the units could not claim,
    or the note is reporting a number that does not add up."""
    for budget in range(1, 33):
        for units in (1, 2, 5, 12, 40):
            alloc = CpuBudget.allocate(budget=budget, num_units=units)
            spare = alloc.budget - alloc.unit_workers
            tool_reachable = alloc.unit_workers * (alloc.tool_threads - 1)
            assert tool_reachable + alloc.unallocated_cores == spare, (budget, units)


def test_no_shortfall_note_when_the_units_claim_the_whole_budget():
    alloc = CpuBudget.allocate(budget=4, num_units=4)
    assert alloc.unit_workers == 4
    assert CpuBudget.budget_shortfall_note(alloc) is None


# --------------------------------------------------------------- component forking


def test_the_component_ceiling_reserves_nothing():
    """The ceiling is what a unit MAY ask for, never cores set aside for it.

    Proof that nothing is reserved: the unit-worker count and the tool-thread count are
    identical whether the ceiling is the whole budget or zero.
    """
    with_ceiling = CpuBudget.allocate(budget=16, num_units=2)
    without = CpuBudget.allocate(budget=16, num_units=2, max_component_workers=0)
    assert with_ceiling.component_workers == 16
    assert without.component_workers == 0
    assert with_ceiling.unit_workers == without.unit_workers
    assert with_ceiling.tool_threads == without.tool_threads
    assert with_ceiling.unallocated_cores == without.unallocated_cores


def test_a_ceiling_that_cannot_reach_two_is_reported_as_off():
    """One fork is not parallelism: it builds the manager and its result transport to run
    one component serially."""
    assert CpuBudget.allocate(16, 2, max_component_workers=1).component_workers == 0
    assert CpuBudget.allocate(16, 2, max_component_workers=0).component_workers == 0


def test_a_grant_is_declined_when_only_one_component_is_eligible():
    """Forking runs DIFFERENT components concurrently; it cannot subdivide one.

    So a tail unit whose remaining work is one enormous component gains nothing from a
    grant of any width, and charging the budget for it would be pure loss.
    """
    lease = CpuBudget.CoreLease(16)
    granted, release = CpuBudget.component_worker_grant(lease, ceiling=8, eligible_components=1)
    assert granted == 0
    release()
    # the refusal cost the budget nothing
    assert lease.take(16) == 16


def test_a_grant_is_capped_by_the_outstanding_eligible_component_count():
    """More workers than components is process creation, serialisation and shard fsync
    charged to the parent whose latency the grant was meant to reduce."""
    lease = CpuBudget.CoreLease(16)
    granted, release = CpuBudget.component_worker_grant(lease, ceiling=16, eligible_components=3)
    assert granted == 3
    release()


def test_a_grant_takes_only_what_is_free_and_gives_it_back():
    """Zero when the queue holds the budget, the spare once units finish, and the permits
    return so a later unit can be granted them in turn."""
    lease = CpuBudget.CoreLease(4)
    assert lease.take(4) == 4  # four unit workers hold the whole budget
    granted, release = CpuBudget.component_worker_grant(lease, ceiling=4, eligible_components=4)
    assert granted == 0, "nothing may be granted while the queue holds every core"

    lease.give_back(3)  # three units finish
    granted, release = CpuBudget.component_worker_grant(lease, ceiling=4, eligible_components=4)
    assert granted == 3, "the tail unit may now use the drained cores"
    release()
    assert lease.take(3) == 3, "the grant returned its permits"


def test_a_grant_never_exceeds_the_budget_even_across_several_units():
    """Component workers are children of a unit worker, so if they did not count against
    the budget the hard cap would be a fiction."""
    lease = CpuBudget.CoreLease(8)
    held = lease.take(2)  # two unit workers
    first, release_first = CpuBudget.component_worker_grant(lease, 8, 8)
    second, release_second = CpuBudget.component_worker_grant(lease, 8, 8)
    assert held + first + second <= 8
    release_first()
    release_second()
    assert lease.take(6) == 6


def test_releasing_a_grant_twice_does_not_inflate_the_budget():
    """The release runs from a finally on paths that can also reach the normal exit, so it
    has to be idempotent -- otherwise a double release would manufacture cores."""
    lease = CpuBudget.CoreLease(4)
    granted, release = CpuBudget.component_worker_grant(lease, 4, 4)
    assert granted == 4
    release()
    release()
    assert lease.take(4) == 4
    assert lease.take(1) == 0


def test_no_lease_means_no_forking():
    """A caller that never built a lease -- a direct library user, a test -- must not fork."""
    granted, release = CpuBudget.component_worker_grant(None, ceiling=8, eligible_components=8)
    assert granted == 0
    release()


# ------------------------------------------------------------- the prep phase only


def test_prep_phase_distributes_the_remainder_because_threads_are_consumable_there():
    """Prep is native tool invocations -- strand split, index, read count -- so a second
    thread genuinely helps and the remainder is handed out rather than floored.

    This is the one place remainder distribution is valid, and it is why prep and the
    main per-unit phase are sized by different functions.
    """
    prep = CpuBudget.allocate_prep(budget=16, num_steps=12)
    assert sum(prep.threads) == 16
    assert prep.threads.count(2) == 4
    assert prep.threads.count(1) == 8


def test_prep_phase_gives_a_single_serial_step_the_whole_budget():
    prep = CpuBudget.allocate_prep(budget=16, num_steps=1)
    assert prep.concurrency == 1
    assert prep.threads == [16]


def test_prep_phase_never_exceeds_the_budget():
    for budget in range(1, 33):
        for steps in (1, 2, 3, 12, 40):
            prep = CpuBudget.allocate_prep(budget, steps)
            assert sum(prep.threads) == budget
            assert prep.concurrency == min(budget, steps)
            assert all(t >= 1 for t in prep.threads)


# --------------------------------------------------------------- the startup line


def test_startup_line_states_units_workers_and_threads_on_one_line():
    """A user who gets unexpected wall time must see the division without instrumenting
    anything, so the whole allocation is one greppable line."""
    alloc = CpuBudget.allocate(budget=16, num_units=2)
    line = CpuBudget.format_allocation(alloc, phase="main")
    assert "\n" not in line
    assert line == (
        "CPU allocation: budget=16 cores; phase=main; units=2; unit_workers=2; "
        "tool_threads_per_worker=8; component_worker_ceiling=16; unallocated_cores=0"
    )


def test_startup_line_carries_every_derived_number():
    alloc = CpuBudget.allocate(budget=16, num_units=12, max_component_workers=2)
    line = CpuBudget.format_allocation(alloc, phase="chunks")
    for field, value in (
        ("budget", alloc.budget),
        ("units", alloc.units),
        ("unit_workers", alloc.unit_workers),
        ("tool_threads_per_worker", alloc.tool_threads),
        ("component_worker_ceiling", alloc.component_workers),
        ("unallocated_cores", alloc.unallocated_cores),
    ):
        assert re.search(r"\b{}={}\b".format(field, value), line), (field, line)
    assert "phase=chunks" in line


# ------------------------------------------------------- the flat queue and ordering


def test_a_whole_contig_strand_unit_keeps_the_historical_artifact_stem():
    """Resume from a pre-flat-queue run has to keep matching its own checkpoints."""
    unit = CpuBudget.WorkUnit("chr20", "+")
    assert unit.artifact_token == "chr20.+"


def test_chunk_units_of_one_contig_strand_get_distinct_artifact_stems():
    """The collision this prevents is total: quant.expr, quant.tracking, gtf, the .ok
    resume checkpoint, the error log and the resource trace are all keyed on this stem,
    so two chunks of one contig-strand would overwrite every one of them."""
    first = CpuBudget.WorkUnit("chr20", "+", chunk_index=0, region=(1, 10_000_000))
    second = CpuBudget.WorkUnit("chr20", "+", chunk_index=1, region=(10_000_001, 20_000_000))
    assert first.artifact_token != second.artifact_token
    assert first.artifact_token == "chr20.+.chunk0_1-10000000"
    assert second.artifact_token == "chr20.+.chunk1_10000001-20000000"


def test_two_regions_of_one_contig_strand_cannot_share_a_stem():
    """This also closes a live --region collision: two --region runs sharing an
    --output_prefix previously resumed from each other's checkpoint."""
    a = CpuBudget.WorkUnit("chr20", "+", region=(1, 1000))
    b = CpuBudget.WorkUnit("chr20", "+", region=(2000, 3000))
    assert a.artifact_token != b.artifact_token


def test_an_inverted_region_is_rejected_rather_than_silently_keyed():
    with pytest.raises(ValueError):
        CpuBudget.WorkUnit("chr20", "+", region=(1000, 1000))
    with pytest.raises(ValueError):
        CpuBudget.WorkUnit("chr20", "+", region=(2000, 1000))


def test_longest_first_puts_the_most_expensive_unit_first():
    """Unit cost varies enormously across the flat queue and span does not bound read
    count, so a small highly expressed contig can outweigh a large quiet one."""
    chrM = CpuBudget.WorkUnit("chrM", "+", cost=500_000)
    chr1_chunk = CpuBudget.WorkUnit("chr1", "+", chunk_index=3,
                                    region=(30_000_001, 40_000_000), cost=120_000)
    sirv = CpuBudget.WorkUnit("SIRV1", "-", cost=800)
    ordered = CpuBudget.order_longest_first([sirv, chr1_chunk, chrM])
    assert [u.contig_acc for u in ordered] == ["chrM", "chr1", "SIRV1"]


def test_longest_first_is_deterministic_across_equal_costs():
    """It replaces an unseeded random.shuffle, so ties must not reintroduce run-to-run
    variation in which worker wins a race on any shared artifact."""
    units = [
        CpuBudget.WorkUnit("chr2", "-", cost=100),
        CpuBudget.WorkUnit("chr1", "+", cost=100),
        CpuBudget.WorkUnit("chr1", "-", cost=100),
        CpuBudget.WorkUnit("chr2", "+", cost=100),
    ]
    once = CpuBudget.order_longest_first(units)
    twice = CpuBudget.order_longest_first(list(reversed(units)))
    assert once == twice
    assert [(u.contig_acc, u.contig_strand) for u in once] == [
        ("chr1", "+"), ("chr1", "-"), ("chr2", "+"), ("chr2", "-")
    ]


def test_ordering_preserves_the_queue_it_was_given():
    units = [CpuBudget.WorkUnit("chr{}".format(i), "+", cost=i) for i in range(10)]
    assert sorted(CpuBudget.order_longest_first(units)) == sorted(units)


def test_units_are_hashable_so_a_scheduler_can_map_them_back_to_its_own_records():
    """Both callers order the units and then look up the job or chunk each one came
    from, which needs the unit to be a usable dict key."""
    unit = CpuBudget.WorkUnit("chr20", "+", region=(1, 1000), cost=5)
    assert {unit: "job"}[CpuBudget.WorkUnit("chr20", "+", region=(1, 1000), cost=5)] == "job"


# ------------------------------------------------------- levels divide, never multiply


def test_the_orchestrator_division_cannot_exceed_the_total():
    """The orchestrator runs `unit_workers` chunks at once and hands each chunk's LRAA
    `tool_threads` as its own --cpu_budget. A chunk is single-contig and strand-pure, so
    LRAA's queue inside it is one unit and its own pool clamps to 1 -- the share is not
    double-counted against that. Their product is what must not exceed the budget, and
    that is precisely what --concurrency times --num_parallel_contigs times
    --num_threads_per_worker used to be free to do.
    """
    for budget in range(1, 33):
        for chunks in (1, 2, 6, 12, 25, 60):
            alloc = CpuBudget.allocate(budget=budget, num_units=chunks)
            per_chunk_budget = alloc.tool_threads
            inner = CpuBudget.allocate(budget=per_chunk_budget, num_units=1)
            assert inner.unit_workers == 1
            assert alloc.unit_workers * per_chunk_budget <= budget, (budget, chunks)
            assert alloc.unit_workers * inner.unit_workers * inner.tool_threads <= budget
