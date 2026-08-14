#!/usr/bin/env python3

"""The reconciliation of shard publication with the CPU-budget core lease.

Two changes landed on the same code path and each carries an invariant the other
must not weaken.  Shard publication says a component worker's result exists only
as a file renamed into place, and a worker that exits 0 without publishing is a
lost result and must fail the run.  The CPU budget says a component worker is a
core charged against a hard cap, so the permit it holds must come back on every
exit path or the node quietly shrinks for the rest of a long run.

The interleaving is what makes this worth a test.  ``LRAA.reconstruct_isoforms``
grants permits before it forks and returns them after it joins, and the shard
accounting is enforced inside ``wait_for_remaining_processes`` -- which RAISES.
Draining outside the handler that releases the grant reads as correct and leaks a
permit every time a unit exits 0 without publishing.

Every test here drives the real ``reconstruct_isoforms`` against a stand-in graph
and a real ``CpuBudget.CoreLease``, and every test has a deadline: a parent that
blocks forever is one of the two defects shard publication exists to retire, so a
suite that could hang would hide the regression it is here to catch.
"""

import logging
import multiprocessing
import os
import signal
import sys
import time
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import CpuBudget
import LRAA_Globals
from LRAA import LRAA
from MultiProcessManager import WorkUnitAccountingError

# Generous next to the work each test does -- the longest sleeps 0.6s -- and far
# below any plausible hang.
TEST_DEADLINE_SEC = 60.0


@pytest.fixture(autouse=True)
def fail_rather_than_hang():
    """Bound every test, and leave no component worker behind if one fails."""

    def _expired(signum, frame):
        raise AssertionError(
            "test exceeded its {}s deadline -- a wedged parent is one of the "
            "failure modes under repair".format(TEST_DEADLINE_SEC)
        )

    previous_handler = signal.signal(signal.SIGALRM, _expired)
    signal.setitimer(signal.ITIMER_REAL, TEST_DEADLINE_SEC)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_handler)
        for child in multiprocessing.active_children():
            child.terminate()
            child.join(timeout=5)


# ---------------------------------------------------------------------------
# the smallest graph reconstruct_isoforms will accept
# ---------------------------------------------------------------------------


class _Plan:
    """What the stand-in worker for one component should do, and what to record.

    Carried on the component itself rather than in a module global so the
    component travels intact under either start method.
    """

    def __init__(self, role, lease, budget, sleep_sec=0.0, observed_free=None):
        self.role = role
        self.lease = lease
        self.budget = budget
        self.sleep_sec = sleep_sec
        self.observed_free = observed_free


class _FakeMpgn:
    """The only thing reconstruct_isoforms asks of a node before it forks."""

    def __init__(self, lend, rend, plan=None):
        self._lend = lend
        self._rend = rend
        self.plan = plan

    def get_coords(self):
        return [self._lend, self._rend]


class _FakeMultipathGraph:
    def __init__(self, components):
        self._components = components

    def define_disjoint_graph_components_via_shared_splice_graph_vertex(self):
        return list(self._components)

    def remove_small_components(self, components, min_transcript_length):
        return list(components)

    def write_mp_graph_nodes_to_gtf(self, filename):  # pragma: no cover
        raise AssertionError("these tests do not run with DEBUG artifacts on")


# ---------------------------------------------------------------------------
# the stand-in component worker, module scope so it survives either start method
# ---------------------------------------------------------------------------


def _component_worker(
    shard_store, mpg_component, component_counter, mpg_token, single_best_only
):
    plan = mpg_component[0].plan

    if plan.observed_free is not None:
        # Read the lease from INSIDE a live component worker.  Without this the
        # post-run permit count proves nothing: a grant that was never taken is
        # also a grant that never leaked.
        free = plan.lease.take(plan.budget)
        plan.lease.give_back(free)
        plan.observed_free.value = free

    if plan.sleep_sec:
        time.sleep(plan.sleep_sec)

    if plan.role == "publish":
        shard_store.publish(
            component_counter,
            ["tx-{}".format(component_counter)],
            token=mpg_token,
        )

    elif plan.role == "exit_clean_without_publishing":
        # The exactly-once violation: a payload that was never delivered, from a
        # process the OS reports as a clean success.
        return

    elif plan.role == "sigkill_mid_write":
        # Real bytes in the scratch file, then a signal that cannot be caught, so
        # nothing tidies up and the .tmp genuinely survives.  On the old queue
        # this left a truncated frame and the parent blocked in _recv_bytes for
        # the life of the run.
        tmp_path = shard_store.tmp_path(component_counter)
        with open(tmp_path, "wb") as fh:
            fh.write(b"\x80\x05\x95\xff\xff\xff\xff")
            fh.flush()
            os.fsync(fh.fileno())
            os.kill(os.getpid(), signal.SIGKILL)

    elif plan.role == "sigkill":
        os.kill(os.getpid(), signal.SIGKILL)

    else:  # pragma: no cover
        raise AssertionError("unknown role {!r}".format(plan.role))


# ---------------------------------------------------------------------------
# harness
# ---------------------------------------------------------------------------


def _free_permits(lease, budget):
    """How many of the budget's permits are free right now."""
    taken = lease.take(budget)
    lease.give_back(taken)
    return taken


def _shard_dirs(tmp_path):
    return sorted(p.name for p in Path(tmp_path).glob("__lraa_shards.*"))


def _build_round(
    tmp_path, monkeypatch, roles, budget, ceiling, sleeps=None, measure_index=None
):
    """A real LRAA instance whose components fork into `_component_worker`."""

    monkeypatch.setenv("LRAA_TMP_DIR", str(tmp_path))
    monkeypatch.setattr(LRAA_Globals, "DEBUG", False)
    monkeypatch.setitem(LRAA_Globals.config, "min_transcript_length", 0)
    # Every component is eligible, so every component forks and the in-process
    # branch stays out of the way of what is being measured.
    monkeypatch.setitem(LRAA_Globals.config, "min_mpgn_component_size_for_spawn", 1)

    lease = CpuBudget.CoreLease(budget)
    observed_free = multiprocessing.Value("i", -1)
    sleeps = sleeps or [0.0] * len(roles)

    components = []
    for index, role in enumerate(roles):
        plan = _Plan(
            role=role,
            lease=lease,
            budget=budget,
            sleep_sec=sleeps[index],
            observed_free=observed_free if index == measure_index else None,
        )
        lend = 1000 * index + 1
        components.append(
            [_FakeMpgn(lend, lend + 99, plan), _FakeMpgn(lend + 200, lend + 299)]
        )

    lraa = LRAA(splice_graph=None, component_workers=ceiling, core_lease=lease)
    lraa._multipath_graph = _FakeMultipathGraph(components)
    lraa._contig_acc = "chr1"
    lraa._contig_strand = "+"
    lraa._restrict_splice_type = "ME"
    lraa._reconstruct_isoforms_single_component = _component_worker

    return lraa, lease, observed_free


# ---------------------------------------------------------------------------
# both invariants on the success path
# ---------------------------------------------------------------------------


def test_a_clean_round_publishes_shards_holds_permits_and_returns_them(
    tmp_path, monkeypatch
):
    """The positive control: neither invariant is satisfied by doing nothing.

    Results come back in unit order although the workers finish in reverse, the
    grant was genuinely held while they ran, the permits are all free afterwards,
    and the shard directory is gone because the round succeeded.
    """

    lraa, lease, observed_free = _build_round(
        tmp_path,
        monkeypatch,
        roles=["publish", "publish", "publish"],
        budget=4,
        ceiling=3,
        sleeps=[0.6, 0.3, 0.0],
        measure_index=2,
    )

    transcripts = lraa.reconstruct_isoforms()

    assert transcripts == ["tx-1", "tx-2", "tx-3"]
    # Three of four permits were out while the workers ran.
    assert observed_free.value == 1
    assert _free_permits(lease, 4) == 4
    assert _shard_dirs(tmp_path) == []


# ---------------------------------------------------------------------------
# invariant 1: a worker killed mid-write fails loudly instead of wedging
# ---------------------------------------------------------------------------


def test_a_component_worker_sigkilled_mid_write_fails_loudly(
    tmp_path, monkeypatch, caplog
):
    """SIGKILL with bytes in the scratch file must fail the round, not hang it.

    The deadline on this test is the assertion that matters most: under the
    result queue the parent blocked in _recv_bytes on the truncated frame and
    never returned at all.
    """

    lraa, lease, _observed = _build_round(
        tmp_path,
        monkeypatch,
        roles=["publish", "sigkill_mid_write", "publish"],
        budget=4,
        ceiling=3,
    )

    with caplog.at_level(logging.ERROR):
        with pytest.raises(RuntimeError, match="1 component failures encountered"):
            lraa.reconstruct_isoforms()

    # The scratch file is surveyed, named and removed, and the survey is what
    # tells the operator a writer was killed rather than merely slow.
    assert "orphaned shard scratch file" in caplog.text
    assert "unit.00000002" in caplog.text

    # The shard directory is deliberately retained on a failed round: it sits
    # with the failing worker's evidence.
    assert _shard_dirs(tmp_path), "a failed round must not discard its shards"

    assert _free_permits(lease, 4) == 4


# ---------------------------------------------------------------------------
# invariant 2: exit 0 without publishing is still a failure, by name
# ---------------------------------------------------------------------------


def test_a_component_worker_exiting_zero_without_publishing_still_raises(
    tmp_path, monkeypatch
):
    """The exactly-once check survived the rebase, and still names the unit."""

    lraa, lease, _observed = _build_round(
        tmp_path,
        monkeypatch,
        roles=["publish", "exit_clean_without_publishing"],
        budget=4,
        ceiling=2,
    )

    with pytest.raises(WorkUnitAccountingError) as excinfo:
        lraa.reconstruct_isoforms()

    message = str(excinfo.value)
    assert "exited 0 but published no result shard" in message
    # The unit, its exit code and the absolute path that is missing.
    assert "2" in message
    assert "unit.00000002.done" in message
    assert "exit 0" in message


# ---------------------------------------------------------------------------
# invariant 3: the permit comes back on every abnormal exit path
# ---------------------------------------------------------------------------


def test_a_core_lease_permit_is_released_when_a_component_worker_dies(
    tmp_path, monkeypatch
):
    """A worker that dies by signal must not take its core with it."""

    lraa, lease, observed_free = _build_round(
        tmp_path,
        monkeypatch,
        roles=["publish", "sigkill"],
        budget=4,
        ceiling=2,
        sleeps=[0.3, 0.0],
        measure_index=0,
    )

    with pytest.raises(RuntimeError, match="1 component failures encountered"):
        lraa.reconstruct_isoforms()

    # Non-vacuity: two permits really were out while the workers ran, so the
    # count below is a release rather than a grant that never happened.
    assert observed_free.value == 2
    assert _free_permits(lease, 4) == 4


def test_a_core_lease_permit_is_released_when_the_shard_accounting_raises(
    tmp_path, monkeypatch
):
    """The ordering this test exists to pin: drain INSIDE the release handler.

    ``wait_for_remaining_processes`` is where the shard accounting is enforced,
    and it RAISES -- so the drain is an exception path, not just a collection
    step.  In ``reconstruct_isoforms`` it must therefore sit inside the ``try``
    whose handler calls ``release_component_workers()``.  Moving that call out
    to after the handler reads as a tidy separation of collection from
    teardown and silently reintroduces the leak: every unit that exits 0
    without publishing takes its whole grant with it, and the node's usable
    width shrinks for the rest of the run with nothing in the log to say so.

    Measured on the mutant, with the drain moved back outside the handler:
    this test reports 2 free permits of 4 instead of 4 of 4.  If you are moving
    that call, that number is what you are trading away.
    """

    lraa, lease, observed_free = _build_round(
        tmp_path,
        monkeypatch,
        roles=["publish", "exit_clean_without_publishing"],
        budget=4,
        ceiling=2,
        sleeps=[0.3, 0.0],
        measure_index=0,
    )

    with pytest.raises(WorkUnitAccountingError):
        lraa.reconstruct_isoforms()

    assert observed_free.value == 2
    assert _free_permits(lease, 4) == 4
