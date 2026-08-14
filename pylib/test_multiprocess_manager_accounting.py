#!/usr/bin/env python3

"""Work-unit accounting and reaping order in MultiProcessManager.

The manager's job is not just to keep a pool full; it is to guarantee that every
unit submitted comes back exactly once.  Results travel as per-unit shard files
published by atomic rename (see test_shard_publication.py for the transport
itself); what is pinned here is the bookkeeping around it:

  * reaping order -- a result must be collected after the process that produced
    it has been joined, never on a schedule of the parent's own;
  * tagged, audited results -- a result that never arrives is otherwise
    indistinguishable from a unit that had nothing to say.

Every test has a deadline.  The failure being guarded against is a parent that
never finishes, so a suite that could hang would hide it.
"""

import multiprocessing
import os
import signal
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

from MultiProcessManager import MultiProcessManager, ShardStore, WorkUnitAccountingError

TEST_DEADLINE_SEC = 60.0


@pytest.fixture(autouse=True)
def fail_rather_than_hang():
    def _expired(signum, frame):
        raise AssertionError(
            "test exceeded its {}s deadline".format(TEST_DEADLINE_SEC)
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


@pytest.fixture
def store(tmp_path):
    shard_store = ShardStore.create("chr1.fwd.ME.round1", base_dir=str(tmp_path))
    yield shard_store
    shard_store.discard()


# ---------------------------------------------------------------------------
# workers -- module scope so they survive either start method
# ---------------------------------------------------------------------------


def _worker_publish(shard_store, unit_id):
    shard_store.publish(unit_id, "payload-{}".format(unit_id), token="unit")


def _worker_exit_clean_without_publishing(shard_store, unit_id):
    return


def _worker_die(shard_store, unit_id):
    raise RuntimeError("unit {} is meant to fail".format(unit_id))


def _worker_noop():
    return


def _worker_gated_publish(shard_store, unit_id, gate):
    gate.wait(60)
    shard_store.publish(unit_id, "payload-{}".format(unit_id), token="unit")


def _worker_chained_publish(
    shard_store, unit_id, my_gate, next_gate, completion_counter
):
    """Publish only when released, then release the next unit in the chain.

    Chaining the gates makes completion order a property of the test rather than
    of the scheduler, so the ordering assertion below has nothing to be flaky
    about.
    """
    my_gate.wait(60)
    with completion_counter.get_lock():
        completion_counter.value += 1
        completion_rank = completion_counter.value
    shard_store.publish(unit_id, completion_rank, token="unit")
    if next_gate is not None:
        next_gate.set()


def _spawn(shard_store, target, unit_id, name, extra_args=()):
    return multiprocessing.Process(
        target=target, name=name, args=(shard_store, unit_id) + tuple(extra_args)
    )


# ---------------------------------------------------------------------------
# reaping order
# ---------------------------------------------------------------------------


class _ProcessPublishingAsItIsJoined:
    """Process stub that publishes its shard from inside join().

    A worker renames its shard and then exits, so the last possible instant at
    which the result can appear is the instant the process goes away.  Collecting
    after join() therefore cannot miss it, and collecting anywhere else can: the
    result queue this replaced was drained on the parent's own schedule, and a
    child that flushed and exited in the gap was joined as a success with its
    payload unread.
    """

    def __init__(self, name, shard_store, unit_id, payload):
        self.name = name
        self.pid = -1
        self.exitcode = 0
        self._shard_store = shard_store
        self._unit_id = unit_id
        self._payload = payload
        self.published = False

    def start(self):
        pass

    def is_alive(self):
        return False

    def join(self, timeout=None):
        if not self.published:
            self._shard_store.publish(self._unit_id, self._payload, token="unit")
            self.published = True


def test_a_shard_published_at_the_moment_of_death_is_collected(store):
    """Collection must happen after the join, not on the parent's own schedule."""

    mgr = MultiProcessManager(1, store)

    # Drive the screening pass directly: the interleaving under test is between
    # two statements inside it, and there is no way to schedule a real worker
    # into that gap on demand.
    process = _ProcessPublishingAsItIsJoined(
        "chr1+:100-200:1", store, 1, "payload-1"
    )
    mgr.submitted_unit_ids = {1}
    mgr.process_to_unit_id[process] = 1
    mgr.process_list.add(process)
    mgr.num_running = 1

    mgr._screen_running_processes()

    assert mgr.retrieve_unit_results() == [(1, "payload-1")]
    assert mgr.num_running == 0
    assert mgr.num_successes == 1
    mgr.verify_units_accounted()


# ---------------------------------------------------------------------------
# completeness
# ---------------------------------------------------------------------------


def test_clean_exit_without_a_result_is_not_a_silent_success(store):
    mgr = MultiProcessManager(2, store)

    mgr.launch_process(_spawn(store, _worker_publish, 1, "chr1+"), unit_id=1)
    mgr.launch_process(
        _spawn(store, _worker_exit_clean_without_publishing, 2, "chr1+"), unit_id=2
    )

    with pytest.raises(WorkUnitAccountingError) as excinfo:
        mgr.wait_for_remaining_processes()

    message = str(excinfo.value)
    assert "published no result shard" in message
    assert "2 (exit 0" in message, message
    # The unit really did exit cleanly; that is precisely what made the loss
    # invisible before.
    assert mgr.num_errors == 0
    assert mgr.num_successes == 2


def test_two_unit_ids_that_share_a_shard_name_are_rejected_at_launch(store):
    """Distinct units must be distinct on disk.

    A unit publishes to a file named after its id, so two ids that reduce to one
    name -- 1 and "1" -- would have the second silently overwrite the first.  The
    queue keyed a dict and could not care; a filesystem can, so this is checked
    at launch, before either worker starts.
    """

    mgr = MultiProcessManager(2, store)

    mgr.launch_process(_spawn(store, _worker_publish, 1, "chr1+"), unit_id=1)
    with pytest.raises(WorkUnitAccountingError) as excinfo:
        mgr.launch_process(_spawn(store, _worker_publish, "00000001", "chr1+"), unit_id="00000001")

    message = str(excinfo.value)
    assert "publish to shard" in message, message
    assert "unit.00000001" in message, message

    assert mgr.wait_for_remaining_processes() == 0


def test_resubmitting_a_unit_id_is_rejected_at_launch(store):
    mgr = MultiProcessManager(2, store)

    mgr.launch_process(_spawn(store, _worker_publish, 3, "chr1+"), unit_id=3)
    with pytest.raises(WorkUnitAccountingError) as excinfo:
        mgr.launch_process(_spawn(store, _worker_publish, 3, "chr1+"), unit_id=3)
    assert "submitted more than once" in str(excinfo.value)

    assert mgr.wait_for_remaining_processes() == 0


def test_a_shard_store_requires_a_unit_id(store):
    mgr = MultiProcessManager(1, store)

    with pytest.raises(ValueError):
        mgr.launch_process(_spawn(store, _worker_publish, 1, "chr1+"))


def test_failed_unit_owes_no_result_but_is_named_in_the_failure_record(store):
    mgr = MultiProcessManager(2, store)

    mgr.launch_process(_spawn(store, _worker_publish, 1, "chr1+"), unit_id=1)
    mgr.launch_process(_spawn(store, _worker_die, 2, "chr1+"), unit_id=2)

    assert mgr.wait_for_remaining_processes() == 1

    (failure,) = mgr.get_failures()
    assert failure["unit_id"] == 2
    assert failure["name"] == "chr1+:2"
    assert failure["exitcode"] != 0
    assert failure["duration_sec"] is not None
    assert mgr.retrieve_unit_results() == [(1, "payload-1")]


# ---------------------------------------------------------------------------
# identity and ordering
# ---------------------------------------------------------------------------


def test_process_names_are_unique_across_concurrently_submitted_units(store):
    gate = multiprocessing.Event()
    mgr = MultiProcessManager(4, store)

    # Every component of a contig/strand is named after the same mpg token in
    # LRAA, so the shared base name here is the realistic case.
    processes = []
    for unit_id in range(4):
        process = _spawn(
            store, _worker_gated_publish, unit_id, "chr1+:100-200", (gate,)
        )
        mgr.launch_process(process, unit_id=unit_id)
        processes.append(process)

    # Nothing has been released, so all four units are genuinely concurrent.
    assert len(mgr.process_list) == 4

    names = [process.name for process in processes]
    assert names == [
        "chr1+:100-200:0",
        "chr1+:100-200:1",
        "chr1+:100-200:2",
        "chr1+:100-200:3",
    ]
    assert len(set(names)) == len(names)

    gate.set()
    assert mgr.wait_for_remaining_processes() == 0


def test_results_are_ordered_by_unit_id_not_by_completion(store):
    unit_ids = [0, 1, 2, 3]
    gates = [multiprocessing.Event() for _ in unit_ids]
    completion_counter = multiprocessing.Value("i", 0)
    mgr = MultiProcessManager(len(unit_ids), store)

    # Unit i is released by unit i+1, so results are published 3, 2, 1, 0.
    for unit_id in unit_ids:
        next_gate = gates[unit_id - 1] if unit_id > 0 else None
        mgr.launch_process(
            _spawn(
                store,
                _worker_chained_publish,
                unit_id,
                "chr1+:100-200",
                (gates[unit_id], next_gate, completion_counter),
            ),
            unit_id=unit_id,
        )

    gates[unit_ids[-1]].set()
    assert mgr.wait_for_remaining_processes() == 0

    results = mgr.retrieve_unit_results()
    assert [unit_id for unit_id, _ in results] == unit_ids
    # Payload is the completion rank: descending, so submission order was
    # recovered rather than merely echoed back.
    assert [completion_rank for _, completion_rank in results] == [4, 3, 2, 1]


def test_retrieving_results_clears_only_the_buffer_not_the_audit(store):
    mgr = MultiProcessManager(2, store)

    for unit_id in (1, 2):
        mgr.launch_process(
            _spawn(store, _worker_publish, unit_id, "chr1+"), unit_id=unit_id
        )

    assert mgr.wait_for_remaining_processes() == 0

    assert mgr.retrieve_unit_results(clear=True) == [
        (1, "payload-1"),
        (2, "payload-2"),
    ]
    assert mgr.retrieve_unit_results() == []
    # Collecting incrementally must not make the completeness check forget what
    # already came back.
    mgr.verify_units_accounted()


def test_per_unit_start_times_are_not_overwritten_by_later_launches(store):
    gate = multiprocessing.Event()
    mgr = MultiProcessManager(2, store)

    slow = _spawn(store, _worker_gated_publish, 1, "chr1+:100-200", (gate,))
    mgr.launch_process(slow, unit_id=1)
    first_start = mgr.process_start_time[slow]

    fast = _spawn(store, _worker_gated_publish, 2, "chr1+:100-200", (gate,))
    mgr.launch_process(fast, unit_id=2)

    # Both workers share a base name; keying start times by name collapsed them
    # onto one entry and reported the newest launch as every unit's start.
    assert mgr.process_start_time[slow] == first_start
    assert mgr.process_start_time[fast] >= first_start

    gate.set()
    assert mgr.wait_for_remaining_processes() == 0


# ---------------------------------------------------------------------------
# pools with no shard store (the contig/strand scheduler in the LRAA driver,
# which publishes its own per-unit files and needs no results back)
# ---------------------------------------------------------------------------


def test_storeless_pool_needs_no_results():
    mgr = MultiProcessManager(2)

    for unit_id in (0, 1):
        mgr.launch_process(
            multiprocessing.Process(target=_worker_noop, name="chr17+"),
            unit_id=unit_id,
        )

    assert mgr.wait_for_remaining_processes() == 0
    assert mgr.num_successes == 2


def test_a_torn_down_run_is_exempt_from_the_completeness_check(store):
    """terminate_all_processes means the run is being abandoned, so the units it
    killed owe nothing.  Without the exemption every teardown would raise on the
    way out and bury the error that caused it."""

    gate = multiprocessing.Event()
    mgr = MultiProcessManager(2, store)

    for unit_id in (1, 2):
        mgr.launch_process(
            _spawn(store, _worker_gated_publish, unit_id, "chr1+:100-200", (gate,)),
            unit_id=unit_id,
        )

    mgr.terminate_all_processes()
    mgr.wait_for_remaining_processes()

    assert mgr.returned_unit_ids == set()
    assert not os.path.exists(store.done_path(1))
