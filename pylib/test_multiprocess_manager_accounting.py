#!/usr/bin/env python3

"""Work-unit accounting and reaping order in MultiProcessManager.

The manager's job is not just to keep a pool full; it is to guarantee that every
unit submitted comes back exactly once.  Two things break that guarantee and
both are silent, which is why they are pinned down here:

  * reaping order -- draining the result queue before checking which children
    have died lets a child flush its payload and exit in the gap, so the parent
    joins it as a success with its bytes still unread;
  * untagged, unaudited results -- a payload that never arrives is
    indistinguishable from a unit that had nothing to say.
"""

import multiprocessing
import sys
from pathlib import Path
from queue import Empty as QueueEmpty

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

from MultiProcessManager import MultiProcessManager, WorkUnitAccountingError


# ---------------------------------------------------------------------------
# workers -- module scope so they survive either start method
# ---------------------------------------------------------------------------


def _worker_put(q, unit_id):
    q.put((unit_id, "payload-{}".format(unit_id)))


def _worker_exit_clean_without_payload(q, unit_id):
    return


def _worker_put_twice(q, unit_id):
    q.put((unit_id, "payload-{}".format(unit_id)))
    q.put((unit_id, "payload-{}-again".format(unit_id)))


def _worker_put_foreign_id(q, unit_id, foreign_unit_id):
    q.put((foreign_unit_id, "payload-{}".format(foreign_unit_id)))


def _worker_die(q, unit_id):
    raise RuntimeError("unit {} is meant to fail".format(unit_id))


def _worker_noop():
    return


def _worker_gated_put(q, unit_id, gate):
    gate.wait(60)
    q.put((unit_id, "payload-{}".format(unit_id)))


def _worker_chained_put(q, unit_id, my_gate, next_gate, arrival_counter):
    """Put only when released, then release the next unit in the chain.

    Chaining the gates makes arrival order a property of the test rather than of
    the scheduler, so the ordering assertion below has nothing to be flaky about.
    """
    my_gate.wait(60)
    with arrival_counter.get_lock():
        arrival_counter.value += 1
        arrival_rank = arrival_counter.value
    q.put((unit_id, arrival_rank))
    if next_gate is not None:
        next_gate.set()


def _spawn(q, target, unit_id, name, extra_args=()):
    return multiprocessing.Process(
        target=target, name=name, args=(q, unit_id) + tuple(extra_args)
    )


# ---------------------------------------------------------------------------
# reaping order
# ---------------------------------------------------------------------------


class _FlushOnDeathQueue:
    """Queue stub whose payload becomes readable exactly when the process dies.

    A worker hands its payload to the queue's feeder thread and then exits, so
    the bytes can reach the pipe at any point up to the instant the parent
    observes the process gone.  This stub pins that instant down: nothing is
    readable until `flush`, which the process stub below calls from is_alive().

    A manager that drains before taking the liveness snapshot therefore sees an
    empty queue, then sees a dead process, and joins it as a success with the
    payload still unread.  A manager that snapshots first cannot lose it.
    """

    def __init__(self, pending):
        self._pending = list(pending)
        self._readable = []

    def flush(self):
        self._readable.extend(self._pending)
        self._pending = []

    def empty(self):
        return not self._readable

    def get(self, block=True, timeout=None):
        if not self._readable:
            raise QueueEmpty
        return self._readable.pop(0)

    def get_nowait(self):
        return self.get(block=False)


class _ProcessDeadOnFirstLook:
    """Process stub that completes its queue flush as it is checked for liveness."""

    def __init__(self, name, flush_on_death_queue):
        self.name = name
        self.pid = -1
        self.exitcode = 0
        self._queue = flush_on_death_queue

    def start(self):
        pass

    def is_alive(self):
        self._queue.flush()
        return False

    def join(self, timeout=None):
        pass


def test_payload_flushed_at_the_moment_of_death_is_not_reaped_past():
    """The liveness snapshot must be taken before the drain, not after."""

    q = _FlushOnDeathQueue([(1, "payload-1")])
    mgr = MultiProcessManager(1, q)

    # Drive the screening pass directly: the interleaving under test is between
    # two statements inside it, and there is no way to schedule a real worker
    # into that gap on demand.
    process = _ProcessDeadOnFirstLook("chr1+:100-200:1", q)
    mgr.submitted_unit_ids = {1}
    mgr.process_list.add(process)
    mgr.num_running = 1

    mgr._screen_running_processes()

    assert mgr.retrieve_queue_contents() == [(1, "payload-1")]
    assert mgr.num_running == 0
    assert mgr.num_successes == 1


class _EmptyIsATripwireQueue:
    """Delegating queue proxy that fails the test if empty() is consulted."""

    def __init__(self, wrapped):
        self._wrapped = wrapped

    def empty(self):
        raise AssertionError(
            "Queue.empty() lags put() by the feeder thread and must never gate draining"
        )

    def put(self, *args, **kwargs):
        return self._wrapped.put(*args, **kwargs)

    def get(self, *args, **kwargs):
        return self._wrapped.get(*args, **kwargs)

    def get_nowait(self):
        return self._wrapped.get_nowait()


def test_draining_never_consults_queue_empty():
    real_q = multiprocessing.Queue()
    mgr = MultiProcessManager(2, _EmptyIsATripwireQueue(real_q))

    for unit_id in (1, 2):
        mgr.launch_process(
            _spawn(real_q, _worker_put, unit_id, "chr1+"), unit_id=unit_id
        )

    assert mgr.wait_for_remaining_processes() == 0
    assert mgr.retrieve_queue_contents() == [(1, "payload-1"), (2, "payload-2")]


# ---------------------------------------------------------------------------
# completeness
# ---------------------------------------------------------------------------


def test_clean_exit_without_a_payload_is_not_a_silent_success():
    q = multiprocessing.Queue()
    mgr = MultiProcessManager(2, q)

    mgr.launch_process(_spawn(q, _worker_put, 1, "chr1+"), unit_id=1)
    mgr.launch_process(
        _spawn(q, _worker_exit_clean_without_payload, 2, "chr1+"), unit_id=2
    )

    with pytest.raises(WorkUnitAccountingError) as excinfo:
        mgr.wait_for_remaining_processes()

    message = str(excinfo.value)
    assert "exited 0 but returned no payload" in message
    assert "[2]" in message, message
    # The unit really did exit cleanly; that is precisely what made the loss
    # invisible before.
    assert mgr.num_errors == 0
    assert mgr.num_successes == 2


def test_a_unit_returning_two_payloads_is_rejected():
    q = multiprocessing.Queue()
    mgr = MultiProcessManager(2, q)

    mgr.launch_process(_spawn(q, _worker_put, 1, "chr1+"), unit_id=1)
    mgr.launch_process(_spawn(q, _worker_put_twice, 7, "chr1+"), unit_id=7)

    with pytest.raises(WorkUnitAccountingError) as excinfo:
        mgr.wait_for_remaining_processes()

    message = str(excinfo.value)
    assert "returned more than one payload" in message
    assert "[7]" in message, message


def test_payload_tagged_with_an_unsubmitted_id_is_rejected():
    q = multiprocessing.Queue()
    mgr = MultiProcessManager(1, q)

    mgr.launch_process(
        _spawn(q, _worker_put_foreign_id, 1, "chr1+", extra_args=(99,)), unit_id=1
    )

    with pytest.raises(WorkUnitAccountingError) as excinfo:
        mgr.wait_for_remaining_processes()

    message = str(excinfo.value)
    assert "never-submitted unit id" in message
    assert "[99]" in message, message


def test_resubmitting_a_unit_id_is_rejected_at_launch():
    q = multiprocessing.Queue()
    mgr = MultiProcessManager(2, q)

    mgr.launch_process(_spawn(q, _worker_put, 3, "chr1+"), unit_id=3)
    with pytest.raises(WorkUnitAccountingError) as excinfo:
        mgr.launch_process(_spawn(q, _worker_put, 3, "chr1+"), unit_id=3)
    assert "submitted more than once" in str(excinfo.value)

    assert mgr.wait_for_remaining_processes() == 0


def test_a_result_queue_requires_a_unit_id():
    q = multiprocessing.Queue()
    mgr = MultiProcessManager(1, q)

    with pytest.raises(ValueError):
        mgr.launch_process(_spawn(q, _worker_put, 1, "chr1+"))


def test_failed_unit_owes_no_payload_but_is_named_in_the_failure_record():
    q = multiprocessing.Queue()
    mgr = MultiProcessManager(2, q)

    mgr.launch_process(_spawn(q, _worker_put, 1, "chr1+"), unit_id=1)
    mgr.launch_process(_spawn(q, _worker_die, 2, "chr1+"), unit_id=2)

    assert mgr.wait_for_remaining_processes() == 1

    (failure,) = mgr.get_failures()
    assert failure["unit_id"] == 2
    assert failure["name"] == "chr1+:2"
    assert failure["exitcode"] != 0
    assert failure["duration_sec"] is not None
    assert mgr.retrieve_queue_contents() == [(1, "payload-1")]


# ---------------------------------------------------------------------------
# identity and ordering
# ---------------------------------------------------------------------------


def test_process_names_are_unique_across_concurrently_submitted_units():
    q = multiprocessing.Queue()
    gate = multiprocessing.Event()
    mgr = MultiProcessManager(4, q)

    # Every component of a contig/strand is named after the same mpg token in
    # LRAA, so the shared base name here is the realistic case.
    processes = []
    for unit_id in range(4):
        process = _spawn(q, _worker_gated_put, unit_id, "chr1+:100-200", (gate,))
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


def test_results_are_ordered_by_unit_id_not_by_arrival():
    unit_ids = [0, 1, 2, 3]
    q = multiprocessing.Queue()
    gates = [multiprocessing.Event() for _ in unit_ids]
    arrival_counter = multiprocessing.Value("i", 0)
    mgr = MultiProcessManager(len(unit_ids), q)

    # Unit i is released by unit i+1, so payloads arrive 3, 2, 1, 0.
    for unit_id in unit_ids:
        next_gate = gates[unit_id - 1] if unit_id > 0 else None
        mgr.launch_process(
            _spawn(
                q,
                _worker_chained_put,
                unit_id,
                "chr1+:100-200",
                (gates[unit_id], next_gate, arrival_counter),
            ),
            unit_id=unit_id,
        )

    gates[unit_ids[-1]].set()
    assert mgr.wait_for_remaining_processes() == 0

    results = mgr.retrieve_queue_contents()
    assert [unit_id for unit_id, _ in results] == unit_ids
    # Payload is the arrival rank: descending, so submission order was recovered
    # rather than merely echoed back.
    assert [arrival_rank for _, arrival_rank in results] == [4, 3, 2, 1]


def test_retrieve_queue_contents_clears_only_the_buffer_not_the_audit():
    q = multiprocessing.Queue()
    mgr = MultiProcessManager(2, q)

    for unit_id in (1, 2):
        mgr.launch_process(_spawn(q, _worker_put, unit_id, "chr1+"), unit_id=unit_id)

    assert mgr.wait_for_remaining_processes() == 0

    assert mgr.retrieve_queue_contents(clear=True) == [
        (1, "payload-1"),
        (2, "payload-2"),
    ]
    assert mgr.retrieve_queue_contents() == []
    # Draining incrementally must not make the completeness check forget what
    # already came back.
    mgr.verify_units_accounted()


def test_per_unit_start_times_are_not_overwritten_by_later_launches():
    q = multiprocessing.Queue()
    gate = multiprocessing.Event()
    mgr = MultiProcessManager(2, q)

    slow = _spawn(q, _worker_gated_put, 1, "chr1+:100-200", (gate,))
    mgr.launch_process(slow, unit_id=1)
    first_start = mgr.process_start_time[slow]

    fast = _spawn(q, _worker_gated_put, 2, "chr1+:100-200", (gate,))
    mgr.launch_process(fast, unit_id=2)

    # Both workers share a base name; keying start times by name collapsed them
    # onto one entry and reported the newest launch as every unit's start.
    assert mgr.process_start_time[slow] == first_start
    assert mgr.process_start_time[fast] >= first_start

    gate.set()
    assert mgr.wait_for_remaining_processes() == 0


# ---------------------------------------------------------------------------
# queue-less pools (the contig/strand scheduler in the LRAA driver)
# ---------------------------------------------------------------------------


def test_queueless_pool_needs_no_payloads():
    mgr = MultiProcessManager(2)

    for unit_id in (0, 1):
        mgr.launch_process(
            multiprocessing.Process(target=_worker_noop, name="chr17+"),
            unit_id=unit_id,
        )

    assert mgr.wait_for_remaining_processes() == 0
    assert mgr.num_successes == 2
