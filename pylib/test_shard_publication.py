#!/usr/bin/env python3

"""Per-unit shard publication: the two ways a worker used to lose its result.

Workers returned results through a multiprocessing.Queue, which failed in two
ways that no amount of care at the call site could fix:

  R2, payload drop.  `Queue.put()` pickles on a background feeder thread.  When
  that pickling raised -- an unpicklable member, or a RecursionError on a deep
  splice-graph object graph -- CPython dropped the item and the worker still
  exited 0.  Here the pickling happens in the worker's own stack, so the failure
  is the worker's failure: it raises, it exits nonzero, and it is reported as a
  failed unit rather than a clean exit that delivered nothing.

  R3, truncated-frame hang.  A worker killed mid-write left a partial frame in
  the pipe and the parent blocked in `_recv_bytes` forever -- no EOFError while
  other writers held the pipe open, and no timeout reaching the read.  Here a
  worker killed mid-write leaves a `.tmp` that is never renamed.  The parent
  looks for a `.done`, does not find one, names the unit, and returns.

Every test in this module has a deadline.  A parent that blocks forever is the
defect under repair, so a suite that could hang would hide precisely the
regression it exists to catch.
"""

import logging
import multiprocessing
import os
import pickle
import signal
import sys
import threading
import time
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

from MultiProcessManager import (
    MultiProcessManager,
    ShardPublicationError,
    ShardStore,
    WorkUnitAccountingError,
)

# Generous next to the work each test does (all of it sub-second apart from the
# staggered-sleep test) and far below any plausible hang.
TEST_DEADLINE_SEC = 60.0

# What "the parent did not wedge" is asserted against.  The scenario finishes in
# well under a second; a blocked read would never finish at all.
PARENT_DEADLINE_SEC = 30.0


@pytest.fixture(autouse=True)
def fail_rather_than_hang():
    """Bound every test, and leave no worker behind if one fails."""

    def _expired(signum, frame):
        raise AssertionError(
            "test exceeded its {}s deadline -- a blocked parent is the failure "
            "mode under repair".format(TEST_DEADLINE_SEC)
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


def _payload_for(unit_id):
    return "payload-{}".format(unit_id)


def _worker_publish(shard_store, unit_id):
    shard_store.publish(
        unit_id, _payload_for(unit_id), token="chr1+:{}-{}".format(unit_id, unit_id + 9)
    )


def _worker_exit_clean_without_publishing(shard_store, unit_id):
    return


def _worker_publish_foreign_id(shard_store, unit_id, foreign_unit_id):
    shard_store.publish(foreign_unit_id, _payload_for(foreign_unit_id), token="foreign")


def _bind_stderr(stderr_path):
    """Reproduce LRAA's per-worker stderr redirection.

    The worker's stderr is a file inside the per-contig temp directory, which is
    retained when a run fails.  Both levels are bound because both carry part of
    the report: the fd, which multiprocessing's traceback goes to, and
    sys.stderr, which logging.StreamHandler() resolves against.
    """
    stream = open(stderr_path, "w")
    os.dup2(stream.fileno(), 2)
    sys.stderr = stream

    root = logging.getLogger()
    for handler in list(root.handlers):
        root.removeHandler(handler)
    root.addHandler(logging.StreamHandler(stream))
    root.setLevel(logging.INFO)
    return stream


def _worker_publish_unpicklable(shard_store, unit_id, stderr_path):
    stream = _bind_stderr(stderr_path)
    try:
        # A lock cannot be pickled.  On the queue this was dropped by the feeder
        # thread and the worker exited 0 regardless.
        shard_store.publish(
            unit_id, {"transcripts": threading.Lock()}, token="chr1+:1-100"
        )
    finally:
        logging.shutdown()
        stream.flush()


def _worker_publish_recursive(shard_store, unit_id, depth, stderr_path):
    stream = _bind_stderr(stderr_path)
    try:
        # A deep object graph: pickling recurses per link and dies on the
        # recursion limit, which is how a large splice-graph component failed.
        deep = None
        for _ in range(depth):
            deep = [deep]
        shard_store.publish(unit_id, deep, token="chr1+:1-100")
    finally:
        logging.shutdown()
        stream.flush()


def _worker_killed_mid_write(shard_store, unit_id):
    """Die by SIGKILL with a partial shard on disk.

    The bytes are flushed and fsynced first, so the scratch file genuinely
    survives the kill; SIGKILL cannot be caught, so nothing tidies up.  This is
    the shard-file image of the truncated pipe frame that used to wedge the
    parent forever.
    """
    tmp_path = shard_store.tmp_path(unit_id)
    with open(tmp_path, "wb") as fh:
        fh.write(b"\x80\x05\x95\xff\xff\xff\xff")  # opening bytes of a pickle frame
        fh.flush()
        os.fsync(fh.fileno())
        os.kill(os.getpid(), signal.SIGKILL)


def _worker_sleep_then_publish(shard_store, unit_id, sleep_sec, completion_counter):
    time.sleep(sleep_sec)
    with completion_counter.get_lock():
        completion_counter.value += 1
        completion_rank = completion_counter.value
    shard_store.publish(unit_id, completion_rank, token="chr1+:{}".format(unit_id))


def _spawn(shard_store, target, unit_id, name, extra_args=()):
    return multiprocessing.Process(
        target=target, name=name, args=(shard_store, unit_id) + tuple(extra_args)
    )


# ---------------------------------------------------------------------------
# publication mechanics
# ---------------------------------------------------------------------------


def test_publication_is_a_rename_and_leaves_no_scratch(store):
    done_path = store.publish(4, ["transcript-a"], token="chr1+:1-100")

    assert done_path == store.done_path(4)
    assert os.path.exists(done_path)
    assert not os.path.exists(store.tmp_path(4))
    assert store.load(4) == ["transcript-a"]

    survey = store.survey([4])
    assert survey == {"orphan_tmp": [], "unexpected_done": []}


class _ObservesTheDirectoryWhileBeingPickled:
    """Records the shard directory's contents from inside the serialisation.

    Pickling calls __reduce__, so this runs at the one moment that matters: after
    bytes have started to land and before anything has been published.
    """

    def __init__(self, directory):
        self.directory = directory
        self.seen_mid_write = []

    def __reduce__(self):
        self.seen_mid_write.extend(sorted(os.listdir(self.directory)))
        return (str, ("transcripts",))


def test_nothing_is_visible_as_a_result_until_the_rename(store):
    """Publication has to be write-then-rename, not a write into the final name.

    A parent that can see a partially written result is the whole defect; the
    rename is what makes the result appear all at once or not at all.
    """

    payload = _ObservesTheDirectoryWhileBeingPickled(store.directory)
    store.publish(6, payload, token="chr1+:1-100")

    mid_write = payload.seen_mid_write
    assert mid_write, "__reduce__ never ran, so nothing was observed"
    assert any(entry.endswith(".tmp") for entry in mid_write), mid_write
    assert not any(entry.endswith(".done") for entry in mid_write), mid_write

    assert store.load(6) == "transcripts"


def test_a_shard_records_the_unit_it_was_published_for(store):
    store.publish(4, ["transcript-a"], token="chr1+:1-100")

    with open(store.done_path(4), "rb") as fh:
        frame = pickle.load(fh)

    # The name of a file is not evidence about its contents; the contents say
    # which unit and which region they describe, and load() checks.
    assert frame["unit_id"] == 4
    assert frame["token"] == "chr1+:1-100"


def test_a_shard_holding_another_units_result_is_refused(store):
    store.publish(3, "payload-3", token="chr1+:1-100")
    os.rename(store.done_path(3), store.done_path(4))

    with pytest.raises(WorkUnitAccountingError) as excinfo:
        store.load(4)

    message = str(excinfo.value)
    assert "published for unit 3" in message, message
    assert "read as unit 4" in message, message


def test_a_truncated_shard_fails_loudly_instead_of_blocking(store):
    store.publish(5, ["transcript-a", "transcript-b"], token="chr1+:1-100")

    with open(store.done_path(5), "rb") as fh:
        intact = fh.read()
    with open(store.done_path(5), "wb") as fh:
        fh.write(intact[: len(intact) // 2])

    # The pipe equivalent of this was an unbounded read.  A short file just ends.
    started = time.time()
    with pytest.raises(WorkUnitAccountingError) as excinfo:
        store.load(5)
    assert time.time() - started < 5.0

    assert "could not be read" in str(excinfo.value)
    assert store.done_path(5) in str(excinfo.value)


def test_each_invocation_gets_its_own_directory(tmp_path):
    """The ME round and the SE round call reconstruct_isoforms with the same
    contig and strand.  A directory keyed on those alone would have the second
    round reading the first round's shards, which is a misattributed result
    rather than a crash."""

    first = ShardStore.create("chr1.fwd.ME.round1", base_dir=str(tmp_path))
    second = ShardStore.create("chr1.fwd.ME.round1", base_dir=str(tmp_path))
    try:
        assert first.directory != second.directory

        first.publish(1, "from-the-first-round", token="chr1+:1-100")
        assert second.survey([1]) == {"orphan_tmp": [], "unexpected_done": []}
        assert not os.path.exists(second.done_path(1))
    finally:
        first.discard()
        second.discard()


# ---------------------------------------------------------------------------
# R2: a payload that cannot be serialised
# ---------------------------------------------------------------------------


def test_an_unserialisable_payload_raises_in_the_worker_naming_the_unit(store, caplog):
    with caplog.at_level(logging.ERROR, logger="MultiProcessManager"):
        with pytest.raises(ShardPublicationError) as excinfo:
            store.publish(7, {"transcripts": threading.Lock()}, token="chr1+:1-100")

    message = str(excinfo.value)
    assert "work unit 7" in message, message
    assert "chr1+:1-100" in message, message
    assert "cannot pickle" in message, message

    logged = "\n".join(record.getMessage() for record in caplog.records)
    assert "work unit 7" in logged, logged

    # Nothing is left that could be read as a result, and nothing was published.
    assert not os.path.exists(store.done_path(7))
    assert not os.path.exists(store.tmp_path(7))


def test_a_deep_object_graph_raises_rather_than_being_dropped(store):
    deep = None
    for _ in range(sys.getrecursionlimit() * 20):
        deep = [deep]

    with pytest.raises(ShardPublicationError) as excinfo:
        store.publish(8, deep, token="chr1+:1-100")

    assert "work unit 8" in str(excinfo.value)
    assert "RecursionError" in str(excinfo.value)
    assert not os.path.exists(store.done_path(8))


@pytest.mark.parametrize(
    "worker,extra,expected_reason",
    [
        (_worker_publish_unpicklable, (), "cannot pickle"),
        (_worker_publish_recursive, (sys.getrecursionlimit() * 20,), "RecursionError"),
    ],
    ids=["unpicklable-member", "deep-object-graph"],
)
def test_a_unit_that_cannot_serialise_fails_loudly_and_alone(
    store, tmp_path, worker, extra, expected_reason
):
    """R2 closed: the failure is the unit's own, and the other units are intact.

    Under the queue this unit exited 0 having delivered nothing, and the loss
    surfaced only as an accounting failure that could not say why.  It now shows
    up as a failed unit with a reason, in the worker's own stderr log -- the file
    LRAA retains for a contig whose run failed.
    """

    stderr_path = tmp_path / "worker.err.log"
    mgr = MultiProcessManager(3, store)

    for unit_id in (1, 3):
        mgr.launch_process(
            _spawn(store, _worker_publish, unit_id, "chr1+:1-100"), unit_id=unit_id
        )
    mgr.launch_process(
        _spawn(store, worker, 2, "chr1+:1-100", tuple(extra) + (str(stderr_path),)),
        unit_id=2,
    )

    # One failed unit, and no accounting error: a unit that died owes no shard.
    assert mgr.wait_for_remaining_processes() == 1

    (failure,) = mgr.get_failures()
    assert failure["unit_id"] == 2
    assert failure["exitcode"] != 0

    # No other unit's result was lost.
    assert mgr.retrieve_unit_results() == [(1, "payload-1"), (3, "payload-3")]

    captured = stderr_path.read_text()
    assert "work unit 2" in captured, captured
    assert expected_reason in captured, captured
    assert "ShardPublicationError" in captured, captured


# ---------------------------------------------------------------------------
# R3: a worker killed mid-write
# ---------------------------------------------------------------------------


def test_a_worker_sigkilled_mid_write_does_not_wedge_the_parent(store):
    """R3 closed, and this is the central test of the change.

    A worker with a partial write on disk is SIGKILLed.  Under the queue the
    equivalent -- a frame half-written into the pipe -- blocked the parent in
    `_recv_bytes` for the life of the run: no EOFError arrives while other
    writers hold the pipe open, and no timeout reaches that read, so the job
    simply looked slow forever.  Here the parent reads files it names itself, so
    there is nothing to block on: the killed unit has no `.done`, and the parent
    says which unit and with what exit code.

    Asserted against a wall clock.  Were the parent to block, this fails on the
    elapsed-time assertion, or on the module deadline if it never returns at all.
    """

    mgr = MultiProcessManager(3, store)

    for unit_id in (1, 3):
        mgr.launch_process(
            _spawn(store, _worker_publish, unit_id, "chr1+:1-100"), unit_id=unit_id
        )
    mgr.launch_process(
        _spawn(store, _worker_killed_mid_write, 2, "chr1+:1-100"), unit_id=2
    )

    started = time.time()
    num_errors = mgr.wait_for_remaining_processes()
    elapsed = time.time() - started

    assert elapsed < PARENT_DEADLINE_SEC, (
        "the parent took {:.1f}s to collect three units; a wedged read is exactly "
        "what this change removes".format(elapsed)
    )

    # The killed unit is reported, by id and by signal.
    assert num_errors == 1
    (failure,) = mgr.get_failures()
    assert failure["unit_id"] == 2
    assert failure["exitcode"] == -signal.SIGKILL

    # Its partial write was never mistaken for a result, and it was reported.
    assert not os.path.exists(store.done_path(2))
    assert len(mgr.orphan_shard_names) == 1
    assert mgr.orphan_shard_names[0].startswith("unit.00000002.")
    assert mgr.orphan_shard_names[0].endswith(".tmp")

    # The surviving units are complete and in order.
    assert mgr.retrieve_unit_results() == [(1, "payload-1"), (3, "payload-3")]


def test_a_partial_shard_is_invisible_to_collection(store):
    """The same property at the store level: only a rename publishes."""

    with open(store.tmp_path(2), "wb") as fh:
        fh.write(b"\x80\x05\x95")

    assert not os.path.exists(store.done_path(2))

    survey = store.survey([2], remove_orphans=False)
    assert len(survey["orphan_tmp"]) == 1
    assert survey["unexpected_done"] == []


# ---------------------------------------------------------------------------
# exactly-once
# ---------------------------------------------------------------------------


def test_n_units_in_n_results_out(store):
    unit_ids = list(range(8))
    mgr = MultiProcessManager(4, store)

    for unit_id in unit_ids:
        mgr.launch_process(
            _spawn(store, _worker_publish, unit_id, "chr1+:1-100"), unit_id=unit_id
        )

    assert mgr.wait_for_remaining_processes() == 0
    assert mgr.retrieve_unit_results() == [
        (unit_id, _payload_for(unit_id)) for unit_id in unit_ids
    ]
    assert mgr.returned_unit_ids == set(unit_ids)


def test_a_unit_that_exits_0_without_publishing_still_raises(store):
    """The exactly-once guarantee, unchanged in strength and cheaper to hold: a
    missing result is now a missing file the parent can name a path for."""

    mgr = MultiProcessManager(2, store)

    mgr.launch_process(_spawn(store, _worker_publish, 1, "chr1+:1-100"), unit_id=1)
    mgr.launch_process(
        _spawn(store, _worker_exit_clean_without_publishing, 2, "chr1+:1-100"),
        unit_id=2,
    )

    with pytest.raises(WorkUnitAccountingError) as excinfo:
        mgr.wait_for_remaining_processes()

    message = str(excinfo.value)
    assert "published no result shard" in message, message
    assert "exit 0" in message, message
    assert store.done_path(2) in message, message

    # The unit really did exit cleanly; that is what made the loss invisible.
    assert mgr.num_errors == 0
    assert mgr.num_successes == 2


def test_a_shard_published_for_an_unsubmitted_unit_is_rejected(store):
    mgr = MultiProcessManager(1, store)

    mgr.launch_process(
        _spawn(store, _worker_publish_foreign_id, 1, "chr1+:1-100", extra_args=(99,)),
        unit_id=1,
    )

    with pytest.raises(WorkUnitAccountingError) as excinfo:
        mgr.wait_for_remaining_processes()

    message = str(excinfo.value)
    assert "never-submitted unit id" in message, message
    assert "unit.00000099.done" in message, message


def test_payloads_are_not_read_until_the_caller_asks(store):
    """The parent's copy of a result exists only once it has been requested,
    which is what keeps a contig's transcripts out of the parent's memory."""

    mgr = MultiProcessManager(2, store)
    for unit_id in (1, 2):
        mgr.launch_process(
            _spawn(store, _worker_publish, unit_id, "chr1+:1-100"), unit_id=unit_id
        )

    assert mgr.wait_for_remaining_processes() == 0

    # Collected means "the file is there", not "the payload is in hand".
    assert set(mgr.published_shard_paths) == {1, 2}
    assert all(os.path.exists(path) for path in mgr.published_shard_paths.values())

    assert mgr.retrieve_unit_results(clear=True) == [
        (1, "payload-1"),
        (2, "payload-2"),
    ]
    # Read once, and the shard is gone rather than accumulating on disk.
    assert not os.path.exists(store.done_path(1))
    assert mgr.retrieve_unit_results() == []
    mgr.verify_units_accounted()


# ---------------------------------------------------------------------------
# order independence
# ---------------------------------------------------------------------------


def test_results_concatenate_in_unit_order_whatever_the_completion_order(store):
    """Staggered so the units finish in reverse.  The payload is the completion
    rank, so this shows submission order was recovered rather than echoed."""

    unit_ids = [0, 1, 2, 3]
    completion_counter = multiprocessing.Value("i", 0)
    mgr = MultiProcessManager(len(unit_ids), store)

    for unit_id in unit_ids:
        sleep_sec = 0.2 * (len(unit_ids) - unit_id)
        mgr.launch_process(
            _spawn(
                store,
                _worker_sleep_then_publish,
                unit_id,
                "chr1+:1-100",
                (sleep_sec, completion_counter),
            ),
            unit_id=unit_id,
        )

    assert mgr.wait_for_remaining_processes() == 0

    results = mgr.retrieve_unit_results()
    assert [unit_id for unit_id, _ in results] == unit_ids
    assert [rank for _, rank in results] == [4, 3, 2, 1]


# ---------------------------------------------------------------------------
# orphan handling
# ---------------------------------------------------------------------------


def test_a_stray_scratch_file_is_reported_and_never_read_as_a_result(store, caplog):
    """A `.tmp` left by a previous run is a failed unit, not a result.

    It cannot be collected -- only `.done` files are ever read -- and it is
    reported by name and count, because a nonzero count means a writer was
    killed.  It is not fatal on its own: the unit that died has already failed on
    its exit code, and a stray file must not fail an otherwise complete run.
    """

    stray = os.path.join(store.directory, "unit.00000042.99999.tmp")
    with open(stray, "wb") as fh:
        fh.write(b"half a pickle")

    mgr = MultiProcessManager(2, store)
    for unit_id in (1, 2):
        mgr.launch_process(
            _spawn(store, _worker_publish, unit_id, "chr1+:1-100"), unit_id=unit_id
        )

    with caplog.at_level(logging.ERROR, logger="MultiProcessManager"):
        assert mgr.wait_for_remaining_processes() == 0

    assert mgr.orphan_shard_names == ["unit.00000042.99999.tmp"]
    assert not os.path.exists(stray), "the orphan is cleaned up once reported"

    logged = "\n".join(record.getMessage() for record in caplog.records)
    assert "1 orphaned shard scratch file" in logged, logged
    assert "unit.00000042.99999.tmp" in logged, logged

    # The stray contributed nothing and cost nothing.
    assert mgr.retrieve_unit_results() == [(1, "payload-1"), (2, "payload-2")]
