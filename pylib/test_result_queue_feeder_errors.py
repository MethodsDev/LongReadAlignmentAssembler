#!/usr/bin/env python3

"""A payload that cannot be pickled is dropped by the queue's feeder thread while
the worker still exits 0. CPython's default response is `traceback.print_exc()`,
which in this pipeline lands on a worker stderr inside a per-contig temp directory
that is removed on success -- so the only record of the loss is discarded.

`ResultQueue` redirects that hook into the logger, tagged with the unit id, so the
cause survives beside the accounting failure the drop provokes. These tests pin the
hook, not the drop: the item is still lost, which is R2 and is out of scope here.
"""

import logging
import multiprocessing
import threading
import time

import pytest

from MultiProcessManager import ResultQueue


def _wait_for(predicate, timeout=10.0):
    deadline = time.time() + timeout
    while time.time() < deadline:
        if predicate():
            return True
        time.sleep(0.02)
    return False


def test_unpicklable_payload_is_reported_with_its_unit_id(caplog):
    """The message must name the unit, or it cannot be matched to the accounting
    failure that follows."""
    q = ResultQueue()
    try:
        with caplog.at_level(logging.ERROR, logger="MultiProcessManager"):
            # a lock is not picklable; the feeder thread raises while serialising
            q.put((7, threading.Lock()))
            assert _wait_for(lambda: any(r.levelno >= logging.ERROR for r in caplog.records)), (
                "feeder error was never reported"
            )
        text = "\n".join(r.getMessage() for r in caplog.records)
        assert "unit 7" in text
        assert "LOST" in text
    finally:
        q.close()


def test_the_drop_is_described_as_a_loss_not_a_warning(caplog):
    """Wording matters here: a worker that loses its payload still exits 0, so a
    message that reads as advisory invites someone to ignore it."""
    q = ResultQueue()
    try:
        with caplog.at_level(logging.ERROR, logger="MultiProcessManager"):
            q.put((3, threading.Lock()))
            assert _wait_for(lambda: caplog.records)
        record = next(r for r in caplog.records if r.levelno >= logging.ERROR)
        assert record.levelno == logging.ERROR
        assert record.exc_info is not None, "the underlying exception must be attached"
    finally:
        q.close()


def test_a_payload_without_the_expected_shape_still_reports(caplog):
    """The hook must not raise while reporting; a bare unpicklable object has no
    unit id to recover."""
    q = ResultQueue()
    try:
        with caplog.at_level(logging.ERROR, logger="MultiProcessManager"):
            q.put(threading.Lock())
            assert _wait_for(lambda: caplog.records), "feeder error was never reported"
        text = "\n".join(r.getMessage() for r in caplog.records)
        assert "unit unknown" in text
    finally:
        q.close()


def test_picklable_payloads_are_unaffected():
    """The override must not disturb the normal path."""
    q = ResultQueue()
    try:
        q.put((1, ["transcript-a", "transcript-b"]))
        assert _wait_for(lambda: not q.empty())
        assert q.get(timeout=5) == (1, ["transcript-a", "transcript-b"])
    finally:
        q.close()


def test_the_feeder_survives_a_drop_and_delivers_the_next_payload(caplog):
    """CPython drops the offending item and keeps the feeder running. Pinned because
    the surrounding accounting relies on it: if the feeder died, every later unit on
    the same queue would be reported missing too, and the diagnosis would point at
    the wrong units."""
    q = ResultQueue()
    try:
        with caplog.at_level(logging.ERROR, logger="MultiProcessManager"):
            q.put((11, threading.Lock()))
            assert _wait_for(lambda: caplog.records)
            q.put((12, "survivor"))
            assert _wait_for(lambda: not q.empty()), "feeder stopped after the drop"
        assert q.get(timeout=5) == (12, "survivor")
    finally:
        q.close()


# ---------------------------------------------------------- the production channel


def _child_puts_unpicklable(q, stderr_path, unit_id):
    """Runs in a forked child, reproducing what a component worker does: configure
    logging to its own stderr, put a payload, flush, return.
    """
    import os
    import sys

    # Bind the child's stderr to the file at both levels: the fd, as LRAA's
    # per-contig redirection does, and sys.stderr, which is what
    # logging.StreamHandler() resolves against.
    stream = open(stderr_path, "w")
    os.dup2(stream.fileno(), 2)
    sys.stderr = stream

    root = logging.getLogger()
    for handler in list(root.handlers):
        root.removeHandler(handler)
    root.addHandler(logging.StreamHandler(stream))
    root.setLevel(logging.INFO)

    q.put((unit_id, threading.Lock()))

    # The flush the worker performs. Without it the pickling happens during
    # interpreter shutdown, where Queue._feed skips the error hook.
    q.close()
    q.join_thread()

    logging.shutdown()
    stream.flush()


def test_the_report_reaches_a_real_child_s_stderr(tmp_path):
    """The point of the override, tested through the channel that carries it.

    A worker logs to its own stderr, which LRAA redirects to a per-contig file. The
    hook must fire there -- and it only does if the payload is flushed while the
    child is still alive, because `Queue._feed` takes an `if is_exiting(): return`
    branch during shutdown and never calls the hook at all. An in-process test
    cannot catch that: the main thread is still running, so is_exiting() is False
    and the hook fires whether or not the flush exists.
    """
    ctx = multiprocessing.get_context("fork")
    stderr_path = tmp_path / "worker.err.log"

    q = ResultQueue(ctx=ctx)
    try:
        child = ctx.Process(
            target=_child_puts_unpicklable, args=(q, str(stderr_path), 42)
        )
        child.start()
        child.join(timeout=30)
        assert child.exitcode == 0, "the worker exits cleanly despite losing its payload"
    finally:
        q.close()

    captured = stderr_path.read_text()
    assert "failed to serialise the payload for unit 42" in captured, captured
    assert "LOST" in captured
    assert "cannot pickle" in captured
