#!/usr/bin/env python3

import time
import multiprocessing
import multiprocessing.queues
import random
import logging
import re

from queue import Empty as QueueEmpty

logger = logging.getLogger(__name__)

SLEEPTIME = 0.1

MPM_DEBUG = False


class ResultQueue(multiprocessing.queues.Queue):
    """A result queue that says so when a payload fails to serialise.

    `Queue.put()` hands the object to a background feeder thread which pickles it
    and writes it to the pipe.  When that pickling raises -- an unpicklable
    member, or a RecursionError on a deep object graph -- CPython drops that one
    item, calls `_on_queue_feeder_error`, and lets the feeder carry on.  The
    worker is never told, so it exits 0 having silently delivered nothing.  The
    default handler is `traceback.print_exc()`, which in this pipeline lands on a
    worker's stderr inside a per-contig temp directory that is deleted when the
    run succeeds -- so the one message explaining the loss is written and then
    thrown away.

    Overriding the hook puts the failure in the logger, tagged with the unit id
    carried by the payload, so it survives next to the WorkUnitAccountingError
    that the drop will trigger.  Note this only makes the loss legible: the item
    is still gone, and preventing that needs the payload off the queue entirely.
    """

    def __init__(self, *args, ctx=None, **kwargs):
        if ctx is None:
            ctx = multiprocessing.get_context()
        super().__init__(*args, ctx=ctx, **kwargs)

    @staticmethod
    def _on_queue_feeder_error(exception, obj):
        # Payloads are (unit_id, result); recover the id when the shape holds so
        # the message names the unit that lost its result.
        unit_id = None
        if isinstance(obj, tuple) and len(obj) >= 1:
            unit_id = obj[0]

        logger.error(
            "result queue failed to serialise the payload for unit %s: %s: %s. "
            "That result is LOST -- the worker will still exit 0, and the loss "
            "surfaces as a work-unit accounting failure.",
            "unknown" if unit_id is None else unit_id,
            type(exception).__name__,
            exception,
            exc_info=exception,
        )


class WorkUnitAccountingError(RuntimeError):
    """Raised when the payloads collected off the queue do not match the work
    units that were submitted: a unit that exited cleanly but returned nothing,
    a unit that returned twice, or a payload tagged with an id nobody submitted.

    A lost payload is otherwise invisible -- the child exits 0, the parent joins
    it, and the run reports success with a short result set -- so this is
    deliberately fatal rather than a warning.
    """


def _sorted_ids(unit_ids):
    """Order unit ids for a diagnostic message, tolerating exotic id types."""
    try:
        return sorted(unit_ids)
    except TypeError:
        return sorted(unit_ids, key=repr)


class MultiProcessManager:
    """Rolling scheduler for a fixed-width pool of multiprocessing.Process workers.

    Each submitted process is a *work unit* identified by a caller-supplied
    ``unit_id``.  When the manager was given a result queue, workers put
    ``(unit_id, payload)`` on it, the manager keys the payloads by ``unit_id``,
    and ``wait_for_remaining_processes`` refuses to return unless the set of
    payloads exactly matches the set of units that exited cleanly.
    """

    def __init__(self, num_parallel_processes, queue=None):

        self.num_parallel_processes = num_parallel_processes
        self.num_running = 0
        self.num_successes = 0
        self.num_errors = 0
        self.process_list = set()
        self.queue = queue

        # unit_id -> payload, for units drained but not yet handed to the caller
        self.captured_unit_results = dict()

        # Accounting sets.  returned_unit_ids is the full history and is *not*
        # cleared by retrieve_queue_contents(clear=True), which callers use to
        # drain incrementally.
        self.submitted_unit_ids = set()
        self.returned_unit_ids = set()
        self.succeeded_unit_ids = set()
        self.failed_unit_ids = set()
        self.duplicate_unit_ids = list()
        self.unknown_unit_ids = list()

        # Per-process bookkeeping.  Keyed by the process object rather than by
        # process.name: LRAA names every component worker after its mpg token,
        # and a name-keyed start time is overwritten by the next launch, which
        # made duration_sec the age of the most recent launch instead of the
        # age of the unit being reported.
        self.process_start_time = dict()
        self.process_to_unit_id = dict()
        self.launched_process_names = set()

        # Set by terminate_all_processes: the run is being torn down, so missing
        # payloads are expected and the completeness check is suppressed.
        self.aborted = False

        # Track failure details for downstream diagnostics
        self.failures = []  # list of dicts: {unit_id, name, pid, exitcode, duration_sec}

    def launch_process(self, process, unit_id=None):
        """Start `process` as the work unit identified by `unit_id`.

        `unit_id` is required whenever the manager owns a result queue, because
        it is the only thing tying a payload back to the unit that produced it.
        It must be unique within the manager and orderable against the other
        ids, since it also fixes the order results are returned in.
        """

        time_start = time.time()
        if MPM_DEBUG:
            logger.info("-launching process")

        if unit_id is None:
            if self.queue is not None:
                raise ValueError(
                    "unit_id is required when the manager collects results through a queue"
                )
        else:
            if unit_id in self.submitted_unit_ids:
                raise WorkUnitAccountingError(
                    "work unit {!r} was submitted more than once".format(unit_id)
                )
            # Callers name workers after the region they cover, which is not
            # unique per unit; appending the unit id keeps process names
            # distinguishable in logs and in the failure records.
            process.name = "{}:{}".format(process.name, unit_id)
            if process.name in self.launched_process_names:
                raise WorkUnitAccountingError(
                    "process name {!r} is already in use by another work unit".format(
                        process.name
                    )
                )

        self._screen_running_processes()

        if self.num_running >= self.num_parallel_processes:
            self.wait_for_open_slot()

        process.start()
        self.process_start_time[process] = time.time()
        if unit_id is not None:
            self.submitted_unit_ids.add(unit_id)
            self.process_to_unit_id[process] = unit_id
            self.launched_process_names.add(process.name)
        self.process_list.add(process)
        self.num_running += 1
        time_end = time.time()
        # Include contig/strand context when available (process.name like "chr1+:123-456")
        elapsed = time_end - time_start
        ctx = None
        try:
            name = getattr(process, "name", "") or ""
            m = re.match(r"^(.+?)([+-])(?::|$)", name)
            if m:
                ctx = f"{m.group(1)}{m.group(2)}"
        except Exception:
            ctx = None
        if ctx:
            logger.info("[%s] time to launch process: %.3fs", ctx, elapsed)
        else:
            logger.info("time to launch process: %.3fs", elapsed)

    def wait_for_open_slot(self):
        if MPM_DEBUG:
            logger.info("-waiting for open slot")

        counter = 0
        while self.num_running >= self.num_parallel_processes:
            if MPM_DEBUG:
                logger.info("\twaiting for open slot round({})".format(counter))
            counter += 1
            self._screen_running_processes()
            time.sleep(SLEEPTIME)

    def _drain_queue(self):
        """Read everything currently readable off the result queue.

        Never gated on Queue.empty(): empty() reports whether bytes have reached
        the pipe, which lags the worker's put() by however long the feeder
        thread takes to be scheduled, so a False-y empty() is not evidence that
        a queued payload does not exist.  get_nowait() until Empty is the only
        honest way to ask.
        """

        if self.queue is None:
            return

        while True:
            try:
                entry = self.queue.get_nowait()
            except QueueEmpty:
                break

            if MPM_DEBUG:
                logger.info("\t\t\t-reaped entry from queue")

            self._record_queue_entry(entry)

    def _record_queue_entry(self, entry):

        try:
            unit_id, payload = entry
        except (TypeError, ValueError):
            raise WorkUnitAccountingError(
                "queue payload is not a (unit_id, payload) pair: {!r}".format(entry)
            )

        if unit_id in self.returned_unit_ids:
            # Keep the first payload; the run is going to abort at the
            # completeness check anyway, and dropping the first would make the
            # surviving result depend on arrival order.
            self.duplicate_unit_ids.append(unit_id)
            return

        if unit_id not in self.submitted_unit_ids:
            self.unknown_unit_ids.append(unit_id)

        self.returned_unit_ids.add(unit_id)
        self.captured_unit_results[unit_id] = payload

    def _screen_running_processes(self):

        if MPM_DEBUG:
            logger.info(
                "-screening {} running processes".format(len(self.process_list))
            )

        # Order is load-bearing.  A child flushes its queue payload and then
        # exits, so a process already observed dead has necessarily handed its
        # bytes to the pipe, while a process observed alive may flush at any
        # instant.  Snapshotting the dead set BEFORE draining therefore
        # guarantees that every payload belonging to a process reaped below is
        # read here.  Draining first left a window in which a child could flush
        # and exit between the drain and the liveness check, and be joined and
        # counted a success with its payload still unread -- silently lost when
        # it was the last unit, and merely late when it was not.
        #
        # Only processes observed not-alive are joined.  Joining a live producer
        # deadlocks: a payload larger than the pipe buffer blocks the child's
        # feeder thread until the parent reads, and LRAA component payloads are
        # routinely larger than that.
        completed_processes = [
            process for process in self.process_list if not process.is_alive()
        ]

        if MPM_DEBUG:
            for process in self.process_list:
                logger.info(
                    "\t-process {} is {}.".format(
                        process.name,
                        (
                            "finished"
                            if process in completed_processes
                            else "alive"
                        ),
                    )
                )

        self._drain_queue()

        if not completed_processes:
            if MPM_DEBUG:
                logger.info("\t-done screening running processes.")
            return

        if MPM_DEBUG:
            logger.info(
                "\t-processing {} completed processes.".format(len(completed_processes))
            )

        # Reap in name order rather than set-iteration order so the success /
        # failure bookkeeping, and the failure log it feeds, do not vary with
        # object hashing.
        completed_processes.sort(key=lambda process: process.name)

        for completed_process in completed_processes:
            self._reap_completed_process(completed_process)

        if MPM_DEBUG:
            logger.info("\t-done screening running processes.")

    def _reap_completed_process(self, completed_process):

        if MPM_DEBUG:
            logger.info("\t\t\t<joining process {}>".format(completed_process.name))
        completed_process.join()
        if MPM_DEBUG:
            logger.info(
                "\t\t\t<joined process {} having exit code {}>".format(
                    completed_process.name, completed_process.exitcode
                )
            )

        unit_id = self.process_to_unit_id.pop(completed_process, None)
        start_ts = self.process_start_time.pop(completed_process, None)
        duration = None
        if start_ts is not None:
            duration = max(0.0, time.time() - start_ts)

        if completed_process.exitcode == 0:
            self.num_successes += 1
            if unit_id is not None:
                self.succeeded_unit_ids.add(unit_id)
        else:
            self.num_errors += 1
            if unit_id is not None:
                self.failed_unit_ids.add(unit_id)
            # record failure meta for later reporting
            self.failures.append(
                {
                    "unit_id": unit_id,
                    "name": completed_process.name,
                    "pid": getattr(completed_process, "pid", None),
                    "exitcode": completed_process.exitcode,
                    "duration_sec": duration,
                }
            )
            if MPM_DEBUG:
                logger.info("-captured a failed process")

        # remove completed process from process list
        self.process_list.remove(completed_process)
        self.num_running -= 1

    def wait_for_remaining_processes(self):

        if MPM_DEBUG:
            logger.info("-waiting for remaining processes")

        while self.num_running > 0:
            if MPM_DEBUG:
                logger.info("-waiting on {} processes".format(self.num_running))
            self._screen_running_processes()
            time.sleep(SLEEPTIME)

        # Final drain, after the last unit has been reaped.  The pass that
        # reaped it drained after observing it dead, so this is belt-and-braces
        # rather than the primary guarantee, but it costs one non-blocking poll
        # and it is the only drain that is unconditionally last.
        self._drain_queue()

        if MPM_DEBUG:
            logger.info("-done waiting. All processes are completed")

        self.verify_units_accounted()

        return self.num_errors

    def verify_units_accounted(self):
        """Raise unless every work unit is accounted for exactly once.

        Call after the final drain.  A unit that exited 0 without producing a
        payload is a failure here, not a success: that is the whole point of the
        check, since it is exactly the shape a lost payload takes.
        """

        if self.queue is None or self.aborted:
            return

        problems = []

        unreaped = _sorted_ids(
            self.submitted_unit_ids - self.succeeded_unit_ids - self.failed_unit_ids
        )
        if unreaped:
            problems.append(
                "{} unit(s) never reaped: {}".format(len(unreaped), unreaped)
            )

        missing = _sorted_ids(self.succeeded_unit_ids - self.returned_unit_ids)
        if missing:
            problems.append(
                "{} unit(s) exited 0 but returned no payload: {}".format(
                    len(missing), missing
                )
            )

        if self.duplicate_unit_ids:
            duplicates = _sorted_ids(set(self.duplicate_unit_ids))
            problems.append(
                "{} unit(s) returned more than one payload: {}".format(
                    len(duplicates), duplicates
                )
            )

        if self.unknown_unit_ids:
            unknown = _sorted_ids(set(self.unknown_unit_ids))
            problems.append(
                "{} payload(s) tagged with never-submitted unit id(s): {}".format(
                    len(unknown), unknown
                )
            )

        if problems:
            raise WorkUnitAccountingError(
                "work unit accounting failed ({} submitted, {} returned): {}".format(
                    len(self.submitted_unit_ids),
                    len(self.returned_unit_ids),
                    "; ".join(problems),
                )
            )

    def summarize_status(self):
        return "{} jobs succeeded & {} jobs failed".format(
            self.num_successes, self.num_errors
        )

    def retrieve_queue_contents(self, clear=False):
        """Return captured payloads as (unit_id, payload) pairs ordered by unit_id.

        Ordering is by unit id, never by arrival: arrival order is a function of
        worker scheduling, and letting it reach the caller made the result
        sequence -- and everything order-sensitive downstream of it -- vary
        between runs on identical input.

        If clear is True, the internal buffer is emptied after returning
        the collected objects so callers can drain incrementally without
        reprocessing prior results.  The accounting history is retained.
        """
        contents = [
            (unit_id, self.captured_unit_results[unit_id])
            for unit_id in sorted(self.captured_unit_results)
        ]
        if clear:
            self.captured_unit_results.clear()
        return contents

    def terminate_all_processes(self):
        """Best-effort termination of all tracked child processes.

        Strategy:
        - send terminate() (SIGTERM on POSIX)
        - join with short timeout
        - if still alive, escalate to SIGKILL when possible
        """
        # Torn-down workers owe no payloads, so stop enforcing completeness.
        self.aborted = True

        # First pass: request termination
        for process in list(self.process_list):
            try:
                process.terminate()
            except Exception:
                pass

        # Second pass: wait briefly and escalate if needed
        import time as _time
        import os as _os
        import signal as _signal

        deadline = _time.time() + 5.0  # up to ~5s total waiting across children
        for process in list(self.process_list):
            try:
                # Join remaining time budget but at least a short slice
                remaining = max(0.1, deadline - _time.time())
                process.join(timeout=remaining)
            except Exception:
                pass
            # Escalate if still alive
            try:
                if process.is_alive():
                    pid = getattr(process, "pid", None)
                    if pid:
                        try:
                            _os.kill(pid, _signal.SIGKILL)
                        except Exception:
                            pass
            except Exception:
                pass

        logger.info("-terminated all processes")

    def get_failures(self):
        """Return a list of failure dicts collected during execution.
        Each dict contains: unit_id, name, pid, exitcode, duration_sec (may be None).
        """
        try:
            return list(self.failures)
        except Exception:
            return []


def set_debug():
    logger.setLevel(logging.DEBUG)
    global MPM_DEBUG
    MPM_DEBUG = True


def test_mpm(num_parallel_processes=8, num_total_processes=100):

    def runner(id, q):
        logger.info("running id:{}".format(id))
        x = id / (id % 10)  # should error as div-by-zero on occasion
        time.sleep(random.randint(0, 10))
        q.put((id, id))

    q = multiprocessing.Queue()
    mpm = MultiProcessManager(num_parallel_processes, q)

    set_debug()

    for i in range(num_total_processes):

        p = multiprocessing.Process(target=runner, args=(i, q))

        mpm.launch_process(p, unit_id=i)

    # Every tenth id divides by zero and its worker dies; those units are
    # counted as failures, and the completeness check only demands payloads
    # from units that exited cleanly, so this is expected to return normally.
    mpm.wait_for_remaining_processes()

    logger.info("Test job completed.")
    logger.info("Captured queue contents are: {}".format(mpm.retrieve_queue_contents()))
    logger.info(mpm.summarize_status())


if __name__ == "__main__":

    # run test
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s : %(levelname)s : %(message)s",
        datefmt="%H:%M:%S",
    )

    logger.setLevel(logging.DEBUG)
    test_mpm()
