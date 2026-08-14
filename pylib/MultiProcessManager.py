#!/usr/bin/env python3

import errno
import logging
import multiprocessing
import os
import pickle
import random
import re
import shutil
import tempfile
import time

from multiprocessing.reduction import ForkingPickler

logger = logging.getLogger(__name__)

SLEEPTIME = 0.1

MPM_DEBUG = False

SHARD_TMP_SUFFIX = ".tmp"
SHARD_DONE_SUFFIX = ".done"


def _safe_token(token):
    """Reduce a token to characters that are safe in a path component."""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(token))


class ShardPublicationError(RuntimeError):
    """Raised in the worker when a result cannot be written to its shard file.

    Serialisation happens in the worker's own call stack now, so a payload that
    cannot be pickled raises here instead of being dropped by a queue's feeder
    thread.  The worker dies with the unit id in the message, which makes the
    unit a reported failure rather than a clean exit that delivered nothing.
    """


class WorkUnitAccountingError(RuntimeError):
    """Raised when the shards collected do not match the work units submitted: a
    unit that exited cleanly but published nothing, a shard whose contents
    disagree with the unit it was published for, or a shard for an id nobody
    submitted.

    A lost result is otherwise invisible -- the child exits 0, the parent joins
    it, and the run reports success with a short result set -- so this is
    deliberately fatal rather than a warning.
    """


def _sorted_ids(unit_ids):
    """Order unit ids for a diagnostic message, tolerating exotic id types."""
    try:
        return sorted(unit_ids)
    except TypeError:
        return sorted(unit_ids, key=repr)


class ShardStore:
    """Per-unit result files, published by atomic rename.

    A worker writes its result to ``unit.<id>.<pid>.tmp``, fsyncs it, closes it,
    and renames it to ``unit.<id>.done``.  Rename within a directory is atomic,
    so a ``.done`` file is always a complete result and a surviving ``.tmp`` is
    always a failed unit.  That is the point, and it replaces a
    multiprocessing.Queue for two reasons the queue could not be fixed for:

      * ``Queue.put()`` pickles on a background feeder thread.  A payload that
        fails to serialise there is dropped, the error hook is not even called
        during interpreter shutdown, and the worker still exits 0.  Writing a
        shard pickles inline, so the failure lands in the worker's own traceback
        and the unit becomes a failure with an exit code.

      * a worker killed mid-write leaves a truncated frame in the pipe, and the
        parent blocks in ``_recv_bytes`` forever: no EOFError arrives while other
        writers hold the pipe open, and the read is covered by no timeout.  A
        worker killed mid-write here leaves a ``.tmp`` that is never renamed; the
        parent looks for a ``.done``, does not find one, and says so.  There is
        no shared pipe left to wedge.

    The directory is created per invocation by ``create()``, so nothing from an
    earlier round, an earlier run, or a concurrent job can be found in it.  Each
    shard also carries its unit id and the caller's token inside the file, and
    ``load()`` refuses a shard whose recorded id is not the one being asked for:
    a file name is not evidence about contents, and this project has twice been
    bitten by a name that did not describe what it held.
    """

    def __init__(self, directory):
        self.directory = str(directory)

    @classmethod
    def create(cls, token, base_dir=None):
        """Make a fresh shard directory whose name records what it holds.

        `token` should name everything that determines the contents -- for LRAA
        the contig, the strand, the splice-type round and the invocation -- so a
        directory left behind by a failed run can be identified on sight.  The
        random suffix mkdtemp appends is what makes two concurrent jobs, or two
        rounds inside one process, unable to collide.
        """
        if base_dir is None:
            base_dir = os.environ.get("LRAA_TMP_DIR") or os.getcwd()
        os.makedirs(base_dir, exist_ok=True)
        directory = tempfile.mkdtemp(
            prefix="__lraa_shards.{}.".format(_safe_token(token)), dir=base_dir
        )
        return cls(directory)

    def unit_basename(self, unit_id):
        """Shard name for `unit_id`.

        Integer ids are zero padded so that a directory listing sorts the way the
        results are merged; the merge itself orders by unit id, never by name.
        """
        if isinstance(unit_id, int) and not isinstance(unit_id, bool):
            return "unit.{:08d}".format(unit_id)
        return "unit.{}".format(_safe_token(unit_id))

    def done_path(self, unit_id):
        return os.path.join(
            self.directory, self.unit_basename(unit_id) + SHARD_DONE_SUFFIX
        )

    def tmp_path(self, unit_id):
        # The writer's pid is in the name so a surviving .tmp names the process
        # that died holding it, and so two writers can never share one scratch
        # file.
        return os.path.join(
            self.directory,
            "{}.{}{}".format(
                self.unit_basename(unit_id), os.getpid(), SHARD_TMP_SUFFIX
            ),
        )

    def publish(self, unit_id, payload, token=None):
        """Write, fsync, close, then rename.  Returns the published path.

        Called in the worker.  Everything before the rename is invisible to the
        parent; the rename is the single instant at which the result exists.
        """

        tmp_path = self.tmp_path(unit_id)
        done_path = self.done_path(unit_id)
        frame = {"unit_id": unit_id, "token": token, "payload": payload}

        try:
            with open(tmp_path, "wb") as fh:
                # The serialiser the queue used: ForkingPickler at its default
                # protocol is what multiprocessing.queues hands its feeder
                # thread.  This changes where pickling happens and when it is
                # reported, not how the payload is encoded.
                ForkingPickler(fh).dump(frame)
                fh.flush()
                os.fsync(fh.fileno())
        except Exception as exc:
            self._unlink_quietly(tmp_path)
            logger.error(
                "work unit %s (%s) could not be written to its result shard %s: %s: %s. "
                "The unit fails here, in its own stack, instead of exiting 0 having "
                "delivered nothing.",
                unit_id,
                "no token" if token is None else token,
                tmp_path,
                type(exc).__name__,
                exc,
                exc_info=exc,
            )
            raise ShardPublicationError(
                "work unit {} ({}) could not serialise its result to {}: {}: {}".format(
                    unit_id,
                    "no token" if token is None else token,
                    tmp_path,
                    type(exc).__name__,
                    exc,
                )
            ) from exc

        os.rename(tmp_path, done_path)
        self._fsync_directory()

        if MPM_DEBUG:
            logger.info("\t\t\t-published shard %s", done_path)

        return done_path

    def load(self, unit_id, remove=False):
        """Read the payload published for `unit_id`.

        Raises WorkUnitAccountingError rather than a bare pickle error, so a
        shard that is unreadable, or that holds someone else's result, names the
        unit it was asked about.
        """

        path = self.done_path(unit_id)

        try:
            with open(path, "rb") as fh:
                frame = pickle.load(fh)
        except Exception as exc:
            raise WorkUnitAccountingError(
                "result shard {} for unit {!r} could not be read: {}: {}".format(
                    path, unit_id, type(exc).__name__, exc
                )
            ) from exc

        if (
            not isinstance(frame, dict)
            or "unit_id" not in frame
            or "payload" not in frame
        ):
            raise WorkUnitAccountingError(
                "result shard {} for unit {!r} is not a result frame: {}".format(
                    path, unit_id, type(frame).__name__
                )
            )

        if frame["unit_id"] != unit_id:
            raise WorkUnitAccountingError(
                "result shard {} was published for unit {!r} but is being read as "
                "unit {!r}; the shard's contents disagree with its name".format(
                    path, frame["unit_id"], unit_id
                )
            )

        payload = frame["payload"]

        if remove:
            self._unlink_quietly(path)

        return payload

    def survey(self, submitted_unit_ids=(), remove_orphans=False):
        """Report what is in the shard directory that should not be.

        Returns a dict of `orphan_tmp` -- scratch files never renamed, i.e. units
        that died mid-write -- and `unexpected_done` -- published shards for ids
        nobody submitted.  Orphans are counted and named before removal, because
        a nonzero count is the signature of a writer that was killed.
        """

        expected = {
            self.unit_basename(unit_id) + SHARD_DONE_SUFFIX
            for unit_id in submitted_unit_ids
        }

        orphan_tmp = []
        unexpected_done = []

        try:
            entries = os.listdir(self.directory)
        except FileNotFoundError:
            entries = []

        for entry in sorted(entries):
            if entry.endswith(SHARD_TMP_SUFFIX):
                orphan_tmp.append(entry)
            elif entry.endswith(SHARD_DONE_SUFFIX) and entry not in expected:
                unexpected_done.append(entry)

        if remove_orphans:
            for entry in orphan_tmp:
                self._unlink_quietly(os.path.join(self.directory, entry))

        return {"orphan_tmp": orphan_tmp, "unexpected_done": unexpected_done}

    def discard(self):
        """Remove the shard directory.

        Called on the success path only: a directory left behind sits next to the
        failing worker's stderr log, which is where this pipeline already keeps
        the evidence for a run that did not finish.
        """
        shutil.rmtree(self.directory, ignore_errors=True)

    def _fsync_directory(self):
        # The rename is already atomic with respect to other processes; syncing
        # the directory is what makes it survive a machine that dies, which is
        # the same class of event that produces the SIGKILLs this transport
        # exists to tolerate.  Some filesystems refuse the call, so that refusal
        # is logged rather than raised.
        try:
            fd = os.open(self.directory, os.O_RDONLY)
        except OSError as exc:
            logger.debug("could not open %s to fsync: %s", self.directory, exc)
            return
        try:
            os.fsync(fd)
        except OSError as exc:
            if exc.errno not in (errno.EINVAL, errno.ENOTSUP, errno.EPERM):
                raise
            logger.debug("directory fsync unsupported on %s: %s", self.directory, exc)
        finally:
            os.close(fd)

    @staticmethod
    def _unlink_quietly(path):
        try:
            os.unlink(path)
        except FileNotFoundError:
            pass
        except OSError as exc:
            logger.warning("could not remove %s: %s", path, exc)


class MultiProcessManager:
    """Rolling scheduler for a fixed-width pool of multiprocessing.Process workers.

    Each submitted process is a *work unit* identified by a caller-supplied
    ``unit_id``.  When the manager was given a ShardStore, workers publish
    ``unit.<unit_id>.done`` into it, the manager requires that file of every unit
    that exits 0, and ``wait_for_remaining_processes`` refuses to return unless
    the set of published shards exactly matches the set of units that exited
    cleanly.
    """

    def __init__(self, num_parallel_processes, shard_store=None):

        self.num_parallel_processes = num_parallel_processes
        self.num_running = 0
        self.num_successes = 0
        self.num_errors = 0
        self.process_list = set()
        self.shard_store = shard_store

        # unit_id -> published shard path, for units collected but not yet handed
        # to the caller.  The payload stays on disk until the caller asks for it,
        # so the parent holds no copy of a result it has not been asked for.
        self.published_shard_paths = dict()

        # Accounting sets.  returned_unit_ids is the full history and is *not*
        # cleared by retrieve_unit_results(clear=True), which callers use to
        # collect incrementally.
        self.submitted_unit_ids = set()
        self.returned_unit_ids = set()
        self.succeeded_unit_ids = set()
        self.failed_unit_ids = set()

        # Units that exited 0 having published nothing: dicts of unit_id, name,
        # pid, exitcode, expected_path.
        self.missing_shards = list()

        # Filled by the final survey: scratch files from units that died
        # mid-write, and published shards for ids nobody submitted.
        self.orphan_shard_names = list()
        self.unexpected_shard_names = list()

        # Guards two distinct unit ids that would share one shard file name (1
        # and "1", say).  Distinct units are required to be distinct on disk.
        self.shard_name_to_unit_id = dict()

        # Per-process bookkeeping.  Keyed by the process object rather than by
        # process.name: LRAA names every component worker after its mpg token,
        # and a name-keyed start time is overwritten by the next launch, which
        # made duration_sec the age of the most recent launch instead of the
        # age of the unit being reported.
        self.process_start_time = dict()
        self.process_to_unit_id = dict()
        self.launched_process_names = set()

        # Set by terminate_all_processes: the run is being torn down, so missing
        # shards are expected and the completeness check is suppressed.
        self.aborted = False

        # Track failure details for downstream diagnostics
        self.failures = []  # list of dicts: {unit_id, name, pid, exitcode, duration_sec}

    def launch_process(self, process, unit_id=None):
        """Start `process` as the work unit identified by `unit_id`.

        `unit_id` is required whenever the manager owns a shard store, because it
        is the only thing tying a published shard back to the unit that produced
        it.  It must be unique within the manager and orderable against the other
        ids, since it also fixes the order results are returned in.
        """

        time_start = time.time()
        if MPM_DEBUG:
            logger.info("-launching process")

        if unit_id is None:
            if self.shard_store is not None:
                raise ValueError(
                    "unit_id is required when the manager collects results from shard files"
                )
        else:
            if unit_id in self.submitted_unit_ids:
                raise WorkUnitAccountingError(
                    "work unit {!r} was submitted more than once".format(unit_id)
                )
            if self.shard_store is not None:
                shard_name = self.shard_store.unit_basename(unit_id)
                collides_with = self.shard_name_to_unit_id.get(shard_name)
                if collides_with is not None:
                    raise WorkUnitAccountingError(
                        "work units {!r} and {!r} both publish to shard {!r}; distinct "
                        "units must be distinct on disk".format(
                            collides_with, unit_id, shard_name
                        )
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
            if self.shard_store is not None:
                self.shard_name_to_unit_id[
                    self.shard_store.unit_basename(unit_id)
                ] = unit_id
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

    def _screen_running_processes(self):

        if MPM_DEBUG:
            logger.info(
                "-screening {} running processes".format(len(self.process_list))
            )

        # Only processes observed not-alive are joined: joining a live worker
        # would block this loop until that worker finished, which is what the
        # rolling scheduler exists to avoid.
        #
        # Results are collected in _reap_completed_process, after join() -- not
        # here, and not on a schedule.  A worker publishes its shard by rename
        # before it exits, so a process that has been joined has necessarily
        # published, and no interleaving exists in which a result is present but
        # not yet visible.  Under the result queue one did: bytes reached the
        # pipe some time after put() returned, so the liveness snapshot had to be
        # taken before the drain or a child could flush and exit in the gap and
        # be joined as a success with its payload unread.
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
                self._collect_shard(unit_id, completed_process)
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

    def _collect_shard(self, unit_id, completed_process):
        """Require the shard of a unit that has exited 0.

        The payload is not read here.  Only its existence is established, which
        is all the accounting needs, and which keeps a result out of the parent's
        memory until the caller asks for it.
        """

        if self.shard_store is None:
            return

        done_path = self.shard_store.done_path(unit_id)

        if os.path.exists(done_path):
            self.returned_unit_ids.add(unit_id)
            self.published_shard_paths[unit_id] = done_path
            if MPM_DEBUG:
                logger.info("\t\t\t-collected shard for unit %s", unit_id)
            return

        # Exit 0 with no shard is the exactly-once violation.  It is reported
        # here as well as raised at the end of the run, so the message sits next
        # to the worker's own output in the log.
        self.missing_shards.append(
            {
                "unit_id": unit_id,
                "name": completed_process.name,
                "pid": getattr(completed_process, "pid", None),
                "exitcode": completed_process.exitcode,
                "expected_path": done_path,
            }
        )
        logger.error(
            "work unit %s (%s, pid %s) exited %s but published no result shard; "
            "expected %s",
            unit_id,
            completed_process.name,
            getattr(completed_process, "pid", None),
            completed_process.exitcode,
            done_path,
        )

    def wait_for_remaining_processes(self):

        if MPM_DEBUG:
            logger.info("-waiting for remaining processes")

        while self.num_running > 0:
            if MPM_DEBUG:
                logger.info("-waiting on {} processes".format(self.num_running))
            self._screen_running_processes()
            time.sleep(SLEEPTIME)

        if MPM_DEBUG:
            logger.info("-done waiting. All processes are completed")

        self.verify_units_accounted()

        return self.num_errors

    def verify_units_accounted(self):
        """Raise unless every work unit is accounted for exactly once.

        Call after the last process has been reaped.  A unit that exited 0
        without publishing a shard is a failure here, not a success: that is the
        whole point of the check, since it is exactly the shape a lost result
        takes.
        """

        if self.shard_store is None or self.aborted:
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
            detail_by_id = {record["unit_id"]: record for record in self.missing_shards}
            details = []
            for unit_id in missing:
                record = detail_by_id.get(unit_id)
                if record is None:
                    details.append("{!r}".format(unit_id))
                else:
                    details.append(
                        "{!r} (exit {}, expected {})".format(
                            unit_id, record["exitcode"], record["expected_path"]
                        )
                    )
            problems.append(
                "{} unit(s) exited 0 but published no result shard: {}".format(
                    len(missing), "; ".join(details)
                )
            )

        survey = self.shard_store.survey(self.submitted_unit_ids, remove_orphans=True)
        self.orphan_shard_names = survey["orphan_tmp"]
        self.unexpected_shard_names = survey["unexpected_done"]

        if self.orphan_shard_names:
            # Not fatal by itself.  A unit that died mid-write has already failed
            # on its exit code, and a scratch file left by something else must
            # not fail an otherwise complete run.  It is never mistaken for a
            # result -- only .done files are ever read -- but it is evidence that
            # a writer was killed, so it is reported by name and by count.
            logger.error(
                "%d orphaned shard scratch file(s) in %s -- a unit died mid-write and "
                "its result was never published: %s",
                len(self.orphan_shard_names),
                self.shard_store.directory,
                self.orphan_shard_names,
            )

        if self.unexpected_shard_names:
            problems.append(
                "{} published shard(s) for never-submitted unit id(s): {}".format(
                    len(self.unexpected_shard_names), self.unexpected_shard_names
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

    def retrieve_unit_results(self, clear=False):
        """Return collected payloads as (unit_id, payload) pairs ordered by unit_id.

        Ordering is by unit id, never by completion: completion order is a
        function of worker scheduling, and letting it reach the caller made the
        result sequence -- and everything order-sensitive downstream of it --
        vary between runs on identical input.

        If clear is True, each shard is deleted once its payload has been read,
        so callers can collect incrementally without reprocessing prior results
        and without the shard directory holding a whole contig at once.  The
        accounting history is retained.
        """
        results = []
        for unit_id in _sorted_ids(self.published_shard_paths):
            results.append((unit_id, self.shard_store.load(unit_id, remove=clear)))
        if clear:
            self.published_shard_paths.clear()
        return results

    def terminate_all_processes(self):
        """Best-effort termination of all tracked child processes.

        Strategy:
        - send terminate() (SIGTERM on POSIX)
        - join with short timeout
        - if still alive, escalate to SIGKILL when possible
        """
        # Torn-down workers owe no results, so stop enforcing completeness.
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

    def runner(id, shard_store):
        logger.info("running id:{}".format(id))
        x = id / (id % 10)  # should error as div-by-zero on occasion
        time.sleep(random.randint(0, 10))
        shard_store.publish(id, id, token="unit-{}".format(id))

    shard_store = ShardStore.create("mpm-demo")
    mpm = MultiProcessManager(num_parallel_processes, shard_store)

    set_debug()

    for i in range(num_total_processes):

        p = multiprocessing.Process(target=runner, args=(i, shard_store))

        mpm.launch_process(p, unit_id=i)

    # Every tenth id divides by zero and its worker dies; those units are
    # counted as failures, and the completeness check only demands shards from
    # units that exited cleanly, so this is expected to return normally.
    mpm.wait_for_remaining_processes()

    logger.info("Test job completed.")
    logger.info("Collected shard contents are: {}".format(mpm.retrieve_unit_results()))
    logger.info(mpm.summarize_status())
    shard_store.discard()


if __name__ == "__main__":

    # run test
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s : %(levelname)s : %(message)s",
        datefmt="%H:%M:%S",
    )

    logger.setLevel(logging.DEBUG)
    test_mpm()
