#!/usr/bin/env python3

import time
import multiprocessing
import random
import logging
import re

logger = logging.getLogger(__name__)

SLEEPTIME = 0.1

MPM_DEBUG = False


class MultiProcessManager:

    def __init__(self, num_parallel_processes, queue=None):

        self.num_parallel_processes = num_parallel_processes
        self.num_running = 0
        self.num_successes = 0
        self.num_errors = 0
        self.process_list = set()
        self.queue = queue
        self.captured_queue_contents = list()
        self.process_name_to_start_time = dict()
        # Track failure details for downstream diagnostics
        self.failures = []  # list of dicts: {name, pid, exitcode, duration_sec}

    def launch_process(self, process):

        time_start = time.time()
        if MPM_DEBUG:
            logger.info("-launching process")

        self._screen_running_processes()

        if self.num_running >= self.num_parallel_processes:
            self.wait_for_open_slot()

        process.start()
        self.process_name_to_start_time[process.name] = time.time()
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

        if self.queue is not None and not self.queue.empty():
            if MPM_DEBUG:
                logger.info("\t-reaping queue")
            while not self.queue.empty():
                if MPM_DEBUG:
                    logger.info("\t\t-try reaping entry from queue")

                entry = self.queue.get()

                if MPM_DEBUG:
                    logger.info("\t\t\t-reaped entry")

                self.captured_queue_contents.append(entry)

            if MPM_DEBUG:
                logger.info("\t\t-done reaping queue")

        completed_processes = set()

        for process in self.process_list:
            if process.is_alive():
                if MPM_DEBUG:
                    logger.info("\t-process {} is alive.".format(process.name))
            else:
                if MPM_DEBUG:
                    logger.info("\t-process {} is finished.".format(process.name))
                completed_processes.add(process)

        if completed_processes:
            if MPM_DEBUG:
                logger.info(
                    "\t-processing {} completed processes.".format(
                        len(completed_processes)
                    )
                )

            for completed_process in completed_processes:

                if MPM_DEBUG:
                    logger.info(
                        "\t\t\t<joining process {}>".format(completed_process.name)
                    )
                completed_process.join()
                if MPM_DEBUG:
                    logger.info(
                        "\t\t\t<joined process {} having exit code {}>".format(
                            completed_process.name, completed_process.exitcode
                        )
                    )
                if completed_process.exitcode == 0:
                    self.num_successes += 1
                else:
                    self.num_errors += 1
                    # record failure meta for later reporting
                    try:
                        start_ts = self.process_name_to_start_time.get(
                            completed_process.name
                        )
                        duration = None
                        if start_ts is not None:
                            duration = max(0.0, time.time() - start_ts)
                        self.failures.append(
                            {
                                "name": completed_process.name,
                                "pid": getattr(completed_process, "pid", None),
                                "exitcode": completed_process.exitcode,
                                "duration_sec": duration,
                            }
                        )
                    except Exception:
                        pass
                    if MPM_DEBUG:
                        logger.info("-captured a failed process")

                # remove completed process from process list
                self.process_list.remove(completed_process)
                self.num_running -= 1

        if MPM_DEBUG:
            logger.info("\t-done screening running processes.")

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

        return self.num_errors

    def summarize_status(self):
        return "{} jobs succeeded & {} jobs failed".format(
            self.num_successes, self.num_errors
        )

    def retrieve_queue_contents(self, clear=False):
        """Return captured queue payloads.

        If clear is True, the internal buffer is emptied after returning
        the collected objects so callers can drain incrementally without
        reprocessing prior results.
        """
        contents = list(self.captured_queue_contents)
        if clear:
            self.captured_queue_contents.clear()
        return contents

    def terminate_all_processes(self):
        """Best-effort termination of all tracked child processes.

        Strategy:
        - send terminate() (SIGTERM on POSIX)
        - join with short timeout
        - if still alive, escalate to SIGKILL when possible
        """
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
        Each dict contains: name, pid, exitcode, duration_sec (may be None).
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
        q.put(id)

    q = multiprocessing.Queue()
    mpm = MultiProcessManager(num_parallel_processes, q)

    set_debug()

    for i in range(num_total_processes):

        p = multiprocessing.Process(target=runner, args=(i, q))

        mpm.launch_process(p)

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
