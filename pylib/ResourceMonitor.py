#!/usr/bin/env python3
# Lightweight process resource monitor

import os
import time
import threading
import logging
from typing import Optional

try:
    import psutil
except Exception:  # pragma: no cover
    psutil = None

logger = logging.getLogger(__name__)


def _tree_cpu_seconds(proc):
    """Cumulative CPU seconds burned by ``proc`` and every descendant it has had.

    ``cpu_percent()`` cannot be used for the children. psutil measures it against
    a reference point stored ON THE Process OBJECT, and the children have to be
    re-enumerated every sample, so every child object is on its first call and
    every first call returns 0.0. That is an object-lifetime bug, not a sampling
    one: a shorter interval samples 0.0 more often.

    ``cpu_times`` avoids the reference point entirely. ``children_user`` and
    ``children_system`` come from cutime/cstime, which accumulate a descendant's
    whole subtree at the moment it is REAPED, so a worker that has already exited
    still counts. Live descendants are added directly, and they cannot be double
    counted: a process appears in the parent's cutime only once it is gone from
    the parent's children list.
    """

    def _own(p):
        t = p.cpu_times()
        return (
            t.user
            + t.system
            + getattr(t, "children_user", 0.0)
            + getattr(t, "children_system", 0.0)
        )

    total = _own(proc)
    for child in proc.children(recursive=True):
        try:
            total += _own(child)
        except Exception:
            continue
    return total


def _tree_rss_bytes(proc, include_children):
    """RSS of the process and, separately, of its live descendants."""

    own = proc.memory_info().rss if proc.is_running() else 0
    children = 0
    if include_children:
        for child in proc.children(recursive=True):
            try:
                if child.is_running():
                    children += getattr(child.memory_info(), "rss", 0)
            except Exception:
                continue
    return own, children


class ResourceMonitor:
    """
    Periodically samples CPU%, RSS, and children aggregate for the current process
    and writes a simple TSV log. Intended to be low-overhead and optional.

    Columns: epoch_ts, elapsed_sec, rss_mb, cpu_percent, rss_mb_children,
    cpu_percent_children, cpu_seconds_total, note

    A row summarises the interval that ended at it rather than the instant it was
    taken. RSS is a high-water mark over that interval, sampled at
    ``sample_interval_sec``, and the two RSS columns are the split observed at
    the instant the tree peaked, so they still sum to a total that really
    occurred. CPU is the interval's mean, from cumulative CPU-seconds.

    This matters because the median LRAA work unit is shorter than the default
    write interval: on a spot reading, most units reported whatever happened to
    be true at one arbitrary moment, and a unit shorter than one interval
    reported a peak of nothing.
    """

    def __init__(
        self,
        output_tsv: str,
        interval_sec: float = 15.0,
        include_children: bool = True,
        note: Optional[str] = None,
        sample_interval_sec: Optional[float] = None,
    ):
        self.output_tsv = output_tsv
        self.interval_sec = max(0.5, float(interval_sec))
        # Sampling is decoupled from writing so a long write interval keeps a
        # small log without turning the peak back into a spot reading.
        self.sample_interval_sec = min(
            self.interval_sec,
            1.0 if sample_interval_sec is None else max(0.05, float(sample_interval_sec)),
        )
        self.include_children = include_children
        self._stop_evt = threading.Event()
        self._thread = None
        self._start_time = None
        self._note = note or ""

        if psutil is None:
            logger.warning(
                "psutil is not available; resource monitoring will be disabled for this run"
            )

    def start(self):
        if psutil is None:
            return
        if self._thread is not None:
            return
        self._start_time = time.time()

        try:
            with open(self.output_tsv, "wt") as fh:
                fh.write(
                    "\t".join(
                        [
                            "epoch_ts",
                            "elapsed_sec",
                            "rss_mb",
                            "cpu_percent",
                            "rss_mb_children",
                            "cpu_percent_children",
                            "cpu_seconds_total",
                            "note",
                        ]
                    )
                    + "\n"
                )
        except Exception as e:
            logger.warning(f"Failed opening resource monitor log {self.output_tsv}: {e}")

        self._thread = threading.Thread(target=self._run_loop, name="ResourceMonitor", daemon=True)
        self._thread.start()

    def stop(self, timeout: Optional[float] = 5.0):
        if self._thread is None:
            return
        self._stop_evt.set()
        try:
            self._thread.join(timeout=timeout)
        except Exception:
            pass
        self._thread = None

    def _run_loop(self):
        if psutil is None:
            return
        proc = psutil.Process(os.getpid())

        try:
            last_cpu_s = _tree_cpu_seconds(proc)
            last_own_cpu_s = sum(proc.cpu_times()[:2])
        except Exception:
            last_cpu_s = last_own_cpu_s = 0.0
        last_write = self._start_time or time.time()

        # peak of the TOTAL tree, with the split observed at that same instant,
        # so the two columns still sum to a footprint that really occurred
        peak_total = -1
        peak_own = peak_children = 0

        while True:
            stopping = self._stop_evt.is_set()
            now = time.time()
            try:
                own_rss, children_rss = _tree_rss_bytes(proc, self.include_children)
                if own_rss + children_rss > peak_total:
                    peak_total = own_rss + children_rss
                    peak_own, peak_children = own_rss, children_rss
            except Exception:
                pass

            # A final row is always written, so a unit shorter than one interval
            # still records its peak and its elapsed time instead of nothing.
            if stopping or (now - last_write) >= self.interval_sec:
                try:
                    cpu_s = _tree_cpu_seconds(proc)
                    own_cpu_s = sum(proc.cpu_times()[:2])
                    window = max(now - last_write, 1e-9)
                    cpu_pct = 100.0 * (own_cpu_s - last_own_cpu_s) / window
                    cpu_pct_children = (
                        100.0 * ((cpu_s - last_cpu_s) - (own_cpu_s - last_own_cpu_s)) / window
                    )
                    with open(self.output_tsv, "at") as fh:
                        fh.write(
                            "\t".join(
                                [
                                    f"{now:.3f}",
                                    f"{now - (self._start_time or now):.3f}",
                                    f"{max(peak_own, 0) / (1024 * 1024):.3f}",
                                    f"{max(cpu_pct, 0.0):.2f}",
                                    f"{max(peak_children, 0) / (1024 * 1024):.3f}",
                                    f"{max(cpu_pct_children, 0.0):.2f}",
                                    f"{cpu_s:.3f}",
                                    self._note,
                                ]
                            )
                            + "\n"
                        )
                    last_cpu_s, last_own_cpu_s, last_write = cpu_s, own_cpu_s, now
                    peak_total = -1
                    peak_own = peak_children = 0
                except Exception:
                    # swallow all to keep the monitor from crashing the app
                    pass

            if stopping:
                return
            self._stop_evt.wait(self.sample_interval_sec)



def summarize_resource_log(path):
    """Peak and mean over one resources.tsv. None if it cannot be read.

    One implementation, because there used to be two: LRAA carried an inline
    copy for the per-worker summary and util/misc/summarize_resource_logs.py
    had its own. They agreed until one of them was fixed, at which point the
    per-worker summary kept reporting a 53 second work unit as lasting 0.0 s.
    """

    peak_rss = peak_cpu = sum_rss = sum_cpu = 0.0
    n = 0
    first_ts = last_ts = None
    max_elapsed = 0.0
    try:
        with open(path, "rt") as fh:
            header = next(fh).rstrip("\n").split("\t")
            idx = {k: i for i, k in enumerate(header)}
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                try:
                    ts = float(parts[idx.get("epoch_ts", 0)])
                    max_elapsed = max(max_elapsed, float(parts[idx.get("elapsed_sec", 1)]))
                    rss = float(parts[idx.get("rss_mb", 2)])
                    cpu = float(parts[idx.get("cpu_percent", 3)])
                    rss_ch = float(parts[idx.get("rss_mb_children", 4)])
                    cpu_ch = float(parts[idx.get("cpu_percent_children", 5)])
                except Exception:
                    continue
                peak_rss = max(peak_rss, rss + rss_ch)
                peak_cpu = max(peak_cpu, cpu + cpu_ch)
                sum_rss += rss + rss_ch
                sum_cpu += cpu + cpu_ch
                n += 1
                if first_ts is None:
                    first_ts = ts
                last_ts = ts
    except Exception:
        return None

    # elapsed_sec of the last row, not last_ts - first_ts. A unit shorter than
    # one write interval produces a single row, and differencing one timestamp
    # against itself reports every such unit as lasting no time at all.
    duration = max_elapsed
    if not duration and first_ts is not None and last_ts is not None:
        duration = last_ts - first_ts

    return {
        "file": path,
        "samples": n,
        "duration_sec": duration,
        "peak_rss_mb_total": peak_rss,
        "avg_rss_mb_total": (sum_rss / n) if n else 0.0,
        "peak_cpu_percent_total": peak_cpu,
        "avg_cpu_percent_total": (sum_cpu / n) if n else 0.0,
    }

__all__ = ["ResourceMonitor", "summarize_resource_log"]
