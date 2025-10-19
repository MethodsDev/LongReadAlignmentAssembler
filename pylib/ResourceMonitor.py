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


class ResourceMonitor:
    """
    Periodically samples CPU%, RSS, and children aggregate for the current process
    and writes a simple TSV log. Intended to be low-overhead and optional.

    Columns: epoch_ts, elapsed_sec, rss_mb, cpu_percent, rss_mb_children, cpu_percent_children, note
    """

    def __init__(
        self,
        output_tsv: str,
        interval_sec: float = 2.0,
        include_children: bool = True,
        note: Optional[str] = None,
    ):
        self.output_tsv = output_tsv
        self.interval_sec = max(0.5, float(interval_sec))
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
        # initialize CPU percent calculations
        try:
            p = psutil.Process(os.getpid())
            p.cpu_percent(interval=None)
            if self.include_children:
                for c in p.children(recursive=True):
                    try:
                        c.cpu_percent(interval=None)
                    except Exception:
                        pass
        except Exception:
            pass

        # write header
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
        while not self._stop_evt.is_set():
            try:
                epoch_ts = time.time()
                elapsed = epoch_ts - (self._start_time or epoch_ts)

                rss = proc.memory_info().rss if proc.is_running() else 0
                cpu = proc.cpu_percent(interval=None)

                rss_children = 0
                cpu_children = 0.0
                if self.include_children:
                    for c in proc.children(recursive=True):
                        try:
                            if not c.is_running():
                                continue
                            mi = c.memory_info()
                            rss_children += getattr(mi, "rss", 0)
                            cpu_children += c.cpu_percent(interval=None)
                        except Exception:
                            continue

                with open(self.output_tsv, "at") as fh:
                    fh.write(
                        "\t".join(
                            [
                                f"{epoch_ts:.3f}",
                                f"{elapsed:.3f}",
                                f"{rss/ (1024*1024):.3f}",
                                f"{cpu:.2f}",
                                f"{rss_children/ (1024*1024):.3f}",
                                f"{cpu_children:.2f}",
                                self._note,
                            ]
                        )
                        + "\n"
                    )
            except Exception:
                # swallow all to keep monitor from crashing the app
                pass

            # wait with ability to exit early
            self._stop_evt.wait(self.interval_sec)


__all__ = ["ResourceMonitor"]
