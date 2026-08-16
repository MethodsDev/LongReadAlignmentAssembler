#!/usr/bin/env python3

"""Tests for pylib/ResourceMonitor.py.

The monitor is what a user reads to size a machine, so a wrong number here is
worse than no number: it is a confident wrong number. Two defects made it one.

CPU attributed nothing to children. psutil measures ``cpu_percent()`` against a
reference point held on the Process OBJECT, the children were re-enumerated on
every sample, and a first call on a fresh object always returns 0.0 -- so the
column read 0.0 for every child of every run. Sizing from it under-provisioned
by the whole worker fleet.

RSS was a spot reading. The write interval was also the sampling rate, and at
60 s against a ~53 s median work unit most units recorded whatever happened to
be true at one arbitrary instant, while a unit shorter than one interval
recorded a duration of zero.

These tests run real load rather than mocking psutil, because both defects were
in what psutil actually does rather than in what the calling code looks like.
"""

import importlib.util
import os
import subprocess
import sys
import time
from importlib.machinery import SourceFileLoader
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ResourceMonitor import ResourceMonitor  # noqa: E402

psutil = pytest.importorskip("psutil")


def _load_summarizer():
    path = (
        Path(__file__).resolve().parents[1] / "util" / "misc" / "summarize_resource_logs.py"
    )
    loader = SourceFileLoader("summarize_resource_logs_under_test", str(path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


summarizer = _load_summarizer()

BURNER = "import time\nt=time.time()\nwhile time.time()-t < {}: pass\n"


def _read_rows(path):
    with open(path, "rt") as fh:
        header = next(fh).rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(header)}
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    return idx, rows


def _burn_in_children(count, seconds):
    kids = [
        subprocess.Popen([sys.executable, "-c", BURNER.format(seconds)])
        for _ in range(count)
    ]
    for kid in kids:
        kid.wait()


def test_cpu_burned_by_children_is_attributed(tmp_path):
    """The defect: this column read 0.0 for every child of every run.

    Two children spinning for two seconds is four CPU-seconds that the process
    tree definitely consumed. Anything at all in the children column fails the
    old implementation, which returned exactly zero.
    """

    out = str(tmp_path / "resources.tsv")
    monitor = ResourceMonitor(out, interval_sec=1.0, include_children=True)
    monitor.start()
    _burn_in_children(2, 2.0)
    monitor.stop()

    idx, rows = _read_rows(out)
    assert rows, "monitor wrote no samples"

    total_cpu_s = float(rows[-1][idx["cpu_seconds_total"]])
    peak_children = max(float(r[idx["cpu_percent_children"]]) for r in rows)

    # 2 x 2 s = 4 CPU-seconds; allow wide slack for a loaded box, but zero is
    # not slack, it is the bug
    assert total_cpu_s > 2.0, "cpu_seconds_total {} is not counting children".format(
        total_cpu_s
    )
    assert peak_children > 100.0, (
        "children CPU peaked at {}%, so two concurrent spinners were credited "
        "to nobody".format(peak_children)
    )


def test_cpu_of_a_child_that_already_exited_still_counts(tmp_path):
    """A worker that finishes inside one interval must not vanish from the total.

    Enumerating live children cannot see it. The reaped subtree is recovered
    from cutime/cstime, which is why the total is built from cpu_times rather
    than from a walk of who happens to be running at the sampling instant.
    """

    out = str(tmp_path / "resources.tsv")
    monitor = ResourceMonitor(out, interval_sec=30.0, include_children=True)
    monitor.start()
    _burn_in_children(1, 1.5)
    time.sleep(0.3)  # child is gone well before the first scheduled write
    monitor.stop()

    idx, rows = _read_rows(out)
    assert rows
    assert float(rows[-1][idx["cpu_seconds_total"]]) > 1.0


def test_a_unit_shorter_than_the_write_interval_still_reports(tmp_path):
    """The 53 second unit that logged duration 0.0.

    With the write interval longer than the run, the old loop wrote its first
    row and then waited, so a short unit produced one row or none and the
    summary differenced a single timestamp against itself.
    """

    out = str(tmp_path / "resources.tsv")
    monitor = ResourceMonitor(out, interval_sec=60.0, include_children=True)
    monitor.start()
    time.sleep(2.0)
    monitor.stop()

    idx, rows = _read_rows(out)
    assert rows, "a run shorter than one interval recorded nothing at all"
    assert float(rows[-1][idx["elapsed_sec"]]) >= 1.5
    assert float(rows[-1][idx["rss_mb"]]) > 0.0

    summary = summarizer.summarize_file(out)
    assert summary["duration_sec"] >= 1.5, (
        "a {:.1f} s run summarised as {:.1f} s".format(2.0, summary["duration_sec"])
    )


def test_rss_peak_survives_between_writes(tmp_path):
    """RSS is a high-water mark over the interval, not a reading at its end.

    A transient allocation that is gone by the time the row is written is
    exactly the peak a machine has to be sized for.
    """

    out = str(tmp_path / "resources.tsv")
    monitor = ResourceMonitor(out, interval_sec=30.0, include_children=False,
                              sample_interval_sec=0.1)
    monitor.start()
    time.sleep(0.4)
    baseline = psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024)
    hog = bytearray(300 * 1024 * 1024)
    hog[::4096] = b"\x01" * len(hog[::4096])  # fault the pages in
    time.sleep(0.6)
    del hog
    time.sleep(0.6)
    monitor.stop()

    idx, rows = _read_rows(out)
    assert rows
    peak = max(float(r[idx["rss_mb"]]) for r in rows)
    assert peak > baseline + 200, (
        "peak {:.0f} MB missed a 300 MB transient over a {:.0f} MB baseline; "
        "the row is reporting an instant, not the interval".format(peak, baseline)
    )


def test_columns_the_summarizer_reads_are_all_present(tmp_path):
    """The summary tool indexes by header name; renaming one silently zeroes it."""

    out = str(tmp_path / "resources.tsv")
    monitor = ResourceMonitor(out, interval_sec=0.5)
    monitor.start()
    time.sleep(1.0)
    monitor.stop()

    idx, rows = _read_rows(out)
    for column in (
        "epoch_ts",
        "elapsed_sec",
        "rss_mb",
        "cpu_percent",
        "rss_mb_children",
        "cpu_percent_children",
        "note",
    ):
        assert column in idx, column
    assert summarizer.summarize_file(out) is not None


def test_note_is_carried_through(tmp_path):
    out = str(tmp_path / "resources.tsv")
    monitor = ResourceMonitor(out, interval_sec=0.5, note="chr21+")
    monitor.start()
    time.sleep(0.8)
    monitor.stop()

    idx, rows = _read_rows(out)
    assert rows
    assert all(r[idx["note"]] == "chr21+" for r in rows)


def test_sampling_is_never_slower_than_writing(tmp_path):
    """A high-water mark taken once per row is a spot reading with extra steps."""

    monitor = ResourceMonitor(str(tmp_path / "r.tsv"), interval_sec=2.0)
    assert monitor.sample_interval_sec <= monitor.interval_sec

    tight = ResourceMonitor(str(tmp_path / "r2.tsv"), interval_sec=0.5)
    assert tight.sample_interval_sec <= tight.interval_sec


def test_stop_without_start_is_harmless(tmp_path):
    ResourceMonitor(str(tmp_path / "r.tsv")).stop()
