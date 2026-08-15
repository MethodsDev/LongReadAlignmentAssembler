#!/usr/bin/env python3

"""Tests for util/misc/compare_read_compatibility.py.

The tool exists to separate two things a quant metric cannot separate: a config that
merely REWEIGHTS reads across a fixed candidate set, and a config that changes WHICH
transcripts a read is compatible with. The second reaches through the splice graph and
has a much larger blast radius, so conflating them would let a compatibility change be
reported as a harmless reweighting.

Its known-answer case comes from real output: 3'-agreement weighting on vs off on
CL_HG002_E2_sirv gives 0 of 250,129 changed compatible sets alongside 277 reads whose
fractional assignment moved. Turning the weighting off cannot alter compatibility -- the
flag is read only where weights are applied -- so any nonzero set count there would have
meant the instrument was miscounting. These tests pin that signature on fixtures small
enough to check by eye, plus the opposite case the instrument must not miss.
"""

import gzip
import importlib.util
from importlib.machinery import SourceFileLoader
from pathlib import Path

import pytest


def _load_tool():
    path = (
        Path(__file__).resolve().parents[1]
        / "util"
        / "misc"
        / "compare_read_compatibility.py"
    )
    loader = SourceFileLoader("compare_read_compatibility_under_test", str(path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


tool = _load_tool()

HEADER = ("gene_id\ttranscript_id\ttranscript_splice_hash_code\tnum_exons\t"
          "mp_id\tread_name\tfrac_assigned\tread_weight\n")


def _write_tracking(path, rows, gzipped=False):
    """rows: (transcript_id, mp_id, read_name, frac_assigned, read_weight)."""
    body = HEADER + "".join(
        f"G1\t{tid}\thash\t2\t{mp}\t{read}\t{frac:.6f}\t{w:.6f}\n"
        for tid, mp, read, frac, w in rows
    )
    opener = gzip.open if gzipped else open
    with opener(path, "wt") as fh:
        fh.write("# LRAA version test\n")
        fh.write(body)
    return str(path)


def test_pure_reweighting_changes_no_compatible_sets(tmp_path):
    """The real known-answer case, in miniature.

    r1 is ambiguous between t1 and t2 and its split moves; r2 is unambiguous. No
    candidate set changes, so changed_compatible_set must be 0 while changed_total_frac
    picks up nothing either -- a read's fractions still sum to 1 however they are split.
    """
    a = _write_tracking(tmp_path / "a.tracking", [
        ("t1", "MP1", "r1", 0.50, 0.50),
        ("t2", "MP1", "r1", 0.50, 0.50),
        ("t1", "MP2", "r2", 1.00, 1.00),
    ])
    b = _write_tracking(tmp_path / "b.tracking", [
        ("t1", "MP1", "r1", 0.80, 1.00),
        ("t2", "MP1", "r1", 0.20, 1.00),
        ("t1", "MP2", "r2", 1.00, 1.00),
    ])
    r = tool.compare(a, b)
    assert r["changed_compatible_set"] == 0
    assert r["changed_set_size"] == 0
    assert r["shared_reads"] == 2
    assert r["reads_only_in_A"] == 0 and r["reads_only_in_B"] == 0


def test_detects_a_changed_candidate_set(tmp_path):
    """A read gaining a candidate is a compatibility change, not a reweighting.

    This is the case the tool exists to catch, and the one a MARD delta would hide.
    """
    a = _write_tracking(tmp_path / "a.tracking", [
        ("t1", "MP1", "r1", 0.50, 0.50),
        ("t2", "MP1", "r1", 0.50, 0.50),
    ])
    b = _write_tracking(tmp_path / "b.tracking", [
        ("t1", "MP1", "r1", 0.34, 0.33),
        ("t2", "MP1", "r1", 0.33, 0.33),
        ("t3", "MP1", "r1", 0.33, 0.33),
    ])
    r = tool.compare(a, b)
    assert r["changed_compatible_set"] == 1
    assert r["changed_set_size"] == 1


def test_same_set_size_but_different_members_still_counts(tmp_path):
    """A swap keeps N fixed, so set-size alone would miss it. The set comparison must not.

    Substituting one candidate for another is exactly what a shifted polyA window can do
    without changing how ambiguous the read looks.
    """
    a = _write_tracking(tmp_path / "a.tracking", [
        ("t1", "MP1", "r1", 0.50, 0.50),
        ("t2", "MP1", "r1", 0.50, 0.50),
    ])
    b = _write_tracking(tmp_path / "b.tracking", [
        ("t1", "MP1", "r1", 0.50, 0.50),
        ("t3", "MP1", "r1", 0.50, 0.50),
    ])
    r = tool.compare(a, b)
    assert r["changed_compatible_set"] == 1
    assert r["changed_set_size"] == 0


def test_reports_reads_present_in_only_one_arm(tmp_path):
    """Reads appearing or vanishing are counted apart from set changes.

    They are not comparable per-read, so folding them into changed_compatible_set would
    inflate a compatibility signal with what is really a retention difference.
    """
    a = _write_tracking(tmp_path / "a.tracking", [
        ("t1", "MP1", "r1", 1.00, 1.00),
        ("t1", "MP2", "r2", 1.00, 1.00),
    ])
    b = _write_tracking(tmp_path / "b.tracking", [
        ("t1", "MP1", "r1", 1.00, 1.00),
        ("t1", "MP3", "r3", 1.00, 1.00),
    ])
    r = tool.compare(a, b)
    assert r["reads_only_in_A"] == 1
    assert r["reads_only_in_B"] == 1
    assert r["shared_reads"] == 1
    assert r["changed_compatible_set"] == 0


def test_changed_total_frac_catches_a_partially_assigned_read(tmp_path):
    """Total assigned fraction moving is the reweighting channel, reported separately."""
    a = _write_tracking(tmp_path / "a.tracking", [("t1", "MP1", "r1", 1.00, 1.00)])
    b = _write_tracking(tmp_path / "b.tracking", [("t1", "MP1", "r1", 0.40, 1.00)])
    r = tool.compare(a, b)
    assert r["changed_total_frac"] == 1
    assert r["changed_compatible_set"] == 0


def test_reads_gzipped_tracking(tmp_path):
    """LRAA writes quant.tracking.gz by default, so the gz path is the normal path."""
    a = _write_tracking(tmp_path / "a.tracking.gz", [
        ("t1", "MP1", "r1", 0.50, 0.50),
        ("t2", "MP1", "r1", 0.50, 0.50),
    ], gzipped=True)
    b = _write_tracking(tmp_path / "b.tracking.gz", [
        ("t1", "MP1", "r1", 0.50, 1.00),
        ("t2", "MP1", "r1", 0.50, 1.00),
    ], gzipped=True)
    r = tool.compare(a, b)
    assert r["shared_reads"] == 1
    assert r["changed_compatible_set"] == 0


def test_empty_input_is_an_error_not_a_silent_zero(tmp_path):
    """An unreadable arm must not look like 'nothing changed'."""
    a = _write_tracking(tmp_path / "a.tracking", [("t1", "MP1", "r1", 1.0, 1.0)])
    empty = tmp_path / "empty.tracking"
    empty.write_text("# LRAA version test\n" + HEADER)
    with pytest.raises(ValueError):
        tool.compare(a, str(empty))
