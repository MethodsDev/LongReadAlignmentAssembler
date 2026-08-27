#!/usr/bin/env python3

"""What ``ChunkedRun.merge_read_assignment_summaries`` puts in ``chunk_id``, and which
rows it adds up.

Two contracts, both about telling a SENTINEL apart from a VALUE.

``chunk_id`` is empty on exactly two kinds of row: the merged ``TOTAL``, which spans
every chunk so no single chunk id is true of it, and a worker row whose unit record
carried no chunk id at all (a manifest-driven stage 6 over a manifest predating the
field). Empty therefore means "no chunk applies". Rendering the column with
``str(unit.get("chunk_id") or "")`` broke that: any FALSY id -- integer 0 most
obviously -- came out empty and became indistinguishable from those two cases. Ids
are strings like ``chr19_00`` today, so nothing reached it; a provenance column that
cannot tell a real id from "not applicable" is a trap regardless.

The grand ``TOTAL`` is a sum of the ``worker`` rows. Each per-unit file is itself the
output of a complete LRAA run and already carries its own ``TOTAL``, so that row is
skipped; but skipping only ``TOTAL`` would have admitted any future third row type
into the run's read total, and dropping every non-``worker`` row would have made one
vanish from it. Both are silent, and this table IS the run's read accounting, so an
unrecognized ``row_type`` is refused instead.

Inputs are built by LRAA's own writer and its own single-unit merge, never hand-typed,
so the schema under test is the one a real unit emits.
"""

import csv
import importlib.util
import os
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "pylib"))

import ChunkedRun  # noqa: E402


def _load_lraa():
    loader = SourceFileLoader("lraa_driver_chunk_id_fixture", str(REPO / "LRAA"))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


lraa = _load_lraa()

MERGED_NAME = "chunked.read_assignment.summary.tsv"


def _unit(root, stem, reads, chunk_id="__absent__"):
    """One quant unit's directory contents plus the unit record stage 6 consumes.

    ``chunk_id="__absent__"`` omits the key entirely, which is the manifest-driven
    case; anything else is stored verbatim so a falsy value stays falsy.
    """

    prefix = str(root / stem)
    worker = prefix + ".read_assignment.summary.worker.tsv"
    lraa._write_read_assignment_summary(
        worker, "chrT", "+", {"reads_total": reads, "reads_kept_genome": reads}
    )
    lraa._merge_read_assignment_summary_files(
        [worker], prefix + ".read_assignment.summary.tsv"
    )
    os.remove(worker)

    record = {"unit_id": stem + "_plus", "offset": 0, "quant_prefix": prefix}
    if chunk_id != "__absent__":
        record["chunk_id"] = chunk_id
    return record


def _merge(tmp_path, units):
    merged_dir = tmp_path / "merged"
    merged_dir.mkdir(exist_ok=True)
    result = ChunkedRun.merge_read_assignment_summaries(str(merged_dir), units)
    assert result is not None
    with open(result["path"], "rt", newline="") as fh:
        rows = list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))
    return rows


def _workers(rows):
    return [r for r in rows if r["row_type"] == "worker"]


def _total(rows):
    picked = [r for r in rows if r["row_type"] == "TOTAL"]
    assert len(picked) == 1, rows
    return picked[0]


# ------------------------------------------------------------- the chunk_id column


def test_a_falsy_but_real_chunk_id_renders_as_itself(tmp_path):
    """Integer 0 is a chunk id, not a missing one.

    This is the whole defect: ``or ""`` cannot distinguish 0 from absent, so the row
    that names cut 0 came out looking like the row that names no cut. Both other
    falsy spellings are pinned beside it, because a scheme that numbers cuts is
    exactly the scheme that would produce them.
    """

    units = [
        _unit(tmp_path, "u_zero_int", 10, chunk_id=0),
        _unit(tmp_path, "u_zero_str", 20, chunk_id="0"),
    ]

    rows = _workers(_merge(tmp_path, units))
    assert [r["chunk_id"] for r in rows] == ["0", "0"]


def test_an_absent_chunk_id_leaves_the_column_empty(tmp_path):
    """Deliberate and load-bearing: a unit record without one must not be guessed at."""

    rows = _workers(_merge(tmp_path, [_unit(tmp_path, "u_absent", 7)]))
    assert [r["chunk_id"] for r in rows] == [""]


def test_an_explicit_none_chunk_id_leaves_the_column_empty(tmp_path):
    """None present is the same statement as the key being absent."""

    rows = _workers(_merge(tmp_path, [_unit(tmp_path, "u_none", 7, chunk_id=None)]))
    assert [r["chunk_id"] for r in rows] == [""]


def test_a_string_chunk_id_is_carried_through(tmp_path):
    """Today's shape, so the fix above cannot have been made by breaking it."""

    rows = _workers(_merge(tmp_path, [_unit(tmp_path, "u_real", 7, chunk_id="chr19_00")]))
    assert [r["chunk_id"] for r in rows] == ["chr19_00"]


def test_the_total_row_carries_no_chunk_id(tmp_path):
    """It spans every chunk, so no single id is true of it -- and that is why a real
    falsy id rendering empty was a collision and not a cosmetic problem."""

    rows = _merge(tmp_path, [_unit(tmp_path, "u_t", 5, chunk_id=0)])
    assert _total(rows)["chunk_id"] == ""


# ------------------------------------------------------- which rows enter the TOTAL


def test_the_total_sums_the_worker_rows_only(tmp_path):
    """Each unit file carries its own TOTAL; counting it would double every unit."""

    units = [
        _unit(tmp_path, "u_a", 11, chunk_id="chrT_00"),
        _unit(tmp_path, "u_b", 31, chunk_id="chrT_01"),
    ]
    rows = _merge(tmp_path, units)
    assert sum(int(r["reads_total"]) for r in _workers(rows)) == 42
    assert int(_total(rows)["reads_total"]) == 42


def test_an_unrecognized_row_type_is_refused(tmp_path):
    """Neither summing it nor dropping it is a defensible guess.

    LRAA writes exactly ``worker`` and ``TOTAL``, so this is unreachable from any
    in-tree writer. It is refused rather than skipped because both alternatives move
    the run's reported read total with nothing in the table contradicting the new
    number: summed, the reads are counted twice; dropped, they are gone.
    """

    unit = _unit(tmp_path, "u_odd", 9, chunk_id="chrT_00")
    path = unit["quant_prefix"] + ".read_assignment.summary.tsv"
    with open(path, "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames)
        rows = list(reader)
    intruder = dict(rows[0])
    intruder["row_type"] = "cluster"
    with open(path, "wt", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
        writer.writerow(intruder)

    merged_dir = tmp_path / "merged_refused"
    merged_dir.mkdir()
    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.merge_read_assignment_summaries(str(merged_dir), [unit])

    message = str(err.value)
    assert "'cluster'" in message
    assert "worker" in message
    assert not os.path.exists(str(merged_dir / MERGED_NAME))
