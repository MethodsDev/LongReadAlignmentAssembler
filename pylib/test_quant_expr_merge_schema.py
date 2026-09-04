#!/usr/bin/env python3

"""The quant.expr merge must refuse a shard whose rows do not match its own header.

`_merge_quant_expr_files` reads shards with `csv.DictReader`, which hides schema damage
in two different places:

  SHORT row (fewer fields than the header) -- every missing column becomes a None
      VALUE. This is the shape a writer produces when it omits a conditional column,
      and it is exactly how the `--oversimplify` writers came to emit `RPM_total_reads`
      into a containment field and leave `uniq_FSM_reads` blank in the merged output.
  LONG row (more fields than the header) -- the surplus lands under the None KEY
      (DictReader's restkey).

The obvious guard, `None in row`, tests KEYS only. It catches the long row and misses
the short one -- missing the very shape it would be written for. Both are asserted here
so that guard cannot be reintroduced.

A malformed shard must FAIL the merge rather than propagate: a silently shifted column
is indistinguishable downstream from a real value, and the merge is the last place the
header is still known.
"""

import csv
import io
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

HEADER = ["gene_id", "transcript_id", "uniq_reads", "all_reads", "RPM_total_reads"]


def _row_dict(header_fields, data_fields):
    """What DictReader hands the merge for a given header/row pair."""
    text = "\t".join(header_fields) + "\n" + "\t".join(data_fields) + "\n"
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))[0]


def _guard_rejects(row, fieldnames):
    """The guard as implemented in _merge_quant_expr_files."""
    return any(value is None for value in row.values()) or None in row


def test_short_row_is_rejected():
    """The OVSIMP shape: a writer omitting conditional columns.

    Fails if the guard is written as `None in row`, which tests keys and therefore
    sees nothing wrong with a row whose missing columns are None VALUES.
    """
    row = _row_dict(HEADER, ["g1", "t1", "5", "5.0"])  # one short
    assert None in row.values(), "fixture is not short; DictReader changed behaviour"
    assert _guard_rejects(row, HEADER), (
        "a row with fewer fields than its header was accepted: the guard is testing "
        "keys instead of values, so a shifted-column shard would merge silently"
    )


def test_long_row_is_rejected():
    """Surplus fields land under the None key (restkey)."""
    row = _row_dict(HEADER, ["g1", "t1", "5", "5.0", "1.0", "surplus"])
    assert None in row, "fixture is not long; DictReader changed behaviour"
    assert _guard_rejects(row, HEADER)


def test_well_formed_row_is_accepted():
    """Converse, so the two above are not passing because everything is rejected."""
    row = _row_dict(HEADER, ["g1", "t1", "5", "5.0", "1.0"])
    assert not _guard_rejects(row, HEADER)


def test_empty_string_fields_are_not_treated_as_damage():
    """Empty is a legitimate value and must not trip the guard.

    `splice_compat_contained`/`splice_contained_by` are emitted as empty strings by the
    oversimplify writers under --debug, and a guard rejecting falsy values rather than
    None would refuse every one of those shards.
    """
    header = HEADER[:4] + ["splice_compat_contained", "splice_contained_by", "RPM_total_reads"]
    row = _row_dict(header, ["g1", "t1", "5", "5.0", "", "", "1.0"])
    assert not _guard_rejects(row, header), (
        "a row with legitimately empty containment fields was rejected: the guard is "
        "testing falsiness rather than None"
    )


def test_the_guard_is_actually_wired_into_the_merge():
    """Guard logic asserted above is worthless if the merge does not call it.

    Checks the source rather than the behaviour, deliberately: driving the real merge
    needs a full run's shard layout, while what can regress here is someone
    simplifying the condition back to `None in row`.
    """
    lraa = REPO_ROOT / "LRAA"
    src = lraa.read_text()
    assert "malformed quant.expr shard" in src, (
        "the merge no longer refuses malformed shards"
    )
    assert "any(value is None for value in row.values())" in src, (
        "the merge's guard no longer inspects row VALUES, so short rows -- the "
        "shifted-column shape -- would be accepted again"
    )
