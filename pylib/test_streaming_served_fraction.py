#!/usr/bin/env python3

"""The streaming pass reports how much its table served, and refuses nothing for it.

There used to be a threshold here: a unit whose reads landed on unseen paths more than 25%
of the time raised, and it raised from `_report`, which runs after the streaming loop has
mapped, looked up, written and dropped every read. So the refusal could only discard a
complete and correct tracking file for having been slow. These tests pin the two halves of
its removal -- that such a unit now COMPLETES, and that the measurement it gated on is
still computed and logged.

The served fraction is worth reporting because it varies enormously between units of one
healthy run: MEASURED over the 42 PBMC units that streamed, 5.3% to 93.9%, median 79.1%,
with the low unit's nearest neighbour at 43.7%. That 8x gap is a pointer at one contig's
splice graph, which is what the abort never gave anybody.
"""

import logging
import os
import sys

import pysam
import pytest

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import StreamingQuant

CONTIG = "chr1"
CONTIG_LEN = 5000

# Two reads, both mapping to a path the table does not carry, so every accounted read is
# resolved in-stream: unseen fraction 1.0, four times the retired 0.25 limit.
START_A = 100
START_B = 300


class _FakeSpliceGraph:
    def canonical_simple_path(self, simple_path):
        return "canon:" + ",".join(simple_path)

    def get_contig_acc(self):
        return CONTIG

    def get_contig_strand(self):
        return "+"


class _FakeLRAA:
    def __init__(self):
        self._splice_graph = _FakeSpliceGraph()

    def _map_read_to_graph(self, segments, **kwargs):
        return ["E:1"]


class _Totals:
    """Only the fields served_read_fraction reads, so the arithmetic is pinned alone."""

    def __init__(self, streamed, resolved):
        self.reads_streamed = streamed
        self.reads_on_stream_resolved_paths = resolved


def _write_bam(path):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"},
         "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}]}
    )
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for name, start in (("readA", START_A), ("readB", START_B)):
            a = pysam.AlignedSegment(header)
            a.query_name = name
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigarstring = "100M"
            a.query_sequence = "ACGT" * 25
            a.query_qualities = pysam.qualitystring_to_array("I" * 100)
            a.set_tag("NM", 1, value_type="i")
            fh.write(a)
    pysam.index(str(path))
    return path


@pytest.fixture
def bam(tmp_path):
    return _write_bam(tmp_path / "reads.bam")


def _run(bam, tmp_path, table):
    tracking = tmp_path / "tracking.tsv"
    with open(tracking, "wt") as fh:
        totals = StreamingQuant.stream_assign(
            str(bam), CONTIG, "+", "A" * CONTIG_LEN, _FakeLRAA(), table, fh,
            # The resolver answers, as it does in a real run: an unseen path is resolved
            # once and cached, never dropped.
            resolver=lambda path: [("g1", "t1", "h1", 2, "mp1", 1.0, 1.0, 0)],
            try_correct_alignments=False,
        )
    return totals, tracking.read_text().splitlines()


def test_an_entirely_unserved_unit_completes(bam, tmp_path):
    """100% unseen -- what used to raise -- returns totals and a complete tracking file.

    An empty table is the extreme of the condition the old guard fired on. It is also the
    case where the guard was most clearly wrong: both reads are assigned, both have their
    tracking rows written, and the refusal came after both were on disk.
    """
    totals, rows = _run(bam, tmp_path, StreamingQuant.AssignmentTable())

    assert totals.reads_streamed == 2
    assert totals.reads_assigned == 2
    # One resolve served both reads, and both are charged to it.
    assert totals.paths_resolved_in_stream == 1
    assert totals.reads_on_stream_resolved_paths == 2
    assert StreamingQuant.served_read_fraction(totals) == 0.0
    assert [r.split("\t")[5] for r in rows] == ["readA", "readB"]


def test_the_report_logs_the_served_fraction_and_what_it_means(bam, tmp_path, caplog):
    """The measurement the gate gated on survives the gate, as an actionable line.

    Asserted on the rendered message rather than on the format string, because a %-arg in
    the wrong position renders a wrong number with a healthy-looking template.
    """
    with caplog.at_level(logging.INFO, logger="StreamingQuant"):
        _run(bam, tmp_path, StreamingQuant.AssignmentTable())

    served_lines = [
        r.getMessage() for r in caplog.records if "first-pass table served" in r.getMessage()
    ]
    assert len(served_lines) == 1
    message = served_lines[0]
    assert message.startswith(f"[{CONTIG}+] first-pass table served 0.0% of this unit's 2")
    # Reads per resolve and table size travel with it: two reads on one resolve, and the
    # table the stream ended holding is the empty one plus that resolution.
    assert "2.0 reads per resolve" in message
    assert "against a table of 1 paths" in message
    # A number with no reading attached is what made the old guard's threshold look like
    # the actionable part. The line has to name what to inspect.
    assert "splice graph" in message
    assert "--normalize_max_cov_level" in message
    assert "costs time, not correctness" in message


def test_a_fully_served_unit_reports_the_other_end(bam, tmp_path, caplog):
    """The control: a table carrying the path resolves nothing and serves everything."""
    table = StreamingQuant.AssignmentTable()
    table.add("canon:E:1", [("g1", "t1", "h1", 2, "mp1", 1.0, 1.0, 0)])

    with caplog.at_level(logging.INFO, logger="StreamingQuant"):
        totals, _rows = _run(bam, tmp_path, table)

    assert totals.paths_resolved_in_stream == 0
    assert StreamingQuant.served_read_fraction(totals) == 1.0
    assert any(
        "first-pass table served 100.0% of this unit's 2" in r.getMessage()
        for r in caplog.records
    )


@pytest.mark.parametrize(
    "streamed,resolved,expected",
    [
        # The two ends of the range MEASURED across the 42 PBMC units, 5.3% and 93.9%.
        (1000, 947, 0.053),
        (1000, 61, 0.939),
        (100, 0, 1.0),
        (100, 100, 0.0),
        # Real-run magnitudes, from that run's aggregate counters: 7,591,317 reads landed
        # on unseen paths out of 36,510,028. Included so the identity is checked at a
        # scale where float division actually has to hold up.
        (36510028, 7591317, 0.79208),
    ],
)
def test_served_fraction_is_streamed_minus_resolved_over_streamed(
    streamed, resolved, expected
):
    """The identity, across the whole range a real run spans."""
    frac = StreamingQuant.served_read_fraction(_Totals(streamed, resolved))
    assert frac == pytest.approx((streamed - resolved) / streamed)
    assert frac == pytest.approx(expected, abs=1e-3)


def test_zero_reads_is_undefined_rather_than_zero():
    """A unit that streamed nothing served nothing of nothing.

    Reporting 0.0 there would read as "the table answered none of this unit's work", which
    is the alarming end of the scale, on a unit that asked it nothing. None is the only
    honest answer and `_report` skips the line rather than printing a fabricated rate.
    """
    assert StreamingQuant.served_read_fraction(_Totals(0, 0)) is None


def test_the_retired_guard_leaves_no_trace():
    """No accepted-and-ignored remnant of either name anywhere in the tree.

    Commit b15430d is the cautionary case: a flag declared as a no-op survived releases
    while one path still acted on it. The needles are assembled from fragments so that this
    file is not itself an exception to the assertion it makes.
    """
    needles = ["max_fallback" + "_path_frac", "stream_reads_max_unseen" + "_path_read_frac"]
    repo = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    # CHANGELOG.txt is excluded on purpose: it is the record that these were removed, and a
    # release history that cannot name what it removed is useless.
    skip_dirs = {".git", ".venv", "__pycache__", "node_modules"}
    hits = []
    for dirpath, dirnames, filenames in os.walk(repo):
        dirnames[:] = [d for d in dirnames if d not in skip_dirs]
        for filename in filenames:
            if filename == "CHANGELOG.txt":
                continue
            full = os.path.join(dirpath, filename)
            try:
                with open(full, "rt", encoding="utf-8", errors="strict") as fh:
                    text = fh.read()
            except (OSError, UnicodeDecodeError):
                continue  # binary fixtures and unreadable paths carry no source reference
            for needle in needles:
                if needle in text:
                    hits.append(f"{os.path.relpath(full, repo)}: {needle}")
    assert hits == []
