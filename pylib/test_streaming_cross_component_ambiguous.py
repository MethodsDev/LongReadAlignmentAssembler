#!/usr/bin/env python3

"""A resolved-in-stream path bridging two independently-normalized components is
collateral, not a fatal error -- the same shape of loss chunk boundaries already accept
for a severed alignment: dropped, counted, and never silent.

`rows_for_multipath` (StreamingQuant.py) still refuses to guess at a split for such a
read and raises `CrossComponentAmbiguousPath`. What changed is who is allowed to see that
exception: `stream_assign`'s `assign_path` now catches it, so the run completes and the
read is folded into `totals.reads_cross_component_ambiguous` (and, like any other read
with zero tracking rows, into `totals.reads_unassignable`) instead of aborting the unit.

Fixture mirrors `test_streaming_served_fraction.py`: two real bam records mapped by a
fake LRAA object onto the SAME graph path, so both reads share one canonical path and
the resolver is asked for it only once -- exercising both the first-resolution catch and
the cached-verdict branch a second read on the same path takes.
"""

import os
import sys

import pysam
import pytest

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import LRAA_Globals
import StreamingQuant

CONTIG = "chr1"
CONTIG_LEN = 5000
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
        # Same path for every read, so readA and readB canonicalize identically and the
        # second read exercises the table's cached verdict rather than a second resolve.
        return ["E:1"]


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


def _run(bam, tmp_path, resolver, rescuer=None):
    tracking = tmp_path / "tracking.tsv"
    with open(tracking, "wt") as fh:
        totals = StreamingQuant.stream_assign(
            str(bam), CONTIG, "+", "A" * CONTIG_LEN, _FakeLRAA(),
            StreamingQuant.AssignmentTable(), fh,
            resolver=resolver,
            try_correct_alignments=False,
            rescuer=rescuer,
        )
    return totals, tracking.read_text().splitlines()


def _raising_resolver(path):
    raise StreamingQuant.CrossComponentAmbiguousPath(
        "path {} is compatible with isoforms from 2 separately quantified components "
        "(1,2)".format(path)
    )


def test_the_unit_completes_instead_of_raising(bam, tmp_path):
    """What used to abort the whole unit now just drops the two affected reads."""
    totals, rows = _run(bam, tmp_path, _raising_resolver)

    assert totals.reads_streamed == 2
    assert rows == []


def test_every_read_on_the_ambiguous_path_is_counted_not_just_the_first(bam, tmp_path):
    """The resolver is consulted once (`table` caches the negative verdict), but both
    reads landing on that path must still be individually counted: `stream_resolved`'s
    own `reads_on_stream_resolved_paths` counter uses the identical two-site pattern
    (checked at cache-hit AND incremented at first resolution) for exactly this reason.
    """
    totals, _rows = _run(bam, tmp_path, _raising_resolver)

    assert totals.reads_cross_component_ambiguous == 2
    # The resolver's compatibility cascade itself only ran once for the shared path.
    assert totals.paths_resolved_in_stream == 1


def test_dropped_reads_are_also_unassignable_not_a_separate_bucket(bam, tmp_path):
    """reads_cross_component_ambiguous is a strict subset/breakdown of reads_unassignable,
    the same relationship reads_zero_fraction has to reads_assigned -- not a parallel
    accounting track a naive reader of reads_assigned + reads_unassignable could miss.
    """
    totals, _rows = _run(bam, tmp_path, _raising_resolver)

    assert totals.reads_unassignable == 2
    assert totals.reads_assigned == 0


def test_a_real_isoform_match_on_a_different_path_is_unaffected(bam, tmp_path):
    """The collateral is scoped to the one ambiguous path; a table hit for another read
    on another path must still write its row normally."""
    tracking = tmp_path / "tracking.tsv"
    with open(tracking, "wt") as fh:
        table = StreamingQuant.AssignmentTable()
        table.add("canon:E:1", [("g1", "t1", "h1", 2, "mp1", 1.0, 1.0)])
        totals = StreamingQuant.stream_assign(
            str(_write_bam(tmp_path / "reads2.bam")), CONTIG, "+", "A" * CONTIG_LEN,
            _FakeLRAA(), table, fh,
            resolver=_raising_resolver,
            try_correct_alignments=False,
        )
    rows = tracking.read_text().splitlines()

    # Both reads land on the SAME path here too, but the table already carries a real
    # answer for it, so the resolver (and therefore the exception) is never reached.
    assert totals.reads_cross_component_ambiguous == 0
    assert totals.reads_assigned == 2
    assert [r.split("\t")[5] for r in rows] == ["readA", "readB"]


def test_the_warning_names_the_count_and_the_lever(bam, tmp_path, caplog):
    import logging

    with caplog.at_level(logging.WARNING, logger="StreamingQuant"):
        _run(bam, tmp_path, _raising_resolver)

    lines = [
        r.getMessage() for r in caplog.records
        if "independently normalized component" in r.getMessage()
    ]
    assert len(lines) == 1
    assert lines[0].startswith(f"[{CONTIG}+] 2 reads")
    assert "--normalize_max_cov_level" in lines[0]


class _AcceptsEverything:
    """A rescuer that offers a real, resolvable, single-component path for anything it
    is asked about -- the adversarial case: if the ambiguous-path exclusion in
    `stream_assign` were missing, this is what would let the dropped read back in.
    """

    def __init__(self, path=("E:2",)):
        self._path = tuple(path)
        self.offered = []

    def offer(self, read, category):
        self.offered.append((read.query_name, category))
        return (self._path,)


def test_a_cross_component_ambiguous_read_is_never_offered_to_rescue(bam, tmp_path, monkeypatch):
    """Rescue must not be able to hand the SAME read a second, different, non-ambiguous
    path once its own path has been ruled out as collateral -- that would silently
    reassign a read this pass deliberately refused to attribute anywhere, via a route
    with no knowledge of why the first assignment was refused.

    `stream_reads_rescue_unassigned_to_targets` is the one flag that would route a
    dropped-as-unassignable read to rescue at all, so it has to be on for this to be a
    real adversarial test rather than one where rescue is never consulted anyway.
    """
    monkeypatch.setitem(
        LRAA_Globals.config, "stream_reads_rescue_unassigned_to_targets", True
    )
    rescuer = _AcceptsEverything()

    totals, rows = _run(bam, tmp_path, _raising_resolver, rescuer=rescuer)

    assert rescuer.offered == []
    assert totals.reads_cross_component_ambiguous == 2
    assert totals.rescue_rows_written == 0
    assert rows == []


