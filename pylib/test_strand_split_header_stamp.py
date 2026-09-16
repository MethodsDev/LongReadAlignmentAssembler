#!/usr/bin/env python3

"""A strand-split bam says which orientation it holds, and LRAA believes it.

The splitter writes one bam per orientation. Until the file SAYS so, the only way to
learn which one it is was to find reads in it -- and that fails on exactly the case
that matters: a slice with nothing for its orientation. A reader then cannot tell
"this file is not about that strand" from "that strand is empty here".

That ambiguity is not academic. LRAA enumerates a job per (contig, strand) and the
chunked pipeline invokes it once per strand-split unit, so the opposite-orientation
job runs against a bam that cannot contain it. Anything emitting per orientation --
the de novo oversimplify aggregate, which exists whether or not a read reached it --
had to resolve that from read presence, and emitting unconditionally put a '-'
aggregate in the '+' unit's gtf.

Asserted at both ends: the splitter stamps, and the resolution honours the stamp.
"""

import os
import sys

import pysam
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if os.path.join(REPO_ROOT, "pylib") not in sys.path:
    sys.path.insert(0, os.path.join(REPO_ROOT, "pylib"))

import Util_funcs


def _bam(path, header, reads=()):
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for read in reads:
            fh.write(read)
    return str(path)


BASE_HEADER = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chrM", "LN": 16569}],
}


def test_a_stamped_bam_declares_its_orientation(tmp_path):
    for strand in ("+", "-"):
        header = Util_funcs.stamp_strand_split_header(BASE_HEADER, strand)
        path = _bam(tmp_path / "strand{}.bam".format("plus" if strand == "+" else "minus"), header)
        assert Util_funcs.declared_strand_of_bam(path) == strand


def test_an_empty_stamped_bam_still_declares_it(tmp_path):
    """The whole point: no reads, and the orientation is still knowable."""
    header = Util_funcs.stamp_strand_split_header(BASE_HEADER, "-")
    path = _bam(tmp_path / "empty.bam", header)

    with pysam.AlignmentFile(path, "rb") as reader:
        assert sum(1 for _ in reader.fetch(until_eof=True)) == 0
    assert Util_funcs.declared_strand_of_bam(path) == "-"


def test_an_unstamped_bam_claims_nothing(tmp_path):
    """Every bam that did not come from the splitter, which must stay unrestricted."""
    path = _bam(tmp_path / "plain.bam", BASE_HEADER)
    assert Util_funcs.declared_strand_of_bam(path) is None


def test_restamping_replaces_rather_than_accumulates():
    """Splitting an already-split bam must not leave it claiming both orientations."""
    once = Util_funcs.stamp_strand_split_header(BASE_HEADER, "+")
    twice = Util_funcs.stamp_strand_split_header(once, "-")

    stamps = [
        c
        for c in twice["CO"]
        if str(c).startswith(Util_funcs.STRAND_SPLIT_HEADER_PREFIX)
    ]
    assert stamps == [Util_funcs.strand_split_header_comment("-")]


def test_a_bam_claiming_both_orientations_is_treated_as_claiming_neither(tmp_path):
    """Refusing to pick. Acting on either would silently drop half the reads.

    Reachable only by hand-editing a header, which is why it answers None rather
    than raising: the caller's safe default is "no restriction".
    """
    header = dict(BASE_HEADER)
    header["CO"] = [
        Util_funcs.strand_split_header_comment("+"),
        Util_funcs.strand_split_header_comment("-"),
    ]
    path = _bam(tmp_path / "both.bam", header)
    assert Util_funcs.declared_strand_of_bam(path) is None


def test_other_comments_are_preserved(tmp_path):
    """The stamp is additive: a header carrying provenance keeps it."""
    header = dict(BASE_HEADER)
    header["CO"] = ["some upstream note"]
    stamped = Util_funcs.stamp_strand_split_header(header, "+")
    path = _bam(tmp_path / "noted.bam", stamped)

    with pysam.AlignmentFile(path, "rb") as reader:
        comments = reader.header.to_dict()["CO"]
    assert "some upstream note" in comments
    assert Util_funcs.declared_strand_of_bam(path) == "+"


def test_a_missing_file_claims_nothing():
    assert Util_funcs.declared_strand_of_bam("/nonexistent/none.bam") is None
    assert Util_funcs.declared_strand_of_bam(None) is None


def test_the_splitter_stamps_both_outputs(tmp_path):
    """End to end through the real splitter, not its header helper.

    Guards the wiring: the helper being right is worth nothing if the writers are
    still constructed from the unmodified template.
    """
    import importlib.machinery
    import importlib.util

    path = os.path.join(REPO_ROOT, "util", "separate_bam_by_strand.py")
    loader = importlib.machinery.SourceFileLoader("separate_bam_under_test", path)
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)

    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrM", "LN": 1000}]}
    source = tmp_path / "in.bam"
    with pysam.AlignmentFile(str(source), "wb", header=header) as fh:
        for i, reverse in enumerate((False, True)):
            read = pysam.AlignedSegment()
            read.query_name = "r{}".format(i)
            read.query_sequence = "A" * 50
            read.flag = 16 if reverse else 0
            read.reference_id = 0
            read.reference_start = 100 + i
            read.mapping_quality = 60
            read.cigar = [(0, 50)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 50)
            fh.write(read)
    pysam.index(str(source))

    plus = str(tmp_path / "out.plus.bam")
    minus = str(tmp_path / "out.minus.bam")
    module.split_bam_by_strand(str(source), plus, minus, 100000)

    assert Util_funcs.declared_strand_of_bam(plus) == "+"
    assert Util_funcs.declared_strand_of_bam(minus) == "-"
