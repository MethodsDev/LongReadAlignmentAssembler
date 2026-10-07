#!/usr/bin/env python3

"""Transient bams are written at a low BGZF level; kept bams never are.

``normalize_bam_by_strand.py`` writes strand bams and per-contig parts that are read once
and deleted, so it asks for a fast level (measured: 545 s -> 377 s on a 6.6 GB merged bam
together with the threaded PG scan, records and header identical). The risk this guards
is the level leaking into something that is KEPT: a whole-file unit's part IS the final
output, and a single-strand input concatenates its per-contig parts into the final file.
"""

import os
import sys
from pathlib import Path

import pysam
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "util"))

import separate_bam_by_strand as sep  # noqa: E402


def _write_bam(path, level):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 100000}]}
    options = {} if level is None else {"format_options": ["level={}".format(level).encode()]}
    with pysam.AlignmentFile(str(path), "wb", header=header, **options) as out:
        for i in range(4000):
            a = pysam.AlignedSegment()
            a.query_name = "read{:05d}".format(i)
            a.query_sequence = "ACGT" * 25
            a.flag = 0
            a.reference_id = 0
            a.reference_start = i * 20
            a.mapping_quality = 60
            a.cigarstring = "100M"
            a.query_qualities = pysam.qualitystring_to_array("I" * 100)
            out.write(a)
    return path


def test_the_default_leaves_the_level_to_htslib(monkeypatch):
    monkeypatch.setattr(sep, "_INTERMEDIATE_COMPRESSION_LEVEL", None)
    assert sep._writer_options() == {}


def test_a_level_reaches_the_writer(monkeypatch, tmp_path):
    monkeypatch.setattr(sep, "_INTERMEDIATE_COMPRESSION_LEVEL", 1)
    assert sep._writer_options() == {"format_options": [b"level=1"]}

    fast = _write_bam(tmp_path / "fast.bam", 1)
    best = _write_bam(tmp_path / "best.bam", 9)
    # the option is applied, not accepted and ignored: a higher level is smaller
    assert os.path.getsize(best) < os.path.getsize(fast)
    # and either is an ordinary readable bam with the same records
    with pysam.AlignmentFile(str(fast)) as a, pysam.AlignmentFile(str(best)) as b:
        assert [r.query_name for r in a.fetch(until_eof=True)] == [
            r.query_name for r in b.fetch(until_eof=True)
        ]


def test_the_level_never_reaches_a_kept_output():
    """Pinned on the source: both exclusions are what keep the final bam at the default."""

    src = (REPO / "util" / "normalize_bam_by_strand.py").read_text()
    # single-strand input: parts are concatenated into the FINAL output
    assert "None if args.input_is_single_strand else args.intermediate_compression_level" in src
    # a whole-file unit (scope None) writes the final output directly
    assert 'compression_level=level if unit["scope"] is not None else None' in src
    # and the final merge is not given a level at all
    merge_cmd = src[src.index('f"samtools merge --no-PG'):]
    assert "--output-fmt-option" not in merge_cmd[:200]


def test_the_collapse_scan_is_threaded():
    """The PG:Z: tag scan is a full decompression of the bam; one thread was 125 s of 6.6 GB."""

    src = (REPO / "util" / "normalize_bam_by_strand.py").read_text()
    assert "--threads {max(samtools_threads, 1)}" in src


def test_the_normalize_wdl_passes_the_level_and_it_is_an_intermediate_only_knob():
    wdl = (REPO / "WDL" / "subwdls" / "Normalize_bam.wdl").read_text()
    assert "Int intermediate_compression_level = 1" in wdl
    assert "--intermediate_compression_level ~{intermediate_compression_level}" in wdl
