#!/usr/bin/env python3

import sys
from pathlib import Path

import pysam

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA_Globals
from Bam_alignment_extractor import Bam_alignment_extractor


def _alignment(read_name, flag, ref_id, start):
    aln = pysam.AlignedSegment()
    aln.query_name = read_name
    aln.flag = flag
    aln.reference_id = ref_id
    aln.reference_start = start
    aln.mapping_quality = 60
    aln.cigar = [(0, 50)]
    aln.query_sequence = "A" * 50
    aln.query_qualities = pysam.qualitystring_to_array("I" * 50)
    aln.set_tag("NM", 0)
    aln.set_tag("AS", 50)
    return aln


def _write_test_bam(path):
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(path, "wb", header=header) as outf:
        outf.write(_alignment("read1", 0, 0, 100))
        outf.write(_alignment("read1", 256, 0, 150))
        outf.write(_alignment("read1", 2048, 0, 200))
    pysam.index(str(path))


def test_supplementary_alignments_are_always_filtered(tmp_path):
    bam_path = tmp_path / "input.bam"
    _write_test_bam(bam_path)

    previous_allow_secondary = LRAA_Globals.config.get("allow_secondary_alignments")
    try:
        LRAA_Globals.config["allow_secondary_alignments"] = True
        read_flags = [
            read.flag
            for read in Bam_alignment_extractor(str(bam_path)).get_read_alignments("chr1")
        ]
        assert read_flags == [0, 256]

        LRAA_Globals.config["allow_secondary_alignments"] = False
        read_flags = [
            read.flag
            for read in Bam_alignment_extractor(str(bam_path)).get_read_alignments("chr1")
        ]
        assert read_flags == [0]
    finally:
        LRAA_Globals.config["allow_secondary_alignments"] = previous_allow_secondary
