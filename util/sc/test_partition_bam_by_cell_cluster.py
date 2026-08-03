#!/usr/bin/env python3

import subprocess
import sys
from pathlib import Path

import pysam
import pytest


def _make_alignment(name, start, tag, barcode):
    alignment = pysam.AlignedSegment()
    alignment.query_name = name
    alignment.query_sequence = "A" * 20
    alignment.flag = 0
    alignment.reference_id = 0
    alignment.reference_start = start
    alignment.mapping_quality = 60
    alignment.cigar = ((0, 20),)
    alignment.query_qualities = pysam.qualitystring_to_array("I" * 20)
    alignment.set_tag(tag, barcode)
    return alignment


@pytest.mark.parametrize("cell_barcode_tag", ["CB", "XC"])
def test_partitions_bam_using_configured_cell_barcode_tag(tmp_path, cell_barcode_tag):
    input_bam = tmp_path / "input.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
    }
    other_tag = "XC" if cell_barcode_tag == "CB" else "CB"

    with pysam.AlignmentFile(input_bam, "wb", header=header) as bam_writer:
        bam_writer.write(_make_alignment("selected", 100, cell_barcode_tag, "barcode1"))
        bam_writer.write(
            _make_alignment("unrecognized", 200, cell_barcode_tag, "barcode2")
        )
        bam_writer.write(_make_alignment("other-tag", 300, other_tag, "barcode1"))

    clusters = tmp_path / "clusters.tsv"
    clusters.write_text("barcode1\tcluster1\n")
    output_prefix = tmp_path / "sample"
    script = Path(__file__).with_name("partition_bam_by_cell_cluster.py")

    result = subprocess.run(
        [
            sys.executable,
            str(script),
            "--bam",
            str(input_bam),
            "--cell_clusters",
            str(clusters),
            "--output_prefix",
            str(output_prefix),
            "--cell_barcode_tag",
            cell_barcode_tag,
            "--threads",
            "1",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    with pysam.AlignmentFile(tmp_path / "sample.cluster1.bam", "rb") as bam_reader:
        assert [read.query_name for read in bam_reader] == ["selected"]

    assert f"Reads without cell barcode tag {cell_barcode_tag}: 1" in result.stderr
    assert "Reads with unrecognized cluster: 1" in result.stderr
