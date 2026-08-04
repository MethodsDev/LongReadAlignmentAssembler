#!/usr/bin/env python3

import importlib.util
from importlib.machinery import SourceFileLoader
from io import StringIO
from pathlib import Path

import pysam
import LRAA_Globals
from Transcript import Transcript


def _load_lraa_module():
    lraa_path = Path(__file__).resolve().parents[1] / "LRAA"
    loader = SourceFileLoader("lraa_oversimplify_composite_test", str(lraa_path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


class _FakePrettyAlignment:
    def get_pretty_alignment_segments(self):
        return [[1, 100]]

    def get_read_name(self):
        return "cell-1^umi-1"


class _FakePrettyAlignmentManager:
    def __init__(self, splice_graph=None):
        pass

    def retrieve_pretty_alignments(self, *args, **kwargs):
        return [_FakePrettyAlignment()]


def _write_repeated_alignment_bam(bam_path):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chrM", "LN": 500}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam_writer:
        for flag, reference_start, cigar_tuples in (
            (0, 0, [(0, 30), (4, 20)]),
            (256, 200, [(0, 50)]),
        ):
            read = pysam.AlignedSegment()
            read.query_name = "read1"
            read.query_sequence = "A" * 50
            read.flag = flag
            read.reference_id = 0
            read.reference_start = reference_start
            read.mapping_quality = 60
            read.cigartuples = cigar_tuples
            read.query_qualities = pysam.qualitystring_to_array("I" * 50)
            read.set_tag("NM", 0)
            bam_writer.write(read)
    pysam.index(str(bam_path))
    return str(bam_path)


def _build_overlap_transcripts():
    primary_tx = Transcript("chrM", [[1, 30]], "+")
    primary_tx.set_gene_id("g1")
    primary_tx.set_transcript_id("t1")

    secondary_tx = Transcript("chrM", [[201, 250]], "+")
    secondary_tx.set_gene_id("g1")
    secondary_tx.set_transcript_id("t2")
    return primary_tx, secondary_tx


def test_oversimplify_reports_composite_gene_identifier(monkeypatch):
    lraa = _load_lraa_module()
    monkeypatch.setattr(lraa, "Pretty_alignment_manager", _FakePrettyAlignmentManager)

    transcript = Transcript("chrM", [[1, 100]], "+")
    transcript.set_gene_id("ENSG00000198888")
    transcript.set_transcript_id("ENST00000361390")
    transcript.add_meta("gene_name", "MT-ND1")

    old_num_total_reads = LRAA_Globals.config["num_total_reads"]
    LRAA_Globals.config["num_total_reads"] = 1
    try:
        quant_output = StringIO()
        tracking_output = StringIO()
        lraa._run_oversimplify_best_overlap(
            "chrM",
            "+",
            "A" * 100,
            "unused.bam",
            None,
            None,
            [transcript],
            quant_output,
            tracking_output,
        )
    finally:
        LRAA_Globals.config["num_total_reads"] = old_num_total_reads

    assert quant_output.getvalue().split("\t", 1)[0] == "MT-ND1^ENSG00000198888"
    assert tracking_output.getvalue().split("\t", 1)[0] == "MT-ND1^ENSG00000198888"


def test_oversimplify_best_overlap_uses_primary_alignment_only(tmp_path, monkeypatch):
    lraa = _load_lraa_module()
    bam_path = _write_repeated_alignment_bam(tmp_path / "repeated.bam")
    primary_tx, secondary_tx = _build_overlap_transcripts()
    monkeypatch.setenv("LRAA_TMP_DIR", str(tmp_path))
    monkeypatch.setitem(LRAA_Globals.config, "oversimplify_enabled", True)
    monkeypatch.setitem(LRAA_Globals.config, "oversimplify_contigs", ["chrM"])
    monkeypatch.setitem(LRAA_Globals.config, "allow_secondary_alignments", True)
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 1)

    quant_output = StringIO()
    tracking_output = StringIO()
    lraa._run_oversimplify_best_overlap(
        "chrM",
        "+",
        "A" * 500,
        bam_path,
        None,
        None,
        [primary_tx, secondary_tx],
        quant_output,
        tracking_output,
    )

    quant_by_transcript = {
        row[1]: row
        for row in map(
            lambda line: line.split("\t"), quant_output.getvalue().splitlines()
        )
    }
    assert quant_by_transcript["t1"][2:4] == ["1", "1.0"]
    assert quant_by_transcript["t2"][2:4] == ["0", "0.0"]

    tracking_rows = [
        line.split("\t") for line in tracking_output.getvalue().splitlines()
    ]
    assert len(tracking_rows) == 1
    assert tracking_rows[0][1] == "t1"
    assert tracking_rows[0][5] == "read1"


def test_oversimplify_aggregate_uses_primary_alignment_only(tmp_path, monkeypatch):
    lraa = _load_lraa_module()
    bam_path = _write_repeated_alignment_bam(tmp_path / "repeated.bam")
    monkeypatch.setenv("LRAA_TMP_DIR", str(tmp_path))
    monkeypatch.setitem(LRAA_Globals.config, "oversimplify_enabled", True)
    monkeypatch.setitem(LRAA_Globals.config, "oversimplify_contigs", ["chrM"])
    monkeypatch.setitem(LRAA_Globals.config, "allow_secondary_alignments", True)
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 1)

    gtf_output = StringIO()
    quant_output = StringIO()
    tracking_output = StringIO()
    lraa._run_oversimplify_aggregate(
        "chrM",
        "+",
        "A" * 500,
        bam_path,
        None,
        None,
        gtf_output,
        quant_output,
        tracking_output,
    )

    quant_row = quant_output.getvalue().rstrip().split("\t")
    assert quant_row[2:4] == ["1", "1.0"]
    tracking_rows = tracking_output.getvalue().splitlines()
    assert len(tracking_rows) == 1
    assert tracking_rows[0].split("\t")[5] == "read1"
