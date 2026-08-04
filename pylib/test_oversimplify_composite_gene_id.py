#!/usr/bin/env python3

import importlib.util
from importlib.machinery import SourceFileLoader
from io import StringIO
from pathlib import Path

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
