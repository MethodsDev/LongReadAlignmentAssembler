#!/usr/bin/env python3
"""Gene-body rows from the two LRAA paths outside the streaming pass.

  1. --oversimplify (chrM by default): each best-overlap read gets an "isoform" row for
     its chosen model's gene, so MT- genes reach gene_full; a read with no exonic
     overlap goes to GeneBodyIndex.
  2. run_quant_only's early exits (no splice graph / no model mapped): a geometry-only
     pass gives every retained read its gene-body row.
"""

import importlib.util
from importlib.machinery import SourceFileLoader
from io import StringIO
from pathlib import Path

import pysam

import LRAA_Globals
from Transcript import Transcript


def _load_lraa_module():
    lraa_path = Path(__file__).resolve().parents[1] / "LRAA"
    loader = SourceFileLoader("lraa_gene_body_paths_test", str(lraa_path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def _tx(contig, exons, gene, tid, strand="+"):
    t = Transcript(contig, [list(e) for e in exons], strand)
    t.set_gene_id(gene)
    t.set_transcript_id(tid)
    return t


def _rows(fh):
    return [l.split("\t") for l in fh.getvalue().splitlines()]


def test_oversimplify_writes_gene_body_rows(monkeypatch, tmp_path):
    lraa = _load_lraa_module()
    monkeypatch.setenv("LRAA_TMP_DIR", str(tmp_path))
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 1)

    class _Alignment:
        def __init__(self, name, segs):
            self._name, self._segs = name, segs

        def get_pretty_alignment_segments(self):
            return self._segs

        def get_read_name(self):
            return self._name

    def _alignments(*args, **kwargs):
        yield _Alignment("bc^u^exonic", [[10, 90]])
        yield _Alignment("bc^u^intronic", [[150, 250]])
        yield _Alignment("bc^u^outside", [[2000, 2100]])

    monkeypatch.setattr(lraa, "_oversimplify_alignments", _alignments)
    t = _tx("chrM", [(1, 100), (300, 400)], "MT-X", "t1")
    gene_body = StringIO()
    lraa._run_oversimplify_best_overlap(
        "chrM", "+", "A" * 3000, "unused.bam", None, None, [t], StringIO(), StringIO(),
        gene_body_fh=gene_body,
    )
    rows = _rows(gene_body)
    assert [(r[0], r[1], r[2], r[3], r[4]) for r in rows] == [
        ("MT-X", "bc^u^exonic", "1.000000", "A", "isoform"),
        ("MT-X", "bc^u^intronic", "1.000000", "U", "gene_body:oversimplify_no_exon_overlap"),
    ]


def _write_bam(bam_path):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 20000}]}
    records = [  # (name, flag, start0, cigar, mapq)
        ("exonic", 0, 1000, [(0, 500)], 60),
        ("intronic", 0, 3000, [(0, 800)], 60),
        ("antisense", 16, 3000, [(0, 800)], 60),
        ("secondary", 256, 3000, [(0, 800)], 60),
        ("intergenic", 0, 15000, [(0, 500)], 60),
    ]
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as w:
        for name, flag, start, cigar, mapq in records:
            n = sum(l for op, l in cigar if op in (0, 4))
            r = pysam.AlignedSegment()
            r.query_name = name
            r.query_sequence = "C" * n
            r.flag = flag
            r.reference_id = 0
            r.reference_start = start
            r.mapping_quality = mapq
            r.cigartuples = cigar
            r.query_qualities = pysam.qualitystring_to_array("I" * n)
            r.set_tag("NM", 0)
            w.write(r)
    pysam.index(str(bam_path))
    return str(bam_path)


def test_geometry_only_pass_for_units_without_isoform_quant(monkeypatch, tmp_path):
    lraa = _load_lraa_module()
    monkeypatch.setenv("LRAA_TMP_DIR", str(tmp_path))
    bam = _write_bam(tmp_path / "r.bam")
    t = _tx("chr1", [(1001, 2000), (5001, 6000)], "G", "t1")
    fh = StringIO()
    lraa._gene_body_geometry_only_pass(fh, [t], "chr1", "+", bam, None, None, "no_quant_unit")
    rows = _rows(fh)
    assert [(r[0], r[1], r[2], r[3], r[4]) for r in rows] == [
        ("G", "exonic", "1.000000", "A", "gene_body:no_quant_unit"),
        ("G", "intronic", "1.000000", "U", "gene_body:no_quant_unit"),
    ]
    # no annotation on the strand: nothing written, nothing read
    fh = StringIO()
    lraa._gene_body_geometry_only_pass(fh, [], "chr1", "+", bam, None, None, "no_quant_unit")
    assert fh.getvalue() == ""
