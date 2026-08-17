#!/usr/bin/env python3

"""Weighted quantification must refuse the outputs it cannot make correct.

XW weighting fixes the aggregate expression file: its counts become estimates on the
scale of the original bam. It cannot fix the tracking file, which carries one row per
*retained* read -- the reads normalization discarded have no row to hold a weight. So
any consumer that derives counts by summing tracking rows sees roughly the retained
fraction of the library, unevenly, worst at exactly the deep loci weighting exists to
correct.

Two such consumers ship in this repo, and both must be refused rather than served
output that looks finished:

  - util/sc/singlecell_tracking_to_sparse_matrix.py sums frac_assigned per cell and
    never reads XW, so every cell is undercounted;
  - util/annotate_bam_with_read_tracking_info.py tags reads from tracking, leaving
    thinned reads untagged and indistinguishable from unassigned ones.

Single-cell input cannot be detected reliably -- any bounded scan can miss a bam whose
early records are untagged, and a false negative emits the silently wrong matrix this
guard exists to prevent -- so the declaration is required and detection is used only to
contradict it.
"""

import os
import subprocess
import sys

import pysam
import pytest

LRAA = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "LRAA")
CONTIG = "chr1"


def _bam(path, *, single_cell):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"},
         "SQ": [{"SN": CONTIG, "LN": 10000}]}
    )
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for i in range(10):
            a = pysam.AlignedSegment(header)
            a.query_name = f"r{i}"
            a.reference_id = 0
            a.reference_start = 100 + i
            a.mapping_quality = 60
            a.cigarstring = "100M"
            a.query_sequence = "A" * 100
            a.query_qualities = pysam.qualitystring_to_array("I" * 100)
            if single_cell:
                a.set_tag("CB", "BARCODE-1")
                a.set_tag("XM", f"UMI{i}")
            fh.write(a)
    pysam.index(str(path))
    return path


def _run(bam, *extra, gtf=None, genome=None, tmp=None):
    cmd = [sys.executable, LRAA, "--quant_only", "--bam", str(bam),
           "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp / "out"),
           "--use_XW_read_weights_for_quant"] + list(extra)
    return subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp))


@pytest.fixture
def inputs(tmp_path):
    gtf = tmp_path / "a.gtf"
    gtf.write_text(
        f'{CONTIG}\ttest\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    )
    genome = tmp_path / "g.fa"
    genome.write_text(f">{CONTIG}\n" + "A" * 10000 + "\n")
    return gtf, genome


def test_weighted_quant_demands_an_explicit_library_type(tmp_path, inputs):
    """Refusing to guess is the point: a wrong guess emits a wrong matrix silently."""
    gtf, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    r = _run(bam, gtf=gtf, genome=genome, tmp=tmp_path)
    assert r.returncode != 0
    assert "requires an explicit --library_type" in (r.stdout + r.stderr)


def test_weighted_quant_refuses_single_cell(tmp_path, inputs):
    gtf, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    r = _run(bam, "--library_type", "single_cell", gtf=gtf, genome=genome, tmp=tmp_path)
    assert r.returncode != 0
    assert "single-cell quantification is not supported" in (r.stdout + r.stderr)


def test_weighted_quant_refuses_a_cell_list(tmp_path, inputs):
    gtf, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    cells = tmp_path / "cells.txt"
    cells.write_text("BARCODE-1\n")
    r = _run(bam, "--library_type", "bulk", "--cell_list", str(cells),
             gtf=gtf, genome=genome, tmp=tmp_path)
    assert r.returncode != 0
    assert "single-cell quantification is not supported" in (r.stdout + r.stderr)


def test_a_bulk_declaration_is_overruled_by_tagged_reads(tmp_path, inputs):
    """Detection is only trusted to contradict, never to clear -- a positive is sound."""
    gtf, genome = inputs
    bam = _bam(tmp_path / "sc.bam", single_cell=True)
    r = _run(bam, "--library_type", "bulk", gtf=gtf, genome=genome, tmp=tmp_path)
    assert r.returncode != 0
    out = r.stdout + r.stderr
    assert "--library_type bulk was declared" in out
    assert "CB" in out and "XM" in out


def test_weighted_quant_refuses_bam_tagging(tmp_path, inputs):
    gtf, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    r = _run(bam, "--library_type", "bulk", "--tag_bam",
             gtf=gtf, genome=genome, tmp=tmp_path)
    assert r.returncode != 0
    assert "--tag_bam is not supported" in (r.stdout + r.stderr)


def test_weighted_quant_refuses_discovery_mode(tmp_path, inputs):
    """Selecting isoforms and filtering them on different scales is not a mode.

    Discovery's absolute support gates -- min_reads_novel_isoform,
    min_unique_reads_novel_isoform, the FSM gates -- still count retained reads. Weighted
    EM alongside them would build isoform sets under one notion of support and filter
    them under another, which nothing has evaluated.
    """
    _, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    cmd = [sys.executable, LRAA, "--bam", str(bam), "--genome", str(genome),
           "--output_prefix", str(tmp_path / "out"),
           "--use_XW_read_weights_for_quant", "--library_type", "bulk"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode != 0
    assert "supported only with --quant_only" in (r.stdout + r.stderr)


def test_stream_reads_requires_rescue_turned_off(tmp_path, inputs):
    """The streaming pass maps genomically and cannot reproduce transcriptome rescue.

    Rescue is on by default, and its first pass rescues the reads that failed genomic
    assignment in the NORMALIZED bam. Reads present only in the full bam reach the second
    pass and are assigned genomically or resolved, never rescued -- so the run's
    assignments would differ from the default path's while the rescue summary described
    only the first pass.
    """
    gtf, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    cmd = [sys.executable, LRAA, "--quant_only", "--bam", str(bam),
           "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp_path / "out"),
           "--stream_reads"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode != 0
    assert "cannot reproduce transcriptome rescue" in (r.stdout + r.stderr)
    assert not list(tmp_path.glob("out*quant.expr")), "and must not emit results"


def test_stream_reads_needs_a_thinner_first_pass_bam(tmp_path, inputs):
    """With nothing to thin pass 1, both passes read the same bam and the mode only adds work."""
    gtf, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    cmd = [sys.executable, LRAA, "--quant_only", "--bam", str(bam),
           "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp_path / "out"),
           "--stream_reads",
           "--no_rescue_unassigned_reads_via_transcriptome_alignment",
           "--normalize_max_cov_level", "0"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode != 0
    assert "needs a first-pass bam thinner than the one it streams" in (r.stdout + r.stderr)


def test_stream_reads_accepts_an_externally_normalized_bam_for_sg(tmp_path, inputs):
    """A pre-normalized --bam_for_sg satisfies the gate, which is how chunking composes.

    The chunked pipeline normalizes each chunk itself and then runs its quant stage with
    --bam_for_sg <normalized> --no_norm. Keying the gate on --normalize_max_cov_level alone
    rejected every one of those calls, because the thinning had already happened in another
    process. What the mode actually needs is that pass 1 read a thinner bam than pass 2, and
    a distinct --bam_for_sg is exactly that.

    Asserts only that the gate lets the run proceed -- it must not fail on THIS check.
    """
    gtf, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    sg_bam = _bam(tmp_path / "thinned.bam", single_cell=False)
    cmd = [sys.executable, LRAA, "--quant_only", "--bam", str(bam),
           "--bam_for_sg", str(sg_bam), "--no_norm",
           "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp_path / "out"),
           "--stream_reads",
           "--no_rescue_unassigned_reads_via_transcriptome_alignment"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert "needs a first-pass bam thinner than the one it streams" not in (
        r.stdout + r.stderr
    )


def test_stream_reads_refuses_when_bam_for_sg_is_the_streamed_bam(tmp_path, inputs):
    """The same bam under both flags is the case the gate exists for, spelled differently."""
    gtf, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    cmd = [sys.executable, LRAA, "--quant_only", "--bam", str(bam),
           "--bam_for_sg", str(bam), "--no_norm",
           "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp_path / "out"),
           "--stream_reads",
           "--no_rescue_unassigned_reads_via_transcriptome_alignment"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode != 0
    assert "needs a first-pass bam thinner than the one it streams" in (r.stdout + r.stderr)


def test_stream_reads_refuses_tag_bam(tmp_path, inputs):
    """Tagging needs per-read assignments the streaming pass deliberately does not keep.

    Unguarded, the run would tag from the first pass alone -- that is, from the normalized
    bam -- and emit a bam whose tagged subset silently omits every read normalization removed.
    """
    gtf, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    cmd = [sys.executable, LRAA, "--quant_only", "--bam", str(bam),
           "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp_path / "out"),
           "--stream_reads",
           "--no_rescue_unassigned_reads_via_transcriptome_alignment", "--tag_bam"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode != 0
    assert "cannot be combined with --tag_bam" in (r.stdout + r.stderr)


def test_none_of_these_guards_fire_without_the_weighting_flag(tmp_path, inputs):
    """The default path must be untouched: these refusals belong to the flag alone."""
    gtf, genome = inputs
    bam = _bam(tmp_path / "sc.bam", single_cell=True)
    cmd = [sys.executable, LRAA, "--quant_only", "--bam", str(bam),
           "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp_path / "out"), "--tag_bam"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    out = r.stdout + r.stderr
    assert "not supported with --use_XW_read_weights_for_quant" not in out
    assert "requires an explicit --library_type" not in out
