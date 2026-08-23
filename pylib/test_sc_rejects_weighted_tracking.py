#!/usr/bin/env python3

"""An incomplete tracking file must stay recognizable, and be refused where it is unsafe.

Per-cell counts are sums of `frac_assigned`, one row per read, and `XW` is never
consulted. A tracking file produced with XW weighting applied to a pre-thinned bam has
rows only for the reads coverage normalization retained, so every cell comes out short
-- unevenly, and worst at the high-coverage loci weighting exists to correct.

LRAA no longer produces such a file. `--use_XW_read_weights_for_quant` and its negation
are gone, weighting is unconditional, and the input that made a tracking file incomplete
-- an already-thinned `--bam` -- is now refused at startup by
`LRAA._require_no_thinning_weights`. The marker is therefore no longer emitted by
anything in this tree.

The two READERS keep recognizing it anyway, deliberately, because tracking files written
by earlier versions carry it legitimately and are still safe to hand to these tools. So
two things have to hold. The single-cell converter must refuse such a file: it configures
pandas with `comment="#"`, so the marker is invisible to exactly the consumer that is
unsafe, and has to be read before pandas sees the file. And the marker has to survive
merging: `lraa_merge_header.py` rebuilds the comment block rather than copying it, so a
marker that is not recognized explicitly is dropped and the merged file understates
counts with nothing about it to say so.

Everything below the merge section therefore works on HAND-BUILT tracking files, since
no LRAA invocation can produce one any more. The single live-LRAA test at the bottom
pins the other half: that a stock run writes no marker.
"""

import gzip
import os
import subprocess
import sys

import pytest

UTIL = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "util")
CONVERTER = os.path.join(UTIL, "sc", "singlecell_tracking_to_sparse_matrix.py")
MERGE_HEADER = os.path.join(UTIL, "lraa_merge_header.py")
LRAA = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "LRAA")

HEADER = (
    "gene_id\ttranscript_id\ttranscript_splice_hash_code\tnum_exons\tmp_id\t"
    "read_name\tfrac_assigned\tread_weight\n"
)
ROW = "g1\tt1\th1\t2\tMP1\tBC1^UMI1^r1\t1.000000\t1.000000\n"
MARKER = (
    "# WARNING: produced with --use_XW_read_weights_for_quant. Rows cover only reads "
    "retained by coverage normalization, and frac_assigned is NOT weighted by XW.\n"
)


def _convert(tracking, prefix):
    return subprocess.run(
        [sys.executable, CONVERTER, "--tracking", str(tracking),
         "--output_prefix", str(prefix)],
        capture_output=True, text=True, cwd=str(tracking.parent),
    )


# ---------------------------------------------------------------- converter refusal


def test_weighted_tracking_is_refused(tmp_path):
    t = tmp_path / "weighted.tracking"
    t.write_text("# LRAA version test\n" + MARKER + HEADER + ROW)
    r = _convert(t, tmp_path / "w")
    assert r.returncode != 0, "a file that would undercount every cell must not convert"
    assert "use_XW_read_weights_for_quant" in (r.stdout + r.stderr)


def test_refusal_happens_before_anything_is_written(tmp_path):
    """Otherwise a partial matrix survives the failure and looks usable."""
    t = tmp_path / "weighted.tracking"
    t.write_text("# LRAA version test\n" + MARKER + HEADER + ROW)
    _convert(t, tmp_path / "w")
    assert not list(tmp_path.glob("w*sparseM*")), "no matrix may be emitted"


def test_ordinary_tracking_still_converts(tmp_path):
    """The guard must not cost the normal path: only the marker triggers it."""
    t = tmp_path / "plain.tracking"
    t.write_text("# LRAA version test\n" + HEADER + ROW)
    r = _convert(t, tmp_path / "p")
    assert r.returncode == 0, r.stdout + r.stderr
    assert list(tmp_path.glob("p*sparseM*")), "the normal path must still emit matrices"


def test_the_marker_is_only_honored_in_the_leading_comment_block(tmp_path):
    """A data row mentioning the flag is data, not a marker.

    The scan stops at the first non-comment line, so an identifier containing the string
    cannot make the converter refuse a legitimate file.
    """
    t = tmp_path / "odd.tracking"
    t.write_text(
        "# LRAA version test\n" + HEADER
        + "use_XW_read_weights_for_quant\tt1\th1\t2\tMP1\tBC1^UMI1^r1\t1.000000\t1.000000\n"
    )
    r = _convert(t, tmp_path / "o")
    assert r.returncode == 0, r.stdout + r.stderr


# ------------------------------------------------------------- marker survives merge


def _merged_header(tmp_path, *inputs):
    out = tmp_path / "merged.header"
    r = subprocess.run(
        [sys.executable, MERGE_HEADER,
         "--version_comment", "# LRAA version test",
         "--output", str(out),
         "--inputs"] + [str(p) for p in inputs],
        capture_output=True, text=True, cwd=str(tmp_path),
    )
    assert r.returncode == 0, r.stdout + r.stderr
    return out.read_text()


def test_merged_header_preserves_the_marker(tmp_path):
    a = tmp_path / "a.tracking"
    a.write_text("# LRAA version test\n# LRAA CMD: LRAA --contig chr1\n" + MARKER + HEADER + ROW)
    b = tmp_path / "b.tracking"
    b.write_text("# LRAA version test\n# LRAA CMD: LRAA --contig chr2\n" + MARKER + HEADER + ROW)
    assert "use_XW_read_weights_for_quant" in _merged_header(tmp_path, a, b)


def test_one_marked_input_marks_the_merge(tmp_path):
    """The merge inherits an input's gaps, so any marked input marks the result."""
    a = tmp_path / "a.tracking"
    a.write_text("# LRAA version test\n# LRAA CMD: LRAA --contig chr1\n" + MARKER + HEADER + ROW)
    b = tmp_path / "b.tracking"
    b.write_text("# LRAA version test\n# LRAA CMD: LRAA --contig chr2\n" + HEADER + ROW)
    assert "use_XW_read_weights_for_quant" in _merged_header(tmp_path, a, b)


def test_unmarked_inputs_produce_an_unmarked_merge(tmp_path):
    a = tmp_path / "a.tracking"
    a.write_text("# LRAA version test\n# LRAA CMD: LRAA --contig chr1\n" + HEADER + ROW)
    b = tmp_path / "b.tracking"
    b.write_text("# LRAA version test\n# LRAA CMD: LRAA --contig chr2\n" + HEADER + ROW)
    assert "use_XW_read_weights_for_quant" not in _merged_header(tmp_path, a, b)


def test_a_merged_marked_file_is_still_refused_by_the_converter(tmp_path):
    """End to end: the marker has to survive the merge *and* be honored afterwards."""
    a = tmp_path / "a.tracking"
    a.write_text("# LRAA version test\n# LRAA CMD: LRAA --contig chr1\n" + MARKER + HEADER + ROW)
    b = tmp_path / "b.tracking"
    b.write_text("# LRAA version test\n# LRAA CMD: LRAA --contig chr2\n" + MARKER + HEADER + ROW)
    header = _merged_header(tmp_path, a, b)

    merged = tmp_path / "merged.tracking"
    with open(merged, "wt") as fh:
        fh.write(header if header.endswith("\n") else header + "\n")
        fh.write(HEADER)
        fh.write(ROW)
        fh.write(ROW)

    r = _convert(merged, tmp_path / "m")
    assert r.returncode != 0, "the merged file is just as incomplete as its inputs"
    assert "use_XW_read_weights_for_quant" in (r.stdout + r.stderr)


# ---------------------------------------------------- LRAA no longer writes the marker


def _lraa_bam(path):
    import pysam

    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"},
         "SQ": [{"SN": "chr1", "LN": 10000}]}
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
            fh.write(a)
    pysam.index(str(path))
    return path


def test_a_stock_lraa_run_writes_no_incompleteness_marker(tmp_path):
    """The removal is pinned, not merely untested.

    This used to be a pair of tests asserting LRAA emitted the marker on a non-streaming
    run and omitted it on a streaming one, both driven through
    `--use_XW_read_weights_for_quant --library_type bulk`. Neither the flag nor the
    emission exists now: weighting is unconditional and the input that made a tracking
    file incomplete, an already-thinned `--bam`, is refused at startup. What remains
    worth holding is that no code path writes the marker back, because re-introducing it
    would make every fresh tracking file unreadable to
    `singlecell_tracking_to_sparse_matrix.py` -- which still refuses it, on purpose, for
    the sake of files written by earlier versions.

    The absence is asserted against a file that demonstrably exists and has both a
    comment block and data rows, so "no marker" cannot pass by the tracking file being
    empty, unwritten, or all comment. Note the tracking file also carries a
    "# LRAA CMD: ..." echo of its own command line, which is why the check is scoped to
    the dedicated "# WARNING:" prefix -- the same scoping both readers use.
    """
    gtf = tmp_path / "a.gtf"
    gtf.write_text(
        'chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    )
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\n" + "A" * 10000 + "\n")
    bam = _lraa_bam(tmp_path / "bulk.bam")
    # --no_chunk because chunk dispatch would reroute this small single-locus fixture
    # through a separate merge path; everything else is left at its default, which is
    # the point -- this is what an ordinary user's run emits.
    cmd = [sys.executable, LRAA, "--quant_only", "--bam", str(bam),
           "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp_path / "out"), "--no_chunk"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode == 0, r.stdout + r.stderr

    tracking = list(tmp_path.glob("out*quant.tracking.gz"))
    assert tracking, "expected a tracking file to be produced"
    with gzip.open(tracking[0], "rt") as fh:
        text = fh.read()

    lines = text.splitlines()
    comment_lines = [line for line in lines if line.startswith("#")]
    data_lines = [line for line in lines if line and not line.startswith("#")]
    assert comment_lines, "the tracking file must still carry its comment header block"
    assert data_lines, "the tracking file must still carry read rows"
    assert not [
        line for line in comment_lines if line.startswith("# WARNING:")
    ], f"LRAA emitted an incompleteness marker: {comment_lines}"
