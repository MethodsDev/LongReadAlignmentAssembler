#!/usr/bin/env python3

"""`--chunk` forwards single-cell and streaming settings to its chunk workers.

Formerly ``test_chunk_refuses_single_cell.py``: ``pylib/ChunkedRun.py`` had no
reference to ``cell_barcode_tag``, ``read_umi_tag``, ``cell_list`` or
``stream_reads``, and ``lraa_cmd`` forwarded nine setting-bearing flags plus a
forced ``--no_norm``. A chunked run given single-cell flags therefore dropped
them and produced BULK output -- complete, internally consistent, and
indistinguishable from a correct single-cell result unless someone inspected
read names, because ``Util_funcs.get_read_name_include_sc_encoding`` composes
``barcode^umi^query_name`` only when both tags are configured and returns the
bare query name otherwise. Rather than let that survive, ``--chunk`` REFUSED
the combination outright (this file's previous incarnation).

C2 (CodeBacklog.chunking-single-cell.md) replaces the refusal with the real
forward: ``lraa_cmd`` now emits the tags, ``--cell_list`` (resolved to an
absolute path -- a chunk worker's cwd is its own chunk directory) and
``--stream_reads`` plus its rescue flags to every chunk worker, on presence
exactly as LRAA's own ``argparse.SUPPRESS`` convention already distinguished a
user-supplied tag from the default every bulk run also carries. These tests
assert the forward reaches the worker's actual argv, read from its chunk log
-- not merely that the run exits 0, which a dropped flag would too.

Uses the default ``--strandless_chunks``, not ``--chunk_by_strand``. It did not
always: this minimal fixture (and, separately,
``testing/single_cells/data/chr19.gtf`` at whole-gene scale) used to trip an
unrelated bug in the strandless split's own annotation accounting (``chunk
chr1_00 annotation split accounting: 0 + 0 = 0 transcript line(s) ... but
extraction emitted 1``), reproducible with no single-cell flag at all -- fixed
in c8c09f2 (``split_chunk_gtf_by_strand`` now counts distinct transcript_id
rather than literal ``transcript`` feature rows). Running the production
default here exercises that path instead of avoiding it.
"""

import gzip
import glob
import os
import subprocess
import sys

import pysam
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
LRAA = os.path.join(REPO, "LRAA")
CONTIG = "chr1"


def _bam(path, *, barcode_tag=None, umi_tag=None):
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": CONTIG, "LN": 10000}],
        }
    )
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for i in range(10):
            a = pysam.AlignedSegment(header)
            a.query_name = "r{}".format(i)
            a.reference_id = 0
            a.reference_start = 100 + i
            a.mapping_quality = 60
            a.cigarstring = "100M"
            a.query_sequence = "A" * 100
            a.query_qualities = pysam.qualitystring_to_array("I" * 100)
            if barcode_tag:
                a.set_tag(barcode_tag, "BARCODE-1")
            if umi_tag:
                a.set_tag(umi_tag, "UMI{}".format(i))
            fh.write(a)
    pysam.index(str(path))
    return path


@pytest.fixture
def inputs(tmp_path):
    gtf = tmp_path / "a.gtf"
    gtf.write_text(
        '{}\ttest\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'.format(
            CONTIG
        )
    )
    genome = tmp_path / "g.fa"
    genome.write_text(">{}\n".format(CONTIG) + "A" * 10000 + "\n")
    bam = _bam(tmp_path / "reads.bam", barcode_tag="CB", umi_tag="XM")
    return bam, gtf, genome


def _run(tmp_path, bam, gtf, genome, *extra):
    cmd = [
        sys.executable,
        LRAA,
        "--bam", str(bam),
        "--genome", str(genome),
        "--gtf", str(gtf),
        "--quant_only",
        "--output_prefix", str(tmp_path / "out"),
    ] + list(extra)
    return subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(tmp_path), timeout=900
    )


def _chunked(tmp_path, bam, gtf, genome, *extra):
    return _run(
        tmp_path,
        bam,
        gtf,
        genome,
        "--chunk",
        "--chunk_work_dir", str(tmp_path / "work"),
        *extra,
    )


def _worker_argv(tmp_path):
    """Every chunk worker's real command line, concatenated.

    Read from the chunk logs rather than timing.json: timing.json records
    stage-level aggregates, not the per-chunk quant argv, so the log's own
    ``cmd: ...`` line -- written by ``run_step`` for every subprocess it
    launches -- is the only place this is observable from outside the process.
    """
    logs = sorted(glob.glob(str(tmp_path / "work" / "logs" / "chunk_*.log")))
    assert logs, "no chunk logs written"
    return "\n".join(open(p).read() for p in logs)


def _assert_succeeds_and_forwards(result, tmp_path, *expected):
    combined = result.stdout + result.stderr
    assert result.returncode == 0, combined[-3000:]
    assert list(tmp_path.glob("out*quant.expr")), "no output written"
    argv = _worker_argv(tmp_path)
    for flag in expected:
        assert flag in argv, argv


def test_chunk_forwards_a_user_specified_cell_barcode_tag(tmp_path, inputs):
    bam, gtf, genome = inputs
    r = _chunked(tmp_path, bam, gtf, genome, "--cell_barcode_tag", "ZB")
    _assert_succeeds_and_forwards(r, tmp_path, "--cell_barcode_tag ZB")


def test_chunk_forwards_a_user_specified_read_umi_tag(tmp_path, inputs):
    bam, gtf, genome = inputs
    r = _chunked(tmp_path, bam, gtf, genome, "--read_umi_tag", "ZM")
    _assert_succeeds_and_forwards(r, tmp_path, "--read_umi_tag ZM")


def test_chunk_forwards_a_tag_specified_at_its_own_default_value(tmp_path, inputs):
    """``--cell_barcode_tag CB`` still forwards, and that is deliberate.

    The distinction that matters is "did the user ask for single cell", not
    "is this value unusual". A run that spells the default out is asking for
    the encoding just as much as one that names a custom tag. This is also
    the assertion that keeps the ``argparse.SUPPRESS`` declaration honest: a
    value comparison against the config default would silently drop this
    exact case while every other test in this file still passed.
    """
    bam, gtf, genome = inputs
    r = _chunked(tmp_path, bam, gtf, genome, "--cell_barcode_tag", "CB")
    _assert_succeeds_and_forwards(r, tmp_path, "--cell_barcode_tag CB")


def test_chunk_forwards_a_cell_list(tmp_path, inputs):
    bam, gtf, genome = inputs
    cells = tmp_path / "cells.txt"
    cells.write_text("BARCODE-1\n")
    r = _chunked(tmp_path, bam, gtf, genome, "--cell_list", str(cells))
    # Forwarded as the ABSOLUTE path: a chunk worker's cwd is its own chunk
    # directory, so the relative path this test passed would not resolve there.
    _assert_succeeds_and_forwards(r, tmp_path, "--cell_list {}".format(cells))


def test_chunk_forwards_stream_reads(tmp_path, inputs):
    bam, gtf, genome = inputs
    r = _chunked(
        tmp_path,
        bam,
        gtf,
        genome,
        "--stream_reads",
        # Required unconditionally of this fixture: LRAA:1536-1552 refuses
        # --stream_reads without it or --no_rescue_unassigned_reads_via_...,
        # and every chunk worker is a full LRAA invocation that enforces its
        # own validations same as the unchunked path would.
        "--stream_reads_rescue_unassigned",
    )
    _assert_succeeds_and_forwards(
        r, tmp_path, "--stream_reads", "--stream_reads_rescue_unassigned"
    )


def test_chunk_forwards_all_four_settings_at_once(tmp_path, inputs):
    """All four together reach the same worker command, not four separate ones."""
    bam, gtf, genome = inputs
    cells = tmp_path / "cells.txt"
    cells.write_text("BARCODE-1\n")
    r = _chunked(
        tmp_path,
        bam,
        gtf,
        genome,
        "--cell_barcode_tag", "ZB",
        "--read_umi_tag", "ZM",
        "--cell_list", str(cells),
        "--stream_reads",
        "--stream_reads_rescue_unassigned",
    )
    _assert_succeeds_and_forwards(
        r,
        tmp_path,
        "--cell_barcode_tag ZB",
        "--read_umi_tag ZM",
        "--cell_list {}".format(cells),
        "--stream_reads",
        "--stream_reads_rescue_unassigned",
    )


def test_plain_chunk_is_untouched(tmp_path, inputs):
    """No single-cell flag, so nothing new is forwarded and the run still succeeds.

    The tag defaults are always populated in the config, so ``hasattr``-gated
    forwarding is what a value-based check would have broken -- every chunked
    run in existence would otherwise carry ``--cell_barcode_tag CB
    --read_umi_tag XM`` whether asked for or not. Asserted directly: an
    unadorned run's worker argv must carry NEITHER tag flag.

    ``--stream_reads`` now defaults on, so a plain run has to say
    ``--no_stream_reads`` itself to mean what this test asks it to mean --
    otherwise the worker argv would legitimately carry ``--stream_reads`` and
    the assertion below would be testing the wrong thing.
    """
    bam, gtf, genome = inputs
    r = _chunked(tmp_path, bam, gtf, genome, "--no_stream_reads")
    combined = r.stdout + r.stderr
    assert r.returncode == 0, combined[-3000:]
    argv = _worker_argv(tmp_path)
    assert "--cell_barcode_tag" not in argv, argv
    assert "--read_umi_tag" not in argv, argv
    assert "--cell_list" not in argv, argv
    assert "--stream_reads" not in argv, argv
    assert "--no_stream_reads" in argv, argv


def test_unchunked_single_cell_is_untouched(tmp_path, inputs):
    """Unchunked single cell was always the supported path, and still is."""
    bam, gtf, genome = inputs
    cells = tmp_path / "cells.txt"
    cells.write_text("BARCODE-1\n")
    r = _run(
        tmp_path,
        bam,
        gtf,
        genome,
        "--cell_barcode_tag", "CB",
        "--read_umi_tag", "XM",
        "--cell_list", str(cells),
    )
    assert r.returncode == 0, (r.stdout + r.stderr)[-3000:]


def test_a_non_default_tag_still_reaches_the_read_name_encoding(tmp_path):
    """The SUPPRESS declaration must not stop the flag from doing its job.

    Declaring the tag flags ``argparse.SUPPRESS`` is what makes a user-supplied
    tag detectable, and it changes how the value reaches ``LRAA_Globals.config``.
    If that hand-off broke, a single-cell run would fall back to CB/XM and emit
    exactly the bulk-shaped read names this file's chunked tests exist to keep
    out of chunked mode too -- silently, and on the unchunked path C2 does not
    touch. So the encoding is asserted end to end, on a bam whose tags are
    deliberately not the defaults.
    """
    gtf = tmp_path / "a.gtf"
    gtf.write_text(
        '{}\ttest\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'.format(
            CONTIG
        )
    )
    genome = tmp_path / "g.fa"
    genome.write_text(">{}\n".format(CONTIG) + "A" * 10000 + "\n")
    bam = _bam(tmp_path / "sc.bam", barcode_tag="ZB", umi_tag="ZM")

    r = _run(
        tmp_path, bam, gtf, genome,
        "--cell_barcode_tag", "ZB",
        "--read_umi_tag", "ZM",
    )
    assert r.returncode == 0, (r.stdout + r.stderr)[-3000:]

    tracking = tmp_path / "out.LRAA.quant-only.quant.tracking.gz"
    assert tracking.exists(), sorted(p.name for p in tmp_path.iterdir())
    with gzip.open(tracking, "rt") as fh:
        rows = [
            line.rstrip("\n").split("\t") for line in fh if not line.startswith("#")
        ]
    # By column NAME, not position: the tracking header has gained columns before now.
    read_col = rows[0].index("read_name")
    read_names = [row[read_col] for row in rows[1:] if len(row) > read_col]
    assert read_names, "no tracking rows to inspect"
    # barcode^umi^query_name, composed only when BOTH configured tags are found on the
    # record. A dropped flag gives the bare "rN" here, which is the bulk shape.
    assert all(
        name.startswith("BARCODE-1^UMI") and name.count("^") == 2
        for name in read_names
    ), read_names[:5]


def test_the_wdl_runner_does_not_hand_every_chunked_run_a_tag_flag():
    """The one caller that emitted both tag flags on every command line.

    ``WDL/subwdls/LRAA_runner.wdl`` declares ``cell_barcode_tag``/``read_umi_tag``
    as non-optional Strings defaulting to CB/XM, so it used to append both flags
    unconditionally -- including under ``chunk = true``, which is how the WDL
    layer was silently getting bulk output for single-cell inputs before the
    refusal this file used to test, and how every chunked run (single cell or
    not) would forward tags it never meant to now that the forward is real. The
    flags are emitted only when they differ from those defaults: a bulk chunked
    run passes nothing and is unaffected, and a real single-cell tag reaches
    LRAA's own presence check.
    """
    with open(os.path.join(REPO, "WDL", "subwdls", "LRAA_runner.wdl")) as fh:
        lines = fh.read().splitlines()
    for flag, default in (("--cell_barcode_tag", "CB"), ("--read_umi_tag", "XM")):
        # An emission is the flag inside a WDL string literal, which is how it gets onto
        # the command line; the prose comments mentioning it are not.
        emitters = [line for line in lines if '"{} '.format(flag) in line]
        assert emitters, "nothing emits {} any more".format(flag)
        for line in emitters:
            assert '!= "{}"'.format(default) in line, line.strip()
