#!/usr/bin/env python3

"""`--chunk` must refuse single-cell and streaming settings instead of dropping them.

``pylib/ChunkedRun.py`` has no reference to ``cell_barcode_tag``, ``read_umi_tag``,
``cell_list`` or ``stream_reads``, and ``lraa_cmd`` forwards nine setting-bearing flags
plus a forced ``--no_norm``. So a chunked run given single-cell flags dropped them and
produced BULK output -- complete, internally consistent, and indistinguishable from a
correct single-cell result unless someone inspected read names, because
``Util_funcs.get_read_name_include_sc_encoding`` composes ``barcode^umi^query_name``
only when both tags are configured and returns the bare query name otherwise. That is
the worst failure shape available: nothing to notice and nothing to audit.

Two properties are held here, and the second is what makes the first safe:

* the combination is REFUSED, before the ``if args.chunk: sys.exit(_run_chunked_mode())``
  dispatch, so the chunked and unchunked paths cannot disagree about it -- the mistake
  ``b15430d`` fixed for ``--LowFi``, where a check sitting past that dispatch made one
  path fatal and the other silent;
* the refusal keys on a USER-SPECIFIED tag rather than on a tag's value.
  ``cell_barcode_tag`` and ``read_umi_tag`` always hold ``CB``/``XM``, so a
  value-based test would refuse every chunked run ever made. The tag flags therefore
  declare ``default=argparse.SUPPRESS`` and the check asks ``hasattr``, which is the
  mechanism ``--min_alt_splice_freq`` already uses to tell an explicit setting from a
  derived one. That declaration change has its own consumer test below: a non-default
  tag must still reach the read-name encoding on an unchunked run.

Making the combination WORK is a separate piece of work (forwarding the flags through
``lraa_cmd``); these tests are about not lying in the meantime.
"""

import gzip
import os
import subprocess
import sys

import pysam
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
LRAA = os.path.join(REPO, "LRAA")
CONTIG = "chr1"

# Text every refusal in this file must carry. Matched as a substring so the message can
# be reworded around it, but the run/rerun advice and the bulk-output reason are what
# make the refusal actionable, so both are asserted.
REFUSAL = "--chunk cannot be combined with"


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


def _assert_refused(result, flag, tmp_path):
    combined = result.stdout + result.stderr
    assert result.returncode != 0, combined[-3000:]
    assert REFUSAL in combined, combined[-3000:]
    # The flag has to be NAMED. A run refused for "some single-cell setting" leaves the
    # user re-reading their own command line to find which of four it was.
    assert flag in combined, combined[-3000:]
    # Why it matters, and what to do instead. Both were absent from the behaviour this
    # replaces, which is why it survived: the run simply succeeded.
    assert "BULK" in combined, combined[-3000:]
    assert "UNCHUNKED" in combined, combined[-3000:]
    # A refusal that still emitted results would leave the bulk output it warns about
    # sitting on disk under the requested prefix.
    assert not list(tmp_path.glob("out*quant.expr")), "refused, so nothing may be written"


def test_chunk_refuses_a_user_specified_cell_barcode_tag(tmp_path, inputs):
    bam, gtf, genome = inputs
    r = _chunked(tmp_path, bam, gtf, genome, "--cell_barcode_tag", "ZB")
    _assert_refused(r, "--cell_barcode_tag", tmp_path)


def test_chunk_refuses_a_user_specified_read_umi_tag(tmp_path, inputs):
    bam, gtf, genome = inputs
    r = _chunked(tmp_path, bam, gtf, genome, "--read_umi_tag", "ZM")
    _assert_refused(r, "--read_umi_tag", tmp_path)


def test_chunk_refuses_a_tag_specified_at_its_own_default_value(tmp_path, inputs):
    """Passing ``--cell_barcode_tag CB`` is refused too, and that is deliberate.

    The distinction the refusal draws is "did the user ask for single cell", not "is
    this value unusual". A run that spells the default out is asking for the encoding
    just as much as one that names a custom tag, and would get bulk output either way.
    This is also the assertion that keeps the ``argparse.SUPPRESS`` declaration honest:
    a value comparison against the config default would pass every other test in this
    file and silently accept this one.
    """
    bam, gtf, genome = inputs
    r = _chunked(tmp_path, bam, gtf, genome, "--cell_barcode_tag", "CB")
    _assert_refused(r, "--cell_barcode_tag", tmp_path)


def test_chunk_refuses_a_cell_list(tmp_path, inputs):
    bam, gtf, genome = inputs
    cells = tmp_path / "cells.txt"
    cells.write_text("BARCODE-1\n")
    r = _chunked(tmp_path, bam, gtf, genome, "--cell_list", str(cells))
    _assert_refused(r, "--cell_list", tmp_path)


def test_chunk_refuses_stream_reads(tmp_path, inputs):
    bam, gtf, genome = inputs
    r = _chunked(tmp_path, bam, gtf, genome, "--stream_reads")
    _assert_refused(r, "--stream_reads", tmp_path)


def test_chunk_names_every_offending_flag_at_once(tmp_path, inputs):
    """All four at once are listed in one message rather than one per rerun."""
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
    )
    combined = r.stdout + r.stderr
    assert r.returncode != 0
    for flag in (
        "--cell_barcode_tag",
        "--read_umi_tag",
        "--cell_list",
        "--stream_reads",
    ):
        assert flag in combined, combined[-3000:]


def test_plain_chunk_is_untouched(tmp_path, inputs):
    """No single-cell flag, so the run reaches the pipeline.

    The tag defaults are always populated in the config, so this is the case a
    value-based check would have broken -- every chunked run in existence. Asserted two
    ways: the refusal is absent, and the run reaches the chunked pipeline's first act
    (the marker ``test_chunked_entry_point.py`` uses for the same purpose) rather than
    dying earlier for some other reason and passing this test by accident.
    """
    bam, gtf, genome = inputs
    r = _chunked(tmp_path, bam, gtf, genome)
    combined = r.stdout + r.stderr
    assert REFUSAL not in combined, combined[-3000:]
    assert "counting genome-mapped reads" in combined, combined[-3000:]


def test_unchunked_single_cell_is_untouched(tmp_path, inputs):
    """The refusal is scoped to --chunk; unchunked single cell is the supported path."""
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
    combined = r.stdout + r.stderr
    assert REFUSAL not in combined, combined[-3000:]
    assert r.returncode == 0, combined[-3000:]


def test_a_non_default_tag_still_reaches_the_read_name_encoding(tmp_path):
    """The SUPPRESS declaration must not stop the flag from doing its job.

    Declaring the tag flags ``argparse.SUPPRESS`` is what makes a user-supplied tag
    detectable, and it changes how the value reaches ``LRAA_Globals.config``. If that
    hand-off broke, a single-cell run would fall back to CB/XM and emit exactly the
    bulk-shaped read names the refusal above exists to prevent -- silently, and on the
    path this item does NOT touch. So the encoding is asserted end to end, on a bam
    whose tags are deliberately not the defaults.
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

    ``WDL/subwdls/LRAA_runner.wdl`` declares ``cell_barcode_tag``/``read_umi_tag`` as
    non-optional Strings defaulting to CB/XM, so it used to append both flags
    unconditionally -- including under ``chunk = true``, which is how the WDL layer was
    silently getting bulk output for single-cell inputs and how it would now be refused
    on EVERY chunked run, single cell or not. The flags are emitted only when they
    differ from those defaults: a bulk chunked run passes nothing and is unaffected, and
    a real single-cell tag reaches LRAA and is refused there rather than in a second
    copy of the rule.
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
