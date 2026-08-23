#!/usr/bin/env python3

"""XW weighting is unconditional now, so the INPUTS carry the rules.

Weighting used to be a mode. A flag turned it on, and a block of guards refused every
combination whose output it could not make correct -- single-cell matrices, --tag_bam,
discovery, an undeclared --library_type. Every one of those refusals rested on the same
fact: a thinned --bam makes the tracking file enumerate fewer reads than the library it
claims to describe, and the consumers that sum tracking rows (
util/sc/singlecell_tracking_to_sparse_matrix.py, util/annotate_bam_with_read_tracking_info.py
) never read XW, so they undercount, unevenly, worst at exactly the deep loci
normalization touched.

The flag is gone and weighting always applies, so those combinations stopped being
special. What is refused instead is the input that made them wrong to begin with. Each
bam has exactly one role and one invariant on it, checked at LRAA:2264 and LRAA:2284:

  --bam         the full library. Reported counts are scaled by it -- --num_total_reads
                counts its records, and the streaming pass sums fractional assignments
                over them -- so it must NOT carry XW. A thinned bam there turns every
                absolute number into a per-retained-read quantity that looks exactly
                like a full-library one.
  --bam_for_sg  the splice-graph evidence, taken verbatim and never re-normalized. The
                splice graph divides each read's acceptance probability back out through
                XW, and an untagged read weighs 1, so it MUST carry XW -- otherwise a
                thinned bam is read as though its survivors were the whole library.

Stated as invariants on the inputs rather than inferred from flags, a weight is present
exactly where thinning happened, which is what lets weighting apply without being asked
for. Both invariants accept an empty bam: there is no support for a missing (or spurious)
weight to distort, and the chunked pipeline produces empty per-orientation bams routinely.

The rest of this module covers the --stream_reads two-pass gates, which are about which
bam each pass reads rather than about weighting, and which share the same setup path.
"""

import os
import subprocess
import sys

import pysam
import pytest

LRAA = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "LRAA")
CONTIG = "chr1"


def _bam(path, *, single_cell, thinning_weight=None):
    """A tiny bam. `thinning_weight` stamps XW on every record, as the normalizer does,
    which is what makes a bam usable as a verbatim --bam_for_sg -- and what makes it
    unusable as --bam.
    """
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
            if thinning_weight is not None:
                a.set_tag("XW", float(thinning_weight), value_type="f")
            fh.write(a)
    pysam.index(str(path))
    return path


def _run(bam, *extra, gtf=None, genome=None, tmp=None):
    cmd = [sys.executable, LRAA, "--quant_only", "--bam", str(bam),
           "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp / "out"),
           # The input-role invariants live in the UNCHUNKED setup path (LRAA:2264,
           # LRAA:2284) -- --chunk now defaults on and dispatches at LRAA:2145, before
           # either of them runs. --no_stream_reads is likewise the base default; a
           # caller that wants streaming passes --stream_reads in `extra`, which
           # (sharing dest with --no_stream_reads) wins by argparse's last-flag-on-the-
           # line rule since extra is appended after this base list.
           "--no_chunk", "--no_stream_reads"] + list(extra)
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


def test_a_pre_thinned_bam_is_refused(tmp_path, inputs):
    """--bam is the population everything is reported against, so it cannot be thinned.

    Reported counts are scaled by --bam: --num_total_reads counts its records, and the
    streaming pass sums fractional assignments over them. Hand it a thinned bam and every
    absolute number becomes a per-retained-read quantity while looking exactly like a
    full-library one, and the tracking file covers a subset of the library while claiming
    to be complete. Weights cannot rescue that -- XW corrects the splice graph and the
    isoform split, not the population being counted -- which is why this is fatal rather
    than a warning, and why thinned input belongs in --bam_for_sg instead.
    """
    gtf, genome = inputs
    thinned = _bam(tmp_path / "thinned.bam", single_cell=False, thinning_weight=4.0)
    r = _run(thinned, "--no_rescue_unassigned_reads_via_transcriptome_alignment",
             gtf=gtf, genome=genome, tmp=tmp_path)
    assert r.returncode != 0
    assert "carries an XW coverage-normalization weight" in (r.stdout + r.stderr)
    assert not list(tmp_path.glob("out*quant.expr")), "and must not emit results"


def test_an_untagged_bam_is_accepted_and_quantified(tmp_path, inputs):
    """The counterweight to the refusal above: an ordinary bam must still get through.

    Without this, a --bam check that rejected every input -- an inverted has_tag, a raise
    on the read itself -- would satisfy the refusal test while breaking every real run.
    Asserts results exist, not merely a zero exit.
    """
    gtf, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    r = _run(bam, "--no_rescue_unassigned_reads_via_transcriptome_alignment",
             gtf=gtf, genome=genome, tmp=tmp_path)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "carries an XW coverage-normalization weight" not in (r.stdout + r.stderr)
    assert list(tmp_path.glob("out*quant.expr")), "must emit results"


def test_stock_discovery_under_stream_reads_writes_its_outputs(tmp_path, inputs):
    """The configuration a stock discovery run resolves to: streaming on, weighting on,
    nothing declared -- the command a plain `LRAA --bam --genome` becomes.

    Discovery used to be refused outright under weighting, on the theory that its
    absolute gates counted retained reads while EM used weighted support. They do not
    disagree: the pre-filter quant that discovery gates on is the one pass that still
    asks for weight_reads=False explicitly (LRAA:6212, under "Do an initial quant"), so
    every gate sees what it always saw and the run must traverse discovery to its output
    stage rather than dying at setup.

    What is asserted is that both output files were WRITTEN, which is the part the old
    refusal prevented. Their isoform content is deliberately not asserted: measured on
    this fixture (10 single-exon 100M reads over one pileup) discovery selects nothing,
    so out.LRAA.ref-free.gtf and .quant.expr arrive carrying their version/CMD header
    and column header only. Getting an isoform out of a synthetic pileup is a
    splice-graph question, tested where the splice graph is; here the files' existence is
    the evidence that setup let the run through.
    """
    _, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    cmd = [sys.executable, LRAA, "--bam", str(bam), "--genome", str(genome),
           "--output_prefix", str(tmp_path / "out"),
           "--no_chunk", "--stream_reads"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    out = r.stdout + r.stderr
    assert r.returncode == 0, out
    assert list(tmp_path.glob("out*.gtf")), "must emit a GTF"
    assert list(tmp_path.glob("out*quant.expr")), "and must reach the final quant"


def test_stream_reads_rescues_by_default_when_rescue_is_enabled(tmp_path, inputs):
    """--stream_reads alone, with rescue left at its own default (on), must SUCCEED.

    Both --stream_reads and --stream_reads_rescue_unassigned now default to True, and
    the latter's unset default tracks whether transcriptome rescue itself is enabled --
    so passing --stream_reads alone, with rescue untouched, resolves the streaming pass
    to rescue reads it cannot reproduce from a batch first pass automatically, rather
    than refusing to run. This used to be the refusal case (rescue on, streaming on,
    the streaming-rescue flag not given); it no longer is, because "not given" no
    longer means "off".
    """
    gtf, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    cmd = [sys.executable, LRAA, "--quant_only", "--bam", str(bam),
           "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp_path / "out"),
           "--stream_reads", "--no_chunk"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode == 0, r.stdout + r.stderr
    assert "cannot reproduce transcriptome rescue" not in (r.stdout + r.stderr)
    assert list(tmp_path.glob("out*quant.expr")), "must emit results"


def test_stream_reads_still_refuses_when_streaming_rescue_is_explicitly_declined(
    tmp_path, inputs
):
    """The refusal survives for the one combination it still describes: a caller who
    explicitly declines rescue INSIDE the stream while leaving transcriptome rescue
    enabled everywhere else. Unlike the default (unstated) case above, an explicit
    --no_stream_reads_rescue_unassigned is not overridden by rescue's own default --
    it is a deliberate request the guard still has to honour.
    """
    gtf, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    cmd = [sys.executable, LRAA, "--quant_only", "--bam", str(bam),
           "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp_path / "out"),
           "--stream_reads", "--no_stream_reads_rescue_unassigned", "--no_chunk"]
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
           "--normalize_max_cov_level", "0",
           "--no_chunk"]
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

    Carries XW, because that is what "pre-normalized" means and what a verbatim
    --bam_for_sg is required to prove. Asserts the run SUCCEEDS and emits results,
    not merely that two specific messages are absent: absence alone would also hold if
    the tagged path died for some unrelated reason, which is exactly how the earlier
    version of this test kept passing while supplying an untagged bam.
    """
    gtf, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    sg_bam = _bam(tmp_path / "thinned.bam", single_cell=False, thinning_weight=4.0)
    cmd = [sys.executable, LRAA, "--quant_only", "--bam", str(bam),
           "--bam_for_sg", str(sg_bam), "--no_norm",
           "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp_path / "out"),
           "--stream_reads",
           "--no_rescue_unassigned_reads_via_transcriptome_alignment",
           "--no_chunk"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    out = r.stdout + r.stderr
    assert "needs a first-pass bam thinner than the one it streams" not in out
    assert "carries no XW coverage-normalization weight" not in out
    assert r.returncode == 0, out
    assert list(tmp_path.glob("out*quant.expr")), "must emit results"


def test_a_verbatim_bam_for_sg_without_weights_is_refused(tmp_path, inputs):
    """The one route by which upstream thinning enters, so the claim has to be checked.

    Taken verbatim as splice-graph evidence, --bam_for_sg is the caller asserting the bam
    is already thinned. The splice graph honours XW unconditionally and an untagged read
    weighs 1, so an unweighted thinned bam is read as though its survivors were the whole
    library -- understating support by the thinning factor at exactly the deep loci
    thinning touched, with nothing downstream able to detect it. Hence fatal, not a
    warning.
    """
    gtf, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    sg_bam = _bam(tmp_path / "untagged.bam", single_cell=False)
    cmd = [sys.executable, LRAA, "--quant_only", "--bam", str(bam),
           "--bam_for_sg", str(sg_bam), "--no_norm",
           "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp_path / "out"),
           "--no_rescue_unassigned_reads_via_transcriptome_alignment",
           "--no_chunk"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode != 0
    assert "carries no XW coverage-normalization weight" in (r.stdout + r.stderr)


def test_a_bam_for_sg_must_carry_weights_even_with_normalization_enabled(tmp_path, inputs):
    """The requirement is unconditional -- it does not depend on --no_norm.

    It used to: with normalization enabled, the run re-normalized the supplied
    --bam_for_sg itself, which is what wrote the tag, so demanding it up front would have
    refused an ordinary invocation. That branch is deleted. --bam_for_sg is now always
    taken verbatim, because thinning an already-thinned bam composes acceptance rates for
    no gain and every artifact it could produce is one the caller could have produced
    directly. Supplying the flag IS the decision about what the splice graph reads --
    so the bam has to arrive already carrying its weights, whatever the normalization
    settings say.

    Omitting --no_norm here is the whole point of the test: the same untagged bam that
    test_a_verbatim_bam_for_sg_without_weights_is_refused supplies must be refused with
    normalization left on.
    """
    gtf, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    sg_bam = _bam(tmp_path / "untagged.bam", single_cell=False)
    cmd = [sys.executable, LRAA, "--quant_only", "--bam", str(bam),
           "--bam_for_sg", str(sg_bam),
           "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp_path / "out"),
           "--no_rescue_unassigned_reads_via_transcriptome_alignment",
           "--no_chunk", "--no_stream_reads"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode != 0
    assert "carries no XW coverage-normalization weight" in (r.stdout + r.stderr)


def test_stream_reads_refuses_when_bam_for_sg_is_the_streamed_bam(tmp_path, inputs):
    """The same bam under both flags is the case the gate exists for, spelled differently."""
    gtf, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    cmd = [sys.executable, LRAA, "--quant_only", "--bam", str(bam),
           "--bam_for_sg", str(bam), "--no_norm",
           "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp_path / "out"),
           "--stream_reads",
           "--no_rescue_unassigned_reads_via_transcriptome_alignment",
           "--no_chunk"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode != 0
    assert "needs a first-pass bam thinner than the one it streams" in (r.stdout + r.stderr)


def test_stream_reads_permits_tag_bam(tmp_path, inputs):
    """--tag_bam is unrestricted, and the merged tracking file under --stream_reads covers
    the full bam, so this combination is the more complete one rather than a refused one.

    Assert real tag content on a known read, not just that a tagged bam exists -- a
    silently mismatched name encoding would produce a real but empty-of-matches tagged
    bam that a returncode/existence-only check would miss.
    """
    gtf, genome = inputs
    bam = _bam(tmp_path / "bulk.bam", single_cell=False)
    cmd = [sys.executable, LRAA, "--quant_only", "--bam", str(bam),
           "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp_path / "out"),
           "--stream_reads",
           "--no_rescue_unassigned_reads_via_transcriptome_alignment", "--tag_bam",
           "--no_chunk"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode == 0, r.stdout + r.stderr
    tagged = tmp_path / "bulk.bam.tagged.bam"
    assert tagged.exists(), "expected a *.tagged.bam to be produced"
    with pysam.AlignmentFile(str(tagged), "rb") as fh:
        reads_by_name = {a.query_name: a for a in fh.fetch(until_eof=True)}
    assert "r0" in reads_by_name
    r0 = reads_by_name["r0"]
    assert r0.has_tag("XG") and r0.get_tag("XG") == "g1"
    assert r0.has_tag("XI") and r0.get_tag("XI") == "t1"


def test_a_stock_invocation_of_single_cell_input_with_tag_bam_succeeds(tmp_path, inputs):
    """Nothing about the combinations weighting used to refuse is refused any more.

    Single-cell input, --tag_bam and no --library_type -- all three of the old refusals
    at once, all three defaults -- must reach results. Weighting is on here without being
    named, since the whole flag is gone, and the run is left entirely on its defaults
    (chunking and streaming both on) so that this covers the command a user actually
    types rather than the unchunked setup path the invariants live in.
    """
    gtf, genome = inputs
    bam = _bam(tmp_path / "sc.bam", single_cell=True)
    cmd = [sys.executable, LRAA, "--quant_only", "--bam", str(bam),
           "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp_path / "out"), "--tag_bam"]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode == 0, r.stdout + r.stderr
    assert list(tmp_path.glob("out*quant.expr")), "must emit results"
