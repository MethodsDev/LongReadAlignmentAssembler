#!/usr/bin/env python3

"""Each caller estimates its first-pass theta from ITS OWN reads.

The defect: cluster-guided single-cell quant hands ~29 per-cluster jobs ONE shared
splice-graph bam (`WDL/LRAA_quant_by_cluster.wdl`), and `_first_pass_assignment_bam`
returned that shared bam under `--stream_reads`. So the theta the final assignment
apportions ambiguous reads by was estimated from POOLED evidence: every cluster's
apportionment became a function of every other cluster's expression, reported as
cluster-local quant. MEASURED before this input existed: all 32 clusters reported an
identical reads_total of 94,908 while assigning 24,083 / 27,616 / 17,414 of them.

`--bam_for_priors` makes the three roles three files, which is what they always
were conceptually:

    splice graph       --bam_for_sg, shared, never reconstructed
    pass-1 theta       --bam_for_priors, THIS caller's own normalized reads
    pass-2 assignment  --bam, the full cluster bam

What these tests hold:

* pass 1 reads the PRIORS bam when there is one and the sg bam when there is not,
  asserted on the read-assignment summary's reads_total -- which is a count over
  whichever bam that pass opened, and is the number the defect made identical
  across clusters;
* the flag REACHES ChunkedRun. `LRAA --chunk` builds its chunk args from an
  explicit allowlist, and an unlisted flag is parsed and silently discarded --
  which for this one means every chunk falling back to the pooled sg slice, with
  no message and plausible output. That omission is the original defect this
  branch started from;
* the priors slice is cut with the SAME geometry as the main and sg slices, so a
  read in more than one of them lands at the same coordinate on the same mini
  contig. Divergence does not fail; it estimates theta over coordinates the reads
  are not at;
* XW survives extraction AND the in-chunk strand split. LRAA:165-192 exits on a
  --bam_for_priors whose first aligned record has no XW, because an untagged read
  weighs 1 and the estimate would silently stop being weighted;
* `--bam_for_priors` naming the same FILE as `--bam_for_sg` is refused at every
  layer that can see both, compared by device and inode so a HARD LINK is caught
  too. A byte-identical COPY is not caught and cannot be without hashing the bams,
  so the refusal NARROWS the pooled configuration rather than eliminating it, and
  the messages say so;
* the `--only_chunk` leaf rebuilds the priors paths rather than dropping them --
  the same omission the sg slice needed fixing for, invisible on the driver's
  machine because no chunk directory there would show it.

The straddle check that replaced the plan's annotation digest is exercised in
pylib/test_shared_cut_plan.py, which owns the plan fixtures.
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

REPO = Path(__file__).resolve().parents[1]
CHUNKED_RUN = REPO / "pylib" / "ChunkedRun.py"
CONTIG = "chr1"
CONTIG_LEN = 10000

sys.path.insert(0, str(REPO / "pylib"))
sys.path.insert(0, str(REPO / "util" / "misc"))

import ChunkedRun  # noqa: E402
import extract_contig_region_inputs as extractor  # noqa: E402
import Util_funcs  # noqa: E402


# --------------------------------------------------------------------- fixtures


def _header():
    return pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}],
        }
    )


def _write_bam(path, reads, xw=None):
    """``reads`` is (name, 0-based start, aligned length, is_reverse), sorted.

    ``xw`` writes that coverage-normalization weight on every record, which is
    what makes a bam admissible as ``--bam_for_sg`` or ``--bam_for_priors``: LRAA
    refuses either one whose first aligned record has none.
    """

    header = _header()
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for name, start, length, reverse in reads:
            aln = pysam.AlignedSegment(header)
            aln.query_name = name
            aln.reference_id = 0
            aln.reference_start = start
            aln.mapping_quality = 60
            aln.cigarstring = "{}M".format(length)
            aln.query_sequence = "A" * length
            aln.query_qualities = pysam.qualitystring_to_array("I" * length)
            aln.is_reverse = reverse
            if xw is not None:
                aln.set_tag("XW", xw)
            fh.write(aln)
    pysam.index(str(path))
    return Path(path)


def _reads(prefix, starts, length=100, reverse=False):
    return [
        ("{}{}".format(prefix, i), start, length, reverse)
        for i, start in enumerate(starts)
    ]


@pytest.fixture
def inputs(tmp_path):
    """Three bams with THREE DIFFERENT read counts, and that is the point.

    16 reads, 4 of evidence, 8 of priors. Equal counts would let a test pass while
    pass 1 read the wrong file, which is exactly the failure being pinned: the
    defect was invisible precisely because every cluster's numbers looked
    plausible.
    """

    genome = tmp_path / "genome.fa"
    genome.write_text(">{}\n{}\n".format(CONTIG, "A" * CONTIG_LEN))
    pysam.faidx(str(genome))

    gtf = tmp_path / "annot.gtf"
    gtf.write_text(
        "\n".join(
            '{}\ttest\texon\t{}\t{}\t.\t{}\t.\tgene_id "{}"; transcript_id "{}";'.format(
                CONTIG, lend, rend, strand, gene, gene + ".t1"
            )
            for gene, lend, rend, strand in (
                ("gplus", 1001, 1400, "+"),
                ("gminus", 3001, 3400, "-"),
            )
        )
        + "\n"
    )

    forward = _reads("fwd", [1000 + 20 * i for i in range(8)])
    reverse = _reads("rev", [3000 + 20 * i for i in range(8)], reverse=True)
    ordered = sorted(forward + reverse, key=lambda r: r[1])

    # The reads being counted must NOT be pre-thinned: LRAA rejects an XW-tagged
    # --bam outright, because a thinned library invalidates the read denominator.
    reads_bam = _write_bam(tmp_path / "reads.bam", ordered)
    # Shared evidence, thinned hardest: in the real shape this is the pool of
    # every cluster's normalized reads, and no cluster's own.
    sg_bam = _write_bam(tmp_path / "evidence.bam", ordered[::4], xw=0.25)
    # THIS caller's own normalized reads. A different subset, so a run that reads
    # the wrong one reports a different number rather than the same one.
    priors_bam = _write_bam(tmp_path / "priors.bam", ordered[::2], xw=0.5)
    return genome, gtf, reads_bam, sg_bam, priors_bam


def _count(bam_path):
    with pysam.AlignmentFile(str(bam_path), "rb") as fh:
        return sum(1 for aln in fh.fetch(until_eof=True) if not aln.is_unmapped)


def _run_chunked(tmp_path, workdir, genome, gtf, bam, *extra):
    cmd = [
        sys.executable,
        str(CHUNKED_RUN),
        "--bam",
        str(bam),
        "--genome_fa",
        str(genome),
        "--gtf",
        str(gtf),
        "--output_dir",
        str(workdir),
        "--cpu_budget",
        "2",
        "--num_total_reads",
        "16",
        # Every chunk directory self-contained, which is also what the scattered
        # workflow runs with. Without it a whole-contig chunk names the source
        # bams and no slice is materialised, so the paths under test would not
        # exist to check.
        "--no_reuse_source_bam",
    ] + [str(x) for x in extra]
    return subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(tmp_path), timeout=1800
    )


def _ok(result):
    combined = result.stdout + result.stderr
    assert result.returncode == 0, combined[-6000:]
    return combined


def _refused(result):
    combined = result.stdout + result.stderr
    assert result.returncode != 0, combined[-6000:]
    return combined


def _manifests(workdir):
    paths = sorted(Path(workdir).glob("chunks/*/chunk.partition.json"))
    assert paths, "no chunk manifests under {}".format(workdir)
    return [json.loads(p.read_text()) for p in paths]


def _chunk_logs(workdir):
    logs = sorted(Path(workdir).glob("logs/chunk_*.log"))
    assert logs, "no chunk logs under {}".format(workdir)
    return "\n".join(p.read_text() for p in logs)


def _first_record(bam_path):
    with pysam.AlignmentFile(str(bam_path), "rb") as fh:
        for aln in fh.fetch(until_eof=True):
            if not aln.is_unmapped:
                return aln
    return None


def _positions(bam_path):
    with pysam.AlignmentFile(str(bam_path), "rb") as fh:
        contigs = set(fh.references)
        return contigs, {
            aln.query_name: aln.reference_start
            for aln in fh.fetch(until_eof=True)
            if not aln.is_unmapped
        }


def _extract(tmp_path, genome, gtf, bam, region, sg_bam=None, priors_bam=None):
    return extractor.extract_partition(
        genome_fa=str(genome),
        bam=str(bam),
        region=region,
        output_prefix=str(tmp_path / "chunk"),
        gtf=str(gtf) if gtf else None,
        max_intron_length=0,
        sg_bam=str(sg_bam) if sg_bam else None,
        priors_bam=str(priors_bam) if priors_bam else None,
    )


# --------------------------------------------------- pass 1 reads which bam


def _lraa_quant(tmp_path, tag, genome, gtf, bam, sg_bam, priors_bam=None):
    """One plain LRAA quant-only run, and the reads_total its pass 1 reported."""

    workdir = tmp_path / tag
    workdir.mkdir()
    cmd = [
        sys.executable,
        str(REPO / "LRAA"),
        "--bam",
        str(bam),
        "--genome",
        str(genome),
        "--gtf",
        str(gtf),
        "--quant_only",
        "--no_chunk",
        "--output_prefix",
        str(workdir / "out"),
        "--num_total_reads",
        "16",
        "--cpu_budget",
        "2",
        "--bam_for_sg",
        str(sg_bam),
        "--no_norm",
        # Rescue would add read names the first pass never saw, which is a second
        # population on top of the one under test. Off, so reads_total is exactly
        # the pass-1 bam's reads.
        "--no_rescue_unassigned_reads_via_transcriptome_alignment",
    ]
    if priors_bam is not None:
        cmd += ["--bam_for_priors", str(priors_bam)]
    result = subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(workdir), timeout=1800
    )
    _ok(result)
    summary = workdir / "out.LRAA.quant-only.read_assignment.summary.tsv"
    assert summary.exists(), result.stdout[-4000:]
    rows = [line.split("\t") for line in summary.read_text().splitlines()]
    header, body = rows[0], rows[1:]
    total = [row for row in body if row[header.index("row_type")] == "TOTAL"]
    assert len(total) == 1, summary.read_text()
    return int(total[0][header.index("reads_total")])


def test_pass_one_reads_the_priors_bam_when_there_is_one(tmp_path, inputs):
    """THE FIX, asserted on which file pass 1 actually opened.

    reads_total is a count over the bam the first assignment pass read (see
    `_first_pass_assignment_bam`), so the two runs below differ by nothing except
    which bam that was. With the pooled sg bam it is 4; with this caller's own
    normalized reads it is 8 -- the caller's own population, which is what a
    cluster-local theta has to be estimated from.
    """

    genome, gtf, reads_bam, sg_bam, priors_bam = inputs
    assert (_count(reads_bam), _count(sg_bam), _count(priors_bam)) == (16, 4, 8)

    without = _lraa_quant(tmp_path, "without", genome, gtf, reads_bam, sg_bam)
    with_priors = _lraa_quant(
        tmp_path, "with", genome, gtf, reads_bam, sg_bam, priors_bam=priors_bam
    )

    assert without == _count(sg_bam)
    assert with_priors == _count(priors_bam)
    assert with_priors != without, "pass 1 read the same bam both times"


def test_lraa_refuses_priors_that_are_the_splice_graph_bam(tmp_path, inputs):
    """The pooled configuration, refused where both bams are first in one hand."""

    genome, gtf, reads_bam, sg_bam, _priors = inputs
    result = subprocess.run(
        [
            sys.executable,
            str(REPO / "LRAA"),
            "--bam",
            str(reads_bam),
            "--genome",
            str(genome),
            "--gtf",
            str(gtf),
            "--quant_only",
            "--no_chunk",
            "--output_prefix",
            str(tmp_path / "out"),
            "--num_total_reads",
            "16",
            "--bam_for_sg",
            str(sg_bam),
            "--bam_for_priors",
            str(sg_bam),
            "--no_norm",
        ],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
        timeout=1800,
    )
    combined = _refused(result)
    assert "are ONE file" in combined
    assert "POOLED configuration" in combined


def test_lraa_refuses_a_hard_link_to_the_splice_graph_bam(tmp_path, inputs):
    """A hard link IS the pooled bam, under a second name.

    ``os.path.realpath`` cannot see it: there is no link to resolve, so the two
    names come back distinct and the equality check passed. The configuration is the
    same one -- one file, two flags, a theta estimated over pooled evidence -- and it
    exits 0 with per-cluster numbers that look local.
    """

    genome, gtf, reads_bam, sg_bam, _priors = inputs
    linked = tmp_path / "linked_priors.bam"
    os.link(str(sg_bam), str(linked))
    os.link(str(sg_bam) + ".bai", str(linked) + ".bai")
    assert os.path.realpath(str(linked)) != os.path.realpath(str(sg_bam))

    result = subprocess.run(
        [
            sys.executable,
            str(REPO / "LRAA"),
            "--bam",
            str(reads_bam),
            "--genome",
            str(genome),
            "--gtf",
            str(gtf),
            "--quant_only",
            "--no_chunk",
            "--output_prefix",
            str(tmp_path / "out_linked"),
            "--num_total_reads",
            "16",
            "--bam_for_sg",
            str(sg_bam),
            "--bam_for_priors",
            str(linked),
            "--no_norm",
        ],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
        timeout=1800,
    )
    combined = _refused(result)
    assert "are ONE file" in combined
    assert "POOLED configuration" in combined


def test_the_refusal_does_not_claim_to_catch_a_copy(tmp_path, inputs):
    """A byte-identical copy PASSES, and the message must not imply otherwise.

    The limit of the check, asserted rather than left implicit: the pooled
    configuration stays reachable by ``cp``. Detecting it needs the contents, and
    these are bams -- the cost argument Util_funcs.file_identity_token already
    records for hashing one on startup, and which that helper could not settle here
    anyway because its token includes the resolved path, so two paths can never
    match however identical their bytes.

    So the requirement is the WORDING: the refusal that does fire must tell the
    operator what it does not cover, or a passing run reads as proof of something it
    is not.
    """

    genome, gtf, reads_bam, sg_bam, _priors = inputs
    copied = tmp_path / "copied_priors.bam"
    shutil.copyfile(str(sg_bam), str(copied))
    with open(str(copied), "rb") as a, open(str(sg_bam), "rb") as b:
        assert a.read() == b.read(), "the fixture is not actually a byte-identical copy"
    assert not Util_funcs.paths_name_one_file(str(copied), str(sg_bam))

    linked = tmp_path / "linked_for_message.bam"
    os.link(str(sg_bam), str(linked))
    os.link(str(sg_bam) + ".bai", str(linked) + ".bai")
    caught = subprocess.run(
        [
            sys.executable,
            str(REPO / "LRAA"),
            "--bam",
            str(reads_bam),
            "--genome",
            str(genome),
            "--gtf",
            str(gtf),
            "--quant_only",
            "--no_chunk",
            "--output_prefix",
            str(tmp_path / "out_msg"),
            "--num_total_reads",
            "16",
            "--bam_for_sg",
            str(sg_bam),
            "--bam_for_priors",
            str(linked),
            "--no_norm",
        ],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
        timeout=1800,
    )
    message = _refused(caught)
    assert "byte-identical COPY is not" in message
    assert "does not rule that shape out" in message
    # And the wording that promised more than the check delivers is gone.
    assert "exists to eliminate" not in message


def test_paths_name_one_file_distinguishes_a_link_from_a_copy(tmp_path):
    """The primitive all three layers share, on its own.

    Four cases, because the point is WHICH of them are one file: the same path and a
    symlink were already caught by realpath equality, a hard link was not and is the
    gap being closed, and a byte-identical copy is genuinely two files.
    """

    original = tmp_path / "a.bin"
    original.write_bytes(b"identical payload")

    linked = tmp_path / "hard.bin"
    os.link(str(original), str(linked))
    symlinked = tmp_path / "soft.bin"
    symlinked.symlink_to(original)
    copied = tmp_path / "copy.bin"
    shutil.copyfile(str(original), str(copied))
    other = tmp_path / "other.bin"
    other.write_bytes(b"different payload")

    assert Util_funcs.paths_name_one_file(str(original), str(original))
    assert Util_funcs.paths_name_one_file(str(original), str(symlinked))
    assert Util_funcs.paths_name_one_file(str(original), str(linked))
    # realpath equality, which this replaced, reports the hard link as two files.
    assert os.path.realpath(str(linked)) != os.path.realpath(str(original))

    assert not Util_funcs.paths_name_one_file(str(original), str(copied))
    assert not Util_funcs.paths_name_one_file(str(original), str(other))


def test_paths_name_one_file_survives_a_missing_path(tmp_path):
    """A guard must not raise before the code that can report the absence usefully."""

    present = tmp_path / "here.bin"
    present.write_bytes(b"x")
    absent = tmp_path / "gone.bin"

    assert not Util_funcs.paths_name_one_file(str(present), str(absent))
    assert Util_funcs.paths_name_one_file(str(absent), str(absent))


def test_lraa_refuses_an_unweighted_priors_bam(tmp_path, inputs):
    """An untagged read weighs 1, so a raw bam here silently unweights the estimate."""

    genome, gtf, reads_bam, sg_bam, _priors = inputs
    raw = _write_bam(tmp_path / "raw_priors.bam", [("only", 1000, 100, False)])
    result = subprocess.run(
        [
            sys.executable,
            str(REPO / "LRAA"),
            "--bam",
            str(reads_bam),
            "--genome",
            str(genome),
            "--gtf",
            str(gtf),
            "--quant_only",
            "--no_chunk",
            "--output_prefix",
            str(tmp_path / "out"),
            "--num_total_reads",
            "16",
            "--bam_for_sg",
            str(sg_bam),
            "--bam_for_priors",
            str(raw),
            "--no_norm",
        ],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
        timeout=1800,
    )
    combined = _refused(result)
    assert "--bam_for_priors" in combined
    assert "coverage-normalization weight" in combined


# ------------------------------------------------- extractor: geometry and tags


def test_the_priors_slice_is_rebased_exactly_like_the_other_two(tmp_path, inputs):
    """Same mini contig, same coordinate, at a NONZERO offset.

    Offset 0 would pass on a slice that never rebased at all, so the region
    deliberately starts inside the contig: every emitted record has to move by
    2,000 bp in all three files or pass 1 estimates over coordinates its own
    chunk's reads are not at.
    """

    genome, gtf, reads_bam, sg_bam, priors_bam = inputs
    manifest = _extract(
        tmp_path,
        genome,
        gtf,
        reads_bam,
        "chr1:2001-7000",
        sg_bam=sg_bam,
        priors_bam=priors_bam,
    )

    assert manifest["offset"] == 2000
    main_contigs, main_pos = _positions(manifest["files"]["bam"])
    sg_contigs, sg_pos = _positions(manifest["files"]["sg_bam"])
    priors_contigs, priors_pos = _positions(manifest["files"]["priors_bam"])

    assert (
        main_contigs == sg_contigs == priors_contigs == {manifest["mini_contig_name"]}
    )
    shared = set(main_pos) & set(priors_pos)
    assert shared, "no read is in both slices, so nothing is being compared"
    assert {name: priors_pos[name] for name in shared} == {
        name: main_pos[name] for name in shared
    }
    shared_sg = set(sg_pos) & set(priors_pos)
    assert shared_sg, "the two auxiliary slices share no read"
    assert {name: priors_pos[name] for name in shared_sg} == {
        name: sg_pos[name] for name in shared_sg
    }
    # And the rebase really happened: source starts are >= 3000, mini ones are not.
    assert max(priors_pos.values()) < 5000


def test_priors_counts_are_separate_from_the_other_accountings(tmp_path, inputs):
    """Three thinning levels, three tallies. Conflating any two hides a loss."""

    genome, gtf, reads_bam, sg_bam, priors_bam = inputs
    manifest = _extract(
        tmp_path,
        genome,
        gtf,
        reads_bam,
        "chr1:1-10000",
        sg_bam=sg_bam,
        priors_bam=priors_bam,
    )
    counts = manifest["counts"]

    assert counts["alignments_emitted"] == 16
    assert counts["sg_alignments_emitted"] == 4
    assert counts["priors_alignments_emitted"] == 8
    assert (
        counts["priors_alignments_emitted_forward"]
        + counts["priors_alignments_emitted_reverse"]
        == counts["priors_alignments_emitted"]
    )
    assert manifest["priors_bam_reused_from_source"] is False
    assert manifest["files"]["priors_bam"].endswith(".priors.bam")


def test_xw_survives_extraction_of_the_priors_slice(tmp_path, inputs):
    """`_rebase_alignment` reaches tags only through `fromstring`, so prove it."""

    genome, gtf, reads_bam, sg_bam, priors_bam = inputs
    manifest = _extract(
        tmp_path,
        genome,
        gtf,
        reads_bam,
        "chr1:1-10000",
        sg_bam=sg_bam,
        priors_bam=priors_bam,
    )

    first = _first_record(manifest["files"]["priors_bam"])
    assert first is not None
    assert first.has_tag("XW")
    assert first.get_tag("XW") == pytest.approx(0.5)
    # And not the sg bam's weight, which would mean the two slices were crossed.
    assert _first_record(manifest["files"]["sg_bam"]).get_tag("XW") == pytest.approx(
        0.25
    )


def test_no_priors_bam_leaves_the_manifest_as_it_was(tmp_path, inputs):
    """Absent, not zero: a measurement never taken is not an empty priors set."""

    genome, gtf, reads_bam, sg_bam, _priors = inputs
    manifest = _extract(
        tmp_path, genome, gtf, reads_bam, "chr1:1-10000", sg_bam=sg_bam
    )

    assert manifest["files"]["priors_bam"] is None
    assert manifest["priors_bam_reused_from_source"] is False
    assert "priors_dropped_read_names" not in manifest
    assert not [key for key in manifest["counts"] if key.startswith("priors_")]
    # The sg accounting is untouched by the new role existing beside it.
    assert manifest["counts"]["sg_alignments_emitted"] == 4


def test_extractor_refuses_priors_that_are_the_sg_bam(tmp_path, inputs):
    """One file for both roles IS the pooled prior, spelled with two flags."""

    genome, gtf, reads_bam, sg_bam, _priors = inputs
    with pytest.raises(extractor.ExtractionError) as excinfo:
        _extract(
            tmp_path,
            genome,
            gtf,
            reads_bam,
            "chr1:1-10000",
            sg_bam=sg_bam,
            priors_bam=sg_bam,
        )
    assert "POOLED" in str(excinfo.value)


def test_extractor_refuses_a_hard_link_to_the_sg_bam(tmp_path, inputs):
    """Same gap, one layer down: the extractor compared realpaths too."""

    genome, gtf, reads_bam, sg_bam, _priors = inputs
    linked = tmp_path / "linked_sg.bam"
    os.link(str(sg_bam), str(linked))
    os.link(str(sg_bam) + ".bai", str(linked) + ".bai")

    with pytest.raises(extractor.ExtractionError) as excinfo:
        _extract(
            tmp_path,
            genome,
            gtf,
            reads_bam,
            "chr1:1-10000",
            sg_bam=sg_bam,
            priors_bam=linked,
        )
    message = str(excinfo.value)
    assert "are ONE file" in message
    assert "byte-identical COPY is not" in message


def test_extractor_refuses_priors_that_are_the_read_bam(tmp_path, inputs):
    """Pass 1 over the population pass 2 assigns is the no-priors configuration."""

    genome, gtf, reads_bam, sg_bam, _priors = inputs
    with pytest.raises(extractor.ExtractionError) as excinfo:
        _extract(
            tmp_path,
            genome,
            gtf,
            reads_bam,
            "chr1:1-10000",
            sg_bam=sg_bam,
            priors_bam=reads_bam,
        )
    assert "resolve to" in str(excinfo.value)


# ------------------------------------------------------------ argument refusals


def test_chunked_run_refuses_priors_that_are_the_sg_bam(tmp_path, inputs):
    """Refused before a single region is read, and named as what it is."""

    genome, gtf, reads_bam, sg_bam, _priors = inputs
    combined = _refused(
        _run_chunked(
            tmp_path,
            tmp_path / "work",
            genome,
            gtf,
            reads_bam,
            "--bam_for_sg",
            sg_bam,
            "--bam_for_priors",
            sg_bam,
            "--bam_for_cut_selection",
            reads_bam,
        )
    )
    assert "--bam_for_sg and --bam_for_priors are ONE file" in combined
    assert "POOLED configuration" in combined


def test_chunked_run_refuses_a_hard_link_to_the_sg_bam(tmp_path, inputs):
    """And at the driver's own pre-flight, before a single region is read."""

    genome, gtf, reads_bam, sg_bam, _priors = inputs
    linked = tmp_path / "linked_sg_chunked.bam"
    os.link(str(sg_bam), str(linked))
    os.link(str(sg_bam) + ".bai", str(linked) + ".bai")

    combined = _refused(
        _run_chunked(
            tmp_path,
            tmp_path / "work_linked",
            genome,
            gtf,
            reads_bam,
            "--bam_for_sg",
            sg_bam,
            "--bam_for_priors",
            linked,
            "--bam_for_cut_selection",
            reads_bam,
        )
    )
    assert "--bam_for_sg and --bam_for_priors are ONE file" in combined
    assert "byte-identical COPY is not" in combined


def test_chunked_run_refuses_priors_with_strand_first_chunking(tmp_path, inputs):
    """Stage 3b, which splits the priors bam, does not run in that mode."""

    genome, gtf, reads_bam, _sg, priors_bam = inputs
    combined = _refused(
        _run_chunked(
            tmp_path,
            tmp_path / "work",
            genome,
            gtf,
            reads_bam,
            "--chunk_by_strand",
            "--bam_for_priors",
            priors_bam,
        )
    )
    assert "--bam_for_priors is strandless-chunking only" in combined


def test_lraa_chunk_forwards_bam_for_priors_instead_of_ignoring_it(tmp_path, inputs):
    """`LRAA --chunk` builds its chunk args from an explicit allowlist.

    An unlisted flag is parsed and silently discarded, which for this one means
    every chunk's pass 1 falling back to the shared sg slice -- the pooled prior,
    restored without a message. Asserted by the run reaching a refusal only
    ChunkedRun produces: `LRAA`'s own validation of this flag lives past the
    chunked dispatch (`main` exits at the `_run_chunked_mode` call), so a message
    naming a missing --bam_for_priors can only have come from the pipeline. If the
    flag were dropped, the run would simply proceed.
    """

    genome, gtf, reads_bam, _sg, _priors = inputs
    missing = tmp_path / "absent_priors.bam"
    result = subprocess.run(
        [
            sys.executable,
            str(REPO / "LRAA"),
            "--bam",
            str(reads_bam),
            "--genome",
            str(genome),
            "--gtf",
            str(gtf),
            "--quant_only",
            "--output_prefix",
            str(tmp_path / "out"),
            "--chunk",
            "--chunk_work_dir",
            str(tmp_path / "work"),
            "--bam_for_priors",
            str(missing),
        ],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
        timeout=1800,
    )
    combined = _refused(result)
    assert "--bam_for_priors {} does not exist".format(missing) in combined


# ------------------------------------------------------------- end to end runs


@pytest.fixture
def priors_run(tmp_path, inputs):
    """All three roles at once, which is the configuration the WDL runs."""

    genome, gtf, reads_bam, sg_bam, priors_bam = inputs
    workdir = tmp_path / "work_priors"
    _ok(
        _run_chunked(
            tmp_path,
            workdir,
            genome,
            gtf,
            reads_bam,
            "--bam_for_sg",
            sg_bam,
            "--bam_for_priors",
            priors_bam,
            "--bam_for_cut_selection",
            reads_bam,
        )
    )
    return workdir


@pytest.fixture
def priors_only_run(tmp_path, inputs):
    """Priors with NO caller-supplied evidence. A distinct path, not a subset.

    Stage 4 is skipped for a unit carrying an sg slice and only for that reason;
    a unit carrying priors alone still has to normalize its own reads, because it
    still needs a splice graph and its priors are not one.
    """

    genome, gtf, reads_bam, _sg, priors_bam = inputs
    workdir = tmp_path / "work_priors_only"
    _ok(
        _run_chunked(
            tmp_path,
            workdir,
            genome,
            gtf,
            reads_bam,
            "--bam_for_priors",
            priors_bam,
        )
    )
    return workdir


def test_priors_without_evidence_still_normalizes_per_unit(priors_only_run):
    """Stage 4 RUNS, stage 5 gets the normalized bam as evidence and the slice as
    priors, and no sg artifact is invented for a run that supplied none."""

    argv = _chunk_logs(priors_only_run)
    assert "normalize_bam_by_strand.py" in argv
    assert "--bam_for_priors" in argv
    assert "chunk.priors.strand." in argv
    assert "chunk.sg.strand" not in argv
    assert list(Path(priors_only_run).glob("chunks/*/chunk.*.norm.bam"))

    for manifest in _manifests(priors_only_run):
        assert manifest["files"]["sg_bam"] is None
        assert manifest["files"]["priors_bam"].endswith(".priors.bam")
        assert not [key for key in manifest["counts"] if key.startswith("sg_")]
        assert "priors_dropped_read_names" in manifest


def _unit_reads_total(workdir):
    """Each quant unit's pass-1 read count, from the summary stage 5 wrote."""

    totals = {}
    for path in sorted(
        Path(workdir).glob("chunks/*/*read_assignment.summary.tsv")
    ):
        rows = [line.split("\t") for line in path.read_text().splitlines()]
        header, body = rows[0], rows[1:]
        total = [r for r in body if r[header.index("row_type")] == "TOTAL"]
        assert len(total) == 1, path.read_text()
        totals[path.name.split(".")[0]] = int(total[0][header.index("reads_total")])
    assert totals, "no per-unit read-assignment summary under {}".format(workdir)
    return totals


def test_a_chunked_units_pass_one_reads_its_own_priors_slice(tmp_path, inputs):
    """THE WHOLE CHAIN, at the level the WDL runs it.

    Two chunked runs differing only in --bam_for_priors. Every unit's reads_total
    is a count over the bam ITS pass 1 opened, so this proves the slice survives
    extraction, the strand split, the stage-5 argv and LRAA's own resolution --
    any break in that chain silently restores the pooled prior, which is what the
    defect looked like in production.
    """

    genome, gtf, reads_bam, sg_bam, priors_bam = inputs
    shared = (
        "--bam_for_sg",
        sg_bam,
        "--bam_for_cut_selection",
        reads_bam,
    )
    _ok(
        _run_chunked(
            tmp_path, tmp_path / "w_off", genome, gtf, reads_bam, *shared
        )
    )
    _ok(
        _run_chunked(
            tmp_path,
            tmp_path / "w_on",
            genome,
            gtf,
            reads_bam,
            *shared,
            "--bam_for_priors",
            priors_bam,
        )
    )

    without = _unit_reads_total(tmp_path / "w_off")
    with_priors = _unit_reads_total(tmp_path / "w_on")

    assert set(without) == set(with_priors)
    # Each orientation holds half of its bam's reads, and the two bams differ 4
    # against 8 -- so the sum over units IS the pass-1 bam's read count.
    assert sum(without.values()) == _count(sg_bam)
    assert sum(with_priors.values()) == _count(priors_bam)
    for unit in without:
        assert with_priors[unit] > without[unit], unit


def test_every_unit_priors_slice_matches_its_own_reads(priors_run):
    """Same mini contig, same coordinates, orientation for orientation.

    The per-unit files, not the chunk-level ones: this is what stage 5 opens, and
    they are produced by a SECOND tool (the strand split) after the extractor, so
    the geometry has two chances to drift before a unit sees it.
    """

    chunks = sorted(Path(priors_run).glob("chunks/*"))
    assert chunks
    for cdir in chunks:
        manifest = json.loads((cdir / "chunk.partition.json").read_text())
        for strand in ("+", "-"):
            main = cdir / "chunk.strand.{}.bam".format(strand)
            priors = cdir / "chunk.priors.strand.{}.bam".format(strand)
            assert priors.exists(), "no priors slice for {} of {}".format(
                strand, cdir.name
            )
            main_contigs, main_pos = _positions(main)
            priors_contigs, priors_pos = _positions(priors)
            assert priors_contigs == main_contigs == {manifest["mini_contig_name"]}
            shared = set(main_pos) & set(priors_pos)
            if not shared:
                # A chunk-orientation this caller's priors do not cover is legal;
                # a MISPLACED read is not, which is what the equality below holds.
                continue
            assert {n: priors_pos[n] for n in shared} == {
                n: main_pos[n] for n in shared
            }


def test_xw_survives_the_strand_split_of_the_priors_slice(priors_run):
    """The split is a separate tool from the extractor and can lose tags on its own.

    LRAA:165-192 exits on a --bam_for_priors whose first aligned record has no XW,
    so this is what stands between a correct run and one that dies in stage 5 --
    and the run below completed, which is half the proof. The weight is checked by
    VALUE too: 0.5 is the priors bam's and 0.25 the sg bam's, so a crossed pair
    would still carry a tag.
    """

    seen = 0
    for priors in sorted(Path(priors_run).glob("chunks/*/chunk.priors.strand.*.bam")):
        first = _first_record(priors)
        if first is None:
            continue
        seen += 1
        assert first.has_tag("XW"), "{} lost XW in the strand split".format(priors)
        assert first.get_tag("XW") == pytest.approx(0.5)
    assert seen, "no non-empty priors slice, so the tag check ran on nothing"


def test_stage_5_is_handed_the_priors_slice_and_stage_4_still_skipped(priors_run):
    """The forwarding, and the skip it must not have changed.

    Stage 4 is skipped because the SG evidence arrived pre-normalized; the priors
    slice arrives pre-normalized too and decides nothing about that skip. Both
    stated here so a future change to one is not mistaken for the other.
    """

    argv = _chunk_logs(priors_run)
    assert "normalize_bam_by_strand.py" not in argv
    for strand in ("+", "-"):
        assert "chunk.priors.strand.{}.bam".format(strand) in argv
        assert "chunk.sg.strand.{}.bam".format(strand) in argv
    assert "--bam_for_priors" in argv
    assert "--bam_for_sg" in argv
    assert "--no_norm" in argv
    assert not list(Path(priors_run).glob("chunks/*/chunk.*.norm.bam"))

    for manifest in _manifests(priors_run):
        assert manifest["files"]["priors_bam"].endswith(".priors.bam")
        assert "priors_dropped_read_names" in manifest
        assert manifest["counts"]["priors_alignments_emitted"] >= 0
        # The split was verified against the PRIORS tally, not the read tally.
        assert manifest["counts"]["priors_alignments_emitted"] <= (
            manifest["counts"]["alignments_emitted"]
        )


def test_the_priors_identity_is_in_the_extraction_sentinel(tmp_path, inputs):
    """A resumed run must not serve chunks extracted against other priors.

    The slice is a FILE IN the chunk directory, so a directory extracted without
    it -- or against a different priors bam -- is not the directory this run
    wants. Asserted on the token rather than through an output, because a wrong
    reuse produces no error at all.
    """

    genome, _gtf, reads_bam, sg_bam, priors_bam = inputs
    common = dict(
        bam=str(reads_bam),
        genome_fa=str(genome),
        output_dir=str(tmp_path / "work"),
        margin=200,
        max_intron_length=100000,
        gtf=None,
        discovery=True,
        strandless_chunks=True,
    )
    planned = {
        "key": "",
        "chunk_id": "chr1_00",
        "region": "chr1:1-10000",
        "spans_whole_contig": True,
        "bam": str(reads_bam),
        "index": 0,
        "order": 0,
        "chrom": CONTIG,
    }

    def token(**extra):
        args = ChunkedRun.default_args(**dict(common, **extra))
        plan = ChunkedRun.extraction_plan(
            args, str(tmp_path / "work"), str(tmp_path / "work" / "chunks"), planned, ""
        )
        return plan["token"], plan["cmd"]

    bare, bare_cmd = token()
    with_priors, priors_cmd = token(bam_for_priors=str(priors_bam))
    with_sg, _sg_cmd = token(bam_for_sg=str(sg_bam))

    assert len({bare, with_priors, with_sg}) == 3
    assert "--priors_bam" not in bare_cmd
    assert priors_cmd[priors_cmd.index("--priors_bam") + 1] == str(
        priors_bam.resolve()
    )


# ----------------------------------------------------- the scattered leaf's view


def _leaf_plan():
    return {
        "version": ChunkedRun.CHUNK_PLAN_VERSION,
        "lraa_suffix": "LRAA.quant-only",
        "num_total_reads": 10,
        "discovery": False,
        "chunks": [
            {
                "chunk_id": "chr1_00",
                "chrom": CONTIG,
                "strand": "",
                "strandless": True,
                "region": "chr1:1-10000",
                "index": 0,
                "order": 0,
                "has_gtf": False,
            }
        ],
    }


def _leaf_manifest(cdir, files):
    (cdir / "chunk.partition.json").write_text(
        json.dumps(
            {
                "offset": 0,
                "window_origin": 0,
                "counts": {"alignments_emitted": 0},
                "files": files,
            }
        )
    )


def test_only_chunk_rebuilds_the_priors_paths_it_was_given(tmp_path):
    """The leaf re-derives every path from ITS OWN outdir, dropping the manifest's.

    A rebuild that omits the priors key does not fail: the leaf falls back to the
    shared sg slice and estimates theta from the pool -- the original defect,
    reintroduced on the scattered path only, where no chunk directory on the
    driver's machine would show it.
    """

    outdir = tmp_path / "leaf"
    cdir = outdir / "chunks" / "chr1_00"
    cdir.mkdir(parents=True)
    prefix = cdir / "chunk"
    _leaf_manifest(
        cdir,
        {
            "fasta": "/elsewhere/chunk.fa",
            "bam": "/elsewhere/chunk.bam",
            "gtf": None,
            "dropped_reads": "/elsewhere/chunk.dropped_reads.txt",
            "sg_bam": "/elsewhere/chunk.sg.bam",
            "priors_bam": "/elsewhere/chunk.priors.bam",
        },
    )

    chunk = ChunkedRun.rebuild_chunk_record(_leaf_plan(), "chr1_00", str(outdir))

    assert chunk["manifest"]["files"]["priors_bam"] == "{}.priors.bam".format(prefix)
    assert [u["priors_bam"] for u in chunk["units"]] == [
        "{}.priors.strand.+.bam".format(prefix),
        "{}.priors.strand.-.bam".format(prefix),
    ]
    # The sg slice is still rebuilt beside it: the two roles are independent.
    assert chunk["manifest"]["files"]["sg_bam"] == "{}.sg.bam".format(prefix)


def test_only_chunk_without_priors_keeps_todays_resolution(tmp_path):
    """The same rebuild must not invent priors for a chunk that has none."""

    outdir = tmp_path / "leaf2"
    cdir = outdir / "chunks" / "chr1_00"
    cdir.mkdir(parents=True)
    _leaf_manifest(
        cdir,
        {
            "fasta": "/elsewhere/chunk.fa",
            "bam": "/elsewhere/chunk.bam",
            "gtf": None,
            "dropped_reads": "/elsewhere/chunk.dropped_reads.txt",
            "sg_bam": "/elsewhere/chunk.sg.bam",
        },
    )

    chunk = ChunkedRun.rebuild_chunk_record(_leaf_plan(), "chr1_00", str(outdir))

    assert chunk["manifest"]["files"]["priors_bam"] is None
    assert [u["priors_bam"] for u in chunk["units"]] == [None, None]


def test_only_chunk_refuses_a_priors_slice_reused_from_the_source(tmp_path):
    """A reused source is not this chunk's slice, and the leaf cannot open it.

    The same refusal the sg slice gets, for the same reason: under reuse the
    manifest names the WHOLE priors bam, which on a scattered leaf is not present
    and, if it were, holds every other chunk's reads.
    """

    outdir = tmp_path / "leaf3"
    cdir = outdir / "chunks" / "chr1_00"
    cdir.mkdir(parents=True)
    (cdir / "chunk.partition.json").write_text(
        json.dumps(
            {
                "offset": 0,
                "window_origin": 0,
                "counts": {"alignments_emitted": 0},
                "priors_bam_reused_from_source": True,
                "files": {
                    "fasta": "/elsewhere/chunk.fa",
                    "bam": "/elsewhere/chunk.bam",
                    "gtf": None,
                    "dropped_reads": "/elsewhere/chunk.dropped_reads.txt",
                    "priors_bam": "/elsewhere/whole_priors.bam",
                },
            }
        )
    )

    with pytest.raises(ChunkedRun.PipelineError) as excinfo:
        ChunkedRun.rebuild_chunk_record(_leaf_plan(), "chr1_00", str(outdir))
    assert "priors-source reuse" in str(excinfo.value)
