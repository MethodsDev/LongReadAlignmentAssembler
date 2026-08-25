#!/usr/bin/env python3

"""A caller-supplied splice-graph bam survives chunking.

The defect: chunked mode ignored `--bam_for_sg` outright and manufactured its own
per-unit evidence by normalizing each unit's own reads. For the single-cell shape
that input exists for -- WDL/LRAA_quant_by_cluster.wdl hands ONE merged
normalized bam to every per-cluster quant job while each job's `--bam` is its own
cluster's raw reads -- that silently replaced one shared splice graph with ~29
different ones, and nothing in any output said so.

What these tests hold:

* the sg slice is cut with the SAME geometry as the main slice, so a read in both
  files lands at the same coordinate on the same mini contig. Divergence here does
  not fail; it builds a graph out of coordinates the reads are not at;
* XW survives extraction AND the in-chunk strand split. `LRAA:165-192` exits on a
  `--bam_for_sg` whose first aligned record has no XW, so a slice that lost its
  tags stops the run -- and `_rebase_alignment` reaches tags only indirectly, via
  `AlignedSegment.fromstring` on the whole record string. That is proven, not
  assumed;
* stage 4 does not run for a unit that has one, and stage 5 is handed it. The
  evidence arrived pre-normalized; normalizing again composes two acceptance
  rates;
* a run without the flag is byte-for-byte the old behaviour;
* the overhang drop applies to both slices identically;
* cut selection is driven by an explicit shared SOURCE, not by `--bam`, and
  refuses rather than guessing when one is not supplied. Cut placement reads the
  bam, so per-caller reads mean per-caller geometry -- and the thinned sg bam
  cannot be that source either, because a cut it shows as unspanned can still be
  spanned by a raw read the extractor then drops with nobody having named it.
"""

import json
import os
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
    what makes a bam admissible as ``--bam_for_sg``: LRAA refuses one whose first
    aligned record has none.
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
    """Reads on both orientations so the in-chunk strand split has work to do."""

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
    # The evidence is thinned, and says so. A subset of the reads by construction,
    # which is what a normalized bam is.
    sg_bam = _write_bam(tmp_path / "evidence.bam", ordered[::2], xw=0.5)
    return genome, gtf, reads_bam, sg_bam


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
        # bams and the sg slice is never materialised, so the paths under test
        # would not exist to check.
        "--no_reuse_source_bam",
    ] + [str(x) for x in extra]
    return subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(tmp_path), timeout=1800
    )


def _ok(result):
    combined = result.stdout + result.stderr
    assert result.returncode == 0, combined[-6000:]
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


# ------------------------------------------------- extractor: geometry and tags


def _extract(tmp_path, genome, gtf, bam, region, sg_bam=None, prefix="chunk"):
    return extractor.extract_partition(
        genome_fa=str(genome),
        bam=str(bam),
        region=region,
        output_prefix=str(tmp_path / prefix),
        gtf=str(gtf) if gtf else None,
        max_intron_length=0,
        sg_bam=str(sg_bam) if sg_bam else None,
    )


def test_sg_slice_is_rebased_exactly_like_the_main_slice(tmp_path, inputs):
    """Same mini contig, same coordinate, at a NONZERO offset.

    Offset 0 would pass on a slice that never rebased at all, so the region
    deliberately starts inside the contig: every emitted record has to move by
    2,000 bp in both files or the graph and the reads disagree.
    """

    genome, gtf, reads_bam, sg_bam = inputs
    manifest = _extract(
        tmp_path, genome, gtf, reads_bam, "chr1:2001-7000", sg_bam=sg_bam
    )

    assert manifest["offset"] == 2000
    main_contigs, main_pos = _positions(manifest["files"]["bam"])
    sg_contigs, sg_pos = _positions(manifest["files"]["sg_bam"])

    assert main_contigs == sg_contigs == {manifest["mini_contig_name"]}
    shared = set(main_pos) & set(sg_pos)
    assert shared, "no read is in both slices, so nothing is being compared"
    assert {name: sg_pos[name] for name in shared} == {
        name: main_pos[name] for name in shared
    }
    # And the rebase really happened: source starts are >= 3000, mini ones are not.
    assert max(sg_pos.values()) < 5000


def test_sg_counts_are_separate_from_the_read_counts(tmp_path, inputs):
    """The sg bam is thinned, so conflating the two accountings hides evidence loss."""

    genome, gtf, reads_bam, sg_bam = inputs
    manifest = _extract(tmp_path, genome, gtf, reads_bam, "chr1:1-10000", sg_bam=sg_bam)
    counts = manifest["counts"]

    assert counts["alignments_emitted"] == 16
    assert counts["sg_alignments_emitted"] == 8
    assert (
        counts["sg_alignments_emitted_forward"]
        + counts["sg_alignments_emitted_reverse"]
        == counts["sg_alignments_emitted"]
    )
    assert manifest["files"]["sg_bam"].endswith(".sg.bam")
    assert os.path.exists(manifest["files"]["sg_bam"] + ".bai")
    assert manifest["sg_bam_reused_from_source"] is False


def test_reuse_names_the_source_sg_bam_rather_than_copying_it(tmp_path, inputs):
    """The reuse optimization is symmetric, and says which file it named.

    A whole-contig strandless chunk at offset 0 would write the source's own
    records unmoved, for both files. Skipping the copy is only safe if the
    manifest states it, because stage 3b then has to restrict the source by
    contig instead of reading a mini bam that covers exactly the chunk.
    """

    genome, gtf, reads_bam, sg_bam = inputs
    manifest = extractor.extract_partition(
        genome_fa=str(genome),
        bam=str(reads_bam),
        region="chr1:1-10000",
        output_prefix=str(tmp_path / "reused"),
        gtf=str(gtf),
        max_intron_length=0,
        reuse_source_bam=True,
        sg_bam=str(sg_bam),
    )

    assert manifest["bam_reused_from_source"] is True
    assert manifest["sg_bam_reused_from_source"] is True
    assert manifest["files"]["sg_bam"] == os.path.abspath(str(sg_bam))
    assert not os.path.exists(str(tmp_path / "reused.sg.bam"))
    # Counted anyway: the manifest states what a full extraction would have
    # stated, and stage 3b's split is checked against these numbers.
    assert manifest["counts"]["sg_alignments_emitted"] == 8


def test_xw_survives_extraction(tmp_path, inputs):
    """`_rebase_alignment` reaches tags only through `fromstring`, so prove it."""

    genome, gtf, reads_bam, sg_bam = inputs
    manifest = _extract(
        tmp_path, genome, gtf, reads_bam, "chr1:2001-7000", sg_bam=sg_bam
    )

    first = _first_record(manifest["files"]["sg_bam"])
    assert first is not None
    assert first.has_tag("XW"), "the extractor dropped the tag LRAA:165-192 requires"
    assert first.get_tag("XW") == pytest.approx(0.5)


def test_overhang_drops_apply_to_both_slices_identically(tmp_path, inputs):
    """A record escaping the boundary is dropped from the sg slice too.

    Truncating it would invent an alignment; keeping it in one file and not the
    other would put a junction in the graph with no read able to support it.
    """

    genome, gtf, reads_bam, sg_bam = inputs
    straddler = ("straddle", 4950, 100, False)
    with_overhang = sorted(
        [("keep", 2500, 100, False), straddler], key=lambda r: r[1]
    )
    reads = _write_bam(tmp_path / "over.reads.bam", with_overhang)
    evidence = _write_bam(tmp_path / "over.sg.bam", with_overhang, xw=0.5)

    manifest = _extract(
        tmp_path, genome, None, reads, "chr1:2001-5000", sg_bam=evidence, prefix="ov"
    )

    assert manifest["counts"]["alignments_dropped_overhang"] == 1
    assert manifest["counts"]["sg_alignments_dropped_overhang"] == 1
    assert manifest["dropped_read_names"] == ["straddle"]
    assert manifest["sg_dropped_read_names"] == ["straddle"]
    _, sg_pos = _positions(manifest["files"]["sg_bam"])
    _, main_pos = _positions(manifest["files"]["bam"])
    assert set(sg_pos) == set(main_pos) == {"keep"}


def test_no_sg_bam_leaves_the_manifest_as_it_was(tmp_path, inputs):
    """Absent, not zero: a measurement never taken is not an empty splice graph."""

    genome, gtf, reads_bam, _sg = inputs
    manifest = _extract(tmp_path, genome, gtf, reads_bam, "chr1:1-10000")

    assert manifest["files"]["sg_bam"] is None
    assert manifest["sg_bam_reused_from_source"] is False
    assert "sg_dropped_read_names" not in manifest
    assert not [key for key in manifest["counts"] if key.startswith("sg_")]


def test_extractor_refuses_an_sg_bam_that_is_the_read_bam(tmp_path, inputs):
    """LRAA refuses the same composition at LRAA:1945-1957, one stage later."""

    genome, gtf, reads_bam, _sg = inputs
    with pytest.raises(extractor.ExtractionError) as excinfo:
        _extract(
            tmp_path, genome, gtf, reads_bam, "chr1:1-10000", sg_bam=reads_bam
        )
    assert "resolve to" in str(excinfo.value)


# ------------------------------------------------------------ argument refusals


def _refusal(tmp_path, genome, gtf, bam, *extra):
    result = _run_chunked(tmp_path, tmp_path / "work", genome, gtf, bam, *extra)
    assert result.returncode != 0
    return result.stdout + result.stderr


def test_bam_for_sg_without_a_shared_cut_source_is_refused(tmp_path, inputs):
    """Both available answers are wrong, so neither is chosen.

    Cuts from --bam give every caller its own geometry; cuts from the thinned sg
    bam lose raw reads that span a position it shows as free. The refusal names
    both, because a user who reaches this will otherwise try the second one.
    """

    genome, gtf, reads_bam, sg_bam = inputs
    combined = _refusal(tmp_path, genome, gtf, reads_bam, "--bam_for_sg", sg_bam)

    assert "--bam_for_sg needs --bam_for_cut_selection" in combined
    assert "its own cut geometry" in combined
    assert "coverage-normalized" in combined


def test_the_sg_bam_may_not_be_the_cut_selection_source(tmp_path, inputs):
    """The trap worth a message: normalization thins, so it is a subset."""

    genome, gtf, reads_bam, sg_bam = inputs
    combined = _refusal(
        tmp_path,
        genome,
        gtf,
        reads_bam,
        "--bam_for_sg",
        sg_bam,
        "--bam_for_cut_selection",
        sg_bam,
    )

    assert "resolve to the same file" in combined
    assert "unthinned superset" in combined


def test_bam_for_sg_is_refused_with_strand_first_chunking(tmp_path, inputs):
    """Stage 3b, which splits the sg bam, does not run in that mode."""

    genome, gtf, reads_bam, sg_bam = inputs
    combined = _refusal(
        tmp_path,
        genome,
        gtf,
        reads_bam,
        "--chunk_by_strand",
        "--bam_for_sg",
        sg_bam,
        "--bam_for_cut_selection",
        reads_bam,
    )

    assert "strandless-chunking only" in combined


def test_lraa_chunk_forwards_bam_for_sg_instead_of_ignoring_it(tmp_path, inputs):
    """`LRAA --chunk` builds its chunk args from an explicit allowlist.

    An unlisted flag is parsed and silently discarded, which for this one means
    every chunk manufacturing its own splice graph. Asserted by the run reaching
    ChunkedRun's own refusal -- a message only the pipeline produces -- rather
    than completing as though the flag had not been passed.
    """

    genome, gtf, reads_bam, sg_bam = inputs
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
            "--bam_for_sg",
            str(sg_bam),
        ],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
        timeout=1800,
    )
    combined = result.stdout + result.stderr
    assert result.returncode != 0
    assert "--bam_for_sg needs --bam_for_cut_selection" in combined


# ------------------------------------------------------------- end to end runs


@pytest.fixture
def sg_run(tmp_path, inputs):
    genome, gtf, reads_bam, sg_bam = inputs
    workdir = tmp_path / "work_sg"
    result = _run_chunked(
        tmp_path,
        workdir,
        genome,
        gtf,
        reads_bam,
        "--bam_for_sg",
        sg_bam,
        "--bam_for_cut_selection",
        reads_bam,
    )
    _ok(result)
    return workdir


@pytest.fixture
def plain_run(tmp_path, inputs):
    genome, gtf, reads_bam, _sg = inputs
    workdir = tmp_path / "work_plain"
    _ok(_run_chunked(tmp_path, workdir, genome, gtf, reads_bam))
    return workdir


def test_every_unit_sg_slice_matches_its_own_reads(sg_run):
    """Same mini contig, same coordinates, orientation for orientation.

    The per-unit files, not the chunk-level ones: this is what stage 5 opens, and
    it is produced by a SECOND tool (the strand split) after the extractor, so
    the geometry has two chances to drift before a unit sees it.
    """

    chunks = sorted(Path(sg_run).glob("chunks/*"))
    assert chunks
    for cdir in chunks:
        manifest = json.loads((cdir / "chunk.partition.json").read_text())
        for strand in ("+", "-"):
            main = cdir / "chunk.strand.{}.bam".format(strand)
            sg = cdir / "chunk.sg.strand.{}.bam".format(strand)
            assert sg.exists(), "no sg slice for {} of {}".format(strand, cdir.name)
            main_contigs, main_pos = _positions(main)
            sg_contigs, sg_pos = _positions(sg)
            assert sg_contigs == main_contigs == {manifest["mini_contig_name"]}
            shared = set(main_pos) & set(sg_pos)
            assert shared, "{} {} shares no read with its evidence".format(
                cdir.name, strand
            )
            assert {n: sg_pos[n] for n in shared} == {
                n: main_pos[n] for n in shared
            }


def test_xw_survives_the_strand_split_too(sg_run):
    """The split is a separate tool from the extractor and can lose tags on its own.

    LRAA:165-192 exits on a --bam_for_sg whose first aligned record has no XW, so
    this is what stands between a correct run and one that dies in stage 5 -- and
    the run below completed, which is half the proof.
    """

    for sg in sorted(Path(sg_run).glob("chunks/*/chunk.sg.strand.*.bam")):
        first = _first_record(sg)
        if first is None:
            continue
        assert first.has_tag("XW"), "{} lost XW in the strand split".format(sg)


def test_stage_4_is_skipped_and_stage_5_gets_the_sg_slice(sg_run):
    """Not merely "stage 5 got something" -- stage 4 must not have run at all.

    Normalizing supplied evidence a second time composes two acceptance rates,
    and the splice graph divides each read's rate back out through XW, so the
    weights would stop matching the reads that produced them.
    """

    argv = _chunk_logs(sg_run)
    assert "normalize_bam_by_strand.py" not in argv
    for strand in ("+", "-"):
        assert "--bam_for_sg" in argv
        assert "chunk.sg.strand.{}.bam".format(strand) in argv
    assert "--no_norm" in argv
    assert "chunk.norm.bam" not in argv
    assert not list(Path(sg_run).glob("chunks/*/chunk.*.norm.bam"))

    # The step record itself, not just the absent subprocess: a skip nobody wrote
    # down is indistinguishable from a stage that was never reached.
    timing = json.loads((Path(sg_run) / "timing.json").read_text())
    skips = [
        step
        for chunk in timing["arms"]["chunked"]["chunks"]
        for step in chunk["steps"]
        if step.get("step") == "stage4_normalize"
    ]
    assert skips, "no stage 4 record at all, so nothing states the skip"
    for step in skips:
        assert "skipped" in step, step
        assert "pre-normalized" in step["skipped"]
        assert step["sg_bam"].endswith(".bam")


def test_a_run_without_the_flag_still_normalizes_per_unit(plain_run):
    """The old path, unchanged. Stage 4 runs and stage 5 reads what it wrote."""

    argv = _chunk_logs(plain_run)
    assert "normalize_bam_by_strand.py" in argv
    assert "chunk.norm.bam" in argv or "norm.bam" in argv
    assert "chunk.sg.strand" not in argv
    for manifest in _manifests(plain_run):
        assert manifest["files"]["sg_bam"] is None
        assert not [k for k in manifest["counts"] if k.startswith("sg_")]
    assert not list(Path(plain_run).glob("chunks/*/chunk.sg.bam"))


def test_the_shared_cut_source_gives_two_different_callers_one_geometry(
    tmp_path, inputs
):
    """The cross-caller guarantee, end to end.

    Two runs whose --bam differ and whose --bam_for_cut_selection and
    --bam_for_sg are the same must partition identically and slice the shared
    evidence identically. Cut placement READS the bam -- blocked positions and
    the spanning-read cost both -- so without an explicit shared source these two
    runs are free to disagree, and every locus near a cut then gets a different
    splice graph per caller.
    """

    genome, gtf, reads_bam, sg_bam = inputs
    # Two callers' own reads: disjoint subsets of the library, as per-cluster bams
    # are. `superset` is what they were split out of, which is the legal source.
    left = _write_bam(
        tmp_path / "left.bam", _reads("fwd", [1000 + 40 * i for i in range(4)])
    )
    right = _write_bam(
        tmp_path / "right.bam",
        _reads("rev", [3000 + 40 * i for i in range(4)], reverse=True),
    )

    geometries = []
    slices = []
    for label, bam in (("left", left), ("right", right)):
        workdir = tmp_path / ("work_" + label)
        _ok(
            _run_chunked(
                tmp_path,
                workdir,
                genome,
                gtf,
                bam,
                "--bam_for_sg",
                sg_bam,
                "--bam_for_cut_selection",
                reads_bam,
                "--num_total_reads",
                "4",
            )
        )
        geometries.append(
            [
                (m["chrom"], m["partition_lend"], m["partition_rend"], m["offset"])
                for m in _manifests(workdir)
            ]
        )
        slices.append(
            [
                _positions(p)
                for p in sorted(Path(workdir).glob("chunks/*/chunk.sg.bam"))
            ]
        )

    assert geometries[0] == geometries[1], geometries
    assert slices[0] == slices[1], "the shared evidence was sliced differently"


def test_cut_selection_reads_the_shared_source_not_the_callers_bam(tmp_path, inputs):
    """The mechanism, asserted where it is decided rather than through the output.

    Cut positions are a deterministic function of the selector's ARGV, so two
    callers get identical geometry exactly when their selection invocations are
    identical. Checked here because the end-to-end assertion above cannot
    distinguish "one shared source" from "these two inputs happened to agree".

    The sentinel is checked separately and for a different property: it must
    NAME the source, so that a directory holding cuts computed one way is not
    reused by a run that wants the other, and must be byte-identical to today's
    when no source is supplied, so an existing output directory still resumes.
    """

    genome, gtf, reads_bam, sg_bam = inputs
    other = _write_bam(tmp_path / "other.bam", _reads("x", [500, 700, 900]))

    def plan(bam, **overrides):
        args = ChunkedRun.default_args(
            bam=str(bam),
            genome_fa=str(genome),
            gtf=str(gtf),
            output_dir=str(tmp_path / "w"),
            **overrides,
        )
        source = ChunkedRun.cut_sources(args, None, "inputs", None)[0]
        return ChunkedRun.cut_selection_plan(
            args, str(tmp_path / "w"), str(tmp_path / "w" / "cuts"), source, CONTIG
        )

    shared = dict(bam_for_cut_selection=str(reads_bam), bam_for_sg=str(sg_bam))
    a = plan(reads_bam, **shared)
    b = plan(other, **shared)
    assert a["cmd"] == b["cmd"]
    assert str(reads_bam) in a["cmd"]
    assert str(other) not in a["cmd"]

    # And without it the two DIFFER, which is what makes the input load-bearing
    # rather than decorative: each caller selects on its own reads.
    bare_a = plan(reads_bam)
    bare_b = plan(other)
    assert bare_a["cmd"] != bare_b["cmd"]
    assert str(other) in bare_b["cmd"]
    # The sentinel names the source only when one was given, so a default run's
    # is unchanged.
    assert ".cutsrc_" not in bare_a["token"]
    assert ".cutsrc_" in a["token"]


# -------------------------------------------------- severed accounting per input


def _chunk_record(chrom, lend, rend, dropped=(), sg_dropped=None):
    manifest = {
        "partition_lend": lend,
        "partition_rend": rend,
        "dropped_read_names": list(dropped),
        "counts": {},
    }
    if sg_dropped is not None:
        manifest["sg_dropped_read_names"] = list(sg_dropped)
    return {
        "chunk_id": "{}_{}".format(chrom, lend),
        "chrom": chrom,
        "strand": "",
        "strandless": True,
        "manifest": manifest,
    }


def test_shared_cut_positions_are_the_interior_boundaries(tmp_path):
    """The last base of every segment but the last, which is what a cut IS.

    `spanning_alignments` calls a read severed when `start <= position < end`, so
    an off-by-one here would make the re-derived set disagree with the extractor
    on every boundary rather than on none.
    """

    chunks = [
        _chunk_record(CONTIG, 1, 4000),
        _chunk_record(CONTIG, 4001, 8000),
        _chunk_record(CONTIG, 8001, 10000),
    ]
    assert ChunkedRun.shared_cut_positions(chunks) == {CONTIG: [4000, 8000]}


def test_per_input_accounting_accepts_a_correct_two_input_partition(tmp_path, inputs):
    """A read the shared selector named but this caller does not hold is fine.

    That is the whole reason the check moved per input: the selection source is a
    superset, so "named but not dropped here" is the expected state and the old
    one-source equality would fail on correct geometry.
    """

    genome, _gtf, reads_bam, sg_bam = inputs
    cut_dir = tmp_path / "cuts"
    cut_dir.mkdir()
    # What a superset selector would have written: a name from another caller.
    (cut_dir / "chr1.strandless.dropped_reads.txt").write_text("somebody_elses_read\n")

    straddler = ("straddle", 3950, 100, False)
    both = sorted([("inside", 1000, 100, False), straddler], key=lambda r: r[1])
    reads = _write_bam(tmp_path / "acct.reads.bam", both)
    evidence = _write_bam(tmp_path / "acct.sg.bam", both, xw=0.5)

    chunks = [
        _chunk_record(CONTIG, 1, 4000, dropped=["straddle"], sg_dropped=["straddle"]),
        _chunk_record(
            CONTIG, 4001, 10000, dropped=["straddle"], sg_dropped=["straddle"]
        ),
    ]
    report = ChunkedRun.verify_severed_accounting(
        str(cut_dir),
        chunks,
        inputs=[
            {"label": "reads", "bam": str(reads), "names_key": "dropped_read_names"},
            {
                "label": "splice_graph_evidence",
                "bam": str(evidence),
                "names_key": "sg_dropped_read_names",
            },
        ],
        max_intron_length=0,
    )

    assert report["shared_cut_selection"] is True
    assert report["named_absent_from_extraction_input"] == 1
    assert report["per_input"]["reads"]["severed_at_shared_cuts"] == 1
    assert report["per_input"]["splice_graph_evidence"]["dropped_by_extraction"] == 1


def test_per_input_accounting_still_fails_on_a_real_geometry_bug(tmp_path, inputs):
    """A check that cannot fail is worse than no check.

    The bug: a chunk boundary that does not sit where the reads were selected
    against, so a chunk holds a record crossing its own edge. Modelled by an
    extraction that dropped nothing while a read straddles the boundary its
    manifests declare -- exactly what a mis-shifted cut produces.
    """

    genome, _gtf, _reads, _sg = inputs
    cut_dir = tmp_path / "cuts2"
    cut_dir.mkdir()
    (cut_dir / "chr1.strandless.dropped_reads.txt").write_text("")

    both = sorted(
        [("inside", 1000, 100, False), ("straddle", 3950, 100, False)],
        key=lambda r: r[1],
    )
    reads = _write_bam(tmp_path / "bug.reads.bam", both)

    chunks = [
        _chunk_record(CONTIG, 1, 4000),
        _chunk_record(CONTIG, 4001, 10000),
    ]
    with pytest.raises(ChunkedRun.PipelineError) as excinfo:
        ChunkedRun.verify_severed_accounting(
            str(cut_dir),
            chunks,
            inputs=[
                {
                    "label": "reads",
                    "bam": str(reads),
                    "names_key": "dropped_read_names",
                }
            ],
            max_intron_length=0,
        )
    message = str(excinfo.value)
    assert "severed-read accounting is inexact for the reads input" in message
    assert "straddle" in message


def test_single_source_accounting_is_untouched(tmp_path):
    """No inputs means the check this has always been, fatal direction included."""

    cut_dir = tmp_path / "cuts3"
    cut_dir.mkdir()
    (cut_dir / "chr1.strandless.dropped_reads.txt").write_text("never_dropped\n")

    with pytest.raises(ChunkedRun.PipelineError) as excinfo:
        ChunkedRun.verify_severed_accounting(
            str(cut_dir), [_chunk_record(CONTIG, 1, 10000)]
        )
    assert "cut selection named 1 read(s) as severed" in str(excinfo.value)


# ----------------------------------------------------- the scattered leaf's view


def test_only_chunk_rebuilds_the_sg_paths_it_was_given(tmp_path):
    """The leaf re-derives every path from ITS OWN outdir, dropping the manifest's.

    A four-key rebuild silently omits the sg slice, and the leaf then builds the
    graph from its own reads -- the original defect, reintroduced on the scattered
    path only, where no chunk directory on the driver's machine would show it.
    """

    outdir = tmp_path / "leaf"
    cdir = outdir / "chunks" / "chr1_00"
    cdir.mkdir(parents=True)
    prefix = cdir / "chunk"
    (cdir / "chunk.partition.json").write_text(
        json.dumps(
            {
                "offset": 0,
                "window_origin": 0,
                "counts": {"alignments_emitted": 0},
                "files": {
                    "fasta": "/elsewhere/chunk.fa",
                    "bam": "/elsewhere/chunk.bam",
                    "gtf": None,
                    "dropped_reads": "/elsewhere/chunk.dropped_reads.txt",
                    "sg_bam": "/elsewhere/chunk.sg.bam",
                },
            }
        )
    )
    plan = {
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

    chunk = ChunkedRun.rebuild_chunk_record(plan, "chr1_00", str(outdir))

    assert chunk["manifest"]["files"]["sg_bam"] == "{}.sg.bam".format(prefix)
    assert [u["sg_bam"] for u in chunk["units"]] == [
        "{}.sg.strand.+.bam".format(prefix),
        "{}.sg.strand.-.bam".format(prefix),
    ]


def test_only_chunk_without_sg_evidence_keeps_normalizing(tmp_path):
    """The same rebuild must not invent evidence for a chunk that has none."""

    outdir = tmp_path / "leaf2"
    cdir = outdir / "chunks" / "chr1_00"
    cdir.mkdir(parents=True)
    (cdir / "chunk.partition.json").write_text(
        json.dumps(
            {
                "offset": 0,
                "window_origin": 0,
                "counts": {"alignments_emitted": 0},
                "files": {
                    "fasta": "/elsewhere/chunk.fa",
                    "bam": "/elsewhere/chunk.bam",
                    "gtf": None,
                    "dropped_reads": "/elsewhere/chunk.dropped_reads.txt",
                },
            }
        )
    )
    plan = {
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

    chunk = ChunkedRun.rebuild_chunk_record(plan, "chr1_00", str(outdir))

    assert chunk["manifest"]["files"]["sg_bam"] is None
    assert [u["sg_bam"] for u in chunk["units"]] == [None, None]
