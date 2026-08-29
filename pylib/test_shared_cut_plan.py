#!/usr/bin/env python3

"""One cut plan, emitted once, applied by every caller.

The defect this closes: cut POSITIONS were per-caller. Cut placement READS a bam
-- blocked positions and the spanning-read cost are both computed from it -- so
the ~29 per-cluster quant jobs of WDL/LRAA_quant_by_cluster.wdl each chose their
own boundaries, and the ONE shared depth-normalized bam they are handed as
`--bam_for_sg` was therefore sliced at a different geometry in every cluster.
Each cluster then built a different splice graph near every cut, which is exactly
the cross-cluster comparability supplying one sg bam exists to provide.

`--bam_for_cut_selection` fixed the geometry by naming a shared superset every
caller re-selects on, which pays selection ~29 times. `--emit_cut_plan` pays it
ONCE, over the pre-partition bam, and `--chunk_plan` applies that geometry to
each caller's own reads.

What these tests hold:

* two runs over DIFFERENT bams sharing one plan produce identical chunk ids,
  regions, mini contig names, offsets and partition bounds. This is the test the
  feature exists for;
* `--emit_cut_plan` extracts NOTHING -- the whole point of splitting it from
  `--stop_after_make_chunks`, which extracts every chunk of the bam it is given;
* a plan emitted from a superset and applied to a subset still ACCOUNTS for every
  read the shared cuts sever in that subset. Shared geometry that lost reads
  quietly would be worse than per-caller geometry, not better;
* `num_total_reads` and the run mode come from the CONSUMING run. The plan carries
  no authoritative denominator, so two callers of one plan get their own TPM
  scale;
* every mismatch between a plan and the run applying it is REFUSED with both
  values named -- geometry parameters, contig length, contig coverage, and the
  contradictory flag combinations. A mismatched plan does not fail downstream;
  it produces plausible output at bounds nobody chose;
* the ANNOTATION is checked DIRECTLY rather than by identity, and the verdict on
  what that finds is SPLIT BY CALL SITE. A plan cut that SEVERS a locus in the
  consuming gtf is refused BY NAME -- neither neighbour contains it, both omit
  it, and nobody quantifies it -- while a cut merely inside `--margin` of a locus
  is a WARNING with a per-run count, because the locus is held whole by one side
  and what it costs is severed reads, which this pipeline already prices as
  collateral. At cut SELECTION the margin stays a hard bar: the position is still
  movable there and the wiggle window exists to move it. One plan is emitted per
  run, before phase 1, so the consuming gtf is never the emitting one -- de novo
  emits with none at all, and ref-guided phase 2 quantifies the consolidated gtf
  -- which is why a margin shortfall at consumption has no remedy behind it;
* a genome-wide plan applied by a run restricted to one contig is ACCEPTED.
  `subwdls/LRAA_runner.wdl` passes `--chunk` unconditionally, so a by_chromosome
  shard is itself a chunked run holding one chromosome, and refusing it would
  break the mode this fix repairs;
* `LRAA --chunk --chunk_plan` reaches ChunkedRun rather than being parsed and
  dropped. LRAA builds its chunk args from an allowlist, which is how
  `--bam_for_sg` came to be silently ignored on every chunked run.
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
LRAA = REPO / "LRAA"

CONTIG = "chr1"
CONTIG_LEN = 20000
OTHER_CONTIG = "chr2"
OTHER_LEN = 8000

# 6 kb nominal spacing on a 20 kb contig: targets at 6000 and 12000, the 18000
# target merged into the tail (minimum_span is half the spacing, floored to a
# depth_window multiple = 3000, and 20000 - 18000 is under it). THREE chunks, so
# the partition has two interior boundaries to agree or disagree about -- one
# chunk would make every geometry assertion below vacuously true.
MB_PER_CUT = 0.006
# 1 kb of search around each target, so the window is [target-500, target+500]
# and every candidate is a multiple of the 100 bp depth window. The 1 Mb default
# would open the whole contig and let a cut land nowhere near its target.
WIGGLE = 0.001

sys.path.insert(0, str(REPO / "pylib"))

import ChunkedRun  # noqa: E402


# --------------------------------------------------------------------- fixtures


def _header(lengths):
    return pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": name, "LN": length} for name, length in lengths],
        }
    )


def _write_bam(path, reads, lengths=((CONTIG, CONTIG_LEN),), xw=None):
    """``reads`` is (name, contig, 0-based start, aligned length, is_reverse)."""

    header = _header(lengths)
    order = {name: i for i, (name, _length) in enumerate(lengths)}
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for name, contig, start, length, reverse in sorted(
            reads, key=lambda r: (order[r[1]], r[2])
        ):
            aln = pysam.AlignedSegment(header)
            aln.query_name = name
            aln.reference_id = order[contig]
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


def _dense_reads(prefix, lend, rend, step=50, length=300, contig=CONTIG):
    """Cover ``[lend, rend]`` so EVERY candidate position in it is spanned.

    Cut selection prices severing rather than forbidding it, so a window with no
    free position still gets a cut -- and that is the case worth testing, because
    a plan whose cuts sever nothing cannot demonstrate that a subset's severed
    reads are accounted for.
    """

    reads = []
    for i, start in enumerate(range(lend, rend, step)):
        reads.append(
            (
                "{}{}".format(prefix, i),
                contig,
                start,
                length,
                bool(i % 2),
            )
        )
    return reads


@pytest.fixture
def inputs(tmp_path):
    """A genome, an annotation clear of both cut windows, and a superset bam.

    The annotation is deliberately far from 6000 and 12000: cuts are inadmissible
    within ``chunk_margin`` of an annotated locus, and a gene sitting in a search
    window would decide the geometry by blocking rather than by depth.
    """

    genome = tmp_path / "genome.fa"
    genome.write_text(
        ">{}\n{}\n>{}\n{}\n".format(
            CONTIG, "A" * CONTIG_LEN, OTHER_CONTIG, "A" * OTHER_LEN
        )
    )
    pysam.faidx(str(genome))

    gtf = tmp_path / "annot.gtf"
    gtf.write_text(
        "\n".join(
            '{}\ttest\texon\t{}\t{}\t.\t{}\t.\tgene_id "{}"; transcript_id "{}";'.format(
                contig, lend, rend, strand, gene, gene + ".t1"
            )
            for gene, contig, lend, rend, strand in (
                ("gplus", CONTIG, 1001, 1400, "+"),
                ("gminus", CONTIG, 16001, 16400, "-"),
                ("gother", OTHER_CONTIG, 2001, 2400, "+"),
            )
        )
        + "\n"
    )

    reads = (
        _dense_reads("near_cut1_", 5000, 7000)
        + _dense_reads("near_cut2_", 11000, 13000)
        + _dense_reads("gene_plus_", 1000, 1400, step=100)
        + _dense_reads("gene_minus_", 16000, 16400, step=100)
        + _dense_reads("other_", 2000, 2400, step=100, contig=OTHER_CONTIG)
    )
    lengths = ((CONTIG, CONTIG_LEN), (OTHER_CONTIG, OTHER_LEN))
    superset = _write_bam(tmp_path / "superset.bam", reads, lengths=lengths)
    return genome, gtf, superset, reads, lengths


# ---------------------------------------------------------------------- helpers


def _run(tmp_path, workdir, genome, gtf, bam, *extra, num_total_reads=None):
    cmd = [
        sys.executable,
        str(CHUNKED_RUN),
        "--bam",
        str(bam),
        "--genome_fa",
        str(genome),
        "--output_dir",
        str(workdir),
        "--cpu_budget",
        "2",
        "--approx_MB_per_cut",
        str(MB_PER_CUT),
        "--approx_MB_per_cut_wiggle_window",
        str(WIGGLE),
        # Every chunk directory self-contained, which is what the scattered
        # workflow runs with: without it a chunk spanning its whole contig names
        # the source bam and writes no mini bam, so there would be nothing on disk
        # to compare between two callers.
        "--no_reuse_source_bam",
    ]
    if gtf is not None:
        cmd += ["--gtf", str(gtf)]
    if num_total_reads is not None:
        cmd += ["--num_total_reads", str(num_total_reads)]
    cmd += [str(x) for x in extra]
    return subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(tmp_path), timeout=3600
    )


def _ok(result):
    combined = result.stdout + result.stderr
    assert result.returncode == 0, combined[-8000:]
    return combined


def _refused(result):
    combined = result.stdout + result.stderr
    assert result.returncode != 0, combined[-8000:]
    return combined


def _emit(tmp_path, genome, gtf, bam, plan_path, label="emit", *extra):
    """Selection only: no -N, because a plan carries no TPM denominator."""

    return _run(
        tmp_path,
        tmp_path / ("work_" + label),
        genome,
        gtf,
        bam,
        "--emit_cut_plan",
        plan_path,
        *extra,
    )


def _apply(tmp_path, genome, gtf, bam, plan_path, label, *extra, num_total_reads=8):
    return _run(
        tmp_path,
        tmp_path / ("work_" + label),
        genome,
        gtf,
        bam,
        "--chunk_plan",
        plan_path,
        *extra,
        num_total_reads=num_total_reads,
    )


def _manifests(workdir):
    paths = sorted(Path(workdir).glob("chunks/*/chunk.partition.json"))
    assert paths, "no chunk manifests under {}".format(workdir)
    return [json.loads(p.read_text()) for p in paths]


def _geometry(workdir):
    """Everything about a chunk that a downstream coordinate depends on."""

    return [
        (
            os.path.basename(os.path.dirname(str(path))),
            m["chrom"],
            m["partition_lend"],
            m["partition_rend"],
            m["offset"],
            m["mini_contig_name"],
            m["mini_length"],
        )
        for path, m in (
            (p, json.loads(p.read_text()))
            for p in sorted(Path(workdir).glob("chunks/*/chunk.partition.json"))
        )
    ]


def _cut_positions(plan):
    return {
        selection["chrom"]: [cut["position"] for cut in selection["cuts"]]
        for selection in plan["geometry"]["by_source"][0]["selections"]
    }


def _timing(workdir):
    return json.loads((Path(workdir) / "timing.json").read_text())


def _chunk_logs(workdir):
    logs = sorted(Path(workdir).glob("logs/chunk_*.log"))
    assert logs, "no chunk logs under {}".format(workdir)
    return "\n".join(p.read_text() for p in logs)


@pytest.fixture
def shared_plan(tmp_path, inputs):
    """ONE plan over the whole pre-partition bam, and the cuts it chose."""

    genome, gtf, superset, reads, lengths = inputs
    plan_path = tmp_path / "shared_cut_plan.json"
    _ok(_emit(tmp_path, genome, gtf, superset, plan_path))
    plan = json.loads(plan_path.read_text())
    cuts = _cut_positions(plan)
    assert cuts[CONTIG], "the fixture placed no cut, so nothing below is tested"
    return plan_path, plan, cuts


def _subset_spanning(tmp_path, name, reads, lengths, position, contig=CONTIG):
    """A caller's own bam: a SUBSET of the superset that spans a shared cut.

    Built from the plan's chosen position rather than from a guess, so the test
    cannot silently degrade into one where nothing is severed.
    """

    kept = [
        r
        for r in reads
        # SEVERED, by the extractor's own containment rule: the left chunk ends AT
        # the cut, so a read whose 0-based half-open interval [start, end) reaches
        # strictly past it fits in neither neighbour. One ending exactly at the cut
        # is contained, not severed -- the earlier form of this predicate included
        # such a read and the test correctly refused to find it among the drops.
        if r[1] == contig and r[2] < position < r[2] + r[3]
    ]
    assert kept, "no read of the superset spans {}:{}".format(contig, position)
    # Plus some reads that span nothing, so the bam is not made entirely of
    # boundary cases.
    kept += [r for r in reads if r[1] == contig and r[2] < 1500]
    return _write_bam(tmp_path / name, kept, lengths=lengths), kept


# ------------------------------------------------- the guarantee this exists for


def test_two_callers_sharing_a_plan_get_identical_geometry(
    tmp_path, inputs, shared_plan
):
    """Different reads, one plan, ONE partition.

    Both bams are subsets of the superset the plan was selected on, and they share
    no read. Without a shared plan each would select on its own reads and get its
    own boundaries; with one, every coordinate a downstream stage depends on is
    identical -- chunk ids, regions, offsets, and the mini contig each unit is
    quantified against.
    """

    genome, gtf, _superset, reads, lengths = inputs
    plan_path, _plan, cuts = shared_plan

    left = _write_bam(
        tmp_path / "left.bam",
        [r for r in reads if r[2] < cuts[CONTIG][0]],
        lengths=lengths,
    )
    right = _write_bam(
        tmp_path / "right.bam",
        [r for r in reads if r[2] >= cuts[CONTIG][0]],
        lengths=lengths,
    )

    geometries = []
    for label, bam in (("left", left), ("right", right)):
        _ok(
            _apply(
                tmp_path,
                genome,
                gtf,
                bam,
                plan_path,
                label,
                "--stop_after_make_chunks",
            )
        )
        geometries.append(_geometry(tmp_path / ("work_" + label)))

    # Three chunks on chr1 plus one on chr2: a single-chunk partition would make
    # the comparison below meaningless.
    assert len(geometries[0]) >= 4, geometries[0]
    assert geometries[0] == geometries[1], geometries


def test_emit_cut_plan_extracts_nothing(tmp_path, inputs, shared_plan):
    """The whole difference from --stop_after_make_chunks.

    That flag extracts every chunk of the bam it is handed; on a whole-genome bam
    that is ~300 extractions of the pre-partition library, and every caller would
    then repeat them for its own reads anyway.
    """

    plan_path, plan, cuts = shared_plan
    workdir = tmp_path / "work_emit"

    assert plan_path.exists()
    assert plan["chunks_extracted"] is False
    assert len(plan["chunks"]) >= 4, plan["chunks"]
    assert not list(workdir.glob("chunks/*/chunk.partition.json"))
    assert not list(workdir.glob("chunks/*/chunk.bam"))
    assert not (workdir / "merged").exists()
    # Selection DID run and is on disk, which is what the plan describes.
    assert list(workdir.glob("cuts/*.cuts.json"))
    assert cuts[CONTIG]


@pytest.mark.parametrize("mechanism", ("chunk_plan", "bam_for_cut_selection"))
def test_shared_geometry_extracts_the_callers_own_reads(
    tmp_path, inputs, shared_plan, mechanism
):
    """Shared cuts, OWN reads. Both mechanisms, because both got this wrong.

    `extraction_plan` took its `--bam` from the planned chunk, and for a
    strandless chunk that field holds the SELECTION source -- which under either
    mechanism is the shared superset. Every per-cluster run would then have
    extracted the pre-partition library's reads under its own cluster's name: 29
    identical quantifications, no error anywhere, and the per-input severed
    accounting still satisfied, because a superset drops a superset of the names
    derived from `--bam`. MEASURED before the fix: both chunks of a run whose
    --bam was cluster.bam were extracted with `--bam superset.bam`.
    """

    genome, gtf, superset, reads, lengths = inputs
    plan_path, _plan, _cuts = shared_plan
    mine = {"gene_plus_0", "gene_plus_1"}
    caller = _write_bam(
        tmp_path / "own_reads.bam",
        [r for r in reads if r[0] in mine],
        lengths=lengths,
    )

    shared_flag = (
        ["--chunk_plan", plan_path]
        if mechanism == "chunk_plan"
        else ["--bam_for_cut_selection", superset]
    )
    workdir = tmp_path / ("work_own_" + mechanism)
    _ok(
        _run(
            tmp_path,
            workdir,
            genome,
            gtf,
            caller,
            *shared_flag,
            "--contigs",
            CONTIG,
            "--stop_after_make_chunks",
            num_total_reads=len(mine),
        )
    )

    emitted = set()
    for bam_path in sorted(workdir.glob("chunks/*/chunk.bam")):
        with pysam.AlignmentFile(str(bam_path), "rb") as fh:
            emitted.update(aln.query_name for aln in fh.fetch(until_eof=True))
    assert emitted == mine, sorted(emitted)


def test_the_plan_records_the_geometry_that_chose_the_cuts(shared_plan):
    """Positions alone cannot be validated; the parameters behind them can."""

    _path, plan, _cuts = shared_plan
    params = plan["geometry"]["params"]
    assert params["approx_MB_per_cut"] == MB_PER_CUT
    assert params["approx_MB_per_cut_wiggle_window"] == WIGGLE
    assert params["strandless"] is True
    # The ANNOTATION is absent from the geometry in every form -- neither a
    # presence boolean nor a per-contig digest. Both were proxies for "no cut
    # falls inside a model this run quantifies", and both refused correct runs:
    # one plan is emitted for a whole run, so de novo single-cell emits with no
    # annotation while every consumer has one, and ref-guided emits against the
    # reference while phase 2 quantifies the consolidated gtf. The property is
    # checked directly instead; see ChunkedRun.cut_blocking_annotation_models and
    # the collision tests below.
    assert "annotation_present" not in params
    assert "annotation_identity" not in params
    assert "annotation" not in plan["geometry"]
    for key in (
        "depth_window",
        "grid_origin",
        "margin",
        "max_intron_length",
        "min_per_id",
        "min_mapping_quality",
        "severed_multiexon_weight",
    ):
        assert key in params, key
    source = plan["geometry"]["by_source"][0]
    assert source["key"] == ""
    assert source["bam_identity"]
    segments = {
        s["chrom"]: [(seg["lend"], seg["rend"]) for seg in s["segments"]]
        for s in source["selections"]
    }
    assert len(segments[CONTIG]) == len(_cut_positions(plan)[CONTIG]) + 1
    assert segments[CONTIG][0][0] == 1
    assert segments[CONTIG][-1][1] == CONTIG_LEN


def test_the_plan_identity_is_in_the_extraction_sentinel(tmp_path, inputs):
    """A resumed run must not serve chunks extracted on another plan's geometry.

    Cut positions are not a function of anything else the stage-3 token carries,
    so two plans over the same bam at the same --margin produce the same sentinel
    without this -- and every offset downstream would then be silently wrong. The
    no-plan token must be byte-identical to today's, so an existing output
    directory still resumes.
    """

    genome, gtf, superset, _reads, _lengths = inputs
    first = tmp_path / "plan_a.json"
    second = tmp_path / "plan_b.json"
    first.write_text("{}")
    second.write_text("{}  ")

    def token(plan_path):
        args = ChunkedRun.default_args(
            bam=str(superset),
            genome_fa=str(genome),
            gtf=str(gtf),
            output_dir=str(tmp_path / "sentinel_work"),
            chunk_plan=None if plan_path is None else str(plan_path),
        )
        source = ChunkedRun.cut_sources(args, None, "inputs", None)[0]
        selection = {
            "chrom": CONTIG,
            "contig_length": CONTIG_LEN,
            "segments": [{"lend": 1, "rend": CONTIG_LEN}],
        }
        planned = next(
            iter(ChunkedRun.planned_chunks_for_selection(source, selection))
        )
        planned["order"] = 0
        return ChunkedRun.extraction_plan(
            args,
            str(tmp_path / "sentinel_work"),
            str(tmp_path / "sentinel_work" / "chunks"),
            planned,
            "upstream",
        )["token"]

    bare, a, b = token(None), token(first), token(second)
    assert ".plan_" not in bare
    assert ".plan_" in a and ".plan_" in b
    assert len({bare, a, b}) == 3


def test_a_dry_run_on_a_plan_reports_the_shared_partition(
    tmp_path, inputs, shared_plan
):
    """Neither selects nor extracts: the plan already answered the first."""

    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, plan, _cuts = shared_plan
    combined = _ok(
        _apply(tmp_path, genome, gtf, superset, plan_path, "dry", "--dry_run")
    )
    workdir = tmp_path / "work_dry"
    assert "PLAN (--dry_run)" in combined
    for chunk in plan["chunks"]:
        assert chunk["chunk_id"] in combined
    assert not list(workdir.glob("chunks/*/chunk.partition.json"))
    assert not list(workdir.glob("cuts/*.cuts.json"))


def test_shared_cuts_that_sever_a_subsets_reads_are_counted(
    tmp_path, inputs, shared_plan
):
    """Shared geometry must not lose a read without saying so.

    The subset holds reads that span a cut chosen on the superset, so extraction
    HAS to refuse them -- an alignment crossing its chunk's boundary cannot be
    placed whole. What matters is that they are counted and named rather than
    quietly absent, and that the per-input check re-derives the severed set from
    THIS bam and demands extraction dropped all of it.
    """

    genome, gtf, _superset, reads, lengths = inputs
    plan_path, _plan, cuts = shared_plan
    position = cuts[CONTIG][0]
    subset, kept = _subset_spanning(
        tmp_path, "spanning.bam", reads, lengths, position
    )
    spanning_names = {
        r[0] for r in kept if r[2] < position < r[2] + r[3]
    }

    workdir = tmp_path / "work_span"
    _ok(
        _apply(
            tmp_path,
            genome,
            gtf,
            subset,
            plan_path,
            "span",
            "--stop_after_make_chunks",
        )
    )

    dropped = set()
    overhang = 0
    for manifest in _manifests(workdir):
        dropped.update(manifest["dropped_read_names"])
        overhang += manifest["counts"]["alignments_dropped_overhang"]
    assert spanning_names <= dropped, sorted(spanning_names - dropped)
    assert overhang >= len(spanning_names)

    accounting = _timing(workdir)["severed_read_accounting"]
    assert accounting["shared_cut_selection"] is True
    per_input = accounting["per_input"]["reads"]
    # Re-derived from this bam by the selector's own spanning test, and every one
    # of them accounted for: `verify_severed_accounting` raises when extraction
    # placed one whole.
    assert per_input["severed_at_shared_cuts"] >= len(spanning_names)
    assert per_input["dropped_by_extraction"] >= per_input["severed_at_shared_cuts"]
    # No selection ran in this output directory, so nothing named these drops in
    # advance; they are persisted where the baseline's own glob finds them.
    assert accounting["dropped_not_named"] == len(dropped)
    named = (workdir / "cuts" / ChunkedRun.EXTRACTION_ONLY_DROPS).read_text().split()
    assert spanning_names <= set(named)


def test_num_total_reads_comes_from_the_consuming_run(tmp_path, inputs, shared_plan):
    """One plan, two callers, two denominators.

    The plan states none: it is written by a run that never quantifies, and a
    denominator copied from it would give every cluster the pre-partition
    library's TPM scale.
    """

    genome, gtf, _superset, reads, lengths = inputs
    plan_path, plan, cuts = shared_plan
    assert plan["num_total_reads"] is None

    left = _write_bam(
        tmp_path / "n_left.bam",
        [r for r in reads if r[1] == CONTIG and r[2] < 1500],
        lengths=lengths,
    )
    right = _write_bam(
        tmp_path / "n_right.bam",
        [r for r in reads if r[1] == CONTIG and 16000 <= r[2]],
        lengths=lengths,
    )

    for label, bam, denominator in (("n7", left, 7), ("n13", right, 13)):
        _ok(
            _apply(
                tmp_path,
                genome,
                gtf,
                bam,
                plan_path,
                label,
                "--contigs",
                CONTIG,
                num_total_reads=denominator,
            )
        )
        workdir = tmp_path / ("work_" + label)
        assert _timing(workdir)["num_total_reads"] == denominator
        argv = _chunk_logs(workdir)
        assert "--num_total_reads {}".format(denominator) in argv, argv[-4000:]
        other = 13 if denominator == 7 else 7
        assert "--num_total_reads {}".format(other) not in argv


@pytest.mark.parametrize("restriction", ("--contig", "--contigs"))
def test_a_genome_wide_plan_applies_to_a_single_contig_run(
    tmp_path, inputs, shared_plan, restriction
):
    """A by_chromosome shard is itself a chunked run, holding one chromosome.

    `subwdls/LRAA_runner.wdl` passes `--chunk` unconditionally, so refusing a run
    that touches only part of a genome-wide plan would break that mode outright.
    The contigs it does not process are not its problem; the ones it does must
    match the plan exactly.
    """

    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, plan, _cuts = shared_plan
    assert {c["chrom"] for c in plan["chunks"]} == {CONTIG, OTHER_CONTIG}

    _ok(
        _apply(
            tmp_path,
            genome,
            gtf,
            superset,
            plan_path,
            "one_contig" + restriction.strip("-"),
            # Both spellings, because a by_chromosome shard uses the singular
            # (LRAA forwards --contig) while a workflow naming a main-chromosome
            # set uses the plural. They reach the same enumeration and must reach
            # the same validation.
            restriction,
            CONTIG,
            "--stop_after_make_chunks",
        )
    )
    geometry = _geometry(tmp_path / ("work_one_contig" + restriction.strip("-")))
    assert {row[1] for row in geometry} == {CONTIG}
    planned = [
        (c["chunk_id"], c["region"])
        for c in plan["chunks"]
        if c["chrom"] == CONTIG
    ]
    assert [(row[0], "{}:{}-{}".format(row[1], row[2], row[3])) for row in geometry] == (
        planned
    )


def test_chunk_plan_with_bam_for_sg_is_accepted_and_slices_the_shared_evidence(
    tmp_path, inputs, shared_plan
):
    """The refusal accepts a plan, because a plan IS the shared geometry.

    `--bam_for_sg` without a shared geometry source is refused: cuts from `--bam`
    give every caller its own bounds, and cuts from the thinned sg bam let a cut
    look unspanned where a raw read spans it. A plan answers both.
    """

    genome, gtf, _superset, reads, lengths = inputs
    plan_path, _plan, _cuts = shared_plan
    caller = _write_bam(
        tmp_path / "sg_caller.bam",
        [r for r in reads if r[1] == CONTIG],
        lengths=lengths,
    )
    evidence = _write_bam(
        tmp_path / "sg_evidence.bam",
        [r for r in reads if r[1] == CONTIG][::3],
        lengths=lengths,
        xw=0.5,
    )

    _ok(
        _apply(
            tmp_path,
            genome,
            gtf,
            caller,
            plan_path,
            "sg",
            "--bam_for_sg",
            evidence,
            "--contigs",
            CONTIG,
            "--stop_after_make_chunks",
        )
    )
    workdir = tmp_path / "work_sg"
    slices = sorted(workdir.glob("chunks/*/chunk.sg.bam"))
    assert len(slices) == len(_geometry(workdir))
    for manifest in _manifests(workdir):
        assert manifest["files"]["sg_bam"]
        assert "sg_dropped_read_names" in manifest
    accounting = _timing(workdir)["severed_read_accounting"]
    assert "splice_graph_evidence" in accounting["per_input"]


# ----------------------------------------------------------------- the refusals


def test_a_contig_the_plan_does_not_partition_is_refused(
    tmp_path, inputs, shared_plan
):
    """Selecting cuts locally for it would give this caller unique bounds."""

    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, plan, _cuts = shared_plan
    trimmed = dict(plan)
    geometry = json.loads(json.dumps(plan["geometry"]))
    geometry["by_source"][0]["selections"] = [
        s for s in geometry["by_source"][0]["selections"] if s["chrom"] != OTHER_CONTIG
    ]
    trimmed["geometry"] = geometry
    partial = tmp_path / "partial_plan.json"
    partial.write_text(json.dumps(trimmed))

    combined = _refused(
        _apply(tmp_path, genome, gtf, superset, partial, "partial")
    )
    assert "processes contig {}".format(OTHER_CONTIG) in combined
    assert "does not partition" in combined


def test_a_contig_length_disagreement_is_refused(tmp_path, inputs, shared_plan):
    """The last segment of every contig runs to its end."""

    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, plan, _cuts = shared_plan
    edited = json.loads(json.dumps(plan))
    for selection in edited["geometry"]["by_source"][0]["selections"]:
        if selection["chrom"] == CONTIG:
            selection["contig_length"] = CONTIG_LEN + 1000
    stale = tmp_path / "stale_length.json"
    stale.write_text(json.dumps(edited))

    combined = _refused(_apply(tmp_path, genome, gtf, superset, stale, "stale"))
    assert str(CONTIG_LEN + 1000) in combined
    assert str(CONTIG_LEN) in combined
    assert "final boundary" in combined


def test_a_geometry_parameter_disagreement_is_refused(tmp_path, inputs, shared_plan):
    """A plan cut at one spacing applied by a run configured for another."""

    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, _plan, _cuts = shared_plan
    result = _run(
        tmp_path,
        tmp_path / "work_spacing",
        genome,
        gtf,
        superset,
        "--chunk_plan",
        plan_path,
        "--margin",
        "777",
        num_total_reads=8,
    )
    combined = _refused(result)
    assert "margin" in combined
    assert "777" in combined
    assert "different partition" in combined


def test_a_model_straddling_a_plan_cut_is_refused(tmp_path, inputs, shared_plan):
    """The property the digest was a proxy for, asked directly.

    A cut at 1-based p splits the contig into [lend, p] and [p+1, rend], and the
    extractor emits an annotated locus whole or not at all, so a locus spanning p
    is contained by neither neighbour and BOTH omit it: quantified by nobody,
    with every chunk reporting success. The novel model below sits across the
    6000 cut this fixture's plan chose, which is the single-cell hazard -- phase
    2 quantifies the CONSOLIDATED gtf, carrying models phase-1 discovery found.

    The refusal must NAME it: the transcript, the cut, and both coordinate pairs.
    A message that says only "the annotation differs" leaves the operator with a
    1.5 GB gtf and no idea which model.
    """

    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, _plan, cuts = shared_plan
    assert 6000 in cuts[CONTIG], cuts
    consolidated = tmp_path / "consolidated.gtf"
    consolidated.write_text(
        gtf.read_text()
        + '{}\ttest\texon\t5800\t6200\t.\t+\t.\tgene_id "novel"; '
        'transcript_id "novel.t1";\n'.format(CONTIG)
    )

    combined = _refused(
        _apply(tmp_path, genome, consolidated, superset, plan_path, "straddle")
    )
    assert "transcript novel.t1 of gene novel" in combined
    assert "places a cut at {}:6000".format(CONTIG) in combined
    assert "{}:5800-6200".format(CONTIG) in combined
    assert "quantified by nobody" in combined


def test_a_different_but_non_straddling_annotation_is_accepted(
    tmp_path, inputs, shared_plan
):
    """THE REASON THE DIGEST WENT, stated as the case it got wrong.

    Same shape as the refusal above -- a plan selected against the reference,
    applied by a run holding a consolidated gtf that the plan has never seen --
    but the added model sits INSIDE a chunk instead of across its boundary. That
    is what phase 2 of the single-cell pipeline always looks like, because phase
    1 discovered its models inside these very chunks and so cannot have produced
    one spanning a cut. A sha256 of the annotation refuses this and the run that
    straddles alike; the straddle check separates them.
    """

    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, plan, cuts = shared_plan
    interior = tmp_path / "consolidated_interior.gtf"
    interior.write_text(
        gtf.read_text()
        + '{}\ttest\texon\t7000\t7400\t.\t+\t.\tgene_id "novel_interior"; '
        'transcript_id "novel_interior.t1";\n'.format(CONTIG)
    )
    assert interior.read_text() != gtf.read_text()
    # Genuinely between two cuts, so acceptance is about the straddle rule rather
    # than about the fixture happening to have one chunk.
    assert 6000 in cuts[CONTIG] and 12000 in cuts[CONTIG], cuts

    _ok(
        _apply(
            tmp_path,
            genome,
            interior,
            superset,
            plan_path,
            "interior_gtf",
            "--stop_after_make_chunks",
        )
    )
    geometry = _geometry(tmp_path / "work_interior_gtf")
    planned = [(c["chunk_id"], c["region"]) for c in plan["chunks"]]
    assert [
        (row[0], "{}:{}-{}".format(row[1], row[2], row[3])) for row in geometry
    ] == planned


# The interior cut this fixture's plan places on CONTIG, and the margin every run
# below uses -- neither is passed on any command line, so both are read from the
# same places the pipeline reads them.
CUT = 6000
MARGIN = ChunkedRun.LRAA_Globals.config["chunk_margin"]

# (label, locus lend, locus rend, kind, clearance) -- ONE single-exon locus placed
# against the cut at CUT. ``kind`` is what the plan check must call it, and what
# the CONSUMER then does with it:
#
#   "straddle"  the cut severs the locus, so both neighbours omit it and nobody
#               quantifies it. FATAL at every call site;
#   "margin"    one chunk holds the locus WHOLE but the boundary sits closer to
#               it than --margin. A hard bar at cut SELECTION, where the position
#               is still movable, and a WARNING at consumption, where it is not;
#   "clear"     the plan is fine and passes silently.
#
# ``clearance`` is the bases strictly between the locus and the cut JUNCTION,
# written as a LITERAL. The rule it samples -- a locus blocks a cut exactly when
# its clearance is under --margin -- is asserted below rather than assumed, on
# both sides of the cut, because ``blocks_cut``'s form is asymmetric
# (``lend <= cut + margin`` against ``rend >= cut - margin + 1``) and only the
# clearance makes the two sides one rule.
CUT_CLEARANCE_CASES = (
    ("straddles", 5900, 6100, "straddle", None),
    # THE DEMONSTRATED GAP: severs nothing, so the straddle test alone passed it
    # and every chunk touching the cut then died at margin 200.
    ("inside_margin_after", 6050, 6500, "margin", 49),
    ("inside_margin_before", 5500, 5950, "margin", 50),
    ("exactly_margin_after", 6200, 6600, "margin", 199),
    ("exactly_margin_before", 5401, 5801, "margin", 199),
    ("margin_plus_one_after", 6201, 6601, "clear", 200),
    ("margin_plus_one_before", 5400, 5800, "clear", 200),
    ("comfortably_clear_after", 7000, 7400, "clear", 999),
    ("comfortably_clear_before", 4600, 5000, "clear", 1000),
)

# THE JUNCTION, base by base. A cut at ``p`` is the boundary BETWEEN bases ``p``
# and ``p + 1``: the left chunk is ``[.., p]``, the right is ``[p + 1, ..]``. So a
# locus straddles exactly when it owns at least one base in each, which is what
# the test below derives by ENUMERATING its bases rather than by re-deriving
# ``lend <= cut < rend``.
#
# These rows exist as their own table because straddle is now the only FATAL
# verdict a plan consumer can reach: a margin shortfall warns, so a one-off here
# is not a spurious refusal any more, it is a locus quantified by nobody while
# every chunk reports success. The 1 bp between ``ends_at_cut`` and
# ``ends_one_past_cut`` is the whole difference between advisory and fatal.
JUNCTION_BOUNDARY_CASES = (
    # last base of the locus IS the last base of the left chunk
    ("ends_at_cut", CUT - 200, CUT, "margin", 0),
    # ... one base further and that base belongs to the right chunk
    ("ends_one_past_cut", CUT - 200, CUT + 1, "straddle", None),
    ("starts_at_cut", CUT, CUT + 200, "straddle", None),
    ("starts_one_before_cut", CUT - 1, CUT + 200, "straddle", None),
    # the mirror of ends_at_cut: the right chunk holds it whole
    ("starts_one_past_cut", CUT + 1, CUT + 300, "margin", 0),
    # one base cannot lie either side of a junction, whichever side it is on
    ("single_base_at_cut", CUT, CUT, "margin", 0),
    ("contains_cut_strictly_inside", CUT - 100, CUT + 100, "straddle", None),
)

# The subset run END TO END. Both verdicts, both sides of the cut, an
# at-the-boundary and a comfortably-clear acceptance, and the 1 bp junction pair
# that separates fatal from advisory. The rest are pinned against the extractor's
# own predicate by the tests below, which call the same function extraction
# calls, so agreement there is by construction.
END_TO_END_LABELS = (
    "straddles",
    "inside_margin_after",
    "exactly_margin_before",
    "margin_plus_one_before",
    "comfortably_clear_after",
    "ends_at_cut",
    "ends_one_past_cut",
)
END_TO_END_CASES = tuple(
    case
    for case in CUT_CLEARANCE_CASES + JUNCTION_BOUNDARY_CASES
    if case[0] in END_TO_END_LABELS
)


def _locus_gtf(path, name, lend, rend, base="", strand="+", contig=CONTIG):
    """``base`` plus ONE extra single-exon locus, named for the assertion."""

    path = Path(path)
    path.write_text(
        base
        + '{}\ttest\texon\t{}\t{}\t.\t{}\t.\tgene_id "{}"; '
        'transcript_id "{}.t1";\n'.format(contig, lend, rend, strand, name, name)
    )
    return path


@pytest.mark.parametrize(
    "label,lend,rend,kind,clearance",
    JUNCTION_BOUNDARY_CASES,
    ids=[case[0] for case in JUNCTION_BOUNDARY_CASES],
)
def test_the_straddle_boundary_is_the_junction_between_bases(
    tmp_path, label, lend, rend, kind, clearance
):
    """PART 1: straddle exactness, derived from the partition, not from the code.

    A cut at ``p`` puts base ``p`` in the left chunk and base ``p + 1`` in the
    right one. So the verdict is settled by asking which side each of the locus'
    bases lands on -- enumerated below -- and a locus straddles exactly when it
    owns at least one base on each side. That derivation is independent of
    ``blocks_cut``'s ``lend <= cut + margin and rend >= cut - margin + 1``, whose
    asymmetry is what encodes the same between-bases fact; asserting the two
    agree is the point.

    Load-bearing as of the verdict split: a margin shortfall is now a warning, so
    straddle detection is the ONLY thing left standing between a plan and a locus
    that no chunk contains and no chunk reports missing. The pair ``ends_at_cut``
    / ``ends_one_past_cut`` differs by one base and by everything else.

    The margin-0 reduction is asserted for every row because it is the identity
    the fatal branch rests on: at margin 0 ``blocks_cut`` IS the straddle test, so
    a row that offends at margin 0 is severed and a row that does not is not.
    """

    extractor = ChunkedRun._extractor_module()
    bases = range(lend, rend + 1)
    in_left = [b for b in bases if b <= CUT]
    in_right = [b for b in bases if b > CUT]
    derived = "straddle" if (in_left and in_right) else "margin"
    assert derived == kind, (label, in_left[:3], in_right[:3])

    if kind == "margin":
        # bases strictly between the locus and the junction, counted on whichever
        # side the locus actually lies
        derived_gap = (
            len(range(rend + 1, CUT + 1)) if not in_right else len(range(CUT + 1, lend))
        )
        assert derived_gap == clearance, (label, derived_gap)

    # margin 0 reduces blocks_cut to the straddle test, and that is the reduction
    # the only-fatal-verdict branch is built on
    assert extractor.blocks_cut(lend, rend, CUT, 0) == (kind == "straddle")
    # at the run's real margin every row here blocks: the straddles by straddling,
    # the rest at clearance 0, which is under any positive margin
    assert extractor.blocks_cut(lend, rend, CUT, MARGIN) is True

    path = _locus_gtf(tmp_path / "junction.gtf", label, lend, rend)
    (offence,) = ChunkedRun.cut_blocking_annotation_models(
        str(path), {CONTIG: [CUT]}, MARGIN
    )
    assert offence["kind"] == kind, offence
    assert offence["gap"] == (None if kind == "straddle" else clearance), offence
    assert offence["gene_id"] == label
    assert (offence["gene_lend"], offence["gene_rend"]) == (lend, rend)

    at_zero = ChunkedRun.cut_blocking_annotation_models(
        str(path), {CONTIG: [CUT]}, 0
    )
    assert [o["kind"] for o in at_zero] == (
        ["straddle"] if kind == "straddle" else []
    )


def test_a_locus_is_defined_per_gene_id_and_strand_as_the_extractor_defines_it(
    tmp_path,
):
    """Two parity defects the verdict split made load-bearing. Both were real.

    ``_GtfIngest`` keys a locus on ``(gene_id, strand)`` -- re-keyed by 2d28609f,
    which fixed strandless extraction raising "gene X appears on both strands" --
    and sets ``gene_id = transcript_id`` for a line that names only a transcript.
    ``cut_blocking_annotation_models`` did neither, so:

    * OVER-REFUSAL. A gene_id present on BOTH strands had its two spans unioned
      into one, and every cut between them read as a straddle. The SIRV reference
      does this by design -- antisense isoforms per locus, gene named after the
      contig -- so ``SIRV1+`` at [1000, 2000] and ``SIRV1-`` at [9000, 9500] read
      as one locus [1000, 9500] straddling a cut at 6000, while the extractor saw
      two clear loci and extracted happily. Harmless while both verdicts refused;
      an unconditional hard failure once straddle is the only fatal one;
    * UNDER-REFUSAL. A line carrying ``transcript_id`` and no ``gene_id`` was
      skipped as "not placeable in a locus", but the extractor forms a locus from
      it and emits it by gene containment. So a genuinely severed model was
      invisible to the one check that is now the sole guard against silent loss.
    """

    both_strands = tmp_path / "both_strands.gtf"
    both_strands.write_text(
        '{c}\ttest\texon\t1000\t2000\t.\t+\t.\tgene_id "SIRV1"; '
        'transcript_id "SIRV1.p";\n'
        '{c}\ttest\texon\t9000\t9500\t.\t-\t.\tgene_id "SIRV1"; '
        'transcript_id "SIRV1.m";\n'.format(c=CONTIG)
    )
    assert (
        ChunkedRun.cut_blocking_annotation_models(
            str(both_strands), {CONTIG: [CUT]}, MARGIN
        )
        == []
    )
    # and the extractor agrees, over the same cut, on both neighbours
    extractor = ChunkedRun._extractor_module()
    annotation = extractor.load_gtf(str(both_strands), CONTIG, "")
    for region in ("{}:1-{}".format(CONTIG, CUT), "{}:{}-{}".format(CONTIG, CUT + 1, CONTIG_LEN)):
        assert (
            extractor.admissibility_offences(
                annotation, extractor.parse_region(region), CONTIG_LEN, MARGIN
            )
            == []
        ), region

    # the same gene_id where the + locus REALLY does straddle is still caught, so
    # the fix is per-strand rather than a blanket exemption
    real = tmp_path / "both_strands_real_straddle.gtf"
    real.write_text(
        '{c}\ttest\texon\t{l}\t{r}\t.\t+\t.\tgene_id "SIRV1"; '
        'transcript_id "SIRV1.p";\n'
        '{c}\ttest\texon\t9000\t9500\t.\t-\t.\tgene_id "SIRV1"; '
        'transcript_id "SIRV1.m";\n'.format(c=CONTIG, l=CUT - 100, r=CUT + 100)
    )
    (offence,) = ChunkedRun.cut_blocking_annotation_models(
        str(real), {CONTIG: [CUT]}, MARGIN
    )
    assert offence["kind"] == "straddle"
    assert (offence["gene_id"], offence["strand"]) == ("SIRV1", "+")
    assert (offence["gene_lend"], offence["gene_rend"]) == (CUT - 100, CUT + 100)

    # a transcript_id with no gene_id is a locus to the extractor, so it must be
    # one here too
    orphan = tmp_path / "transcript_only.gtf"
    orphan.write_text(
        '{}\ttest\texon\t{}\t{}\t.\t+\t.\ttranscript_id "orphan.t1";\n'.format(
            CONTIG, CUT - 100, CUT + 100
        )
    )
    (offence,) = ChunkedRun.cut_blocking_annotation_models(
        str(orphan), {CONTIG: [CUT]}, MARGIN
    )
    assert offence["kind"] == "straddle"
    assert offence["gene_id"] == "orphan.t1"
    orphan_annotation = extractor.load_gtf(str(orphan), CONTIG, "")
    assert [
        o["kind"]
        for o in extractor.admissibility_offences(
            orphan_annotation,
            extractor.parse_region("{}:1-{}".format(CONTIG, CUT)),
            CONTIG_LEN,
            MARGIN,
        )
    ] == ["straddle"]


@pytest.mark.parametrize(
    "label,lend,rend,kind,clearance",
    CUT_CLEARANCE_CASES,
    ids=[case[0] for case in CUT_CLEARANCE_CASES],
)
def test_the_plan_check_reports_exactly_what_extraction_checks(
    tmp_path, label, lend, rend, kind, clearance
):
    """THE AGREEMENT, at the boundary, in both directions.

    The plan check used to test STRICT straddle only, while extraction checks any
    boundary whose clearance from an annotated locus is under ``--margin``. So a
    locus at ``[p + 50, p + 500]`` against a cut at ``p`` was not reported here at
    all, and every chunk touching ``p`` then died on it mid-run -- which is what
    validating the plan against the arguments exists to prevent.

    REPORTING, not refusing: this asserts what the predicate SEES, which is the
    same set extraction sees. What each caller then DOES with a record differs by
    call site -- fatal on a straddle everywhere, hard at selection and advisory at
    consumption on a margin shortfall -- and is asserted by the end-to-end table
    below and by ``test_selection_still_refuses_a_margin_only_position``.

    Both directions matter. Seeing less than extraction is the defect that cost
    the mid-run deaths; seeing MORE would attach a verdict to a boundary
    extraction is happy with, so every acceptance here is asserted against
    ``blocks_cut`` -- the extractor's own predicate, which the check calls rather
    than reimplements.

    The margin-0 assertion at the end is the straddle reduction, stated as data:
    at margin 0 ``blocks_cut`` IS the strict straddle, so only the severed row
    offends. Everything else in this table is what the check was missing.
    """

    extractor = ChunkedRun._extractor_module()
    blocked = extractor.blocks_cut(lend, rend, CUT, MARGIN)
    assert blocked == (kind != "clear")
    if kind != "straddle":
        # ONE rule, both sides: clearance under --margin blocks, at or above it
        # does not. This is what "exactly margin" and "margin + 1" mean.
        assert (clearance < MARGIN) == blocked

    path = _locus_gtf(tmp_path / "one_locus.gtf", label, lend, rend)
    offences = ChunkedRun.cut_blocking_annotation_models(
        str(path), {CONTIG: [CUT]}, MARGIN
    )
    if kind == "clear":
        assert offences == []
    else:
        assert len(offences) == 1, offences
        (offence,) = offences
        assert offence["kind"] == kind
        assert offence["gene_id"] == label
        assert (offence["lend"], offence["rend"]) == (lend, rend)
        assert offence["gap"] == (None if kind == "straddle" else clearance)
        assert offence["margin"] == MARGIN

    at_zero = ChunkedRun.cut_blocking_annotation_models(
        str(path), {CONTIG: [CUT]}, 0
    )
    assert [o["kind"] for o in at_zero] == (
        ["straddle"] if kind == "straddle" else []
    )


@pytest.mark.parametrize(
    "label,lend,rend,kind,clearance",
    END_TO_END_CASES,
    ids=[case[0] for case in END_TO_END_CASES],
)
def test_a_plan_cut_colliding_with_a_locus_is_fatal_only_when_it_severs_it(
    tmp_path, inputs, shared_plan, label, lend, rend, kind, clearance
):
    """WAS ``test_a_locus_colliding_with_a_plan_cut_is_refused_before_extraction``,
    which asserted the ``"margin"`` rows are REFUSED. They are now ACCEPTED with a
    warning, and the straddle rows are still refused; that is the whole change.

    Why the previous verdict was wrong: the margin is satisfiable only while a cut
    is still being PLACED. Once the geometry is fixed, a consumer whose model set
    differs from the emitter's -- which is the normal case, since ref-guided phase
    2 quantifies the consolidated gtf and a de novo plan is emitted against no
    annotation at all -- has no remedy for a margin shortfall. Enforcing it cost
    four multi-hour runs, every one margin-only with zero straddles, clearances of
    95, 3, 23 and 52 bp, one of them 5 h 26 m in on a cut the plan had already
    moved 211 bp to clear the annotated locus.

    What the margin rows must now show is that the run PROCEEDS, extracts at the
    plan's own geometry, and says why in a warning naming the locus, the cut, the
    exact clearance and this run's --margin -- plus the per-run count, which is
    what turns hundreds of individually benign shortfalls into a visible signal.

    A straddle refusal is worth having only if it happens BEFORE the extractor's
    own -- the extractor refuses it too, loudly, but minutes into a run and once
    per chunk. So the straddle assertions are that the message is the PLAN
    check's, that the extractor's is absent, and that no chunk was written.
    """

    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, plan, cuts = shared_plan
    assert CUT in cuts[CONTIG], cuts
    consumer = _locus_gtf(
        tmp_path / "consumer_{}.gtf".format(label),
        label,
        lend,
        rend,
        base=gtf.read_text(),
    )
    workdir = tmp_path / ("work_" + label)
    result = _apply(
        tmp_path,
        genome,
        consumer,
        superset,
        plan_path,
        label,
        "--stop_after_make_chunks",
    )

    if kind == "straddle":
        combined = _refused(result)
        assert "chunk plan" in combined
        assert "places a cut at {}:{}".format(CONTIG, CUT) in combined
        assert "transcript {}.t1 of gene {}".format(label, label) in combined
        assert "{}:{}-{}".format(CONTIG, lend, rend) in combined
        assert "quantified by nobody" in combined
        # the plan check got there first, and nothing was extracted
        assert "REJECTED:" not in combined
        assert not sorted(workdir.glob("chunks/*/chunk.partition.json"))
        return

    # both "clear" and "margin" now run to completion at the plan's geometry
    combined = _ok(result)
    geometry = _geometry(workdir)
    planned = [(c["chunk_id"], c["region"]) for c in plan["chunks"]]
    assert [
        (row[0], "{}:{}-{}".format(row[1], row[2], row[3])) for row in geometry
    ] == planned

    if kind == "clear":
        assert "margin audit" not in combined
        return

    # the advisory, carrying everything the refusal used to carry
    assert "places a cut at {}:{}".format(CONTIG, CUT) in combined
    assert "transcript {}.t1 of gene {}".format(label, label) in combined
    assert "clears by only {} bp".format(clearance) in combined
    assert "--margin of {} bp".format(MARGIN) in combined
    assert "holds the locus WHOLE -- nothing is severed" in combined
    assert "PROCEEDING" in combined
    # and the per-run count, with the straddle tally stated as zero: that is what
    # says the fatal guard ran rather than having been skipped
    assert (
        "margin audit: 1 annotated locus/loci sit within this run's --margin of "
        "{} bp of a plan cut, and 0 straddle one".format(MARGIN) in combined
    )
    # the locus is emitted WHOLE by exactly one chunk, which is the fact that
    # makes the shortfall survivable
    holders = [
        m
        for m in _manifests(workdir)
        if "{}.t1".format(label) in (m.get("emitted_transcript_ids") or [])
    ]
    assert len(holders) == 1, [m["mini_contig_name"] for m in holders]


def _plan_args(genome, gtf, bam, plan_path, outdir):
    return ChunkedRun.default_args(
        bam=str(bam),
        genome_fa=str(genome),
        gtf=str(gtf),
        output_dir=str(outdir),
        chunk_plan=str(plan_path),
        approx_MB_per_cut=MB_PER_CUT,
        approx_MB_per_cut_wiggle_window=WIGGLE,
    )


def test_a_straddling_plan_is_still_refused_and_only_the_straddle_branch_does_it(
    tmp_path, inputs, shared_plan, capsys
):
    """The fatal branch, and a NEGATIVE CONTROL that it is the one doing the work.

    After the split, a straddle is the only thing standing between a plan and a
    locus that no chunk contains and no chunk reports missing -- so a test that
    would still pass with that branch downgraded is worth nothing. The control
    relabels the offence's ``kind`` from ``"straddle"`` to ``"margin"``, which is
    exactly what downgrading the branch does to the verdict, and asserts the
    refusal DISAPPEARS: the run then proceeds and warns. Part 1 above is what
    makes that label trustworthy enough to hang a fatal verdict on.
    """

    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, _plan, cuts = shared_plan
    assert CUT in cuts[CONTIG], cuts
    consumer = _locus_gtf(
        tmp_path / "spanner.gtf",
        "spanner",
        CUT - 100,
        CUT + 100,
        base=gtf.read_text(),
    )
    args = _plan_args(genome, consumer, superset, plan_path, tmp_path / "unused")
    loaded = ChunkedRun.load_chunk_plan(str(plan_path))

    with pytest.raises(ChunkedRun.PipelineError) as excinfo:
        ChunkedRun.selections_from_chunk_plan(
            loaded, str(plan_path), args, [CONTIG]
        )
    message = str(excinfo.value)
    assert "transcript spanner.t1 of gene spanner" in message
    assert "quantified by nobody" in message

    # NEGATIVE CONTROL: the same inputs with the straddle relabelled as a margin
    # shortfall -- the verdict a downgrade of that branch would reach -- are
    # ACCEPTED. So the refusal above comes from the straddle branch and from
    # nothing else, and this test fails if that branch is ever relaxed.
    real = ChunkedRun.cut_blocking_annotation_models

    def downgraded(gtf_path, cuts_by_contig, margin):
        offences = []
        for offence in real(gtf_path, cuts_by_contig, margin):
            if offence["kind"] == "straddle":
                offence = dict(offence, kind="margin", side="before", gap=0)
            offences.append(offence)
        return offences

    ChunkedRun.cut_blocking_annotation_models = downgraded
    try:
        selections = ChunkedRun.selections_from_chunk_plan(
            loaded, str(plan_path), args, [CONTIG]
        )
    finally:
        ChunkedRun.cut_blocking_annotation_models = real
    assert [s["chrom"] for s in selections[""]] == [CONTIG]
    warned = capsys.readouterr().err
    assert "margin audit" in warned
    assert "and 0 straddle one" in warned


def test_selection_still_refuses_a_margin_only_position(tmp_path, inputs):
    """THE OTHER CALL SITE, unchanged, because there the margin can still act.

    The downgrade is a CONSUMPTION verdict and must not leak into placement. At
    selection the position is not yet chosen: the wiggle window exists precisely
    to move it, clearing an annotated boundary costs nothing, and a cut placed
    clear of every locus is a cut no consumer ever has to be warned about. So a
    margin-only position stays inadmissible here -- both in the zone arithmetic
    and in what the selector actually emits.
    """

    extractor = ChunkedRun._extractor_module()
    genome, gtf, superset, _reads, _lengths = inputs

    # a locus 50 bp past the cut target: severs nothing, so it is exactly the
    # case consumption now tolerates
    blocker = _locus_gtf(
        tmp_path / "selection_blocker.gtf",
        "blocker",
        CUT + 50,
        CUT + 500,
        base=gtf.read_text(),
    )

    # ZONE ARITHMETIC: the position is blocked, so it is in no admissible zone
    annotation = extractor.load_gtf(str(blocker), CONTIG, "")
    assert extractor.blocks_cut(CUT + 50, CUT + 500, CUT, MARGIN) is True
    islands = extractor.find_islands(annotation, CONTIG, "", MARGIN)
    zones = extractor.cut_zones(islands, 1, CONTIG_LEN)
    assert not any(lo <= CUT <= hi for lo, hi in zones), zones

    # AND THE SELECTOR HONOURS IT: it moves the cut rather than placing it where
    # a consumer would have to be warned.
    plan_path = tmp_path / "moved_cut_plan.json"
    _ok(_emit(tmp_path, genome, blocker, superset, plan_path, "selection_margin"))
    cuts = _cut_positions(json.loads(plan_path.read_text()))
    assert cuts[CONTIG], cuts
    assert CUT not in cuts[CONTIG], cuts
    for position in cuts[CONTIG]:
        assert not extractor.blocks_cut(CUT + 50, CUT + 500, position, MARGIN), position


def test_the_margin_advisory_counts_every_locus_it_warned_about(
    tmp_path, inputs, shared_plan
):
    """THE PER-RUN COUNT, on a plan carrying SEVERAL margin-only shortfalls.

    One warning per locus is individually benign and says nothing about the plan.
    The count is what does: a locus ending 3 bp from a cut loses most of that
    flank's reads, so dozens of these mean a plan badly matched to the model set
    applying it, and a bare per-locus stream would not surface that. Four loci
    here, two per cut, on both sides of both cuts, and zero straddles -- the
    straddle tally is asserted BECAUSE it is zero: that is what says the fatal
    guard ran rather than having been skipped.
    """

    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, plan, cuts = shared_plan
    assert [6000, 12000] == sorted(cuts[CONTIG]), cuts

    # (gene, lend, rend, cut, clearance): whole-in-left and whole-in-right of each
    # of the plan's two interior cuts
    loci = (
        ("after_first", 6050, 6500, 6000, 49),
        ("before_first", 5500, 5950, 6000, 50),
        ("after_second", 12050, 12400, 12000, 49),
        ("before_second", 11500, 11950, 12000, 50),
    )
    consumer = tmp_path / "four_margin_loci.gtf"
    consumer.write_text(
        gtf.read_text()
        + "".join(
            '{}\ttest\texon\t{}\t{}\t.\t+\t.\tgene_id "{}"; '
            'transcript_id "{}.t1";\n'.format(CONTIG, lend, rend, name, name)
            for name, lend, rend, _cut, _clearance in loci
        )
    )

    combined = _ok(
        _apply(
            tmp_path,
            genome,
            consumer,
            superset,
            plan_path,
            "four_margin",
            "--stop_after_make_chunks",
        )
    )
    workdir = tmp_path / "work_four_margin"

    assert (
        "margin audit: 4 annotated locus/loci sit within this run's --margin of "
        "{} bp of a plan cut, and 0 straddle one, over 2 contig(s)".format(MARGIN)
        in combined
    ), combined[-4000:]
    # every locus is named, and each with its own clearance and cut
    for name, lend, rend, cut, clearance in loci:
        assert name in combined, name
        assert "{} bp clearance".format(clearance) in combined or "clears by only {} bp".format(
            clearance
        ) in combined, name
    assert "{}:{}".format(CONTIG, 6000) in combined
    assert "{}:{}".format(CONTIG, 12000) in combined

    # and it survives into the run record, not just the terminal
    audit = _timing(workdir)["chunk_plan_margin_advisory"]
    assert audit["inside_margin"] == 4
    assert audit["straddling"] == 0
    assert audit["margin"] == MARGIN
    assert sorted(o["gene_id"] for o in audit["loci"]) == sorted(
        name for name, *_rest in loci
    )
    assert sorted((o["position"], o["gap"], o["side"]) for o in audit["loci"]) == [
        (6000, 49, "after"),
        (6000, 50, "before"),
        (12000, 49, "after"),
        (12000, 50, "before"),
    ]
    # the run really did extract, at the plan's geometry
    planned = [(c["chunk_id"], c["region"]) for c in plan["chunks"]]
    assert [
        (row[0], "{}:{}-{}".format(row[1], row[2], row[3])) for row in _geometry(workdir)
    ] == planned


def test_a_per_contig_slice_of_the_plans_annotation_is_accepted(
    tmp_path, inputs, shared_plan
):
    """THE by_chromosome CASE. Its --gtf is a slice, not the whole annotation.

    `WDL/LRAA.wdl` hands each by_chromosome shard `splitByChr.chromosomeGTFs[i]`,
    written by `util/partition_data_by_chromosome.py` as that contig's records
    verbatim. A whole-file digest of the annotation refused every such shard of a
    correct plan. The straddle check never asks what the file is, only whether a
    model in it spans a cut on a contig this run processes, so a slice and the
    whole gtf give the same answer for the contig they share.
    """

    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, plan, _cuts = shared_plan
    slice_gtf = tmp_path / "chr1_only.gtf"
    slice_gtf.write_text(
        "# no header of its own\n"
        + "".join(
            line + "\n"
            for line in gtf.read_text().splitlines()
            if line.startswith(CONTIG + "\t")
        )
    )
    assert slice_gtf.read_text() != gtf.read_text()

    _ok(
        _apply(
            tmp_path,
            genome,
            slice_gtf,
            superset,
            plan_path,
            "gtf_slice",
            "--contig",
            CONTIG,
            "--stop_after_make_chunks",
        )
    )
    geometry = _geometry(tmp_path / "work_gtf_slice")
    planned = [
        (c["chunk_id"], c["region"]) for c in plan["chunks"] if c["chrom"] == CONTIG
    ]
    assert [
        (row[0], "{}:{}-{}".format(row[1], row[2], row[3])) for row in geometry
    ] == planned


def test_an_annotation_differing_off_this_runs_contigs_is_accepted(
    tmp_path, inputs, shared_plan
):
    """Contig N's cuts depend on contig N's records and nothing else.

    The distinction that makes the per-contig scope honest rather than merely
    permissive: a model straddling a cut ON a processed contig is refused
    (above), a difference on one this run never touches is not looked at -- a run
    cannot sever a model it does not quantify.
    """

    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, _plan, _cuts = shared_plan
    edited = tmp_path / "other_contig_edited.gtf"
    edited.write_text(
        gtf.read_text()
        + '{}\ttest\texon\t5000\t5400\t.\t+\t.\tgene_id "novel_other"; '
        'transcript_id "novel_other.t1";\n'.format(OTHER_CONTIG)
    )

    _ok(
        _apply(
            tmp_path,
            genome,
            edited,
            superset,
            plan_path,
            "off_contig",
            "--contig",
            CONTIG,
            "--stop_after_make_chunks",
        )
    )


def test_a_de_novo_run_may_apply_a_plan_selected_around_an_annotation(
    tmp_path, inputs, shared_plan
):
    """Formerly refused on presence alone; now accepted, and correctly.

    `annotation_present` made "no gtf" a state that matched only itself, so this
    pair rejected each other before any contig was compared. Nothing is severed
    here: the plan's cuts avoided the reference's loci, and a run with no
    annotation has no model to lose. The only real hazard is a model spanning a
    cut, and there are no models.
    """

    genome, _gtf, superset, _reads, _lengths = inputs
    plan_path, _plan, _cuts = shared_plan

    _ok(
        _apply(
            tmp_path,
            genome,
            None,
            superset,
            plan_path,
            "no_gtf",
            "--discovery",
            "--stop_after_make_chunks",
        )
    )


@pytest.fixture
def denovo_plan(tmp_path, inputs):
    """ONE plan for a whole de novo run, emitted with NO annotation.

    The shape `WDL/TerraWorkflowConfigs/LRAA_singlecell.cluster_guided.by_chunk.
    de_novo.config.json` runs: geometry is fixed before phase 1, when there is no
    reference to fix it around, while every consumer downstream has a gtf -- the
    init GTF for per-cluster discovery, the consolidated one for final quant.
    """

    genome, _gtf, superset, _reads, _lengths = inputs
    plan_path = tmp_path / "denovo_cut_plan.json"
    _ok(
        _emit(
            tmp_path, genome, None, superset, plan_path, "denovo_emit", "--discovery"
        )
    )
    plan = json.loads(plan_path.read_text())
    cuts = _cut_positions(plan)
    assert cuts[CONTIG], "the de novo fixture placed no cut, so nothing is tested"
    return plan_path, plan, cuts


def test_an_unannotated_plan_is_accepted_by_a_consumer_holding_an_annotation(
    tmp_path, inputs, denovo_plan
):
    """THE DE NOVO SINGLE-CELL SHAPE, which presence equality failed outright.

    The plan carries no annotation because none existed when geometry was fixed;
    the consumer has one because discovery produced it. Accepted, because no
    model in the consumer's gtf spans a cut -- which is the only question that
    matters, and the one a presence boolean could not ask.
    """

    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, plan, _cuts = denovo_plan

    _ok(
        _apply(
            tmp_path,
            genome,
            gtf,
            superset,
            plan_path,
            "denovo_consumer",
            "--stop_after_make_chunks",
        )
    )
    geometry = _geometry(tmp_path / "work_denovo_consumer")
    planned = [(c["chunk_id"], c["region"]) for c in plan["chunks"]]
    assert [
        (row[0], "{}:{}-{}".format(row[1], row[2], row[3])) for row in geometry
    ] == planned


def test_an_unannotated_plan_is_refused_when_the_consumers_gtf_straddles_a_cut(
    tmp_path, inputs, denovo_plan
):
    """A guard that cannot fail is worse than no guard.

    Same shape as the acceptance above, one model moved. Selecting without an
    annotation protects nothing from severing, so a consumer's model CAN land
    across a cut -- and this is the case that proves the de novo path is checked
    rather than waved through. Built from the plan's own chosen position so the
    test cannot degrade into one where nothing straddles.
    """

    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, _plan, cuts = denovo_plan
    position = cuts[CONTIG][0]
    straddling = tmp_path / "denovo_consumer_straddle.gtf"
    straddling.write_text(
        gtf.read_text()
        + '{}\ttest\texon\t{}\t{}\t.\t+\t.\tgene_id "spans"; '
        'transcript_id "spans.t1";\n'.format(CONTIG, position - 200, position + 200)
    )

    combined = _refused(
        _apply(
            tmp_path, genome, straddling, superset, plan_path, "denovo_straddle"
        )
    )
    assert "transcript spans.t1 of gene spans" in combined
    assert "places a cut at {}:{}".format(CONTIG, position) in combined
    assert "{}:{}-{}".format(CONTIG, position - 200, position + 200) in combined


def test_chunk_plan_refuses_a_cut_selection_source(tmp_path, inputs, shared_plan):
    """One says where the cuts are, the other says where to choose them from."""

    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, _plan, _cuts = shared_plan
    combined = _refused(
        _apply(
            tmp_path,
            genome,
            gtf,
            superset,
            plan_path,
            "both_sources",
            "--bam_for_cut_selection",
            superset,
        )
    )
    assert "--chunk_plan supplies the cut positions" in combined


def test_emitting_and_applying_at_once_is_refused(tmp_path, inputs, shared_plan):
    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, _plan, _cuts = shared_plan
    combined = _refused(
        _apply(
            tmp_path,
            genome,
            gtf,
            superset,
            plan_path,
            "emit_and_apply",
            "--emit_cut_plan",
            tmp_path / "again.json",
        )
    )
    assert "--emit_cut_plan produces a partition" in combined


def test_emitting_and_extracting_at_once_is_refused(tmp_path, inputs):
    """--emit_cut_plan exists precisely to NOT extract."""

    genome, gtf, superset, _reads, _lengths = inputs
    combined = _refused(
        _emit(
            tmp_path,
            genome,
            gtf,
            superset,
            tmp_path / "never.json",
            "clash",
            "--stop_after_make_chunks",
        )
    )
    assert "stops BEFORE extraction" in combined


def test_a_plan_may_not_be_the_leaf_plan_path(tmp_path, inputs, shared_plan):
    """--stop_after_make_chunks writes that path; a shared plan is read by all."""

    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, _plan, _cuts = shared_plan
    workdir = tmp_path / "work_clobber"
    workdir.mkdir()
    target = workdir / "chunk_plan.json"
    target.write_text(plan_path.read_text())

    result = _run(
        tmp_path,
        workdir,
        genome,
        gtf,
        superset,
        "--chunk_plan",
        target,
        num_total_reads=8,
    )
    combined = _refused(result)
    assert "leaf plan" in combined


def test_strand_first_may_not_share_a_plan(tmp_path, inputs, shared_plan):
    """Its selection reads the stage-1 split bams, which are per-run files."""

    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, _plan, _cuts = shared_plan
    combined = _refused(
        _apply(
            tmp_path,
            genome,
            gtf,
            superset,
            plan_path,
            "strand_first",
            "--chunk_by_strand",
        )
    )
    assert "strandless-chunking only" in combined


def test_a_geometry_only_plan_is_refused_by_a_leaf(tmp_path, inputs, shared_plan):
    """--only_chunk needs extracted directories, which an emit run never wrote."""

    _genome, _gtf, _superset, _reads, _lengths = inputs
    plan_path, plan, _cuts = shared_plan
    outdir = tmp_path / "leaf"
    (outdir / "logs").mkdir(parents=True)
    (outdir / "chunk_plan.json").write_text(plan_path.read_text())

    args = ChunkedRun.default_args(
        output_dir=str(outdir), only_chunk=plan["chunks"][0]["chunk_id"]
    )
    with pytest.raises(ChunkedRun.PipelineError) as excinfo:
        ChunkedRun.run(args)
    assert "chunks_extracted=false" in str(excinfo.value)


def test_a_contig_absent_from_this_runs_bam_is_refused(tmp_path, inputs, shared_plan):
    """Checked at the validator, because enumeration filters such a contig out.

    A bam cannot hold a record against a reference its own header does not list,
    so a run asked to process one on a plan's geometry has nothing to extract.
    """

    genome, gtf, _superset, reads, _lengths = inputs
    plan_path, plan, _cuts = shared_plan
    one_contig = _write_bam(
        tmp_path / "chr1_only.bam",
        [r for r in reads if r[1] == CONTIG],
        lengths=((CONTIG, CONTIG_LEN),),
    )
    args = ChunkedRun.default_args(
        bam=str(one_contig),
        genome_fa=str(genome),
        gtf=str(gtf),
        output_dir=str(tmp_path / "unused"),
        chunk_plan=str(plan_path),
        approx_MB_per_cut=MB_PER_CUT,
        approx_MB_per_cut_wiggle_window=WIGGLE,
    )
    loaded = ChunkedRun.load_chunk_plan(str(plan_path))
    with pytest.raises(ChunkedRun.PipelineError) as excinfo:
        ChunkedRun.selections_from_chunk_plan(
            loaded, str(plan_path), args, [CONTIG, OTHER_CONTIG]
        )
    message = str(excinfo.value)
    assert OTHER_CONTIG in message
    assert "absent from the bam header" in message


# ------------------------------------------------------------ LRAA's own surface


def test_lraa_chunk_forwards_chunk_plan_instead_of_ignoring_it(tmp_path, inputs):
    """LRAA builds its chunk args from an ALLOWLIST.

    An unlisted flag is parsed and silently discarded, which is how `--bam_for_sg`
    came to be ignored on every chunked run. Dropped here, every caller would go
    back to selecting its own geometry with nothing to say so, since a plan that
    never arrives cannot mismatch. The plan named below does not exist, so the
    proof is that ChunkedRun's own refusal is what comes back.
    """

    genome, gtf, superset, _reads, _lengths = inputs
    missing = tmp_path / "absent_plan.json"
    result = subprocess.run(
        [
            sys.executable,
            str(LRAA),
            "--genome",
            str(genome),
            "--bam",
            str(superset),
            "--gtf",
            str(gtf),
            "--quant_only",
            "--output_prefix",
            str(tmp_path / "lraa_out"),
            "--chunk",
            "--chunk_work_dir",
            str(tmp_path / "lraa_work"),
            "--chunk_plan",
            str(missing),
        ],
        capture_output=True,
        text=True,
        cwd=str(tmp_path),
        timeout=3600,
    )
    combined = result.stdout + result.stderr
    assert result.returncode != 0, combined[-8000:]
    assert "chunk plan {} does not exist".format(missing) in combined


def test_a_bam_header_longer_than_the_fasta_still_takes_the_plan(
    tmp_path, inputs, shared_plan
):
    """Align against the whole assembly, analyse against a slice of it.

    A supported configuration, and the shape the single-cell fixture actually
    has: reads aligned to GRCh38 chr19 carry a 58,617,616 bp header line, while
    the run's genome fasta holds a 2,000,000 bp slice of that contig. Chunk
    geometry is decided by the FASTA -- the last segment runs to the length the
    extractor fetches sequence for, and that fetch is a fasta fetch -- so the
    header's disagreement is not a geometry disagreement.

    Checking the plan against the header instead refused every shard of the
    cluster-guided run for a mismatch that predates the plan and that the
    unchunked path tolerates. The presence check stays: a bam whose header
    omits the contig entirely holds nothing to extract.
    """

    genome, gtf, superset, reads, _lengths = inputs
    plan, _plan_doc, _cuts = shared_plan

    # Same reads, same coordinates, same fasta -- only the header's LN is the
    # full-assembly value, which is what an aligner against GRCh38 writes.
    oversized = _write_bam(
        tmp_path / "full_header.bam",
        reads,
        lengths=((CONTIG, CONTIG_LEN * 2929), (OTHER_CONTIG, OTHER_LEN * 2929)),
    )

    result = _apply(tmp_path, genome, gtf, oversized, plan, "oversized")
    _ok(result)
    # The geometry is the PLAN's, and its last boundary sits at the FASTA length.
    # A run that had re-derived bounds from the header would end its final chunk
    # 2929x further out, and a run that had reselected locally would not match a
    # sibling applying the same plan to the superset bam.
    got = _geometry(tmp_path / "work_oversized")
    control = _apply(tmp_path, genome, gtf, superset, plan, "control")
    _ok(control)
    assert got == _geometry(tmp_path / "work_control")

    ends = [rend for _cid, chrom, _lend, rend, *_ in got if chrom == CONTIG]
    assert max(ends) == CONTIG_LEN, (ends, CONTIG_LEN)
