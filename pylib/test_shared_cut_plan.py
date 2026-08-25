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
  values named -- geometry parameters, annotation identity, contig length, contig
  coverage, and the contradictory flag combinations. A mismatched plan does not
  fail downstream; it produces plausible output at bounds nobody chose;
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
    assert params["annotation_present"] is True
    # WHICH annotation is per contig, not a whole-file hash: by_chromosome hands
    # each shard a slice of the same gtf.
    assert "annotation_identity" not in params
    assert set(plan["geometry"]["annotation"]) == {CONTIG, OTHER_CONTIG}
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


def test_a_different_annotation_on_a_processed_contig_is_refused(
    tmp_path, inputs, shared_plan
):
    """A cut inside an annotated model is inadmissible, so the gtf is geometry.

    The single-cell pipeline selects against the CONSOLIDATED gtf, which carries
    the novel models phase-1 discovery found. A plan selected against the raw
    reference can place a cut straight through one of those, and every cluster
    then quantifies a severed model and reports success -- wrong exactly at the
    novel isoforms discovery exists to find. The novel model below sits across the
    6000 cut this fixture's plan chose.
    """

    genome, gtf, superset, _reads, _lengths = inputs
    plan_path, _plan, _cuts = shared_plan
    consolidated = tmp_path / "consolidated.gtf"
    consolidated.write_text(
        gtf.read_text()
        + '{}\ttest\texon\t5800\t6200\t.\t+\t.\tgene_id "novel"; '
        'transcript_id "novel.t1";\n'.format(CONTIG)
    )

    combined = _refused(
        _apply(tmp_path, genome, consolidated, superset, plan_path, "other_gtf")
    )
    assert "selected contig {}'s cuts against annotation".format(CONTIG) in combined
    assert "severed model" in combined


def test_a_per_contig_slice_of_the_plans_annotation_is_accepted(
    tmp_path, inputs, shared_plan
):
    """THE by_chromosome CASE. Its --gtf is a slice, not the whole annotation.

    `WDL/LRAA.wdl` hands each by_chromosome shard `splitByChr.chromosomeGTFs[i]`,
    written by `util/partition_data_by_chromosome.py` as that contig's records
    verbatim. A whole-file identity refused every such shard of a correct plan,
    which is why the annotation identity is per contig.
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
    permissive: a difference ON a processed contig is refused (above), a
    difference on one this run never touches is not.
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


def test_no_annotation_does_not_match_a_gtf_bearing_plan(
    tmp_path, inputs, shared_plan
):
    """Absent is a DISTINCT state from "some gtf", not a wildcard.

    Refused on presence alone, before any contig is compared, so a de novo run and
    a plan selected around an annotation reject each other even with no contig in
    common.
    """

    genome, _gtf, superset, _reads, _lengths = inputs
    plan_path, _plan, _cuts = shared_plan
    combined = _refused(
        _apply(
            tmp_path,
            genome,
            None,
            superset,
            plan_path,
            "no_gtf",
            "--discovery",
        )
    )
    assert "annotation_present" in combined
    assert "state of its own" in combined


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
