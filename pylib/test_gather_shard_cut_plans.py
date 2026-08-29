#!/usr/bin/env python3

"""A scattered run's own cut geometry, gathered rather than re-selected.

THE GAP THIS CLOSES. `--emit_cut_plan` produces a shareable partition by running a
whole extra cut selection, in a task the run then waits on. That is right for the
cluster-guided shape, where ONE selection serves ~29 sibling runs. It is wrong for
a run that is going to chunk anyway: selection inside `make_chunks` OVERLAPS
extraction -- `run_prep_concurrently` submits a contig's extractions the moment
that contig's selection finishes -- and moving it out into a separate task trades
that overlap for a file. So a basic single-cell run, whose initial pass is the only
LRAA run in the pipeline, emitted no plan at all and its geometry could not be
handed to anything, even though it had chosen a perfectly good one.

Every real chunked run now writes `<work>/cuts/shard_cut_plan.json`
(`ChunkedRun.SHARD_CUT_PLAN_NAME`): the geometry it chose, in the format
`--chunk_plan` consumes, restricted to the contigs THAT run processed.
`util/misc/gather_shard_cut_plans.py` merges the per-shard files of a scatter into
one genome-wide plan.

What these tests hold:

* the shard file carries the selector's own per-contig records VERBATIM -- the same
  objects `<contig>.strandless.cuts.json` holds -- plus the parameters those files
  do not record and `validate_cut_plan_geometry` requires;
* it is written by an ordinary chunked run, not only by the flags whose purpose is
  to produce a plan;
* the gathered envelope has EXACTLY ONE source, keyed with the empty string.
  `ChunkedRun.validate_cut_plan_geometry` refuses anything else, so a gather that
  concatenated the shards' own `by_source` lists would build a file its own
  consumers reject;
* disagreeing geometry is REFUSED with the parameter and both values named, and
  agreeing geometry passes. Taking the first shard's silently would produce a plan
  whose recorded parameters do not match part of what produced it -- which passes
  the consuming run's check and then extracts at bounds selected under different
  parameters, i.e. the failure that check exists to catch, reporting success;
* the ROUND TRIP: a plan gathered from two per-contig shards is accepted by a
  genome-wide run at that geometry and reproduces the shards' own chunk bounds;
* a contig claimed by two shards, and strand-first geometry, are both refused --
  the first because shards are meant to be a partition, the second because
  strand-first selects over per-run stage-1 split bams and is not shareable at all.
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

REPO = Path(__file__).resolve().parents[1]
CHUNKED_RUN = REPO / "pylib" / "ChunkedRun.py"
GATHER = REPO / "util" / "misc" / "gather_shard_cut_plans.py"
WDL_DIR = REPO / "WDL"
RUNNER_WDL = WDL_DIR / "subwdls" / "LRAA_runner.wdl"
LRAA_WDL = WDL_DIR / "LRAA.wdl"
SINGLECELL_WDL = WDL_DIR / "LRAA-singlecell.wdl"

CONTIG = "chr1"
CONTIG_LEN = 20000
OTHER_CONTIG = "chr2"
OTHER_LEN = 12000

# 6 kb nominal spacing on the 20 kb contig and 1 kb of search around each target,
# matching test_shared_cut_plan.py's fixture: three chunks on chr1, so the
# partition has two interior boundaries for a gathered plan to agree or disagree
# about. One chunk per contig would make every geometry assertion vacuous.
MB_PER_CUT = 0.006
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


def _write_bam(path, reads, lengths):
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
            fh.write(aln)
    pysam.index(str(path))
    return Path(path)


def _dense_reads(prefix, lend, rend, contig, step=50, length=300):
    """Cover ``[lend, rend]`` so every candidate position in it is spanned.

    Severing is priced rather than forbidden, so a fully covered window still gets
    a cut -- which is the case worth partitioning on, since a geometry whose cuts
    sever nothing says little about whether two runs agree.
    """

    return [
        ("{}{}".format(prefix, i), contig, start, length, bool(i % 2))
        for i, start in enumerate(range(lend, rend, step))
    ]


LENGTHS = ((CONTIG, CONTIG_LEN), (OTHER_CONTIG, OTHER_LEN))


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    """A two-contig genome, an annotation clear of every cut window, and a bam.

    Module-scoped: the shard runs below are read-only consumers of it, and building
    it per test would pay for the same fasta index and bam index in every one.

    The annotated loci sit far from 6000/12000 on chr1 and 6000 on chr2, because a
    cut is inadmissible within ``chunk_margin`` of an annotated locus -- a gene in
    a search window would decide the geometry by blocking rather than by depth,
    which is a different thing to be testing.
    """

    root = tmp_path_factory.mktemp("gather_inputs")
    genome = root / "genome.fa"
    genome.write_text(
        ">{}\n{}\n>{}\n{}\n".format(
            CONTIG, "A" * CONTIG_LEN, OTHER_CONTIG, "A" * OTHER_LEN
        )
    )
    pysam.faidx(str(genome))

    gtf = root / "annot.gtf"
    gtf.write_text(
        "\n".join(
            '{}\ttest\texon\t{}\t{}\t.\t{}\t.\tgene_id "{}"; transcript_id "{}";'.format(
                contig, lend, rend, strand, gene, gene + ".t1"
            )
            for gene, contig, lend, rend, strand in (
                ("gplus", CONTIG, 1001, 1400, "+"),
                ("gminus", CONTIG, 16001, 16400, "-"),
                ("gother", OTHER_CONTIG, 1001, 1400, "+"),
            )
        )
        + "\n"
    )

    reads = (
        _dense_reads("c1_cut1_", 5000, 7000, CONTIG)
        + _dense_reads("c1_cut2_", 11000, 13000, CONTIG)
        + _dense_reads("c1_gene_plus_", 1000, 1400, CONTIG, step=100)
        + _dense_reads("c1_gene_minus_", 16000, 16400, CONTIG, step=100)
        + _dense_reads("c2_cut1_", 5000, 7000, OTHER_CONTIG)
        + _dense_reads("c2_gene_", 1000, 1400, OTHER_CONTIG, step=100)
    )
    bam = _write_bam(root / "reads.bam", reads, LENGTHS)
    return genome, gtf, bam


# ---------------------------------------------------------------------- helpers


def _chunked_run(workdir, inputs, *extra):
    genome, gtf, bam = inputs
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
        "--approx_MB_per_cut",
        str(MB_PER_CUT),
        "--approx_MB_per_cut_wiggle_window",
        str(WIGGLE),
    ]
    cmd += [str(x) for x in extra]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    combined = result.stdout + result.stderr
    assert result.returncode == 0, combined[-8000:]
    return workdir


def _shard(tmp_path, inputs, contig, label, *extra):
    """One scattered shard: a chunked run restricted to one contig.

    Stopped after make_chunks, which is the phase that chooses and records the
    geometry; everything after it is quantification and would only make these
    tests slower. That the file is ALSO written by a run that does not stop there
    is asserted separately, because that is the run the WDL lifts it out of.
    """

    workdir = _chunked_run(
        tmp_path / ("work_" + label),
        inputs,
        "--contig",
        contig,
        "--num_total_reads",
        "1000",
        "--stop_after_make_chunks",
        *extra,
    )
    plan = workdir / "cuts" / ChunkedRun.SHARD_CUT_PLAN_NAME
    assert plan.exists(), sorted(p.name for p in (workdir / "cuts").iterdir())
    return plan


def _gather(output, *plans):
    result = subprocess.run(
        [
            sys.executable,
            str(GATHER),
            "--shard_plans",
            *[str(p) for p in plans],
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
        timeout=600,
    )
    return result, result.stdout + result.stderr


def _gathered(output, *plans):
    result, combined = _gather(output, *plans)
    assert result.returncode == 0, combined[-8000:]
    return json.loads(Path(output).read_text())


def _refused(output, *plans):
    result, combined = _gather(output, *plans)
    assert result.returncode != 0, combined[-8000:]
    return combined


def _mutated(source, dest, **params):
    """A shard plan with named geometry parameters replaced.

    Mutated on disk rather than produced by a second run at different settings:
    what is under test is the refusal, and editing one key isolates it from
    whether two spacings happen to place different cuts.
    """

    plan = json.loads(Path(source).read_text())
    for name, value in params.items():
        if value is None:
            plan["geometry"]["params"].pop(name, None)
        else:
            plan["geometry"]["params"][name] = value
    Path(dest).write_text(json.dumps(plan, indent=2, sort_keys=True))
    return Path(dest)


@pytest.fixture(scope="module")
def shards(tmp_path_factory, inputs):
    """One shard per contig, at IDENTICAL geometry. The negative control.

    Everything below that tests a refusal starts from these, so the passing case
    and the refused case differ in exactly the thing being refused.
    """

    root = tmp_path_factory.mktemp("gather_shards")
    return (
        _shard(root, inputs, CONTIG, "c1"),
        _shard(root, inputs, OTHER_CONTIG, "c2"),
    )


# ----------------------------------------------------- what a shard actually says


def test_a_shard_carries_the_selector_records_verbatim(shards):
    """The shard plan IS the cut records, plus what makes them applicable.

    ``<contig>.strandless.cuts.json`` is the selector's own output and is what a
    consumer needs the positions from -- but it records only the parameters the
    selector was told about its own search, not ``margin``,
    ``max_intron_length`` or the resolved ``min_per_id`` /
    ``min_mapping_quality``, all of which ``validate_cut_plan_geometry``
    requires. So the shard plan holds the records unchanged and states the rest
    alongside them.
    """

    plan = json.loads(shards[0].read_text())
    selections = plan["geometry"]["by_source"][0]["selections"]
    cuts_json = shards[0].parent / "{}.{}.cuts.json".format(
        CONTIG, ChunkedRun.STRANDLESS_TAG
    )
    assert selections == json.loads(cuts_json.read_text())

    params = plan["geometry"]["params"]
    assert set(params) >= {
        "margin",
        "max_intron_length",
        "min_per_id",
        "min_mapping_quality",
        "strandless",
    }
    assert plan["geometry"]["mode"] == ChunkedRun.STRANDLESS_MODE


def test_a_shard_states_no_denominator_and_no_extracted_chunks(shards):
    """Geometry only, exactly as ``--emit_cut_plan`` records it.

    The TPM denominator belongs to whoever quantifies, and per-shard denominators
    are per-shard by construction. The chunk directories the shard built are its
    own working files at its own paths, so claiming them would offer a consuming
    ``--only_chunk`` leaf a partition it cannot address.
    """

    plan = json.loads(shards[0].read_text())
    assert plan["num_total_reads"] is None
    assert plan["chunks_extracted"] is False


def test_an_ordinary_chunked_run_writes_it_too(tmp_path, inputs):
    """Not gated on a flag. This is the run a WDL task lifts the file out of.

    ``--stop_after_make_chunks`` and ``--emit_cut_plan`` exist to produce a plan;
    the point of this file is that a run whose purpose is to QUANTIFY also records
    what it partitioned into, because otherwise its geometry is unrecoverable.
    """

    workdir = _chunked_run(
        tmp_path / "work_full",
        inputs,
        "--contig",
        OTHER_CONTIG,
        "--num_total_reads",
        "1000",
    )
    plan_path = workdir / "cuts" / ChunkedRun.SHARD_CUT_PLAN_NAME
    assert plan_path.exists()
    plan = json.loads(plan_path.read_text())
    assert [s["chrom"] for s in plan["geometry"]["by_source"][0]["selections"]] == [
        OTHER_CONTIG
    ]
    # ... and the run still did its real job, so recording the geometry did not
    # displace the merge.
    assert (workdir / "merged").is_dir()


def test_a_single_chunk_shard_still_records_its_geometry(tmp_path, inputs):
    """The case that decides whether the WDL output can be non-optional.

    At the default 10 Mb spacing a 12 kb contig gets no cut at all and chunking
    degenerates to one whole-contig chunk. Selection still RUNS, and still records
    zero cuts and one segment -- so the file is present for every chunked run
    however coarse the geometry, and the only case that writes none is a run that
    never chunks at all.
    """

    workdir = _chunked_run(
        tmp_path / "work_one",
        inputs,
        "--contig",
        OTHER_CONTIG,
        "--num_total_reads",
        "1000",
        "--stop_after_make_chunks",
        "--approx_MB_per_cut",
        "10",
        "--approx_MB_per_cut_wiggle_window",
        "1",
    )
    plan = json.loads(
        (workdir / "cuts" / ChunkedRun.SHARD_CUT_PLAN_NAME).read_text()
    )
    (selection,) = plan["geometry"]["by_source"][0]["selections"]
    assert selection["cuts"] == []
    assert len(selection["segments"]) == 1
    assert [c["chunk_id"] for c in plan["chunks"]] == ["{}_00".format(OTHER_CONTIG)]


# ------------------------------------------------------------------- the envelope


def test_the_gathered_envelope_has_one_strandless_source(tmp_path, shards):
    """The constraint the gather exists to respect.

    ``validate_cut_plan_geometry`` refuses a plan that does not have exactly one
    source keyed "", so concatenating N shards' ``by_source`` lists would produce
    a file rejected by its own consumers. Every shard's contigs are merged UNDER
    one source instead.
    """

    plan = _gathered(tmp_path / "gathered.json", *shards)
    by_source = plan["geometry"]["by_source"]
    assert len(by_source) == 1
    assert by_source[0]["key"] == ""
    assert by_source[0]["tag"] == ChunkedRun.STRANDLESS_TAG
    assert [s["chrom"] for s in by_source[0]["selections"]] == sorted(
        [CONTIG, OTHER_CONTIG]
    )


def test_the_gathered_source_claims_no_single_bam(tmp_path, shards):
    """There is no one file this geometry was selected from, so it names none.

    A one-run plan records the bam its cuts were chosen on as provenance. Across
    shards each has its own, so stating one of them for a genome-wide plan would
    be a false claim; the contributing shards are named under ``gathered_from``
    instead.
    """

    plan = _gathered(tmp_path / "gathered.json", *shards)
    source = plan["geometry"]["by_source"][0]
    assert source["bam"] is None
    assert source["bam_identity"] is None
    assert {
        contig
        for entry in plan["gathered_from"]
        for contig in entry["contigs"]
    } == {CONTIG, OTHER_CONTIG}
    assert len(plan["gathered_from"]) == 2


def test_the_gathered_chunks_carry_one_global_order(tmp_path, shards):
    """Ids, regions and ``order`` as a single genome-wide run would assign them.

    Per-shard ``order`` counters each restart at zero, so a concatenation would
    state every value twice. Rebuilt through ``ChunkedRun.planned_chunks`` over the
    merged selections in sorted-contig order -- which is the order
    ``enumerate_prep_contigs`` returns -- so the plan describes a partition LRAA
    would actually build.
    """

    plan = _gathered(tmp_path / "gathered.json", *shards)
    orders = [c["order"] for c in plan["chunks"]]
    assert orders == list(range(len(orders)))
    assert [c["chrom"] for c in plan["chunks"]] == sorted(
        [c["chrom"] for c in plan["chunks"]]
    )
    per_shard = [
        c["chunk_id"]
        for shard in shards
        for c in json.loads(Path(shard).read_text())["chunks"]
    ]
    assert sorted(c["chunk_id"] for c in plan["chunks"]) == sorted(per_shard)


# ------------------------------------------------------------------- the refusals


@pytest.mark.parametrize(
    "param,value",
    [
        ("approx_MB_per_cut", 0.02),
        ("approx_MB_per_cut_wiggle_window", 0.5),
        ("margin", 999),
        ("max_intron_length", 12345),
        ("min_per_id", 97.0),
        ("min_mapping_quality", 30),
        ("severed_multiexon_weight", 3.0),
        ("depth_window", 250),
        # Absence is a disagreement too: the key set is what decided a cut, so a
        # shard that omits one was not selecting on the same rules, and a gathered
        # plan missing it would be refused by the consumer with nothing to check.
        ("margin", None),
    ],
)
def test_the_gather_refuses_shards_whose_geometry_disagrees(
    tmp_path, shards, param, value
):
    """Named, with both values, rather than resolved by taking the first.

    A plan whose recorded parameters do not match what produced part of it PASSES
    ``validate_cut_plan_geometry`` in a consuming run -- the check compares the
    recorded values -- and the run then extracts at bounds selected under
    different ones. That is the exact failure the check exists to catch, arriving
    with the check reporting success, so the only place it can be caught is here.
    """

    odd = _mutated(shards[1], tmp_path / "odd.json", **{param: value})
    combined = _refused(tmp_path / "gathered.json", shards[0], odd)
    assert param in combined
    assert str(shards[0]) in combined and str(odd) in combined


def test_shards_that_agree_are_gathered(tmp_path, shards):
    """The negative control for the refusals above.

    Without it the parametrized cases are satisfied by a gather that refuses
    everything.
    """

    plan = _gathered(tmp_path / "gathered.json", *shards)
    assert len(plan["geometry"]["by_source"][0]["selections"]) == 2
    # The parameters the shards agreed on are the ones recorded, unchanged.
    assert (
        plan["geometry"]["params"]
        == json.loads(shards[0].read_text())["geometry"]["params"]
    )


def test_the_gather_refuses_a_contig_claimed_by_two_shards(tmp_path, shards):
    """Shards are a partition. The same contig in two of them means they are not.

    Deduplicating would pick one shard's partition of that contig arbitrarily, and
    the two need not agree -- cut placement reads a bam, and two shards holding one
    contig are two bams.
    """

    twin = Path(shutil.copy(str(shards[0]), str(tmp_path / "twin.json")))
    combined = _refused(tmp_path / "gathered.json", shards[0], twin)
    assert CONTIG in combined
    assert "one partition" in combined


def test_the_gather_refuses_strand_first_geometry(tmp_path, inputs):
    """Strand-first cuts are selected over per-run stage-1 split bams.

    Not shareable at any number of shards, and ``validate_cut_plan_geometry``
    refuses it downstream anyway. Said here so the message names the shard rather
    than the merged file.
    """

    shard = _shard(tmp_path, inputs, CONTIG, "strandfirst", "--chunk_by_strand")
    combined = _refused(tmp_path / "gathered.json", shard)
    assert ChunkedRun.STRAND_FIRST_MODE in combined
    assert str(shard) in combined


def test_the_gather_refuses_a_plan_with_no_geometry(tmp_path, shards):
    """A file that looks like a plan and states no cut positions.

    ``--only_chunk``'s leaf plans predate the geometry block, and a hand-written
    file is easy to get wrong; either way there is nothing to gather and the
    merged output must not silently omit that shard's contigs.
    """

    stripped = json.loads(shards[0].read_text())
    stripped.pop("geometry")
    path = tmp_path / "no_geometry.json"
    path.write_text(json.dumps(stripped))
    combined = _refused(tmp_path / "gathered.json", path)
    assert "no geometry block" in combined


# ---------------------------------------------------------------- the round trip


def test_a_gathered_plan_is_accepted_by_a_run_at_that_geometry(
    tmp_path, inputs, shards
):
    """THE ACCEPTANCE PROPERTY: what the pipeline emitted, the pipeline consumes.

    A genome-wide run given the gathered plan must pass both refusals a shared
    plan faces -- ``validate_cut_plan_geometry`` on the parameters, and
    ``cut_blocking_annotation_models`` on this run's own gtf -- and must then
    extract at the SAME bounds the two shards chose, contig for contig. A plan that
    cannot be consumed by the pipeline that produced it is worth nothing, and a
    plan that is accepted but shifts the bounds is worse than nothing.
    """

    gathered = tmp_path / "gathered.json"
    _gathered(gathered, *shards)

    workdir = _chunked_run(
        tmp_path / "work_apply",
        inputs,
        "--num_total_reads",
        "1000",
        "--stop_after_make_chunks",
        "--chunk_plan",
        str(gathered),
    )
    applied = json.loads((workdir / "chunk_plan.json").read_text())
    shard_regions = {
        c["chunk_id"]: c["region"]
        for shard in shards
        for c in json.loads(Path(shard).read_text())["chunks"]
    }
    assert {c["chunk_id"]: c["region"] for c in applied["chunks"]} == shard_regions
    # Every chunk directory the plan named was really extracted, so the accepted
    # geometry was applied rather than merely validated.
    for chunk_id in shard_regions:
        assert (workdir / "chunks" / chunk_id / "chunk.partition.json").exists()


def test_a_gathered_plan_is_refused_at_different_geometry(tmp_path, inputs, shards):
    """The other half of the round trip: acceptance above is not unconditional.

    Without this, a validator that accepted everything would satisfy the test
    above.
    """

    gathered = tmp_path / "gathered.json"
    _gathered(gathered, *shards)
    result = subprocess.run(
        [
            sys.executable,
            str(CHUNKED_RUN),
            "--bam",
            str(inputs[2]),
            "--genome_fa",
            str(inputs[0]),
            "--gtf",
            str(inputs[1]),
            "--output_dir",
            str(tmp_path / "work_wrong"),
            "--cpu_budget",
            "2",
            "--num_total_reads",
            "1000",
            "--stop_after_make_chunks",
            "--chunk_plan",
            str(gathered),
            # The one thing changed: a spacing the plan was not selected at.
            "--approx_MB_per_cut",
            "0.012",
            "--approx_MB_per_cut_wiggle_window",
            str(WIGGLE),
        ],
        capture_output=True,
        text=True,
        timeout=3600,
    )
    combined = result.stdout + result.stderr
    assert result.returncode != 0, combined[-8000:]
    assert "approx_MB_per_cut" in combined


# ------------------------------------------------------------------ the wdl path
#
# Static, because none of these has a cheap dynamic check. A WDL output that names
# a file no task writes is legal WDL, passes ``miniwdl check`` and ``womtool
# validate``, and simply resolves to null -- which is indistinguishable from the
# case where there genuinely was no plan. That is the same silent-drop shape
# test_declared_inputs_reach_flags.py exists for, one layer further out.


def test_the_runner_lifts_out_the_file_chunkedrun_actually_writes():
    """The filename is ONE string, and the WDL has to spell the same one.

    ``ChunkedRun.SHARD_CUT_PLAN_NAME`` exists so that the run writing the file and
    the task copying it out cannot drift. A rename on the Python side with the WDL
    left alone produces a task that finds nothing, an optional output that resolves
    to null, and a workflow that reports no geometry -- with no error anywhere.
    """

    text = RUNNER_WDL.read_text()
    assert "cuts/{}".format(ChunkedRun.SHARD_CUT_PLAN_NAME) in text
    assert 'File? LRAA_shard_cut_plan = "' in text
    # ... and the subworkflow forwards the task output, which is the layer the six
    # dropped chunk knobs were lost at.
    assert (
        "File? LRAA_shard_cut_plan = LRAA_runner_task.LRAA_shard_cut_plan" in text
    )


def test_lraa_gathers_every_mode_that_places_cuts():
    """All three scattering modes reach the gather, and it reaches the output.

    by_chromosome is the case the gather exists for; ``off`` is one shard and
    by_chunk's make_chunks already writes this format. Any mode left out would
    surface as a null plan for that mode only, which is exactly the kind of gap
    that gets found in production rather than in a smoke run.
    """

    text = LRAA_WDL.read_text()
    assert "task gather_shard_cut_plans" in text
    assert "gather_shard_cut_plans.py" in text
    for producer in (
        "LRAA_scatter.LRAA_shard_cut_plan",
        "LRAA_direct.LRAA_shard_cut_plan",
        "chunk_scatter.chunkPlan",
    ):
        assert producer in text, producer
    # Gated on there being something to gather, so a no_chunk run does not set the
    # task up at all rather than running it on an empty array.
    assert "if (length(shardCutPlans) > 0)" in text
    assert (
        "File? gatheredChunkPlan = gather_shard_cut_plans.gatheredChunkPlan" in text
    )


def test_singlecell_prefers_a_supplied_plan_then_the_gathered_one():
    """Precedence, and that CLUSTER-GUIDED emission is untouched.

    A supplied ``chunk_plan`` is what every phase applied, so it is what a sibling
    run must be given. The gathered plan is the BASIC-mode fallback only: emitting
    one in basic mode would move cut selection out of the initial pass's own
    make_chunks, where it overlaps extraction, and gathering one in cluster-guided
    mode would report the init pass's geometry as the run's when only the init pass
    applied it.
    """

    text = SINGLECELL_WDL.read_text()
    assert "File? shared_chunk_plan = if defined(chunk_plan) then chunk_plan" in text
    assert "else if run_cluster_guided then emit_run_chunk_plan.chunk_plan" in text
    assert "else LRAA_init.gatheredChunkPlan" in text
    # The emission gate is unchanged: basic mode still sets up no plan task.
    assert "if (want_chunk_plan && !defined(chunk_plan)) {" in text
    assert (
        "Boolean want_chunk_plan = run_cluster_guided\n"
        "      && (scattering_per_cluster != \"off\" "
        "|| scattering_final_quant != \"off\")" in text
    )
