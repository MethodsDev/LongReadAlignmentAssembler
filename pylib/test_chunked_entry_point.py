#!/usr/bin/env python3

"""Tests for `LRAA --chunk`, the single entry point into the chunked pipeline.

Chunked mode used to be reachable only through util/misc/run_chunked_quant_pipeline.py,
which is how it came to be broken for two releases without anything noticing: the
driver was not on any routine path. These tests hold the merged entry point to the
three properties that make it safe to have one:

* chunked mode refuses the configurations it cannot serve, up front;
* a chunk worker cannot re-enter the pipeline;
* NOT passing --chunk leaves LRAA's own argument surface untouched;
* DISCOVERY is reachable, and pays for it with a hard zero-severed rule that a
  quant-only run must not inherit.
"""

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LRAA = REPO / "LRAA"

sys.path.insert(0, str(REPO / "pylib"))

import ChunkedRun  # noqa: E402
import LRAA_Globals  # noqa: E402


def _lraa(*argv, env=None):
    environ = dict(os.environ)
    if env:
        environ.update(env)
    return subprocess.run(
        [sys.executable, str(LRAA), *argv],
        capture_output=True,
        text=True,
        env=environ,
    )


BASE = (
    "--bam", "reads.bam",
    "--genome", "genome.fa",
    "--gtf", "annot.gtf",
    "--output_prefix", "sample",
)


def test_chunked_discovery_is_reachable(tmp_path):
    """The NEW contract: --chunk no longer requires --quant_only.

    It used to refuse outright, because discovery across a boundary can produce
    two partial models of one locus. That mechanism is real but it is caused by
    SEVERING, and severing is now forbidden outright in discovery rather than
    priced -- so the configuration is servable and is offered. See
    ``docs/denovo_chunking.md`` for the measurement the change rests on.

    Asserted by the run reaching the pipeline's FIRST act and dying on the bam
    that does not exist. That act used to be stage 1's strand split; strandless
    chunking is now the default and skips stage 1 entirely, so the marker is the
    library count the chunked path takes before any partitioning. The principle is
    unchanged and is the reason a marker is used at all: a failure only the
    pipeline can produce cannot be confused with the guard still firing under a
    different message.
    """

    result = _lraa(
        "--chunk",
        "--bam", "reads.bam",
        "--genome", "genome.fa",
        "--output_prefix", str(tmp_path / "sample"),
        "--chunk_work_dir", str(tmp_path / "work"),
    )
    combined = result.stdout + result.stderr
    assert result.returncode != 0
    assert "--chunk requires --quant_only" not in combined
    assert "--chunk requires --gtf" not in combined
    assert "counting genome-mapped reads" in combined


def test_chunked_ref_guided_discovery_is_reachable(tmp_path):
    """--gtf WITHOUT --quant_only is ref-guided discovery, and is not refused.

    Separate from the de novo case on purpose: the old guard keyed on
    ``not args.quant_only``, so it rejected this configuration too, and a
    relaxation that only looked at whether a --gtf was present would leave it
    rejected.
    """

    result = _lraa(
        "--chunk",
        "--gtf", "annot.gtf",
        "--bam", "reads.bam",
        "--genome", "genome.fa",
        "--output_prefix", str(tmp_path / "sample"),
        "--chunk_work_dir", str(tmp_path / "work"),
    )
    combined = result.stdout + result.stderr
    assert result.returncode != 0
    assert "--chunk requires --quant_only" not in combined
    assert "counting genome-mapped reads" in combined


def test_chunked_quant_only_still_requires_a_gtf(tmp_path):
    """Quant-only has nothing to quantify without one, and is refused up front.

    The refusal has to happen BEFORE any work runs, because reaching the pipeline
    here would mean the relaxation removed a check instead of narrowing it. The
    negative marker is the library count the chunked path takes first, not stage 1:
    strandless is the default and skips stage 1, so asserting stage 1's absence
    would pass no matter what this code did.
    """

    result = _lraa(
        "--chunk", "--quant_only",
        "--bam", "reads.bam",
        "--genome", "genome.fa",
        "--output_prefix", str(tmp_path / "sample"),
        "--chunk_work_dir", str(tmp_path / "work"),
    )
    combined = result.stdout + result.stderr
    assert result.returncode != 0
    assert "--chunk --quant_only requires --gtf" in combined
    assert "counting genome-mapped reads" not in combined
    assert "stage1_strand_split" not in combined


def test_a_chunk_worker_cannot_re_enter_the_pipeline():
    """The guard that stops a recursive pipeline forking without bound.

    The commands ChunkedRun builds never pass --chunk, so this can only trip if
    someone edits lraa_cmd or invokes LRAA by hand inside a chunk. Both are
    exactly the cases a promise about a command line would not cover.
    """

    result = _lraa(
        "--chunk", "--quant_only", *BASE, env={ChunkedRun.WORKER_ENV: "1"}
    )
    assert result.returncode != 0
    assert "chunk worker" in result.stdout + result.stderr


def test_the_guard_is_off_by_default():
    """A plain environment must not look like a chunk worker."""

    assert os.environ.get(ChunkedRun.WORKER_ENV) is None


def test_not_passing_chunk_leaves_the_normal_path_alone():
    """--chunk now defaults on, so reaching the normal path takes --no_chunk.

    Asserted through a failure that only the normal path produces -- the missing
    --bam check -- so that a stray dispatch into chunked mode would show up as a
    different error rather than the same one.
    """

    result = _lraa("--no_chunk", "--genome", "genome.fa", "--output_prefix", "sample")
    combined = result.stdout + result.stderr
    assert "Must specify --bam" in combined
    assert "chunk" not in combined.lower()


def test_omitting_chunk_and_stream_reads_reaches_chunked_mode_by_default(tmp_path):
    """The mirror image of the test above: naming NEITHER flag must dispatch into
    chunked mode, streaming, by default -- not the pre-v0.25.0 unchunked path.

    Asserted the same way test_chunked_discovery_is_reachable is: through the
    chunked path's own first act (the library count taken before any
    partitioning), on a bam that does not exist so the run dies immediately
    rather than doing real work. The missing --bam check that
    test_not_passing_chunk_leaves_the_normal_path_alone relies on never fires
    here, because chunked mode's own --bam requirement is checked by
    ``count_reads_from_bam`` inside ``_run_chunked_mode``, not by the unchunked
    setup's earlier check -- so reaching "Must specify --bam" instead would
    itself prove a default flip regression.
    """

    result = _lraa(
        "--bam", "reads.bam",
        "--genome", "genome.fa",
        "--output_prefix", str(tmp_path / "sample"),
        "--chunk_work_dir", str(tmp_path / "work"),
    )
    combined = result.stdout + result.stderr
    assert result.returncode != 0
    assert "Must specify --bam" not in combined
    assert "counting genome-mapped reads" in combined


def test_stream_reads_itself_defaults_on_even_without_chunk(tmp_path):
    """Isolates --stream_reads's own default from --chunk's: --no_chunk takes the
    unchunked path, and nothing else says anything about streaming.

    Discriminated through the same-bam-twice guard (LRAA:"needs a first-pass bam
    thinner than the one it streams"), which only exists to fire under
    --stream_reads. --no_norm plus no distinct --bam_for_sg makes both passes read
    the identical bam -- the one combination the guard refuses -- so seeing the
    refusal proves streaming was on despite never being named, and its absence
    would mean the default silently reverted to off.
    """

    result = _lraa(
        "--no_chunk", "--quant_only",
        "--bam", "reads.bam", "--gtf", "annot.gtf", "--genome", "genome.fa",
        "--no_norm",
        "--output_prefix", str(tmp_path / "sample"),
    )
    combined = result.stdout + result.stderr
    assert result.returncode != 0
    assert "needs a first-pass bam thinner than the one it streams" in combined


def test_chunk_appears_in_the_concise_help():
    """A mode nobody can discover is a mode that rots unnoticed."""

    result = _lraa("--help")
    assert result.returncode == 0
    assert "--chunk" in result.stdout


@pytest.mark.parametrize(
    "flag",
    ("--approx_MB_per_cut", "--approx_MB_per_cut_wiggle_window",
     "--chunk_depth_window", "--chunk_margin", "--chunk_work_dir"),
)
def test_chunking_flags_are_accepted(flag):
    """These were accepted, validated, written to config and read by nothing.

    Now they reach the pipeline. The check is that they parse; that they carry
    the canonical values is test_chunking_flags_match_the_pipeline below.
    """

    result = _lraa("--show_full_usage_info")
    assert flag in result.stdout


def test_chunking_flags_match_the_pipeline_defaults():
    """One default per constant, or a cache token can assert one geometry and reuse another.

    Every value here is baked into the stage-2 cut token. Two copies that drift
    do not produce a mismatch -- they produce a stale cache HIT, reusing the old
    cut geometry under a token naming the new parameters.
    """

    pipeline = ChunkedRun.build_parser().parse_args([])
    assert pipeline.approx_MB_per_cut == LRAA_Globals.config["approx_MB_per_cut"]
    assert (
        pipeline.approx_MB_per_cut_wiggle_window
        == LRAA_Globals.config["approx_MB_per_cut_wiggle_window"]
    )
    assert pipeline.depth_window == LRAA_Globals.config["chunk_depth_window"]
    assert pipeline.margin == LRAA_Globals.config["chunk_margin"]
    assert pipeline.random_seed == LRAA_Globals.config["chunk_random_seed"]


def test_every_chunking_constant_has_exactly_one_definition():
    """The duplication that made a silent divergence possible must not come back.

    depth_window existed as four copies, grid_origin three, margin and
    random_seed two each. Each consumer now reads LRAA_Globals; a literal
    reintroduced here is the thing this test exists to catch.
    """

    import importlib.util
    from importlib.machinery import SourceFileLoader

    def _load(name, relpath):
        path = REPO / relpath
        loader = SourceFileLoader(name, str(path))
        spec = importlib.util.spec_from_loader(name, loader)
        module = importlib.util.module_from_spec(spec)
        loader.exec_module(module)
        return module

    extractor = _load("extractor_under_test", "util/misc/extract_contig_region_inputs.py")
    selector = _load("selector_under_test", "util/misc/select_contig_cut_points.py")

    assert extractor.DEFAULT_MARGIN == LRAA_Globals.config["chunk_margin"]
    assert selector.DEFAULT_DEPTH_WINDOW == LRAA_Globals.config["chunk_depth_window"]
    assert selector.DEFAULT_GRID_ORIGIN == LRAA_Globals.config["chunk_grid_origin"]

    # ... and they track the config rather than merely agreeing with it today
    original = LRAA_Globals.config["chunk_margin"]
    try:
        LRAA_Globals.config["chunk_margin"] = original + 17
        reloaded = _load("extractor_reloaded", "util/misc/extract_contig_region_inputs.py")
        assert reloaded.DEFAULT_MARGIN == original + 17
    finally:
        LRAA_Globals.config["chunk_margin"] = original


def test_default_args_rejects_an_unknown_setting():
    """LRAA --chunk fills a namespace by name; a typo must not be silently ignored."""

    with pytest.raises(ChunkedRun.PipelineError):
        ChunkedRun.default_args(approx_MB_per_kut=5)


def test_default_args_carries_every_pipeline_default():
    """Every plain-parser default is carried through, field for field.

    Compared against ``ChunkedRun.parse_args(...)``, not the raw
    ``build_parser().parse_args([])`` -- ``stream_reads_rescue_unassigned``
    defaults to the sentinel ``None`` on the bare parser (unset: neither
    ``--stream_reads_rescue_unassigned`` nor its negation was given) and only
    ``parse_args``/``default_args`` resolve it, via the same
    ``_resolve_stream_reads_rescue_unassigned`` helper, to track
    ``--stream_reads``. The raw parser default would report a mismatch on
    that one field even though both routes agree once resolved.
    """
    built = ChunkedRun.default_args(bam="b", genome_fa="g", gtf="a", output_dir="/tmp/o")
    parsed = ChunkedRun.parse_args(
        ["--bam", "b", "--genome_fa", "g", "--gtf", "a", "--output_dir", "/tmp/o"]
    )
    for key, value in vars(parsed).items():
        if key in ("bam", "genome_fa", "gtf", "output_dir"):
            continue
        assert getattr(built, key) == value, key


def test_run_refuses_an_incomplete_namespace():
    """run() is reachable from two callers, so it validates rather than trusting."""

    with pytest.raises(ChunkedRun.PipelineError):
        ChunkedRun.run(ChunkedRun.default_args(bam="b", genome_fa="g", gtf="a"))


# -- discovery: the gate, and the reporting that makes it visible ---------------


def test_the_cli_requires_a_gtf_only_when_quantifying():
    with pytest.raises(SystemExit):
        ChunkedRun.parse_args(["--bam", "b", "--genome_fa", "g", "--output_dir", "o"])

    args = ChunkedRun.parse_args(
        ["--bam", "b", "--genome_fa", "g", "--output_dir", "o", "--discovery"]
    )
    assert args.discovery is True
    assert args.gtf is None


def test_discovery_refuses_the_baseline_arm(tmp_path):
    """The control arm is a whole-contig QUANT run; it is not a discovery control.

    Refused rather than served, because running it beside a discovery chunked arm
    would produce two artifacts that look comparable and are not.
    """

    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.run(
            ChunkedRun.default_args(
                bam="b",
                genome_fa="g",
                output_dir=str(tmp_path / "out"),
                discovery=True,
                arm="baseline",
            )
        )
    assert "no baseline arm" in str(err.value)


def test_stage_five_drops_quant_only_and_the_gtf_it_has_none_of():
    """The two tokens that separate the modes, and nothing else changes."""

    quant = ChunkedRun.default_args(
        bam="b", genome_fa="g", gtf="a", output_dir="/tmp/o"
    )
    cmd = ChunkedRun.lraa_cmd(
        quant, "q.bam", "sg.bam", "mini.fa", "mini.gtf", "out", 10, 2
    )
    assert "--quant_only" in cmd
    assert cmd[cmd.index("--gtf") + 1] == "mini.gtf"

    denovo = ChunkedRun.default_args(
        bam="b", genome_fa="g", output_dir="/tmp/o", discovery=True
    )
    cmd = ChunkedRun.lraa_cmd(denovo, "q.bam", "sg.bam", "mini.fa", None, "out", 10, 2)
    assert "--quant_only" not in cmd
    assert "--gtf" not in cmd
    # the splice-graph composition is the part that must NOT change
    assert cmd[cmd.index("--bam_for_sg") + 1] == "sg.bam"
    assert "--no_norm" in cmd

    guided = ChunkedRun.default_args(
        bam="b", genome_fa="g", gtf="a", output_dir="/tmp/o", discovery=True
    )
    cmd = ChunkedRun.lraa_cmd(
        guided, "q.bam", "sg.bam", "mini.fa", "mini.gtf", "out", 10, 2
    )
    assert "--quant_only" not in cmd
    assert cmd[cmd.index("--gtf") + 1] == "mini.gtf"


def test_the_output_suffixes_are_the_ones_lraa_actually_writes():
    """These name the files stage 6 then reads; a rename in LRAA must fail here.

    LRAA is a script, so the mapping cannot be imported -- it is restated in
    ChunkedRun and checked against LRAA's own function body, which is the only
    way a divergence shows up as a test failure instead of as paths that nothing
    writes.
    """

    body = (
        LRAA.read_text()
        .split("def _get_lraa_output_mode_suffix(")[1]
        .split("\ndef ")[0]
    )
    for suffix in (
        ChunkedRun.LRAA_QUANT_ONLY_SUFFIX,
        ChunkedRun.LRAA_REF_FREE_SUFFIX,
        ChunkedRun.LRAA_REF_GUIDED_SUFFIX,
    ):
        assert '"{}"'.format(suffix) in body

    assert ChunkedRun.lraa_output_suffix(False, "a.gtf") == "LRAA.quant-only"
    assert ChunkedRun.lraa_output_suffix(True, None) == "LRAA.ref-free"
    assert ChunkedRun.lraa_output_suffix(True, "a.gtf") == "LRAA.ref-guided"


def _cut_plan(args, outdir, contig="chr1"):
    """The stage-2 plan for one contig, which is where the selector command lives.

    Stage 2 is a pool of per-contig ``cut_selection_plan`` calls now rather than one
    loop, and the plan is pure: no subprocess runs and no cuts.json is read, so
    these tests inspect the command and the sentinel directly instead of faking
    ``run_step`` to get at them.
    """

    source = ("", ChunkedRun.STRANDLESS_TAG, "raw.bam", "root")
    return ChunkedRun.cut_selection_plan(
        args, outdir, os.path.join(outdir, "cuts"), source, contig
    )


def test_cut_selection_is_identical_in_both_modes(tmp_path):
    """Discovery must not ask the selector for a different cut rule.

    An earlier revision passed ``--require_zero_severed`` here and only here. That
    contract was rejected -- at depth every base is covered, so refusing to sever
    declines every cut -- so the two modes now differ in what stage 5 runs and in
    NOTHING about placement. Checked on the command line stage 2 builds, because
    that command line is the only thing that reaches the selector.
    """

    outdir = str(tmp_path / "out")
    os.makedirs(os.path.join(outdir, "logs"))

    for discovery, gtf in ((False, "annot.gtf"), (True, None)):
        args = ChunkedRun.default_args(
            bam="b",
            genome_fa="g",
            gtf=gtf,
            output_dir=outdir,
            discovery=discovery,
        )
        cmd = _cut_plan(args, outdir)["cmd"]
        assert "--require_zero_severed" not in cmd
        # the severing cost's shape is passed in both modes, and it is the same
        weight = cmd[cmd.index("--severed_multiexon_weight") + 1]
        assert weight == str(LRAA_Globals.config["chunk_severed_multiexon_weight"])
        assert ("--gtf" in cmd) is (gtf is not None)


def test_the_selector_command_carries_the_weight_it_is_given(tmp_path):
    """The weight decides the cut coordinates, so it must not be a config read.

    A selector left to read its own config would place cuts under one weight while
    the run recorded another.
    """

    outdir = str(tmp_path / "out")
    os.makedirs(os.path.join(outdir, "logs"))
    args = ChunkedRun.default_args(
        bam="b",
        genome_fa="g",
        gtf=None,
        output_dir=outdir,
        discovery=True,
        severed_multiexon_weight=3,
    )

    cmd = _cut_plan(args, outdir)["cmd"]

    assert cmd[cmd.index("--severed_multiexon_weight") + 1] == "3"


def test_the_weight_is_in_the_stage_2_cache_token(tmp_path):
    """It moves the cuts, so a stale hit would serve one geometry as another.

    The values are baked into the sentinel, and a token that omitted the objective
    would turn a changed weight into a cache HIT that reuses the old coordinates
    while the run asserts the new ones. Whether an annotation constrains placement
    at all is the same kind of input and is checked the same way: the per-contig
    sentinel this stage now writes carries both.
    """

    seen = []

    def token_for(**overrides):
        seen.append(overrides)
        outdir = str(tmp_path / "out{}".format(len(seen)))
        os.makedirs(os.path.join(outdir, "logs"))
        params = {
            "bam": "b",
            "genome_fa": "g",
            "gtf": None,
            "output_dir": outdir,
            "discovery": True,
        }
        params.update(overrides)
        args = ChunkedRun.default_args(**params)
        return _cut_plan(args, outdir)["token"]

    assert token_for(severed_multiexon_weight=10) != token_for(
        severed_multiexon_weight=3
    )
    assert token_for(gtf=None) != token_for(gtf="annot.gtf")


def test_a_severed_read_no_longer_fails_the_run(tmp_path):
    """The rejected contract raised here. Severing is now the expected outcome.

    What still has to hold is the IDENTITY: the reads selection NAMED are the reads
    extraction DROPPED. That is a statement about two tools agreeing, and it is
    what the parity comparison's pruned baseline depends on.
    """

    cut_dir = tmp_path / "cuts"
    cut_dir.mkdir()
    (cut_dir / "plus.dropped_reads.txt").write_text("readA\n")
    chunks = [{"manifest": {"dropped_read_names": ["readA"]}}]

    accounting = ChunkedRun.verify_severed_accounting(str(cut_dir), chunks)

    assert accounting["severed_reads"] == 1
    assert accounting["sets_identical"] is True
    assert "zero_severed_required" not in accounting


def test_a_disagreement_between_selection_and_extraction_still_fails(tmp_path):
    """Removing the severing veto must not remove the accounting check with it."""

    cut_dir = tmp_path / "cuts"
    cut_dir.mkdir()
    (cut_dir / "plus.dropped_reads.txt").write_text("readA\n")
    chunks = [{"manifest": {"dropped_read_names": ["readB"]}}]

    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.verify_severed_accounting(str(cut_dir), chunks)
    assert "accounting is inexact" in str(err.value)


def _placement_selection(**overrides):
    selection = {
        "chrom": "chrT",
        "strand": "+",
        "counts": {
            "targets": 3,
            "cuts_placed": 2,
            "targets_unplaced": 1,
            "targets_declined_annotation": 1,
            "targets_tail_merged": 0,
            "segments": 3,
            "alignments_dropped_at_cuts": 4,
            "alignments_dropped_monoexonic": 1,
            "alignments_dropped_multiexon": 3,
        },
        "cuts": [
            {
                "target": 1000,
                "position": 1100,
                "offset_from_target": 100,
                "spanning_alignments_dropped": 4,
                "severed_monoexonic": 1,
                "severed_multiexon": 3,
                "search_radius": 25000,
            },
            {
                "target": 3000,
                "position": 3000,
                "offset_from_target": 0,
                "spanning_alignments_dropped": 0,
                "severed_monoexonic": 0,
                "severed_multiexon": 0,
                "search_radius": 5000,
            },
        ],
        "unplaced_targets": [
            {
                "target": 2000,
                "declined_annotation": True,
                "best_spanning_in_window": None,
                "reason": "DECLINED: no position in the window is both on the "
                "100 bp depth-window grid and outside every annotated locus",
            }
        ],
    }
    selection.update(overrides)
    return {"+": [selection]}


def test_the_placement_report_names_what_each_cut_severed():
    """Severing is expected now, so this report is the only place it is visible.

    Per cut, split monoexonic against multi-exon: a cut that dropped three reads of
    depth and one that dropped three junctions are not the same event, and a bare
    total cannot tell them apart.
    """

    text, summary = ChunkedRun.cut_placement_report(
        _placement_selection(), discovery=True
    )

    assert summary["targets"] == 3
    assert summary["cuts_placed"] == 2
    assert summary["cuts_declined_annotation"] == 1
    assert summary["alignments_severed"] == 4
    assert summary["alignments_severed_monoexonic"] == 1
    assert summary["alignments_severed_multiexon"] == 3
    assert "zero_severed_required" not in summary
    assert summary["per_selection"][0]["cuts"][0]["severed_multiexon"] == 3

    assert "3 requested, 2 placed, 1 declined for annotation" in text
    assert "4 alignment(s) severed (1 monoexonic, 3 multi-exon)" in text
    assert "cut at 1100 (target 1000, +100) severs 4 alignment(s)" in text
    assert "1 monoexonic, 3 multi-exon; searched 25000 bp" in text
    # a clean cut is not listed: only what something cost is worth a line
    assert "cut at 3000" not in text
    assert "DECLINED target 2000" in text
    assert "3 of them" in text or "3 of them carried junctions" in text


def test_the_placement_report_claims_no_refusal_it_did_not_make():
    """Printed in quant-only too, and the wording must match the actual rule."""

    selections = {
        "+": [
            {
                "chrom": "chrT",
                "strand": "+",
                "counts": {
                    "targets": 2,
                    "cuts_placed": 2,
                    "targets_unplaced": 0,
                    "targets_declined_annotation": 0,
                    "targets_tail_merged": 0,
                    "segments": 3,
                    "alignments_dropped_at_cuts": 0,
                    "alignments_dropped_monoexonic": 0,
                    "alignments_dropped_multiexon": 0,
                },
                "cuts": [],
                "unplaced_targets": [],
            }
        ]
    }

    text, summary = ChunkedRun.cut_placement_report(selections, discovery=False)

    assert summary["cuts_declined_annotation"] == 0
    assert summary["alignments_severed"] == 0
    assert "severing is a COST, never a veto" in text
    assert "were DECLINED" not in text
    assert "carried junctions" not in text


def test_the_merged_gtf_is_rebased_and_its_model_ids_cannot_collide(tmp_path):
    """Every chunk names a comp-1, so concatenating unpatched fuses two loci.

    That happened once already, on chr21, and produced 37 spurious
    chromosome-crossing "models" before it was diagnosed. The two rewrites are
    the coordinate shift and the per-unit id prefix, and both are checked on ids
    that are deliberately IDENTICAL between the two units.
    """

    units = []
    for unit_id, offset in (("c0", 0), ("c1", 5000)):
        prefix = str(tmp_path / unit_id)
        attrs = 'gene_id "g:chrT:+:comp-1"; transcript_id "t:chrT:+:comp-1:iso-1";'
        with open(prefix + ".gtf", "wt") as fh:
            print("# LRAA version comment", file=fh)
            print(
                "\t".join(
                    ("chrT", "LRAA", "transcript", "100", "400", ".", "+", ".", attrs)
                ),
                file=fh,
            )
            print(
                "\t".join(
                    ("chrT", "LRAA", "exon", "100", "200", ".", "+", ".", attrs)
                ),
                file=fh,
            )
        units.append(
            {"unit_id": unit_id, "offset": offset, "quant_prefix": prefix}
        )

    merged_dir = str(tmp_path / "merged")
    os.makedirs(merged_dir)
    result = ChunkedRun.merge_discovery_gtf(merged_dir, units)

    with open(result["gtf"], "rt") as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if line[0] != "#"]

    assert result["gtf_transcripts"] == 2
    assert result["gtf_lines"] == 4

    sep = ChunkedRun.NAMESPACE_SEP
    ids = {re.search(r'transcript_id "([^"]*)"', row[8]).group(1) for row in rows}
    assert ids == {
        "c0{}t:chrT:+:comp-1:iso-1".format(sep),
        "c1{}t:chrT:+:comp-1:iso-1".format(sep),
    }
    genes = {re.search(r'gene_id "([^"]*)"', row[8]).group(1) for row in rows}
    assert genes == {
        "c0{}g:chrT:+:comp-1".format(sep),
        "c1{}g:chrT:+:comp-1".format(sep),
    }

    # coordinates land in the whole-contig frame, both columns, every row
    assert [(row[3], row[4]) for row in rows] == [
        ("100", "400"),
        ("100", "200"),
        ("5100", "5400"),
        ("5100", "5200"),
    ]


def test_strandless_is_the_default_chunking_mode(tmp_path):
    """Chunking is strandless unless asked otherwise, on BOTH routes into the pipeline.

    The two routes have separate argparse instances -- LRAA's own and ChunkedRun's --
    and ``default_args`` is documented as the single place both take defaults from. A
    disagreement would make `LRAA --chunk` and the parity driver run different
    pipelines from identical arguments, which is the drift that docstring exists to
    prevent, so both are asserted here rather than one.
    """

    assert ChunkedRun.default_args().strandless_chunks is True
    assert ChunkedRun.default_args(strandless_chunks=False).strandless_chunks is False
    assert LRAA_Globals.config["strandless_chunks"] is True


def test_the_library_count_is_taken_once_before_any_partitioning(tmp_path):
    """The TPM denominator is resolved before chunking, not required from the caller.

    Strandless has no whole-bam serial phase, so nothing downstream can count the
    library. An earlier revision refused the combination and told the caller to pass
    --num_total_reads, which put an extra required argument on the fast path and
    invited callers to compute it with `samtools view -c` -- which counts unmapped and
    supplementary records and inflates the denominator 3.89x on an SG-NEx ONT bam.

    Asserted on ORDER as well as occurrence: the count must precede cut selection, so
    a bam problem surfaces before the partition is computed rather than after.
    """

    result = _lraa(
        "--chunk",
        "--bam", "reads.bam",
        "--genome", "genome.fa",
        "--output_prefix", str(tmp_path / "sample"),
        "--chunk_work_dir", str(tmp_path / "work"),
    )
    combined = result.stdout + result.stderr

    assert "requires --num_total_reads" not in combined, "must not refuse; must count"
    assert "counting genome-mapped reads" in combined
    # the policy, not just that something was counted: -F 0x904 is what makes this
    # agree with the unchunked path's count_reads_from_bam
    assert "0x904" in combined
    # and it happened before any cut selection: the run dies at the count, so no cut
    # stage token appears at all. ("cut" as a bare substring would match "execute".)
    assert "stage2_cuts" not in combined


def test_chunk_by_strand_selects_the_strand_first_ordering():
    """The opt-out is named for what it does, and it has to actually do it.

    Worth its own test because the flag's whole job is to set a boolean FALSE, which
    is the shape that fails silently: an emitted-but-ignored flag, a dest typo, or a
    WDL that omits the flag in its false case all leave the default in place and the
    run looks entirely normal. The WDL did exactly that before this rename -- its
    false case emitted nothing, and nothing means strandless now.
    """

    parser_default = ChunkedRun.build_parser().parse_args([])
    assert parser_default.strandless_chunks is True

    opted_out = ChunkedRun.build_parser().parse_args(["--chunk_by_strand"])
    assert opted_out.strandless_chunks is False

    # and the retired spelling must be REJECTED, not ignored: a script still passing
    # it was asking for strand-first, and silently accepting it would hand back
    # strandless -- the failure this rename exists to remove
    with pytest.raises(SystemExit):
        ChunkedRun.build_parser().parse_args(["--no_strandless_chunks"])


@pytest.mark.parametrize("extra", [[], ["--chunk"]], ids=["unchunked", "chunked"])
def test_the_removed_LowFi_flag_is_rejected_on_every_path(tmp_path, extra):
    """A removed flag has to be removed for the whole program, not for whichever code
    path is reached first.

    --LowFi was declared with argparse.SUPPRESS and a comment calling it a no-op, while
    400 lines later an unchunked run exited on it. The chunked dispatch sits between the
    two, so `--chunk --LowFi` was silently accepted and `--LowFi` alone died in 0.5 s --
    the same flag, two answers, nothing in either output saying which had happened. Found
    by running both arms rather than by reading the declaration, which is how I got it
    wrong: the declaration says no-op and the program does not.

    Parametrized over both paths on purpose. A single-path version of this test is what
    allowed the divergence in the first place.
    """

    result = _lraa(
        "--genome", "genome.fa", "--bam", "reads.bam", "--gtf", "annot.gtf",
        "--quant_only", "--output_prefix", str(tmp_path / "s"), "--LowFi", *extra,
    )
    combined = result.stdout + result.stderr
    assert result.returncode != 0
    assert "--LowFi has been removed" in combined
    # and it must say what to do instead, or the rejection just blocks the user
    assert "--HiFi" in combined
    # rejected BEFORE any work: reaching the pipeline would mean the check moved back
    # behind the dispatch again
    assert "counting genome-mapped reads" not in combined



# -- the scatter seams: --stop_after_make_chunks / --no_reuse_source_bam / --only_chunk


def _write_two_contig_inputs(root):
    """A deterministic two-contig genome and a sorted, indexed bam over it.

    Two contigs because the scatter exists to fan out ACROSS contigs; each carries
    spliced forward reads (a canonical GT..AG intron) and unspliced reverse reads,
    so a strandless chunk yields a non-empty unit in BOTH orientations.
    """

    import random

    import pysam

    rng = random.Random(42)
    genome = {}
    for chrom in ("chrA", "chrB"):
        seq = [rng.choice("ACGT") for _ in range(10000)]
        seq[1500:1502] = ["G", "T"]  # donor of the intron [1500, 2000)
        seq[1998:2000] = ["A", "G"]  # acceptor
        genome[chrom] = "".join(seq)

    fasta = root / "genome.fa"
    with open(fasta, "wt") as fh:
        for chrom, seq in genome.items():
            print(">" + chrom, file=fh)
            print(seq, file=fh)
    pysam.faidx(str(fasta))

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": chrom, "LN": len(seq)} for chrom, seq in genome.items()],
    }
    bam = root / "reads.bam"
    total = 0
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for tid, (chrom, seq) in enumerate(genome.items()):
            for i in range(20):
                aln = pysam.AlignedSegment()
                aln.query_name = "{}_plus_{}".format(chrom, i)
                aln.flag = 0
                aln.reference_id = tid
                aln.reference_start = 1000
                aln.mapping_quality = 60
                aln.cigarstring = "500M500N500M"
                aln.query_sequence = seq[1000:1500] + seq[2000:2500]
                aln.query_qualities = pysam.qualitystring_to_array("I" * 1000)
                out.write(aln)
                total += 1
            for i in range(20):
                aln = pysam.AlignedSegment()
                aln.query_name = "{}_minus_{}".format(chrom, i)
                aln.flag = 16
                aln.reference_id = tid
                aln.reference_start = 5000
                aln.mapping_quality = 60
                aln.cigarstring = "1000M"
                aln.query_sequence = seq[5000:6000]
                aln.query_qualities = pysam.qualitystring_to_array("I" * 1000)
                out.write(aln)
                total += 1
    pysam.index(str(bam))
    return str(fasta), str(bam), total


@pytest.fixture(scope="module")
def make_chunks_run(tmp_path_factory):
    """One --stop_after_make_chunks run, shared by the tests that consume it."""

    root = tmp_path_factory.mktemp("scatter")
    fasta, bam, total = _write_two_contig_inputs(root)
    outdir = root / "work"
    outputs = ChunkedRun.run(
        ChunkedRun.default_args(
            bam=bam,
            genome_fa=fasta,
            output_dir=str(outdir),
            discovery=True,
            num_total_reads=total,
            no_reuse_source_bam=True,
            stop_after_make_chunks=True,
            cpu_budget=2,
        )
    )
    return {
        "root": root,
        "fasta": fasta,
        "bam": bam,
        "total": total,
        "outdir": outdir,
        "outputs": outputs,
    }


def test_stop_after_make_chunks_writes_the_plan_and_runs_nothing_downstream(
    make_chunks_run,
):
    """The prep task's contract: a chunk_plan.json naming exactly the chunk
    directories on disk, and NO merged output -- stages 3b-6 must not have run.
    """

    import json

    outdir = make_chunks_run["outdir"]
    plan_path = outdir / "chunk_plan.json"
    assert plan_path.exists()
    with open(plan_path) as fh:
        plan = json.load(fh)

    assert plan["version"] == ChunkedRun.CHUNK_PLAN_VERSION
    assert plan["num_total_reads"] == make_chunks_run["total"]
    assert plan["discovery"] is True

    planned_ids = {c["chunk_id"] for c in plan["chunks"]}
    on_disk = {p.name for p in (outdir / "chunks").iterdir() if p.is_dir()}
    assert planned_ids == on_disk
    # both contigs are in the fan-out, not one chunk covering everything
    assert {c["chrom"] for c in plan["chunks"]} == {"chrA", "chrB"}

    assert not (outdir / "merged").exists()
    assert make_chunks_run["outputs"]["stopped_after_make_chunks"] is True


def test_no_reuse_source_bam_makes_every_chunk_self_contained(
    make_chunks_run, tmp_path
):
    """WITH the flag every chunk holds a real mini bam; WITHOUT it, on the same
    inputs, at least one chunk manifest names the source instead. The second half
    is what pins the flag as load-bearing rather than cosmetic: these contigs are
    shorter than the segment span, so reuse fires unless suppressed.
    """

    import json

    for cdir in (make_chunks_run["outdir"] / "chunks").iterdir():
        assert (cdir / "chunk.bam").exists(), cdir
        with open(cdir / "chunk.partition.json") as fh:
            manifest = json.load(fh)
        assert manifest["bam_reused_from_source"] is False, cdir

    reuse_outdir = tmp_path / "work_reuse"
    ChunkedRun.run(
        ChunkedRun.default_args(
            bam=make_chunks_run["bam"],
            genome_fa=make_chunks_run["fasta"],
            output_dir=str(reuse_outdir),
            discovery=True,
            num_total_reads=make_chunks_run["total"],
            stop_after_make_chunks=True,
            cpu_budget=2,
        )
    )
    reused = []
    for cdir in (reuse_outdir / "chunks").iterdir():
        with open(cdir / "chunk.partition.json") as fh:
            reused.append(json.load(fh)["bam_reused_from_source"])
    assert any(reused)


def test_only_chunk_processes_one_chunk_and_writes_its_units(make_chunks_run):
    """The leaf task's contract: --only_chunk on a make-chunks directory runs
    stages 3b-5 for that chunk alone and writes units.json whose quant prefixes
    all resolve to real quant.expr files, one unit per orientation.
    """

    import json

    outdir = make_chunks_run["outdir"]
    with open(outdir / "chunk_plan.json") as fh:
        plan = json.load(fh)
    chunk_id = plan["chunks"][0]["chunk_id"]

    result = ChunkedRun.run(
        ChunkedRun.default_args(
            output_dir=str(outdir),
            only_chunk=chunk_id,
            cpu_budget=2,
        )
    )
    assert result["only_chunk"] == chunk_id

    units_path = outdir / "chunks" / chunk_id / "units.json"
    assert units_path.exists()
    with open(units_path) as fh:
        doc = json.load(fh)
    assert doc["chunk_id"] == chunk_id

    # strandless mode: exactly one unit per orientation, ids derived from the chunk's
    assert [u["unit_id"] for u in doc["units"]] == [
        chunk_id + "_plus",
        chunk_id + "_minus",
    ]
    for unit in doc["units"]:
        expr = unit["quant_prefix"] + ".quant.expr"
        assert os.path.exists(expr), expr


def test_only_chunk_refuses_an_id_the_plan_does_not_name(make_chunks_run):
    """No fallback: a leaf that cannot identify its own chunk must fail, not
    quantify something else."""

    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.run(
            ChunkedRun.default_args(
                output_dir=str(make_chunks_run["outdir"]),
                only_chunk="chrZ_99",
                cpu_budget=1,
            )
        )
    assert "not in the plan" in str(err.value)


def test_only_chunk_refuses_a_reuse_extracted_chunk(tmp_path):
    """A manifest naming the SOURCE bam cannot be processed on another machine;
    the refusal must say how to re-extract."""

    import json

    outdir = tmp_path / "work"
    cdir = outdir / "chunks" / "chrA_00"
    cdir.mkdir(parents=True)
    with open(outdir / "chunk_plan.json", "wt") as fh:
        json.dump(
            {
                "version": ChunkedRun.CHUNK_PLAN_VERSION,
                "num_total_reads": 1,
                "discovery": True,
                "lraa_suffix": "LRAA.ref-free",
                "chunks": [
                    {
                        "chunk_id": "chrA_00",
                        "chrom": "chrA",
                        "strand": None,
                        "strandless": True,
                        "region": "chrA:1-10000",
                        "index": 0,
                        "order": 0,
                        "has_gtf": False,
                    }
                ],
            },
            fh,
        )
    with open(cdir / "chunk.partition.json", "wt") as fh:
        json.dump({"bam_reused_from_source": True, "offset": 0}, fh)

    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.run(
            ChunkedRun.default_args(
                output_dir=str(outdir), only_chunk="chrA_00", cpu_budget=1
            )
        )
    assert "--no_reuse_source_bam" in str(err.value)


def test_a_namespaced_model_id_survives_gffcompare_tracking(tmp_path):
    """The namespace separator must not be a character gffcompare delimits with.

    gffcompare writes its query column as
    ``qJ:gene_id|transcript_id|num_exons|FPKM|TPM|cov|len`` and does not escape
    the pipe, and `gene_symbol_utils.parse_gffcompare_mappings` therefore reads
    those subfields POSITIONALLY. While chunked discovery namespaced ids with a
    pipe, every tracking row shattered: on the single-cell fixture 182 rows
    parsed to keys like "chr19_00_plus" (the unit id alone) and
    "g:chr19:+:comp-44" (the model id stripped of its namespace), so no id
    anything downstream held was ever matched and every gene symbol went
    unassigned -- the single-cell pipeline refused the run with "no gene names
    were assigned to feature ids".

    Asserted through the REAL parser on a real tracking record, not by
    inspecting the separator: the contract is that the id round-trips, and only
    the consumer that broke can attest to that.
    """

    from gene_symbol_utils import parse_gffcompare_mappings

    gene = ChunkedRun._namespace_id("chr19_00_plus", "g:chr19:+:comp-1")
    tx = ChunkedRun._namespace_id("chr19_00_plus", "t:chr19:+:comp-1:iso-1")

    tracking = tmp_path / "gffcmp.tracking"
    tracking.write_text(
        "TCONS_00000001\tXLOC_000001\tENSG00000141933.8|ENST00000359315.5\tj\t"
        "q1:{}|{}|1|0.000000|0.000000|0.000000|965\n".format(gene, tx)
    )

    mappings = parse_gffcompare_mappings(str(tracking))

    # the ids the merge actually emits are the ids the parser recovers
    assert gene in mappings, sorted(mappings)
    assert tx in mappings, sorted(mappings)
    assert mappings[gene] == ("ENSG00000141933.8", "ENST00000359315.5")
    assert mappings[tx] == ("ENSG00000141933.8", "ENST00000359315.5")


def test_stage_six_merges_read_assignment_summaries_with_real_counts(
    make_chunks_run,
):
    """A chunked run must publish a read-assignment summary that COUNTS reads.

    The failure this defends against shipped and was invisible: nothing merged
    the per-unit summaries, so the task-level file was never written, the
    workflow's own merge was handed an empty list, and it published a
    correctly-schema'd table whose TOTAL row read 0 reads -- against per-unit
    files reporting 11,768. Asserting the file exists is therefore not enough;
    the assertion has to be on the NUMBER, and on it agreeing with the units it
    came from.

    Also pins the skip: stage 6 returns None when NO unit has a summary, so that
    a chunk directory predating these files still merges. That escape hatch is
    only safe while a real run demonstrably takes the other branch.
    """

    import csv
    import json

    outdir = make_chunks_run["outdir"]
    with open(outdir / "chunk_plan.json") as fh:
        chunk_id = json.load(fh)["chunks"][0]["chunk_id"]

    ChunkedRun.run(
        ChunkedRun.default_args(
            output_dir=str(outdir), only_chunk=chunk_id, cpu_budget=2
        )
    )
    with open(outdir / "chunks" / chunk_id / "units.json") as fh:
        units = json.load(fh)["units"]

    merged_dir = tmp = outdir / "summary_merge"
    os.makedirs(merged_dir, exist_ok=True)
    out = ChunkedRun.merge_read_assignment_summaries(str(merged_dir), units)
    assert out is not None, "a real chunked run must produce a merged summary"

    with open(out, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    workers = [r for r in rows if r["row_type"] == "worker"]
    totals = [r for r in rows if r["row_type"] == "TOTAL"]

    # one worker row per unit, exactly one recomputed TOTAL
    assert len(workers) == len(units), (len(workers), len(units))
    assert len(totals) == 1

    # the number, not merely the file: TOTAL is the sum of the units, and nonzero
    expected = sum(int(r["reads_total"]) for r in workers)
    assert int(totals[0]["reads_total"]) == expected
    assert expected > 0, "the fixture's reads should be counted, not zero"

    # and the per-unit TOTAL rows were NOT folded in a second time
    per_unit_total = 0
    for unit in units:
        with open(unit["quant_prefix"] + ".read_assignment.summary.tsv", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row["row_type"] == "worker":
                    per_unit_total += int(row["reads_total"])
    assert int(totals[0]["reads_total"]) == per_unit_total