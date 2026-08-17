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

import json
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

    Asserted by the run reaching stage 1 and dying on the bam that does not
    exist. A failure only the pipeline can produce cannot be confused with the
    guard still firing under a different message.
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
    assert "stage1_strand_split" in combined


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
    assert "stage1_strand_split" in combined


def test_chunked_quant_only_still_requires_a_gtf(tmp_path):
    """Quant-only has nothing to quantify without one, and is refused up front.

    The refusal has to happen BEFORE any stage runs: reaching stage 1 here would
    mean the relaxation removed a check instead of narrowing it.
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
    """--chunk is opt-in: without it LRAA must reach its usual argument handling.

    Asserted through a failure that only the normal path produces -- the missing
    --bam check -- so that a stray dispatch into chunked mode would show up as a
    different error rather than the same one.
    """

    result = _lraa("--genome", "genome.fa", "--output_prefix", "sample")
    combined = result.stdout + result.stderr
    assert "Must specify --bam" in combined
    assert "chunk" not in combined.lower()


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
    built = ChunkedRun.default_args(bam="b", genome_fa="g", gtf="a", output_dir="/tmp/o")
    parsed = ChunkedRun.build_parser().parse_args([])
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


def test_cut_selection_is_identical_in_both_modes(tmp_path, monkeypatch):
    """Discovery must not ask the selector for a different cut rule.

    An earlier revision passed ``--require_zero_severed`` here and only here. That
    contract was rejected -- at depth every base is covered, so refusing to sever
    declines every cut -- so the two modes now differ in what stage 5 runs and in
    NOTHING about placement. Checked on the command line stage 2 builds, because
    that command line is the only thing that reaches the selector.
    """

    recorded = []

    def fake_run_step(name, cmd, log_path, cwd, rss_interval, append=True):
        recorded.append(cmd)
        prefix = cmd[cmd.index("--output_prefix") + 1]
        with open(prefix + ".cuts.json", "wt") as fh:
            json.dump([], fh)
        return {"step": name, "cmd": cmd}

    monkeypatch.setattr(ChunkedRun, "run_step", fake_run_step)

    outdir = str(tmp_path / "out")
    os.makedirs(os.path.join(outdir, "logs"))
    ckpt = ChunkedRun.Checkpoints(os.path.join(outdir, "__ckpt"))
    sources = [("", ChunkedRun.STRANDLESS_TAG, "raw.bam", "root")]

    for discovery, gtf in ((False, "annot.gtf"), (True, None)):
        recorded.clear()
        args = ChunkedRun.default_args(
            bam="b",
            genome_fa="g",
            gtf=gtf,
            output_dir=outdir,
            discovery=discovery,
        )
        ChunkedRun.stage_select_cuts(args, ckpt, outdir, {}, sources, 0.5)
        assert len(recorded) == 1
        assert "--require_zero_severed" not in recorded[0]
        # the severing cost's shape is passed in both modes, and it is the same
        weight = recorded[0][recorded[0].index("--severed_multiexon_weight") + 1]
        assert weight == str(LRAA_Globals.config["chunk_severed_multiexon_weight"])
        assert ("--gtf" in recorded[0]) is (gtf is not None)


def test_the_selector_command_carries_the_weight_it_is_given(tmp_path, monkeypatch):
    """The weight decides the cut coordinates, so it must not be a config read.

    A selector left to read its own config would place cuts under one weight while
    the run recorded another.
    """

    recorded = []

    def fake_run_step(name, cmd, log_path, cwd, rss_interval, append=True):
        recorded.append(cmd)
        prefix = cmd[cmd.index("--output_prefix") + 1]
        with open(prefix + ".cuts.json", "wt") as fh:
            json.dump([], fh)
        return {"step": name, "cmd": cmd}

    monkeypatch.setattr(ChunkedRun, "run_step", fake_run_step)

    outdir = str(tmp_path / "out")
    os.makedirs(os.path.join(outdir, "logs"))
    ckpt = ChunkedRun.Checkpoints(os.path.join(outdir, "__ckpt"))
    sources = [("", ChunkedRun.STRANDLESS_TAG, "raw.bam", "root")]
    args = ChunkedRun.default_args(
        bam="b",
        genome_fa="g",
        gtf=None,
        output_dir=outdir,
        discovery=True,
        severed_multiexon_weight=3,
    )

    ChunkedRun.stage_select_cuts(args, ckpt, outdir, {}, sources, 0.5)

    cmd = recorded[0]
    assert cmd[cmd.index("--severed_multiexon_weight") + 1] == "3"


def test_the_weight_is_in_the_stage_2_cache_token(tmp_path, monkeypatch):
    """It moves the cuts, so a stale hit would serve one geometry as another.

    The values are baked into the sentinel, and a token that omitted the objective
    would turn a changed weight into a cache HIT that reuses the old coordinates
    while the run asserts the new ones.
    """

    def fake_run_step(name, cmd, log_path, cwd, rss_interval, append=True):
        prefix = cmd[cmd.index("--output_prefix") + 1]
        with open(prefix + ".cuts.json", "wt") as fh:
            json.dump([], fh)
        return {"step": name, "cmd": cmd}

    monkeypatch.setattr(ChunkedRun, "run_step", fake_run_step)
    sources = [("", ChunkedRun.STRANDLESS_TAG, "raw.bam", "root")]

    tokens = []
    for weight in (10, 3):
        outdir = str(tmp_path / "out{}".format(weight))
        os.makedirs(os.path.join(outdir, "logs"))
        ckpt = ChunkedRun.Checkpoints(os.path.join(outdir, "__ckpt"))
        args = ChunkedRun.default_args(
            bam="b",
            genome_fa="g",
            gtf=None,
            output_dir=outdir,
            discovery=True,
            severed_multiexon_weight=weight,
        )
        _, _, cuts_tokens = ChunkedRun.stage_select_cuts(
            args, ckpt, outdir, {}, sources, 0.5
        )
        tokens.append(cuts_tokens[""])

    assert tokens[0] != tokens[1]


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

    ids = {re.search(r'transcript_id "([^"]*)"', row[8]).group(1) for row in rows}
    assert ids == {"c0|t:chrT:+:comp-1:iso-1", "c1|t:chrT:+:comp-1:iso-1"}
    genes = {re.search(r'gene_id "([^"]*)"', row[8]).group(1) for row in rows}
    assert genes == {"c0|g:chrT:+:comp-1", "c1|g:chrT:+:comp-1"}

    # coordinates land in the whole-contig frame, both columns, every row
    assert [(row[3], row[4]) for row in rows] == [
        ("100", "400"),
        ("100", "200"),
        ("5100", "5400"),
        ("5100", "5200"),
    ]
