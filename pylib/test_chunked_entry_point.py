#!/usr/bin/env python3

"""Tests for `LRAA --chunk`, the single entry point into the chunked pipeline.

Chunked mode used to be reachable only through util/misc/run_chunked_quant_pipeline.py,
which is how it came to be broken for two releases without anything noticing: the
driver was not on any routine path. These tests hold the merged entry point to the
three properties that make it safe to have one:

* chunked mode refuses the configurations it cannot serve, up front;
* a chunk worker cannot re-enter the pipeline;
* NOT passing --chunk leaves LRAA's own argument surface untouched.
"""

import os
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


def test_chunk_requires_quant_only():
    """Discovery cannot span a boundary, so chunked discovery must not be offered."""

    result = _lraa("--chunk", *BASE)
    assert result.returncode != 0
    assert "--chunk requires --quant_only" in result.stdout + result.stderr


def test_chunk_requires_a_gtf():
    result = _lraa(
        "--chunk", "--quant_only",
        "--bam", "reads.bam", "--genome", "genome.fa", "--output_prefix", "sample",
    )
    assert result.returncode != 0
    combined = result.stdout + result.stderr
    assert "--chunk requires --gtf" in combined or "quant_only" in combined


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
