#!/usr/bin/env python3

"""A sentinel may not survive a change to anything that decided its contents.

The pipeline is resumable through sentinel files named for a stage and for the
parameters that change that stage's output. That is sound only while every such
parameter is re-listed at every stage it reaches, and it was not:
``baseline_norm`` named the coverage target, window, seed and origin but not
``--max_intron_length``, which it hands the normalizer. Changing the cap reran
the merge into the same ``whole.primary.bam`` and then reused a normalized bam
built under the old cap. Neither quant stage named anything about the bams it
reads, so the same change left both quant results in place as well, and nothing
named the run's own inputs, so pointing one outdir at a different bam reused
everything.

These assertions drive the real stage functions rather than rebuilding their
token strings, because a test that restates the format cannot notice a field
being dropped from it. The stages run no subprocess: ``run_step`` is replaced,
and each function is allowed to fail on its missing outputs after it has written
its sentinels, which is the point at which the tokens are observable.
"""

import argparse
import importlib.util
import os
import sys

import pytest

_HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "..", "pylib"))


def _load_driver():
    """The pipeline now lives in pylib/ChunkedRun.py; the script here is a front end."""
    import ChunkedRun

    return ChunkedRun


DRIVER = _load_driver()


BASE_SETTINGS = dict(
    max_intron_length=200000,
    normalize_max_cov_level=1000,
    depth_window=100,
    random_seed=42,
    margin=1000,
    approx_MB_per_cut=10,
    approx_MB_per_cut_wiggle_window=1,
    HiFi=False,
    contig=None,
    cpu_budget=16,
    # These fixtures are quant-only, which is what the sentinels below encode.
    discovery=False,
    # C2's forwarded settings (see pylib/test_chunked_quant_token_argv.py's
    # identical fixture fix for why cell_barcode_tag/read_umi_tag are absent
    # rather than defaulted).
    cell_list=None,
    stream_reads=False,
    stream_reads_rescue_unassigned=False,
    stream_reads_rescue_unassigned_to_targets=False,
    min_mapping_quality=0,
    min_mapping_quality_for_final_quant=0,
)


def _args(tmp_path, **overrides):
    settings = dict(BASE_SETTINGS)
    settings.update(overrides)
    for name in ("bam", "genome_fa", "gtf"):
        path = tmp_path / name
        if not path.exists():
            path.write_text(name)
        settings[name] = str(path)
    return argparse.Namespace(**settings)


def _baseline_tokens(tmp_path, monkeypatch, split_token="stage1.up_aaaaaaaaaaaa", **kw):
    """Sentinels run_baseline writes, with no subprocess and no outputs."""
    marked = []

    class RecordingCheckpoints(DRIVER.Checkpoints):
        def done(self, token):
            marked.append(token)
            return False

        def mark(self, token, note=""):
            pass

    monkeypatch.setattr(DRIVER, "run_step", lambda *a, **k: {"step": "stub"})

    outdir = tmp_path / "out"
    (outdir / "logs").mkdir(parents=True, exist_ok=True)
    ckpt = RecordingCheckpoints(str(tmp_path / "ckpt"))

    with pytest.raises(DRIVER.PipelineError):
        DRIVER.run_baseline(
            _args(tmp_path, **kw),
            ckpt,
            str(outdir),
            {},
            {"+": str(tmp_path / "plus.bam"), "-": str(tmp_path / "minus.bam")},
            1000,
            0.5,
            split_token,
        )

    assert len(marked) == 3, marked
    return dict(zip(("merge", "norm", "quant"), marked))


def _chunk_tokens(tmp_path, monkeypatch, upstream="stage3.up_aaaaaaaaaaaa", **kw):
    """Sentinels chunk_worker writes, same arrangement."""
    marked = []

    class RecordingCheckpoints(DRIVER.Checkpoints):
        def done(self, token):
            marked.append(token)
            return False

        def mark(self, token, note=""):
            pass

    monkeypatch.setattr(DRIVER, "run_step", lambda *a, **k: {"step": "stub"})

    cdir = tmp_path / "chunk"
    cdir.mkdir(exist_ok=True)
    chunk = {
        "chunk_id": "chr20_plus_00",
        # strand-first: one quant unit, no in-chunk split
        "strandless": False,
        "dir": str(cdir),
        "prefix": str(cdir / "chunk"),
        "log": str(cdir / "chunk.log"),
        "window_origin": 0,
        "upstream_token": upstream,
        "units": DRIVER.chunk_quant_units(
            "chr20_plus_00", str(cdir), str(cdir / "chunk"), "+", 0, 0
        ),
    }

    with pytest.raises(DRIVER.PipelineError):
        DRIVER.chunk_worker(
            _args(tmp_path, **kw),
            RecordingCheckpoints(str(tmp_path / "ckpt")),
            str(tmp_path / "out"),
            chunk,
            1000,
            0.5,
            4,
        )

    assert len(marked) == 2, marked
    return dict(zip(("norm", "quant"), marked))


# ---------------------------------------------------------------- the helper


def test_a_token_moves_when_its_parent_moves():
    """The whole point: a stage inherits its parent rather than re-listing it."""
    a = DRIVER.chain_token("step.k_1", "parent_a")
    b = DRIVER.chain_token("step.k_1", "parent_b")
    assert a != b
    assert a.startswith("step.k_1.") and b.startswith("step.k_1.")


def test_a_token_moves_when_its_own_fields_move():
    assert DRIVER.chain_token("step.k_1", "p") != DRIVER.chain_token("step.k_2", "p")


def test_a_token_is_stable_for_identical_inputs():
    """Sentinels are only useful if an unchanged run hits them."""
    assert DRIVER.chain_token("step.k_1", "p", "x") == DRIVER.chain_token(
        "step.k_1", "p", "x"
    )


def test_opaque_fields_move_the_token_without_entering_the_name():
    """Resolved paths and stat pairs decide contents but cannot go in a filename."""
    bare = DRIVER.chain_token("inputs", None)
    with_id = DRIVER.chain_token("inputs", None, "/some/bam:12:34")
    assert bare != with_id
    assert with_id.startswith("inputs.")
    assert "/" not in with_id


def test_a_token_is_a_usable_filename():
    ckpt_name = DRIVER.chain_token("inputs", None, "/abs/path.bam:1:2")
    assert os.path.basename(ckpt_name) == ckpt_name


# ------------------------------------------------------------- the baseline arm


def test_the_intron_cap_moves_every_baseline_sentinel(tmp_path, monkeypatch):
    """The reported defect. The cap reaches the merge and the normalizer, and the
    quant reads what both produced, so all three have to move with it."""
    wide = _baseline_tokens(tmp_path, monkeypatch, max_intron_length=200000)
    narrow = _baseline_tokens(tmp_path, monkeypatch, max_intron_length=50000)

    for stage in ("merge", "norm", "quant"):
        assert wide[stage] != narrow[stage], stage


def test_normalizer_settings_move_the_baseline_quant_sentinel(tmp_path, monkeypatch):
    """Quant reads the normalized bam, so anything deciding that bam decides quant.

    The merge is upstream of neither, and must not move: a sentinel that changes
    when nothing it describes changed costs a rerun on every parameter sweep.
    """
    base = _baseline_tokens(tmp_path, monkeypatch)
    for setting, value in (
        ("normalize_max_cov_level", 500),
        ("depth_window", 250),
        ("random_seed", 7),
    ):
        moved = _baseline_tokens(tmp_path, monkeypatch, **{setting: value})
        assert moved["norm"] != base["norm"], setting
        assert moved["quant"] != base["quant"], setting
        assert moved["merge"] == base["merge"], setting


def test_an_upstream_change_moves_the_whole_baseline_arm(tmp_path, monkeypatch):
    """Stage 1's own inputs reach here only through its token."""
    base = _baseline_tokens(tmp_path, monkeypatch, split_token="stage1.up_aaaaaaaaaaaa")
    moved = _baseline_tokens(tmp_path, monkeypatch, split_token="stage1.up_bbbbbbbbbbbb")
    for stage in ("merge", "norm", "quant"):
        assert base[stage] != moved[stage], stage


def test_an_unchanged_baseline_run_hits_its_sentinels(tmp_path, monkeypatch):
    assert _baseline_tokens(tmp_path, monkeypatch) == _baseline_tokens(
        tmp_path, monkeypatch
    )


# -------------------------------------------------------------- the chunked arm


def test_the_intron_cap_moves_the_chunk_sentinels(tmp_path, monkeypatch):
    """Stage 4 passes the cap to the normalizer and stage 5 reads what it wrote."""
    wide = _chunk_tokens(tmp_path, monkeypatch, max_intron_length=200000)
    narrow = _chunk_tokens(tmp_path, monkeypatch, max_intron_length=50000)
    assert wide["norm"] != narrow["norm"]
    assert wide["quant"] != narrow["quant"]


def test_an_upstream_change_moves_both_chunk_sentinels(tmp_path, monkeypatch):
    """A re-extracted chunk must not be normalized and quantified from cache."""
    base = _chunk_tokens(tmp_path, monkeypatch, upstream="stage3.up_aaaaaaaaaaaa")
    moved = _chunk_tokens(tmp_path, monkeypatch, upstream="stage3.up_bbbbbbbbbbbb")
    assert base["norm"] != moved["norm"]
    assert base["quant"] != moved["quant"]


def test_an_unchanged_chunk_hits_its_sentinels(tmp_path, monkeypatch):
    assert _chunk_tokens(tmp_path, monkeypatch) == _chunk_tokens(tmp_path, monkeypatch)
