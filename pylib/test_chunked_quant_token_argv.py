#!/usr/bin/env python3

"""The quant sentinels name the COMMAND they guard, not a hand-listed subset of it.

Both quant stages -- the per-chunk one and the whole-contig baseline control --
were keyed on four things: the unit id, the total-read denominator, ``--HiFi``
and the upstream normalization token. Nothing in either key said whether the run
was quant-only or de novo discovery, and nothing said anything about a flag a
caller might later forward. So a resumed or repeated run could answer a
single-cell request, or a discovery request, out of a checkpoint written by a
bulk quant-only run, and the result is complete, internally consistent and
indistinguishable from a correct one without reading the reads.

An enumerated key cannot be fixed by lengthening it, because the list is
maintained by hand and the next forwarded flag is not on it. The keys are
therefore derived from the argv each stage actually issues: ``lraa_cmd`` builds
the command, ``argv_digest`` reduces it, and the sentinel moves for any argument
added, removed or changed -- with no edit at either site.

The tests drive the real ``chunk_worker`` and ``run_baseline`` rather than
rebuilding a token string, because a test that restates the format cannot notice
a field being dropped from it. Where a property depends on a flag the pipeline
does not forward YET (the single-cell tags, ``--stream_reads``), the forward is
simulated by wrapping ``lraa_cmd``: that is the property those forwards need in
place BEFORE they land, and it is what makes an experiment on them trustworthy
instead of answerable from a stale sentinel.

Neither stage runs a subprocess here: ``run_step`` is replaced, and each function
is allowed to fail on its missing outputs once it has produced its sentinels,
which is the point at which the tokens are observable.
"""

import argparse
import os
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "pylib"))

import ChunkedRun  # noqa: E402  (path insert must precede the import)


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
    discovery=False,
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


class _RecordingCheckpoints(ChunkedRun.Checkpoints):
    """Every sentinel the stage asks about, in order, and never a hit."""

    def __init__(self, root, marked):
        super().__init__(root)
        self.marked = marked

    def done(self, token):
        self.marked.append(token)
        return False

    def mark(self, token, note=""):
        pass


def _stage(tmp_path, cmd_hook, run):
    """Run one stage with no subprocesses; return the sentinels it consulted.

    ``cmd_hook`` rewrites what ``lraa_cmd`` hands the stage, which is how a flag
    the pipeline does not forward yet -- or a reordering of the one it does --
    is put in front of the token derivation without inventing a second builder.
    """
    marked = []
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(ChunkedRun, "run_step", lambda *a, **k: {"step": "stub"})
        if cmd_hook is not None:
            real_lraa_cmd = ChunkedRun.lraa_cmd
            mp.setattr(
                ChunkedRun,
                "lraa_cmd",
                lambda *a, **k: cmd_hook(list(real_lraa_cmd(*a, **k))),
            )
        with pytest.raises(ChunkedRun.PipelineError):
            run(_RecordingCheckpoints(str(tmp_path / "ckpt"), marked))
    return marked


CHUNK_ID = "chr20_plus_00"


def chunk_quant_token(tmp_path, cmd_hook=None, **kw):
    """The stage-5 sentinel for one strand-first chunk with one quant unit."""
    cdir = tmp_path / "chunk"
    cdir.mkdir(exist_ok=True)
    chunk = {
        "chunk_id": CHUNK_ID,
        "strandless": False,
        "dir": str(cdir),
        "prefix": str(cdir / "chunk"),
        "log": str(cdir / "chunk.log"),
        "window_origin": 0,
        "upstream_token": "stage3.up_aaaaaaaaaaaa",
        "units": ChunkedRun.chunk_quant_units(
            CHUNK_ID, str(cdir), str(cdir / "chunk"), "+", 0, 0
        ),
    }
    marked = _stage(
        tmp_path,
        cmd_hook,
        lambda ckpt: ChunkedRun.chunk_worker(
            _args(tmp_path, **kw),
            ckpt,
            str(tmp_path / "out"),
            chunk,
            1000,
            0.5,
            4,
        ),
    )
    assert len(marked) == 2, marked
    return marked[1]


def baseline_quant_token(tmp_path, cmd_hook=None, **kw):
    """The baseline arm's quant sentinel -- the second site with the same defect."""
    outdir = tmp_path / "out"
    (outdir / "logs").mkdir(parents=True, exist_ok=True)
    marked = _stage(
        tmp_path,
        cmd_hook,
        lambda ckpt: ChunkedRun.run_baseline(
            _args(tmp_path, **kw),
            ckpt,
            str(outdir),
            {},
            {"+": str(tmp_path / "plus.bam"), "-": str(tmp_path / "minus.bam")},
            1000,
            0.5,
            "stage1.up_aaaaaaaaaaaa",
        ),
    )
    assert len(marked) == 3, marked
    return marked[2]


ARMS = (
    pytest.param(chunk_quant_token, "stage5_quant_" + CHUNK_ID + ".", id="chunk"),
    pytest.param(baseline_quant_token, "baseline_quant.", id="baseline"),
)


def _append(*extra):
    return lambda cmd: cmd + list(extra)


def _mutate_at(pos):
    """Rewrite one argv position, whatever it is -- flag, value or positional."""
    return lambda cmd: cmd[:pos] + ["MUTATED"] + cmd[pos + 1 :]


def _captured_command(tmp_path, token_of):
    """The argv an arm really hands its digest, so a test can enumerate IT.

    The baseline appends ``--contig`` after ``lraa_cmd`` returns, which this does
    not see; that one argument has its own test.
    """
    seen = []
    token_of(tmp_path, cmd_hook=lambda cmd: seen.append(list(cmd)) or cmd)
    assert len(seen) == 1, seen
    return seen[0]


# ------------------------------------------------------- the digest by itself


def test_flag_order_is_inert_but_the_pairing_is_not():
    """Sorting the bare token stream would call these three commands the same."""
    a = ChunkedRun.argv_digest(["prog", "--a", "1", "--b", "2"])
    reordered = ChunkedRun.argv_digest(["prog", "--b", "2", "--a", "1"])
    swapped = ChunkedRun.argv_digest(["prog", "--a", "2", "--b", "1"])
    assert a == reordered
    assert a != swapped


def test_positional_order_is_not_inert():
    """argv[0] and the script are arguments of each other, not of a flag."""
    assert ChunkedRun.argv_digest(["python", "LRAA"]) != ChunkedRun.argv_digest(
        ["LRAA", "python"]
    )


def test_a_bare_flag_is_not_swallowed_by_the_next_flags_value():
    """``--HiFi --gtf x`` must not digest as ``--HiFi=--gtf`` plus a stray token."""
    assert ChunkedRun.argv_digest(["p", "--HiFi", "--gtf", "x"]) == ChunkedRun.argv_digest(
        ["p", "--gtf", "x", "--HiFi"]
    )
    assert ChunkedRun.argv_digest(["p", "--HiFi", "--gtf", "x"]) != ChunkedRun.argv_digest(
        ["p", "--gtf", "x"]
    )


def test_a_negative_value_stays_attached_to_its_flag():
    assert ChunkedRun.argv_digest(["p", "--offset", "-1"]) != ChunkedRun.argv_digest(
        ["p", "--offset", "1"]
    )


# ------------------------------------------------------------- both quant arms


@pytest.mark.parametrize("token_of,prefix", ARMS)
def test_an_unchanged_invocation_keeps_its_sentinel(tmp_path, token_of, prefix):
    """Resumability is the whole point of the sentinel; a digest must not cost it."""
    assert token_of(tmp_path) == token_of(tmp_path)


@pytest.mark.parametrize("token_of,prefix", ARMS)
def test_the_sentinel_is_still_readable_at_a_glance(tmp_path, token_of, prefix):
    """The digest is an ADDITIONAL component: a directory of hashes is its own hazard."""
    token = token_of(tmp_path)
    assert token.startswith(prefix), token
    assert os.path.basename(token) == token


@pytest.mark.parametrize("token_of,prefix", ARMS)
def test_quant_only_and_discovery_do_not_share_a_sentinel(tmp_path, token_of, prefix):
    """The mode distinction, absent from the old key entirely.

    Nothing is added for ``args.discovery`` at either site: it is exactly the
    presence of ``--quant_only`` in the argv, so this fails if that stops being
    true as much as if the argv were dropped from the key.
    """
    assert token_of(tmp_path, discovery=False) != token_of(tmp_path, discovery=True)


@pytest.mark.parametrize("token_of,prefix", ARMS)
def test_a_forwarded_single_cell_tag_moves_the_sentinel(tmp_path, token_of, prefix):
    """The failure this item exists for: bulk output served for a single-cell run.

    The tags are not forwarded yet. This asserts the property that has to hold
    the moment they are -- the sentinel follows whatever ``lraa_cmd`` emits --
    so the experiment that forwards them cannot be answered from a bulk
    checkpoint.
    """
    bulk = token_of(tmp_path)
    single_cell = token_of(
        tmp_path,
        cmd_hook=_append("--cell_barcode_tag", "CB", "--read_umi_tag", "XM"),
    )
    assert bulk != single_cell


@pytest.mark.parametrize("token_of,prefix", ARMS)
def test_a_partial_single_cell_forward_moves_it_too(tmp_path, token_of, prefix):
    """``get_read_name_include_sc_encoding`` needs BOTH tags, so a partial forward
    silently produces bulk. It must still be a distinct sentinel from both the
    bulk run and the complete forward, or one of the three serves the others."""
    bulk = token_of(tmp_path)
    partial = token_of(tmp_path, cmd_hook=_append("--cell_barcode_tag", "CB"))
    both = token_of(
        tmp_path,
        cmd_hook=_append("--cell_barcode_tag", "CB", "--read_umi_tag", "XM"),
    )
    assert len({bulk, partial, both}) == 3


@pytest.mark.parametrize("token_of,prefix", ARMS)
def test_stream_reads_moves_the_sentinel(tmp_path, token_of, prefix):
    assert token_of(tmp_path) != token_of(tmp_path, cmd_hook=_append("--stream_reads"))


@pytest.mark.parametrize("token_of,prefix", ARMS)
def test_reordering_the_command_alone_does_not_move_the_sentinel(
    tmp_path, token_of, prefix
):
    """A cosmetic edit to ``lraa_cmd``'s flag order must not orphan every checkpoint.

    The permutation is deliberate rather than a shuffle: the positionals stay put
    (they are ordered) and the flag/value pairs are emitted last-to-first.
    """

    def reverse_flags(cmd):
        head = [tok for tok in cmd[:2]]
        rest, pairs, i = cmd[2:], [], 0
        while i < len(rest):
            if i + 1 < len(rest) and not rest[i + 1].startswith("--"):
                pairs.append(rest[i : i + 2])
                i += 2
            else:
                pairs.append(rest[i : i + 1])
                i += 1
        return head + [tok for pair in reversed(pairs) for tok in pair]

    assert token_of(tmp_path) == token_of(tmp_path, cmd_hook=reverse_flags)


@pytest.mark.parametrize("token_of,prefix", ARMS)
def test_a_changed_output_prefix_moves_the_sentinel(tmp_path, token_of, prefix):
    """The regression guard: this is what fails if the argv leaves the key.

    ``--output_prefix`` decides which files the step writes and is named by no
    other component of either token -- not by the stage prefix, not by the
    inherited normalization token -- so a key rebuilt from a hand-listed set of
    settings passes most of the file and fails here.
    """

    def rename(cmd):
        out = list(cmd)
        out[out.index("--output_prefix") + 1] = "somewhere_else"
        return out

    assert token_of(tmp_path) != token_of(tmp_path, cmd_hook=rename)


@pytest.mark.parametrize("token_of,prefix", ARMS)
def test_the_cpu_budget_value_is_gated_out_of_the_sentinel(tmp_path, token_of, prefix):
    """Resuming at a lower concurrency must not discard finished work.

    ``--cpu_budget`` is this invocation's SHARE, and a chunk is one quant unit
    whose internal pool clamps to 1 regardless, so it determines nothing about
    the output. A signature naming a non-determining input is as defective as one
    missing a field -- and concretely, a whole-genome run OOM-killed at 16 and
    resumed at 8 would redo every finished unit. This fails if the gating in
    ``quant_command_digest`` is removed.
    """

    def rebudget(cmd):
        out = list(cmd)
        out[out.index("--cpu_budget") + 1] = "99"
        return out

    assert token_of(tmp_path) == token_of(tmp_path, cmd_hook=rebudget)


@pytest.mark.parametrize("token_of,prefix", ARMS)
def test_dropping_cpu_budget_entirely_still_moves_the_sentinel(
    tmp_path, token_of, prefix
):
    """The VALUE is gated, not the flag.

    Gating by deleting the pair would make "budget forwarded" and "budget not
    forwarded" the same token, so a stage that stopped receiving it would be
    invisible. The flag stays in the payload with an inert value instead, which
    is the ``onone`` convention the normalize cache stem uses for an unset field.
    """

    def unbudget(cmd):
        out = list(cmd)
        at = out.index("--cpu_budget")
        del out[at : at + 2]
        return out

    assert token_of(tmp_path) != token_of(tmp_path, cmd_hook=unbudget)


@pytest.mark.parametrize("token_of,prefix", ARMS)
def test_cpu_budget_is_the_only_gated_field(tmp_path, token_of, prefix):
    """Derived from the real command rather than from a list restated here.

    Every argument the stage issues -- interpreter, script, genome, both bams,
    the gtf, the denominator, the output prefix, --quant_only, --no_norm and the
    --cpu_budget FLAG itself -- is mutated in turn, and exactly one position, the
    budget's value, may leave the token where it was. A second gated field, or
    the argv leaving the key altogether, fails this.
    """
    real = _captured_command(tmp_path, token_of)
    unmutated = token_of(tmp_path)
    inert = [
        pos
        for pos in range(len(real))
        if token_of(tmp_path, cmd_hook=_mutate_at(pos)) == unmutated
    ]
    assert inert == [real.index("--cpu_budget") + 1], [
        (pos, real[pos]) for pos in inert
    ]


@pytest.mark.parametrize("token_of,prefix", ARMS)
def test_hifi_still_moves_the_sentinel(tmp_path, token_of, prefix):
    """It was named explicitly before; deriving must not lose what enumerating had."""
    assert token_of(tmp_path, HiFi=False) != token_of(tmp_path, HiFi=True)


@pytest.mark.parametrize("token_of,prefix", ARMS)
def test_the_upstream_normalization_token_still_reaches_it(tmp_path, token_of, prefix):
    """The chain is unchanged: a stage-4 input change still invalidates stage 5.

    ``--depth_window`` reaches the normalizer and NOT the quant argv, so only the
    inherited parent token can carry it.
    """
    assert token_of(tmp_path) != token_of(tmp_path, depth_window=250)


# ---------------------------------------------------- the baseline arm's tail


def test_the_baseline_contig_restriction_is_in_its_sentinel(tmp_path):
    """``--contig`` is appended AFTER lraa_cmd returns, so the token must be built
    after the append. Restricting the control to one contig changes what it
    quantifies, and the old key did not name it at all."""
    assert baseline_quant_token(tmp_path) != baseline_quant_token(
        tmp_path, contig="chr20"
    )


def test_the_two_arms_do_not_share_a_sentinel(tmp_path):
    """Same digest function, different stage prefixes and different commands."""
    assert chunk_quant_token(tmp_path) != baseline_quant_token(tmp_path)
