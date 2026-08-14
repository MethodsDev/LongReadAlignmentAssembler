#!/usr/bin/env python3

"""What util/normalize_bam_by_strand.py must guarantee to be usable per chunk.

Two things, and they are independent.

Checkpoint tokens. Every artifact this script writes is checkpointed by the
Pipeliner and a checkpoint is trusted on sight, so a token has to name everything
that determines the contents it guards. The strand split discards records having
an intron longer than max_intron_length, so both the split checkpoint and the run
token that gates sampling, merge and index depend on that threshold: without it,
a rerun at a different threshold would be answered with the previous run's strand
bams and silently reproduce the previous run's result. --window_origin moves the
depth-window boundaries and so decides which reads survive sampling, and
--input_is_single_strand decides whether the split runs at all, so both belong in
the tokens on the same argument. Neither may perturb the tokens when it is not
supplied, because that would orphan every existing cache.

The window grid. The grid used to be anchored on the first aligned base present
in the input, which makes window boundaries a function of the input's
composition: a chunk of a contig begins at some other read, so the same absolute
locus bins differently than it did in a whole-contig run, the window medians move,
the acceptance probability moves with them, and reads near the threshold flip.
An explicit origin pins the grid to the absolute coordinate grid instead, which
is what lets normalization move inside a chunk. The chunk used here comes from
the real extractor, util/misc/extract_contig_region_inputs.py, and the origin is
the offset that tool reports, because the equivalence rests on the origin being
that offset under the extractor's own coordinate convention -- an off-by-one
there would shift every window by a base and flip reads with nothing erroring.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

UTIL_DIR = Path(__file__).resolve().parent
REPO_ROOT = UTIL_DIR.parent
for _import_dir in (str(UTIL_DIR), str(REPO_ROOT / "pylib")):
    if _import_dir not in sys.path:
        sys.path.insert(0, _import_dir)

import LRAA_Globals
import Util_funcs
import normalize_bam_by_strand as nbs

NORMALIZER = UTIL_DIR / "normalize_bam_by_strand.py"
EXTRACTOR = UTIL_DIR / "misc" / "extract_contig_region_inputs.py"

BASE_ARGS = dict(
    normalize_max_cov_level=1000,
    depth_window=100,
    random_seed=42,
    window_origin=None,
    input_is_single_strand=False,
)


def _tokens(tmp_path, max_intron_length, **overrides):
    kwargs = dict(BASE_ARGS)
    kwargs.update(overrides)

    return nbs.compute_tokens(
        str(tmp_path / "input.bam"),
        str(tmp_path / "output.bam"),
        max_intron_length,
        kwargs["normalize_max_cov_level"],
        kwargs["depth_window"],
        kwargs["random_seed"],
        window_origin=kwargs["window_origin"],
        input_is_single_strand=kwargs["input_is_single_strand"],
    )


def test_tokens_are_stable_for_identical_inputs(tmp_path):
    assert _tokens(tmp_path, 200000) == _tokens(tmp_path, 200000)


def test_changing_max_intron_length_invalidates_both_tokens(tmp_path):
    split_token, run_token = _tokens(tmp_path, 200000)
    other_split_token, other_run_token = _tokens(tmp_path, 500000)

    # the split emits different records, so its checkpoint name must differ ...
    assert split_token != other_split_token
    # ... and so must the token gating everything consuming the split
    assert run_token != other_run_token


def test_disabling_intron_filtering_invalidates_both_tokens(tmp_path):
    filtered = _tokens(tmp_path, 200000)
    disabled = _tokens(tmp_path, 0)

    assert filtered[0] != disabled[0]
    assert filtered[1] != disabled[1]


def test_other_run_inputs_still_invalidate_the_run_token(tmp_path):
    split_token, run_token = _tokens(tmp_path, 200000)

    for override in (
        {"normalize_max_cov_level": 500},
        {"depth_window": 50},
        {"random_seed": 7},
    ):
        other_split_token, other_run_token = _tokens(tmp_path, 200000, **override)
        # these do not change what the split emits, only what consumes it
        assert other_split_token == split_token
        assert other_run_token != run_token


def test_window_origin_invalidates_both_tokens(tmp_path):
    """The origin decides which reads survive, so it decides the contents."""
    default_split, default_run = _tokens(tmp_path, 200000)
    at_zero = _tokens(tmp_path, 200000, window_origin=0)
    at_a_chunk = _tokens(tmp_path, 200000, window_origin=1000)

    assert at_zero[0] != default_split
    assert at_zero[1] != default_run
    assert at_a_chunk[0] != default_split
    assert at_a_chunk[1] != default_run
    # and one origin must not be answered with another origin's artifacts
    assert at_zero != at_a_chunk


def test_single_strand_flag_invalidates_both_tokens(tmp_path):
    """It changes what stands in for the split's output, so the split token moves too."""
    default_split, default_run = _tokens(tmp_path, 200000)
    single_split, single_run = _tokens(tmp_path, 200000, input_is_single_strand=True)

    assert single_split != default_split
    assert single_run != default_run


def test_the_two_new_options_are_distinguished_from_each_other(tmp_path):
    """Neither may be mistaken for the other, nor for both together."""
    origin_only = _tokens(tmp_path, 200000, window_origin=0)
    single_only = _tokens(tmp_path, 200000, input_is_single_strand=True)
    both = _tokens(tmp_path, 200000, window_origin=0, input_is_single_strand=True)

    assert len({origin_only, single_only, both}) == 3


def test_tokens_are_untouched_when_the_new_options_are_absent(tmp_path):
    """A run that supplies neither option must still hit the caches it has.

    The two fields enter the hashed string only when set, so the default
    invocation hashes what it hashed before the options existed. Recomputed here
    from the recipe rather than compared against a checked-in digest, so the
    assertion still names which inputs are supposed to be in there.
    """
    source_token = Util_funcs.file_identity_token(str(tmp_path / "input.bam"))
    expected_split = Util_funcs.get_hash_code("|".join([source_token, "200000"]))[:12]
    expected_run = Util_funcs.get_hash_code(
        "|".join(
            [
                source_token,
                "200000",
                str(BASE_ARGS["normalize_max_cov_level"]),
                str(BASE_ARGS["depth_window"]),
                str(BASE_ARGS["random_seed"]),
                nbs.METHOD,
                os.path.realpath(str(tmp_path / "output.bam")),
            ]
        )
    )[:12]

    assert _tokens(tmp_path, 200000) == (expected_split, expected_run)


def test_help_reports_max_intron_length_with_its_default():
    help_text = subprocess.run(
        [sys.executable, str(NORMALIZER), "--help"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout

    assert "--max_intron_length" in help_text
    assert str(LRAA_Globals.config["max_intron_length"]) in help_text
    assert "disable intron length filtering" in help_text
    assert "--input_is_single_strand" in help_text
    assert "--window_origin" in help_text


# ---------------------------------------------------------------------------
# the window grid, end to end

CONTIG = "chr1"
CONTIG_LEN = 3000
DEPTH_WINDOW = 100
# low enough that the piles below are genuinely thinned: an equivalence between
# two runs that both kept everything would prove nothing about the grid
TARGET_DEPTH = 5
# 1-based last base of the left chunk, a multiple of DEPTH_WINDOW as the cut rule
# requires, so no window straddles it
CUT = 1000


def _bam_header(length=CONTIG_LEN):
    return pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": CONTIG, "LN": length}],
        }
    )


def _aligned_read(header, name, start, length):
    """A forward, primary, unspliced alignment.

    Unspliced on purpose: a read carrying a junction below the target is kept
    unconditionally, which would mask the grid effect being measured.
    """
    read = pysam.AlignedSegment(header)
    read.query_name = name
    read.query_sequence = "A" * length
    read.query_qualities = pysam.qualitystring_to_array("I" * length)
    read.flag = 0
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = 60
    read.cigartuples = [(0, length)]
    return read


def _corpus(header, shift=0):
    """Two piles, one either side of the cut, each above the target depth.

    Both start at a coordinate that is not a multiple of DEPTH_WINDOW, and the
    right-hand pile is 250 bp into its chunk, so the first-read anchor and the
    absolute grid are half a window apart and the two disagree visibly.
    """
    reads = [_aligned_read(header, "L%04d" % i, 350 + shift, 200) for i in range(40)]
    reads += [_aligned_read(header, "R%04d" % i, 1250 + shift, 200) for i in range(60)]
    return reads


def _write_bam(path, reads, header):
    with pysam.AlignmentFile(str(path), "wb", header=header) as writer:
        for read in sorted(reads, key=lambda r: r.reference_start):
            writer.write(read)
    pysam.index(str(path))


def _write_genome(path):
    with open(str(path), "wt") as ofh:
        print(">" + CONTIG, file=ofh)
        sequence = "ACGT" * (CONTIG_LEN // 4)
        for offset in range(0, len(sequence), 60):
            print(sequence[offset : offset + 60], file=ofh)
    pysam.faidx(str(path))


def _retained(path):
    """read name -> sampling weight, for the records that survived"""
    with pysam.AlignmentFile(str(path), "rb") as reader:
        return {
            read.query_name: round(read.get_tag("XW"), 6)
            for read in reader.fetch(until_eof=True)
        }


def _run_normalizer(input_bam, output_bam, cwd, window_origin=None, single_strand=False):
    cmd = [
        sys.executable,
        str(NORMALIZER),
        "--input_bam",
        str(input_bam),
        "--output_bam",
        str(output_bam),
        "--normalize_max_cov_level",
        str(TARGET_DEPTH),
        "--depth_window",
        str(DEPTH_WINDOW),
    ]
    if window_origin is not None:
        cmd += ["--window_origin", str(window_origin)]
    if single_strand:
        cmd += ["--input_is_single_strand"]

    completed = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    return completed


def _extract_chunk(work_dir, genome_fa, whole_bam, region, tag):
    """One rebased chunk, from the tool the partitioned pipeline actually uses."""
    chunk_dir = work_dir / tag
    chunk_dir.mkdir()
    completed = subprocess.run(
        [
            sys.executable,
            str(EXTRACTOR),
            "--genome_fa",
            str(genome_fa),
            "--bam",
            str(whole_bam),
            "--region",
            region,
            # the piles sit well clear of both boundaries; 0 asks only that
            # nothing crosses, which is the condition this equivalence needs
            "--margin",
            "0",
            "--output_prefix",
            str(chunk_dir / tag),
        ],
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    return json.loads((chunk_dir / (tag + ".partition.json")).read_text())


@pytest.fixture(scope="module")
def chunked(tmp_path_factory):
    """A whole-contig run and the same contig cut in two, normalized per chunk.

    The control is run WITH an explicit origin, because that is the only baseline
    a chunk can be compared against: left on its default the control would keep
    the first-aligned-base anchor while the chunks used the absolute grid, and the
    comparison would fail for a reason that has nothing to do with chunking.
    """
    work = tmp_path_factory.mktemp("chunked")
    header = _bam_header()
    whole_bam = work / "whole.bam"
    _write_bam(whole_bam, _corpus(header), header)
    genome_fa = work / "genome.fa"
    _write_genome(genome_fa)

    manifests = {
        "left": _extract_chunk(
            work, genome_fa, whole_bam, "{}+:1-{}".format(CONTIG, CUT), "left"
        ),
        "right": _extract_chunk(
            work,
            genome_fa,
            whole_bam,
            "{}+:{}-{}".format(CONTIG, CUT + 1, CONTIG_LEN),
            "right",
        ),
    }

    baseline_dir = work / "whole_with_origin"
    baseline_dir.mkdir()
    baseline_bam = baseline_dir / "whole.norm.bam"
    _run_normalizer(whole_bam, baseline_bam, baseline_dir, window_origin=0)

    default_dir = work / "whole_default"
    default_dir.mkdir()
    default_bam = default_dir / "whole.norm.bam"
    _run_normalizer(whole_bam, default_bam, default_dir)

    with_origin = dict()
    without_origin = dict()
    for tag, manifest in manifests.items():
        chunk_bam = Path(manifest["files"]["bam"])

        run_dir = work / (tag + "_with_origin")
        run_dir.mkdir()
        out = run_dir / "chunk.norm.bam"
        _run_normalizer(
            chunk_bam,
            out,
            run_dir,
            window_origin=manifest["offset"],
            single_strand=True,
        )
        with_origin.update(_retained(out))

        bare_dir = work / (tag + "_without_origin")
        bare_dir.mkdir()
        bare_out = bare_dir / "chunk.norm.bam"
        _run_normalizer(chunk_bam, bare_out, bare_dir, single_strand=True)
        without_origin.update(_retained(bare_out))

    return {
        "manifests": manifests,
        "whole_bam": whole_bam,
        "baseline": _retained(baseline_bam),
        "whole_default": _retained(default_bam),
        "chunks_with_origin": with_origin,
        "chunks_without_origin": without_origin,
    }


def test_the_origin_to_pass_is_the_extractors_own_offset(chunked):
    """The equivalence rests on this, and an off-by-one here would not raise.

    The extractor rebases a 0-based reference_start by subtracting `offset`, so
    the absolute coordinate of chunk position 0 IS `offset` and that is what
    --window_origin takes. Asserted against the real manifest, not against the
    convention this test would like the extractor to have.
    """
    left, right = chunked["manifests"]["left"], chunked["manifests"]["right"]

    assert left["offset"] == 0
    assert right["offset"] == CUT
    assert right["offset"] % DEPTH_WINDOW == 0, "a window would straddle the cut"

    with pysam.AlignmentFile(str(right["files"]["bam"]), "rb") as reader:
        rebased = {r.query_name: r.reference_start for r in reader.fetch(until_eof=True)}
    with pysam.AlignmentFile(str(chunked["whole_bam"]), "rb") as reader:
        absolute = {r.query_name: r.reference_start for r in reader.fetch(until_eof=True)}

    assert rebased
    for name, start in rebased.items():
        assert start == absolute[name] - right["offset"]


def test_an_explicit_origin_bins_a_chunk_as_the_whole_contig_does(chunked):
    """The point of the option: same locus, same window, same fate.

    Every record of the contig is in exactly one chunk, so the chunk runs put
    together must retain exactly what the whole-contig run retained -- and carry
    the same sampling weights, which is the stronger statement, since the weight
    is 1/p and p is what the window medians decide.
    """
    assert chunked["chunks_with_origin"] == chunked["baseline"]
    # non-trivial: both sides actually thinned something
    assert 0 < len(chunked["baseline"]) < 100


def test_without_an_explicit_origin_a_chunk_bins_differently(chunked):
    """Why the option has to exist, and why the test above has teeth.

    Same chunk bams, same seed, same target: only the grid moves, because the
    anchor becomes the chunk's first read instead of the absolute grid. Reads
    near the threshold flip, and nothing reports it.
    """
    assert chunked["chunks_without_origin"] != chunked["baseline"]
    assert set(chunked["chunks_without_origin"]) != set(chunked["baseline"])


def test_the_whole_contig_run_can_take_the_origin_too(chunked):
    """The control arm needs the same option the chunks use, and it is a real one.

    Supplying it to a whole-contig run changes the answer here, which is exactly
    why the baseline for a chunk comparison has to be the explicit-origin run.
    """
    assert set(chunked["whole_default"]) != set(chunked["baseline"])


def test_the_default_grid_still_follows_the_first_aligned_base(tmp_path):
    """Unchanged behaviour for every caller that supplies no origin.

    The defining property of the old anchor is that it travels with the reads:
    translate a locus and the grid translates too, so the same reads are kept.
    Asserted as the property rather than as a digest, and paired with the same
    two inputs judged on a fixed grid, which must NOT agree -- otherwise this
    would pass on a corpus that is simply insensitive to the grid.
    """
    header = _bam_header()
    at_origin = tmp_path / "at_origin.bam"
    _write_bam(at_origin, _corpus(header), header)
    translated = tmp_path / "translated.bam"
    _write_bam(translated, _corpus(header, shift=37), header)

    outputs = dict()
    for tag, source, origin in (
        ("default_a", at_origin, None),
        ("default_b", translated, None),
        ("fixed_a", at_origin, 0),
        ("fixed_b", translated, 0),
    ):
        out = tmp_path / (tag + ".norm.bam")
        nbs.sift_bam(
            str(source),
            str(out),
            TARGET_DEPTH,
            DEPTH_WINDOW,
            42,
            window_origin=origin,
        )
        outputs[tag] = _retained(out)

    assert outputs["default_a"] == outputs["default_b"]
    assert outputs["fixed_a"] != outputs["fixed_b"]


# ---------------------------------------------------------------------------
# skipping the internal strand split


@pytest.fixture(scope="module")
def single_strand(tmp_path_factory):
    """One orientation-pure bam, normalized both ways."""
    work = tmp_path_factory.mktemp("single_strand")
    header = _bam_header()
    source = work / "plus_only.bam"
    _write_bam(source, _corpus(header), header)

    split_dir = work / "via_split"
    split_dir.mkdir()
    split_out = split_dir / "norm.bam"
    _run_normalizer(source, split_out, split_dir)

    direct_dir = work / "direct"
    direct_dir.mkdir()
    direct_out = direct_dir / "norm.bam"
    _run_normalizer(source, direct_out, direct_dir, single_strand=True)

    return {
        "source": source,
        "split_dir": split_dir,
        "split_out": split_out,
        "direct_dir": direct_dir,
        "direct_out": direct_out,
        "via_split": _retained(split_out),
        "direct": _retained(direct_out),
    }


def test_single_strand_retains_what_the_split_path_retains(single_strand):
    """Skipping the split may not change the answer, only the work done."""
    assert single_strand["direct"] == single_strand["via_split"]
    # thinning happened, so this is not two empty sets agreeing
    assert 0 < len(single_strand["direct"]) < 100


def test_single_strand_writes_no_opposite_strand_artifact(single_strand):
    """The files the split path leaves behind, one of them guaranteed empty."""
    split_names = {p.name for p in single_strand["split_dir"].iterdir()}
    assert any(name.endswith(".SS.-.bam") for name in split_names)
    assert any(name.endswith(".SS.+.bam") for name in split_names)

    direct_names = {p.name for p in single_strand["direct_dir"].iterdir()}
    assert not [name for name in direct_names if ".SS." in name]


def test_single_strand_checkpoints_name_what_produced_them(single_strand):
    """A checkpoint is trusted on sight, so it may not describe a stage that never ran."""
    _, run_token = nbs.compute_tokens(
        str(single_strand["source"]),
        str(single_strand["direct_out"]),
        LRAA_Globals.config["max_intron_length"],
        TARGET_DEPTH,
        DEPTH_WINDOW,
        BASE_ARGS["random_seed"],
        input_is_single_strand=True,
    )

    checkpoints = {p.name for p in (single_strand["direct_dir"] / "__chckpts").iterdir()}
    # no split ran and nothing was merged, so no checkpoint may claim otherwise
    assert not [name for name in checkpoints if "sep_by_strand" in name]
    assert not [name for name in checkpoints if "merge" in name]
    assert "index_single_strand_norm.{}.ok".format(run_token) in checkpoints

    sift_checkpoint = "{}.norm_{}.{}.ok".format(
        single_strand["direct_out"].name, TARGET_DEPTH, run_token
    )
    assert (single_strand["direct_dir"] / sift_checkpoint).exists()


def test_repeating_a_single_strand_run_reuses_its_checkpoints(single_strand):
    """Same inputs, same tokens: the second run must do nothing at all."""
    before = {
        p.name: p.stat().st_mtime_ns
        for p in single_strand["direct_dir"].rglob("*")
        if p.is_file()
    }

    _run_normalizer(
        single_strand["source"],
        single_strand["direct_out"],
        single_strand["direct_dir"],
        single_strand=True,
    )

    after = {
        p.name: p.stat().st_mtime_ns
        for p in single_strand["direct_dir"].rglob("*")
        if p.is_file()
    }

    assert after == before
