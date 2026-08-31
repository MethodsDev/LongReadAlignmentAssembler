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
import re
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
# A real multi-contig bam: 11,763 records over chr1, chr2, chr3 and chrM, of very
# unequal size (5,627 / 4,510 / 15 / 1,611) and carrying alignments that overhang
# their contig by up to 238 kb.  A synthetic corpus can be made to look like this
# but cannot be relied on to; the per-contig fan-out is exactly the change a real
# corpus's asymmetries break.
REAL_MULTI_CONTIG_BAM = (
    REPO_ROOT
    / "testing"
    / "single_cells"
    / "sc_full_pipe_scattered"
    / "data"
    / "chr1_2_3_M.PBMCs.mini.bam"
)

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


def test_intron_length_does_not_key_a_run_whose_split_never_runs(tmp_path):
    """--input_is_single_strand skips the split, so the intron cap decides nothing.

    The split is the only place the threshold is enforced -- depth measurement passes
    max_intron_length=0 explicitly, and on this path removing long-intron alignments is
    the caller's documented responsibility. So two single-strand runs differing only in
    --max_intron_length produce byte-identical output, and naming the threshold would key
    every artifact on an input that cannot change them.

    Over-keying is not the safe direction. It costs rebuilds, and it makes the name
    unauditable: a field that sometimes determines the contents and always appears cannot
    be checked against behaviour. This test fails if the gating is removed, which is the
    only thing keeping the two claims in step.
    """
    filtered = _tokens(tmp_path, 200000, input_is_single_strand=True)
    disabled = _tokens(tmp_path, 0, input_is_single_strand=True)
    wildly_different = _tokens(tmp_path, 1, input_is_single_strand=True)

    assert filtered == disabled == wildly_different

    # and it must still key the mode where the split does enforce it
    assert _tokens(tmp_path, 200000) != _tokens(tmp_path, 0)


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


def _run_normalizer(
    input_bam,
    output_bam,
    cwd,
    window_origin=None,
    single_strand=False,
    num_workers=None,
):
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
    if num_workers is not None:
        cmd += ["--num_workers", str(num_workers)]

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


# ---------------------------------------------------------------------------
# alignments the header does not fit

# far enough past CONTIG_LEN that the windows these reads occupy fall outside
# the array sized from the header (3000 // 100 + 2 = 32 entries)
PAST_THE_END = 6500


def _write_bam_unindexed(path, reads, header, sort=True):
    """Write records pysam can hold but samtools would refuse to index.

    An alignment reaching past the reference length, or an input out of
    coordinate order, is exactly what these tests are about, so the index step
    `_write_bam` performs is not available to them.
    """
    ordered = sorted(reads, key=lambda r: r.reference_start) if sort else reads
    with pysam.AlignmentFile(str(path), "wb", header=header) as writer:
        for read in ordered:
            writer.write(read)


def test_an_alignment_past_the_contig_end_is_measured_rather_than_dropped(tmp_path):
    """The depth an overhanging pile carries has to reach its own acceptance test.

    Coordinate-remapped bams routinely carry alignments running off the contig
    they were rebased onto -- 18% of the records in testing/sep_contigs do. Pass
    1 sized its depth array from the header and wrote to it unchecked, so such a
    record raised IndexError and took the whole run with it.

    Not crashing is the smaller half. Bounding the write by skipping the window
    instead would leave this pile unmeasured, pass 2 would find no depth for it
    and return probability 1, and forty reads over a target of five would be
    retained whole -- a silent normalization failure in place of a loud crash.
    So the assertion is that the pile is *thinned*, and the anchor read, whose
    depth is genuinely below target, is not.
    """
    header = _bam_header()
    anchor_read = _aligned_read(header, "anchor", 350, 200)
    past = [
        _aligned_read(header, "P%04d" % i, PAST_THE_END, 100)
        for i in range(40)
    ]

    source = tmp_path / "overhang.bam"
    _write_bam_unindexed(source, [anchor_read] + past, header)

    completed = _run_normalizer(
        source, tmp_path / "overhang.norm.bam", tmp_path, single_strand=True
    )

    retained = _retained(tmp_path / "overhang.norm.bam")
    survivors = [name for name in retained if name.startswith("P")]

    assert 0 < len(survivors) < len(past)
    assert all(retained[name] > 1.0 for name in survivors)
    assert retained["anchor"] == 1.0

    assert "extend past the reference length" in completed.stderr
    assert str(len(past)) in completed.stderr


def test_a_read_before_the_anchor_is_refused_while_the_anchor_depends_on_order(tmp_path):
    """Refuse where the default anchor is load-bearing, and only there.

    With no origin the grid is anchored on the first aligned base seen, so a
    read starting before it yields a negative window -- which Python resolves
    against the array tail, corrupting an unrelated locus rather than failing.
    Given an absolute origin the anchor sits at or below zero, no aligned
    position can precede it, and the same input bins perfectly well; refusing it
    there would reject every legitimate chunk. Both halves are asserted so that
    removing the gate fails a test.

    The header claims `SO:unsorted` here, but the check reads the records: a bam
    that merely says it is sorted must not get past it.
    """
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "unsorted"},
            "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}],
        }
    )
    reads = [
        _aligned_read(header, "late", 2000, 100),
        _aligned_read(header, "early", 100, 100),
    ]

    source = tmp_path / "unsorted.bam"
    _write_bam_unindexed(source, reads, header, sort=False)

    refused = subprocess.run(
        [
            sys.executable,
            str(NORMALIZER),
            "--input_bam",
            str(source),
            "--output_bam",
            str(tmp_path / "unsorted.norm.bam"),
            "--normalize_max_cov_level",
            str(TARGET_DEPTH),
            "--depth_window",
            str(DEPTH_WINDOW),
            "--input_is_single_strand",
        ],
        cwd=str(tmp_path),
        capture_output=True,
        text=True,
    )

    assert refused.returncode != 0
    assert "coordinate-sorted" in refused.stderr
    assert not (tmp_path / "unsorted.norm.bam").exists()

    # the cli indexes what it writes, which unsorted output cannot support, so
    # the inert half is asserted against the function that owns the invariant
    accepted = tmp_path / "with_origin.norm.bam"
    nbs.sift_bam(
        str(source),
        str(accepted),
        TARGET_DEPTH,
        DEPTH_WINDOW,
        BASE_ARGS["random_seed"],
        window_origin=0,
    )

    assert set(_retained(accepted)) == {"late", "early"}


# ---------------------------------------------------------------------------
# the worker count
#
# This script farms out a strand split that used to be one thread over the whole
# bam -- 15 minutes of a 72-minute step, measured over 48.1 M records with 27 of
# 28 cores idle -- plus two depth normalizations that are independent files. The
# division of the work may not change the answer: this bam is the splice-graph
# input for every cluster's quantification, so a different record set here is a
# different splice graph and a different number in every downstream table.

MULTI_CONTIGS = (("chr1", CONTIG_LEN), ("chr2", CONTIG_LEN), ("chrM", 16569))


def _multi_contig_header():
    return pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": name, "LN": length} for name, length in MULTI_CONTIGS],
        }
    )


def _multi_contig_corpus(header):
    """piles above the target depth on every contig, in both orientations.

    Both orientations because the split sends them to different files, and
    unequal per-contig counts because the fan-out hands the biggest contig to a
    worker first -- so the order the passes finish in is not the header order the
    output has to come out in.
    """

    reads = []
    for reference_id, (name, _length) in enumerate(MULTI_CONTIGS):
        depth = 40 + reference_id * 25
        for index in range(depth):
            for orientation, flag in (("f", 0), ("r", 16)):
                read = pysam.AlignedSegment(header)
                read.query_name = "{}_{}_{:04d}".format(name, orientation, index)
                read.query_sequence = "A" * 200
                read.query_qualities = pysam.qualitystring_to_array("I" * 200)
                read.flag = flag
                read.reference_id = reference_id
                read.reference_start = 350 + (index % 3) * 900
                read.mapping_quality = 60
                read.cigartuples = [(0, 200)]
                reads.append(read)
    return reads


@pytest.fixture(scope="module")
def worker_counts(tmp_path_factory):
    """the same multi-contig bam normalized serially and on four workers"""

    root = tmp_path_factory.mktemp("workers")
    header = _multi_contig_header()
    whole = root / "whole.bam"
    with pysam.AlignmentFile(str(whole), "wb", header=header) as writer:
        for read in sorted(
            _multi_contig_corpus(header),
            key=lambda r: (r.reference_id, r.reference_start),
        ):
            writer.write(read)
    pysam.index(str(whole))

    outputs = dict()
    for label, num_workers in (("serial", 1), ("concurrent", 4)):
        work = root / label
        work.mkdir()
        output = work / "norm.bam"
        _run_normalizer(whole, output, work, num_workers=num_workers)
        outputs[label] = output

    return outputs


def test_the_worker_count_does_not_change_which_reads_survive(worker_counts):
    """The acceptance bar: same records, same XW weights, either way.

    XW as well as the record set, because the weight is what the splice graph
    reads support off; two runs keeping the same reads under different weights
    would still build different graphs.
    """

    serial = _retained(worker_counts["serial"])
    concurrent = _retained(worker_counts["concurrent"])

    # the corpus has to be thinned, or the comparison is between two runs that
    # both kept everything
    assert 0 < len(serial) < 2 * sum(40 + i * 25 for i in range(len(MULTI_CONTIGS)))
    assert concurrent == serial


def test_the_concurrent_run_output_is_coordinate_sorted_and_region_fetchable(
    worker_counts,
):
    """What the merge of per-contig parts must not lose.

    Record-set equality above cannot see ordering; this walks the records, and
    then fetches the LAST contig in the header by region, which is the query a
    concatenation in completion order answers wrongly.
    """

    output = worker_counts["concurrent"]
    with pysam.AlignmentFile(str(output), "rb") as reader:
        positions = [
            (read.reference_id, read.reference_start)
            for read in reader.fetch(until_eof=True)
        ]
    assert positions == sorted(positions)

    assert Path(str(output) + ".bai").exists()
    late_contig = MULTI_CONTIGS[-1][0]
    with pysam.AlignmentFile(str(output), "rb") as reader:
        fetched = [read.query_name for read in reader.fetch(late_contig)]
    assert fetched
    assert all(name.startswith(late_contig) for name in fetched)


def test_changing_the_worker_count_does_not_invalidate_a_warm_run(worker_counts):
    """A cache may not miss over a choice that decides nothing about content.

    The worker count is not part of either checkpoint token, so re-running the
    serial arm's directory on four workers must do no work at all.
    """

    work_dir = worker_counts["serial"].parent
    before = {
        path.relative_to(work_dir): path.stat().st_mtime_ns
        for path in work_dir.rglob("*")
        if path.is_file()
    }
    assert before

    _run_normalizer(
        work_dir.parent / "whole.bam",
        worker_counts["serial"],
        work_dir,
        num_workers=4,
    )

    after = {
        path.relative_to(work_dir): path.stat().st_mtime_ns
        for path in work_dir.rglob("*")
        if path.is_file()
    }

    assert after == before


def test_help_reports_the_worker_count():
    help_text = subprocess.run(
        [sys.executable, str(NORMALIZER), "--help"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout

    assert "--num_workers" in help_text


# ---------------------------------------------------------------------------
# the per-contig fan-out of the depth normalization
#
# The two strand bams are two units, so --num_workers governed the strand split
# and the samtools thread counts but not the normalization itself: on a 28 GB
# 10x bam (179 M mapped-primary records) with --num_workers 8, the split took
# 8.3 min and the normalization 50.5 min of a 65.9 min task, its own log saying
# "normalizing 2 strand bam(s) concurrently" with 6 of 8 cores idle.  The unit is
# now (strand bam, populated reference), so the same flag governs this stage too.
#
# The division of the work may not change the answer.  This bam is the
# splice-graph input for every cluster's quantification, so a different record
# set OR a different XW weight here is a different splice graph and a different
# number in every downstream table -- and nothing downstream would notice, it
# would just be different.


def _units_logged(stderr):
    """(num_bams, num_units, num_workers) from the fan-out's own summary line."""

    for line in stderr.splitlines():
        match = re.search(
            r"normalizing (\d+) bam\(s\) as (\d+) unit\(s\) on (\d+) worker\(s\)", line
        )
        if match:
            return tuple(int(group) for group in match.groups())
    raise AssertionError("no fan-out summary line in:\n" + stderr)


def _identities(bam_path):
    """name, flag, reference, position and XW of every retained record, in order.

    Compared instead of the file's bytes because BGZF framing may legitimately
    differ: the fanned-out arm's output is a concatenation of per-contig parts,
    each compressed independently, so block boundaries fall elsewhere than in a
    single writer's stream while the records are identical.
    """

    with pysam.AlignmentFile(str(bam_path), "rb") as reader:
        return [
            (
                read.query_name,
                read.flag,
                read.reference_name,
                read.reference_start,
                round(float(read.get_tag("XW")), 6) if read.has_tag("XW") else None,
            )
            for read in reader.fetch(until_eof=True)
        ]


@pytest.fixture(scope="module")
def real_bam_arms(tmp_path_factory):
    """the real multi-contig bam normalized as 2 units and as 2 x N units.

    The serial arm is the whole-file pass this script has always run -- one
    ``sift_bam`` over each whole strand bam -- which is what ``--num_workers 1``
    still means.  The concurrent arm divides each strand bam per reference.
    """

    root = tmp_path_factory.mktemp("real_bam")
    arms = dict()
    for label, num_workers in (("whole_file", 1), ("per_contig", 8)):
        work = root / label
        work.mkdir()
        output = work / "norm.bam"
        completed = _run_normalizer(
            REAL_MULTI_CONTIG_BAM, output, work, num_workers=num_workers
        )
        arms[label] = {"output": output, "stderr": completed.stderr}
    return arms


@pytest.mark.skipif(
    not REAL_MULTI_CONTIG_BAM.exists(),
    reason="the real multi-contig fixture is not in this checkout",
)
def test_the_fan_out_retains_the_same_records_with_the_same_weights(real_bam_arms):
    """The acceptance bar, on a real bam: same records, same XW, either way.

    A single-contig bam cannot exercise this at all -- it has one unit either way
    -- which is why the fixture is the four-reference one.
    """

    whole_file = _identities(real_bam_arms["whole_file"]["output"])
    per_contig = _identities(real_bam_arms["per_contig"]["output"])

    # the corpus has to be thinned, or this compares two runs that kept
    # everything, and it has to span every reference, or the fan-out was not
    # exercised across them
    assert 0 < len(whole_file) < 11763
    assert len({identity[2] for identity in whole_file}) == 4

    assert set(per_contig) == set(whole_file)
    assert per_contig == whole_file


@pytest.mark.skipif(
    not REAL_MULTI_CONTIG_BAM.exists(),
    reason="the real multi-contig fixture is not in this checkout",
)
def test_the_fan_out_is_one_unit_per_populated_reference_per_strand(real_bam_arms):
    """What makes --num_workers govern this stage rather than the number 2.

    Four populated references and two strand bams is eight units, which is the
    whole point: the two-way pass could not use more than two workers however
    many were asked for.
    """

    assert _units_logged(real_bam_arms["per_contig"]["stderr"]) == (2, 8, 8)


@pytest.mark.skipif(
    not REAL_MULTI_CONTIG_BAM.exists(),
    reason="the real multi-contig fixture is not in this checkout",
)
def test_one_worker_is_still_the_unchanged_whole_file_pass(real_bam_arms):
    """The documented escape hatch, and it may not quietly become the fan-out.

    --num_workers 1 needs no index and reads each strand bam whole, which is what
    an input the fan-out cannot address falls back to as well.
    """

    stderr = real_bam_arms["whole_file"]["stderr"]
    assert _units_logged(stderr) == (2, 2, 1)
    assert stderr.count("(whole file)") == 2


@pytest.mark.skipif(
    not REAL_MULTI_CONTIG_BAM.exists(),
    reason="the real multi-contig fixture is not in this checkout",
)
def test_the_fan_outs_output_is_sorted_and_fetchable_on_a_late_contig(real_bam_arms):
    """Record-set equality cannot see ordering, and it is what a merge loses.

    Two assertions, because they fail on different mistakes.  Walking
    ``(reference_id, reference_start)`` catches a concatenation in completion
    order.  An indexed fetch of the LAST reference in the header catches an index
    built over a file whose references are out of order -- which samtools does
    without complaint, see the negative control below.
    """

    output = real_bam_arms["per_contig"]["output"]
    with pysam.AlignmentFile(str(output), "rb") as reader:
        positions = [
            (read.reference_id, read.reference_start)
            for read in reader.fetch(until_eof=True)
        ]
        references = list(reader.references)
    assert positions == sorted(positions)

    assert Path(str(output) + ".bai").exists()
    late_contig = references[-1]
    with pysam.AlignmentFile(str(output), "rb") as reader:
        fetched = [
            (read.query_name, read.reference_name)
            for read in reader.fetch(late_contig)
        ]
    assert fetched
    assert all(reference == late_contig for _name, reference in fetched)


def test_the_part_order_guard_is_what_catches_a_reversed_concatenation(tmp_path):
    """The negative control for the ordering guarantee, on the shared guard.

    ``verify_part_order`` is the SAME function the strand split uses, called here
    with the normalization's own part key -- two copies of an ordering assertion
    drift, and the drifting copy is the one that stops catching anything.  The
    second half of this test is why the guard has to exist at all: samtools
    indexes the reversed concatenation with exit 0 and no warning, so nothing
    downstream would say the file is unsorted.
    """

    header = _multi_contig_header()
    parts = list()
    for index, (name, _length) in enumerate(MULTI_CONTIGS):
        part = tmp_path / "part_{:05d}.bam".format(index)
        read = pysam.AlignedSegment(header)
        read.query_name = "{}_only".format(name)
        read.query_sequence = "A" * 50
        read.query_qualities = pysam.qualitystring_to_array("I" * 50)
        read.flag = 0
        read.reference_id = index
        read.reference_start = 100
        read.mapping_quality = 60
        read.cigartuples = [(0, 50)]
        with pysam.AlignmentFile(str(part), "wb", header=header) as writer:
            writer.write(read)
        parts.append(
            {"scope": name, "reference_id": index, "part": str(part)}
        )

    # header order passes
    nbs.verify_part_order(parts, "part", label="norm.bam")

    with pytest.raises(RuntimeError) as raised:
        nbs.verify_part_order(list(reversed(parts)), "part", label="norm.bam")
    assert "not in header order" in str(raised.value)
    assert "norm.bam" in str(raised.value)

    reversed_output = tmp_path / "reversed.bam"
    nbs.concatenate_bam_parts(
        str(reversed_output), [part["part"] for part in reversed(parts)]
    )
    with pysam.AlignmentFile(str(reversed_output), "rb") as reader:
        positions = [
            (read.reference_id, read.reference_start)
            for read in reader.fetch(until_eof=True)
        ]
    assert positions != sorted(positions), "the negative control must be unsorted"

    indexed = subprocess.run(
        ["samtools", "index", str(reversed_output)], capture_output=True, text=True
    )
    assert indexed.returncode == 0, (
        "if samtools started refusing this, the guard could be simplified"
    )


def test_a_part_holding_the_wrong_reference_is_refused_by_the_shared_guard(tmp_path):
    """The other half of the guard, reached through the normalization's part key."""

    header = _multi_contig_header()
    part = tmp_path / "part_00000.bam"
    read = pysam.AlignedSegment(header)
    read.query_name = "misplaced"
    read.query_sequence = "A" * 50
    read.query_qualities = pysam.qualitystring_to_array("I" * 50)
    read.flag = 0
    read.reference_id = 2
    read.reference_start = 100
    read.mapping_quality = 60
    read.cigartuples = [(0, 50)]
    with pysam.AlignmentFile(str(part), "wb", header=header) as writer:
        writer.write(read)

    with pytest.raises(RuntimeError) as raised:
        nbs.verify_part_order(
            [{"scope": "chr1", "reference_id": 0, "part": str(part)}],
            "part",
            label="norm.bam",
        )
    assert "rather than 0" in str(raised.value)
    # the normalize's call shape: an explicit label, because "part" would read as
    # "the part part". The split calls the same function with two positional
    # arguments and gets its strand as the label; util/test_separate_bam_by_strand.py
    # pins that default, so neither caller's spelling is left implicit.
    assert "the norm.bam part for scope" in str(raised.value)


# ---------------------------------------------------------------------------
# the window grid under the fan-out
#
# The default anchor is a function of the input -- the first aligned base seen on
# a contig -- which is why the CONTIG is the unit and not a fragment of one.  A
# caller supplying --window_origin is on an absolute grid instead, and is also
# the caller whose input is a single rebased contig, so the fan-out is not taken
# for it at all.  Both halves are asserted, because "buys nothing" and "changes
# nothing" are different claims.


def test_a_window_origin_caller_gets_the_whole_file_pass_and_todays_answer(tmp_path):
    """A rebased chunk bam is ONE contig, so it takes the unfanned pass.

    Same records and same weights as the same run before this change, which the
    module-scoped `chunked` fixture above already pins against a whole-contig
    control; what is asserted here is that the fan-out is not entered, so those
    equivalences are untouched however many workers are asked for.
    """

    header = _bam_header()
    src = tmp_path / "chunk.bam"
    _write_bam(src, _corpus(header), header)

    output = tmp_path / "chunk.norm.bam"
    completed = _run_normalizer(
        src, output, tmp_path, window_origin=0, single_strand=True, num_workers=8
    )

    num_bams, num_units, _workers = _units_logged(completed.stderr)
    assert (num_bams, num_units) == (1, 1)
    assert "(whole file)" in completed.stderr


def test_an_explicit_origin_survives_the_fan_out_unchanged(tmp_path):
    """Multi-contig plus an origin: the fan-out is entered and is still exact.

    ``grid_anchor`` is then ``-(window_origin % depth_window)`` -- one constant,
    independent of which records or contigs a unit holds -- so the boundaries a
    unit measures over are the boundaries the whole-file pass measures over.  The
    origin is deliberately NOT a multiple of the window, so it moves the
    boundaries off the coordinate grid and a unit silently falling back to the
    default anchor would be visible.
    """

    header = _multi_contig_header()
    src = tmp_path / "multi.bam"
    with pysam.AlignmentFile(str(src), "wb", header=header) as writer:
        for read in sorted(
            _multi_contig_corpus(header),
            key=lambda r: (r.reference_id, r.reference_start),
        ):
            writer.write(read)
    pysam.index(str(src))

    outputs = dict()
    for label, num_workers in (("whole_file", 1), ("per_contig", 4)):
        work = tmp_path / label
        work.mkdir()
        output = work / "norm.bam"
        completed = _run_normalizer(
            src, output, work, window_origin=37, num_workers=num_workers
        )
        outputs[label] = (output, completed.stderr)

    assert _units_logged(outputs["per_contig"][1])[1] == 6  # 3 contigs x 2 strands
    assert "offset 37 within a 100 bp window" in outputs["per_contig"][1]

    whole_file = _identities(outputs["whole_file"][0])
    per_contig = _identities(outputs["per_contig"][0])
    assert 0 < len(whole_file)
    assert per_contig == whole_file


def test_an_input_the_fan_out_cannot_address_falls_back_rather_than_failing(tmp_path):
    """A unit reaches its reference by name, and not every input can be reached.

    Refusing would regress a caller that works today, so the fan-out is declined
    and the whole-file pass runs.  Provoked with a header declaring a sort order
    other than coordinate, which is the case no amount of indexing works around:
    the parts' order is the input's order, and the concatenation is only
    coordinate-sorted for a coordinate-sorted input.

    Run with --input_is_single_strand because the strand split REFUSES such an
    input outright (separate_bam_by_strand.require_coordinate_sorted) and so the
    normalization would never be reached; that refusal is the split's, is correct,
    and is not what this asserts.
    """

    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "queryname"},
            "SQ": [{"SN": name, "LN": length} for name, length in MULTI_CONTIGS],
        }
    )
    src = tmp_path / "queryname.bam"
    # records still in coordinate order, so the pass itself succeeds; only the
    # header's claim stops the parts from being concatenable
    with pysam.AlignmentFile(str(src), "wb", header=header) as writer:
        for read in sorted(
            _multi_contig_corpus(header),
            key=lambda r: (r.reference_id, r.reference_start),
        ):
            writer.write(read)
    pysam.index(str(src))

    output = tmp_path / "norm.bam"
    completed = _run_normalizer(
        src, output, tmp_path, single_strand=True, num_workers=4
    )

    # three populated references, so a fan-out would be three units
    assert _units_logged(completed.stderr) == (1, 1, 1)
    assert "whole-file pass rather than per reference" in completed.stderr
    assert _identities(output)
