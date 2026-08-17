#!/usr/bin/env python3

"""The two arms have to consume the same records, or comparing them proves nothing.

A cut severs some reads, and both neighbouring chunks drop them.  That loss is
accepted -- it is what makes chunking possible -- but the whole-contig control
has no cut to lose them at.  Left in, the control quantifies records the chunked
arm never saw, and every difference between the arms is confounded by exactly
those reads.  Measured on a 3.4 Mb contig cut into seven segments per strand:
2,266 records against 2,261.

That is the whole reason cut selection emits the severed reads.  These tests
pin the subtraction and the two ways it can silently stop happening: the
baseline reading the unpruned bam, and the pruned bam being reused after the
cuts moved.
"""

import argparse
import importlib.util
import os
import sys

import pysam
import pytest

_HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "..", "pylib"))


def _load_driver():
    """The pipeline now lives in pylib/ChunkedRun.py; the script here is a front end."""
    import ChunkedRun

    return ChunkedRun


DRIVER = _load_driver()

SETTINGS = dict(
    max_intron_length=200000,
    normalize_max_cov_level=1000,
    depth_window=100,
    random_seed=42,
    HiFi=False,
    contig=None,
    cpu_budget=16,
    # The parity arms are quant-only; discovery has no baseline arm to compare
    # against, and run() refuses that combination.
    discovery=False,
)


def _args(tmp_path):
    settings = dict(SETTINGS)
    for name in ("bam", "genome_fa", "gtf"):
        path = tmp_path / name
        if not path.exists():
            path.write_text(name)
        settings[name] = str(path)
    return argparse.Namespace(**settings)


def _write_bam(path, names):
    header = pysam.AlignmentHeader.from_references(["chrT"], [10000])
    with pysam.AlignmentFile(str(path), "wb", header=header) as writer:
        for i, name in enumerate(names):
            aln = pysam.AlignedSegment(header)
            aln.query_name = name
            aln.flag = 0
            aln.reference_id = 0
            aln.reference_start = 100 * (i + 1)
            aln.mapping_quality = 60
            aln.cigartuples = [(0, 50)]
            aln.query_sequence = "A" * 50
            aln.query_qualities = pysam.qualitystring_to_array("I" * 50)
            writer.write(aln)
    pysam.index(str(path))


def _run_baseline(tmp_path, severed_names, record):
    """Drive the real run_baseline with no subprocess and no pruning side effect."""

    class RecordingCheckpoints(DRIVER.Checkpoints):
        def done(self, token):
            record.setdefault("tokens", []).append(token)
            return False

        def mark(self, token, note=""):
            pass

    def fake_run_step(name, cmd, *a, **k):
        record.setdefault("cmds", {})[name] = cmd
        return {"step": name}

    def fake_prune(source, names, dest):
        record["pruned"] = (source, sorted(names), dest)
        return len(names)

    DRIVER.run_step = fake_run_step
    DRIVER.write_bam_excluding = fake_prune

    outdir = tmp_path / "out"
    (outdir / "logs").mkdir(parents=True, exist_ok=True)

    with pytest.raises(DRIVER.PipelineError):
        DRIVER.run_baseline(
            _args(tmp_path),
            RecordingCheckpoints(str(tmp_path / "ckpt")),
            str(outdir),
            {},
            {"+": str(tmp_path / "plus.bam"), "-": str(tmp_path / "minus.bam")},
            1000,
            0.5,
            "stage1.up_aaaaaaaaaaaa",
            severed_names=severed_names,
        )
    return record


# ------------------------------------------------------------------ the name set


def test_severed_names_are_unioned_across_orientations(tmp_path):
    """One file per strand, and a read severed on either is severed for the run."""
    (tmp_path / "plus.dropped_reads.txt").write_text("read_a\nread_b\n")
    (tmp_path / "minus.dropped_reads.txt").write_text("read_c\n\n")

    assert DRIVER.severed_read_names(str(tmp_path)) == {"read_a", "read_b", "read_c"}


def test_no_cut_selection_means_no_names(tmp_path):
    """A baseline-only run has no chunked arm to be comparable with."""
    assert DRIVER.severed_read_names(str(tmp_path)) == set()


# -------------------------------------------------------------------- the prune


def test_pruning_removes_exactly_the_named_reads(tmp_path):
    source = tmp_path / "whole.bam"
    _write_bam(source, ["keep_1", "sever_1", "keep_2", "sever_2", "keep_3"])

    dest = tmp_path / "parity.bam"
    kept = DRIVER.write_bam_excluding(str(source), {"sever_1", "sever_2"}, str(dest))

    with pysam.AlignmentFile(str(dest), "rb") as reader:
        remaining = [aln.query_name for aln in reader.fetch(until_eof=True)]

    assert kept == 3
    assert remaining == ["keep_1", "keep_2", "keep_3"]


# ------------------------------------------------------------- what quant reads


def test_the_baseline_quantifies_the_pruned_bam(tmp_path):
    """Both the normalizer and LRAA must read the pruned input, not the merge.

    Normalization too, not only quantification: the chunked arm normalizes chunk
    bams that already exclude these reads, so a control normalizing the full set
    measures depth the other arm never had and thins against it.
    """
    record = _run_baseline(tmp_path, {"sever_1", "sever_2"}, {})

    assert record["pruned"][1] == ["sever_1", "sever_2"]
    parity_bam = record["pruned"][2]
    assert parity_bam.endswith("whole.parity.bam")

    norm_cmd = record["cmds"]["baseline_normalize"]
    assert norm_cmd[norm_cmd.index("--input_bam") + 1] == parity_bam

    quant_cmd = record["cmds"]["baseline_quant"]
    assert quant_cmd[quant_cmd.index("--bam") + 1] == parity_bam


def test_with_nothing_severed_the_baseline_reads_the_merge(tmp_path):
    """No cuts, no subtraction: the whole input is the right input."""
    record = _run_baseline(tmp_path, set(), {})

    assert "pruned" not in record
    norm_cmd = record["cmds"]["baseline_normalize"]
    assert norm_cmd[norm_cmd.index("--input_bam") + 1].endswith("whole.primary.bam")


def test_a_different_severed_set_invalidates_the_pruned_bam(tmp_path):
    """Cuts move, the excluded reads move, and the pruned bam may not be reused.

    Keyed on a digest of the names rather than their count, which would collide
    between two runs severing the same number of different reads.
    """
    one = _run_baseline(tmp_path, {"sever_1", "sever_2"}, {})["tokens"]
    two = _run_baseline(tmp_path, {"sever_1", "sever_3"}, {})["tokens"]

    parity_one = [t for t in one if t.startswith("baseline_parity.")]
    parity_two = [t for t in two if t.startswith("baseline_parity.")]

    assert parity_one and parity_two
    assert parity_one != parity_two
    # and everything downstream of it moves too
    assert one != two
