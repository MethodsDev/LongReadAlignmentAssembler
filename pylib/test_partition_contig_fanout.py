#!/usr/bin/env python3

"""Fanning the partition's contig loop out must not change what it produces.

The script used to fan out over JOB TYPES -- BAM, BAM_FOR_SG, FASTA, GTF -- with a
pool fixed at four, while the contig loop inside ran serially. Three of those four
jobs are trivial or absent, so the task was one serial pass over the largest bam in
the pipeline: REPORTED on a 188 GB library, 27+ minutes with one core busy, ahead of
all shard work on an idle box.

Contigs now fan out. These tests hold the two properties that make that safe -- the
output does not depend on the worker count, and the emitted BAMs are readable and
correctly headed at any worker count -- plus the budget arithmetic, because the
caller's cpu reservation is derived from it and an off-by-one there oversubscribes a
task that has no cgroup to cap it.
"""

import os
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "util" / "partition_data_by_chromosome.py"

# Uneven on purpose: the pool is meant to start the biggest first, and equal-sized
# contigs would hide an ordering mistake behind a symmetric workload.
CONTIGS = [("chrA", 4000, 40), ("chrB", 3000, 25), ("chrC", 2000, 10), ("chrD", 1000, 3)]


def _make_bam(path, contigs):
    """A sorted, indexed BAM with `reads` alignments spread along each contig."""

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": name, "LN": length} for name, length, _ in contigs],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for tid, (name, length, reads) in enumerate(contigs):
            step = max(1, (length - 120) // max(1, reads))
            for i in range(reads):
                a = pysam.AlignedSegment(fh.header)
                a.query_name = "{}_{}".format(name, i)
                a.query_sequence = "A" * 100
                a.flag = 0
                a.reference_id = tid
                a.reference_start = 10 + i * step
                a.mapping_quality = 60
                a.cigartuples = [(0, 100)]
                a.query_qualities = pysam.qualitystring_to_array("I" * 100)
                fh.write(a)
    pysam.index(str(path))
    return path


def _make_fasta(path, contigs):
    with open(path, "wt") as fh:
        for name, length, _ in contigs:
            fh.write(">{}\n".format(name))
            seq = "ACGT" * (length // 4 + 1)
            seq = seq[:length]
            for i in range(0, length, 60):
                fh.write(seq[i : i + 60] + "\n")
    return path


def _run(tmp, bam, fasta, workers, out_root, bam_for_sg=None):
    out_root.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(SCRIPT),
        "--input-bam", str(bam),
        "--genome-fasta", str(fasta),
        "--chromosomes", *[c[0] for c in CONTIGS],
        "--samtools-threads", "1",
        "--num-workers", str(workers),
        "--bam-out-dir", str(out_root / "bams"),
        "--fasta-out-dir", str(out_root / "fa"),
        "--gtf-out-dir", str(out_root / "gtf"),
    ]
    if bam_for_sg is not None:
        cmd += ["--bam-for-sg", str(bam_for_sg), "--bam-for-sg-out-dir", str(out_root / "sg")]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, res.stderr[-3000:]
    return res


def _counts(bam_dir):
    """Per-contig record counts, read back through pysam so an unreadable or
    mis-headed output fails here rather than silently comparing byte sizes."""

    out = {}
    for path in sorted(Path(bam_dir).glob("*.bam")):
        with pysam.AlignmentFile(str(path), "rb") as fh:
            out[path.name] = sum(1 for _ in fh.fetch(until_eof=True))
    return out


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    d = tmp_path_factory.mktemp("partition_inputs")
    bam = _make_bam(d / "reads.bam", CONTIGS)
    sg = _make_bam(d / "sg.bam", CONTIGS)
    fasta = _make_fasta(d / "ref.fa", CONTIGS)
    return d, bam, sg, fasta


@pytest.mark.parametrize("workers", [2, 3, 5])
def test_the_output_does_not_depend_on_the_worker_count(inputs, workers, tmp_path):
    """The property the fan-out has to preserve, stated as equality against serial.

    Worker counts chosen around the work: 2 divides it, 3 does not, and 5 exceeds the
    4 contigs so the pool is wider than there is work for -- the case where a
    min()/floor mistake would either strand a contig or oversubscribe.
    """

    _d, bam, _sg, fasta = inputs
    _run(tmp_path, bam, fasta, 1, tmp_path / "serial")
    _run(tmp_path, bam, fasta, workers, tmp_path / "fanned")

    serial = _counts(tmp_path / "serial" / "bams")
    fanned = _counts(tmp_path / "fanned" / "bams")

    assert serial == fanned, (workers, serial, fanned)
    # and the work was actually there to divide
    assert sum(serial.values()) == sum(c[2] for c in CONTIGS)


def test_both_bam_kinds_are_partitioned_under_one_budget(inputs, tmp_path):
    """bam_for_sg doubles the work, and it must not double the concurrency.

    The budget is shared, so this asserts the OUTPUT of that sharing: every contig of
    both kinds is present and complete. A per-kind pool would still pass this, which
    is why the arithmetic itself is asserted separately below -- but a flattening bug
    that dropped or duplicated a kind's work would fail here.
    """

    _d, bam, sg, fasta = inputs
    _run(tmp_path, bam, fasta, 3, tmp_path / "both", bam_for_sg=sg)

    primary = _counts(tmp_path / "both" / "bams")
    graph = _counts(tmp_path / "both" / "sg")
    expected = {"{}.bam".format(c[0]): c[2] for c in CONTIGS}

    assert primary == expected
    assert graph == expected


def test_emitted_bams_are_readable_and_carry_their_contig(inputs, tmp_path):
    """Index validity and header correctness, not just record counts.

    A fanned write that raced would most likely surface as a truncated BGZF block or a
    header from the wrong contig, neither of which a count comparison alone catches --
    fetch() on a freshly built index is what exercises both.
    """

    _d, bam, _sg, fasta = inputs
    _run(tmp_path, bam, fasta, 3, tmp_path / "idx")

    for name, _length, reads in CONTIGS:
        path = tmp_path / "idx" / "bams" / "{}.bam".format(name)
        pysam.index(str(path))
        with pysam.AlignmentFile(str(path), "rb") as fh:
            assert name in fh.references, (name, fh.references)
            assert sum(1 for _ in fh.fetch(name)) == reads


@pytest.mark.parametrize(
    "workers,extra_threads,expected_cpu",
    [(1, 4, 7), (2, 4, 12), (4, 4, 22), (1, 0, 3)],
)
def test_the_reservation_covers_what_the_script_can_run(workers, extra_threads, expected_cpu):
    """The formula the WDL reserves against, checked against its own definition.

    samtools' -@ is ADDITIONAL threads, so a worker needs extra_threads + 1 runnable
    threads, and the two single-threaded FASTA/GTF jobs run alongside the bam pool.
    Reserving less oversubscribes a task that has no cgroup to cap it; reserving more
    repeats the regression that cut this task's cpu to 5, because cpu nobody uses
    cannot be placed until the box drains.
    """

    assert workers * (extra_threads + 1) + 2 == expected_cpu

    wdl = (REPO / "WDL" / "subwdls" / "Partition_data_by_chromosome.wdl").read_text()
    assert "effective_partition_workers * (samtools_extra_threads + 1) + 2" in wdl
    # clamped, because the reservation is derived from it and the script floors at 1
    assert "if partition_workers < 1 then 1 else partition_workers" in wdl
    assert "--num-workers ~{effective_partition_workers}" in wdl



def test_the_single_cell_workflow_widens_only_its_initial_partition():
    """The fan-out has to REACH the workflow that motivated it, and only there.

    A knob on LRAA.wdl alone ships nothing: the reported case is a single-cell run,
    and LRAA-singlecell.wdl calls LRAA_wf three ways -- once for the initial
    full-library phase, then once per cluster from LRAA-cell_cluster_guided.wdl and
    LRAA_quant_by_cluster.wdl. Only the first should widen. The per-cluster calls run
    14 to 32 times for seconds each, and a 22-core reservation there cannot be placed
    until the box drains, which is the regression that cut this task's cpu to 5.

    Asserted on the wiring rather than on a run, because reaching this in a real
    invocation costs a whole-library partition -- but the failure mode is silent
    (default 1 everywhere, no error, no slower-than-before), so it needs pinning.
    """

    sc = (REPO / "WDL" / "LRAA-singlecell.wdl").read_text()
    top = (REPO / "WDL" / "LRAA.wdl").read_text()

    # the single-cell workflow declares it and forwards it to the initial call
    assert "Int initial_partition_workers = 4" in sc
    assert "partition_workers = initial_partition_workers" in sc

    # LRAA.wdl accepts it and hands it to the partition subworkflow
    assert "Int partition_workers = 1" in top
    assert "partition_workers = partition_workers" in top

    # and the per-cluster paths do NOT pass it, so they keep the default of 1
    for name in ("LRAA-cell_cluster_guided.wdl", "LRAA_quant_by_cluster.wdl"):
        text = (REPO / "WDL" / name).read_text()
        assert "partition_workers" not in text, name