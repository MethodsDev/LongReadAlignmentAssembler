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

# the tag the repo's own build script publishes, so this follows a rebuild rather than
# pinning a revision that goes stale
DOCKER_IMAGE = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:testing"

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

    # the single-cell workflow declares it and forwards it to the initial call. 2, not
    # the measured knee of 4: cpu is workers * (-@ + 1) + 2, so 2 asks 12 and fits a
    # 16-core machine while 4 asks 22 and forces a 32-core instance to be billed for
    # the whole task, buying the 18 s between 1:13.8 and 0:55.8. Pinned because either
    # number is defensible and the reason for choosing between them is not recoverable
    # from the value.
    assert "Int initial_partition_workers = 2" in sc
    assert "partition_workers = initial_partition_workers" in sc

    # LRAA.wdl accepts it and hands it to the partition subworkflow
    assert "Int partition_workers = 1" in top
    assert "partition_workers = partition_workers" in top

    # the per-cluster paths take a SEPARATE knob, because the cost differs: the initial
    # partition is one task widened once, while this one is billed 14 to 32 times over,
    # and locally miniwdl puts all of them on ONE host where the multiplied reservation
    # cannot be placed until the box drains.
    assert "Int cluster_partition_workers = 1" in sc
    assert "cluster_partition_workers = cluster_partition_workers" in sc
    for name in ("LRAA-cell_cluster_guided.wdl", "LRAA_quant_by_cluster.wdl"):
        text = (REPO / "WDL" / name).read_text()
        assert "Int cluster_partition_workers = 1" in text, name
        assert "partition_workers = cluster_partition_workers" in text, name

def test_the_pool_is_held_inside_the_reservation(inputs, tmp_path):
    """A cpu declaration is a promise to the scheduler; the pool must fit inside it.

    Nothing necessarily ENFORCES it. miniwdl adjusts cpu as a scheduling share and
    reports "cpu adjusted to host limit", but the task still sees every core -- so an
    affinity check does not bind there and the argv it already built still says
    --num-workers 4. Passing the reservation is what makes the cap independent of the
    backend. Reserved 7 -- the WDL's own one-worker default -- affords (7 - 2) // 5 = 1,
    so it forces one worker on any host, including one with cores to spare.

    Deliberately not parametrized over a generous reservation: on a box whose visible
    cores afford less than the reservation does, affinity binds first and the expected
    count would be a property of the test machine rather than of the code.
    """

    import subprocess as sp

    reserved, requested = 7, 8
    _d, bam, _sg, fasta = inputs
    out = tmp_path / "reserved{}".format(reserved)
    out.mkdir()
    res = sp.run(
        [
            sys.executable, str(SCRIPT),
            "--input-bam", str(bam),
            "--genome-fasta", str(fasta),
            "--chromosomes", *[c[0] for c in CONTIGS],
            "--samtools-threads", "4",
            "--num-workers", str(requested),
            "--reserved-cpu", str(reserved),
            "--bam-out-dir", str(out / "bams"),
            "--fasta-out-dir", str(out / "fa"),
            "--gtf-out-dir", str(out / "gtf"),
        ],
        capture_output=True,
        text=True,
    )
    assert res.returncode == 0, res.stderr[-2000:]
    # the reservation is the binding cap, and it says which one bound
    assert "reservation(7 core)" in res.stderr, res.stderr[-2000:]
    assert "1 at a time" in res.stderr, res.stderr[-2000:]
    # capped, not failed: the output is still complete
    assert _counts(out / "bams") == {"{}.bam".format(c[0]): c[2] for c in CONTIGS}


def test_the_wdl_passes_its_reservation_to_the_script():
    """The knob and the number it must respect travel together, or the cap is inert."""

    wdl = (REPO / "WDL" / "subwdls" / "Partition_data_by_chromosome.wdl").read_text()
    assert "--num-workers ~{effective_partition_workers}" in wdl
    assert "--reserved-cpu ~{partition_cpu}" in wdl


@pytest.mark.parametrize("reserved", [1, 3, 6])
def test_a_reservation_below_the_floor_is_refused(inputs, tmp_path, reserved):
    """Capping cannot honour a reservation smaller than one worker plus the light jobs.

    max(1, ...) would floor the pool at one worker and still exceed the reservation --
    the same oversubscription, just quieter. The caller has two real fixes (raise the
    reservation, lower --samtools-threads) and the script cannot pick between them, so
    it refuses and names both.
    """

    import subprocess as sp

    _d, bam, _sg, fasta = inputs
    out = tmp_path / "floor{}".format(reserved)
    out.mkdir()
    res = sp.run(
        [
            sys.executable, str(SCRIPT),
            "--input-bam", str(bam),
            "--genome-fasta", str(fasta),
            "--chromosomes", *[c[0] for c in CONTIGS],
            "--samtools-threads", "4",   # floor is 4 + 3 = 7
            "--num-workers", "4",
            "--reserved-cpu", str(reserved),
            "--bam-out-dir", str(out / "bams"),
            "--fasta-out-dir", str(out / "fa"),
            "--gtf-out-dir", str(out / "gtf"),
        ],
        capture_output=True,
        text=True,
    )
    assert res.returncode != 0
    assert "cannot run this task" in res.stderr
    assert "the floor is 7" in res.stderr


def test_the_default_reservation_is_exactly_the_floor():
    """cpu at one worker is 7, which is the smallest this script can run inside.

    Not a coincidence worth losing: if the WDL default and the script's floor ever
    drift apart, the default configuration either refuses to start or overruns.
    """

    threads_each = 4  # samtools_extra_threads at the WDL default of samtools_threads 5
    assert 1 * (threads_each + 1) + 2 == threads_each + 3 == 7


@pytest.mark.parametrize(
    "layout,expected",
    [
        ({"cpu.max": "800000 100000"}, 8),          # v2, docker --cpus=8
        ({"cpu.max": "250000 100000"}, 2),          # v2, fractional 2.5 rounds DOWN
        ({"cpu.max": "50000 100000"}, 1),           # v2, half a core still runs one
        ({"cpu.max": "max 100000"}, None),          # v2, unlimited
        ({"cpu/cpu.cfs_quota_us": "600000",
          "cpu/cpu.cfs_period_us": "100000"}, 6),   # v1
        ({"cpu/cpu.cfs_quota_us": "-1",
          "cpu/cpu.cfs_period_us": "100000"}, None),  # v1, unlimited
        ({}, None),                                 # no cgroup at all
    ],
)
def test_the_granted_cpu_is_read_from_the_cgroup(tmp_path, layout, expected):
    """The grant has to come from what the runtime ENFORCES, not what was requested.

    A WDL cpu declaration is a request. miniwdl applies it as docker --limit-cpu and
    reports "cpu adjusted to host limit" when it trims one, so the task can be granted
    fewer cores than its argv was built for. Measured inside `docker run --cpus=8`:
    sched_getaffinity still reports every host core (16 here), while the cgroup quota
    reads 8 -- so this is the only signal that binds in the case that motivated it.
    """

    sys.path.insert(0, str(REPO / "util"))
    from partition_data_by_chromosome import _cgroup_cpu_quota

    for name, text in layout.items():
        f = tmp_path / name
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(text)

    assert _cgroup_cpu_quota(str(tmp_path)) == expected


@pytest.mark.skipif(
    not __import__("shutil").which("docker"), reason="needs docker to constrain a cgroup"
)
def test_the_pool_is_capped_by_a_real_container_grant(inputs, tmp_path):
    """End-to-end through the mechanism miniwdl itself uses to apply cpu.

    This is the case a reservation argument cannot catch: the caller asked for enough,
    the backend granted less, and nothing in the argv changed. Affinity does not bind
    inside `--cpus`; the quota does.
    """

    import shutil
    import subprocess as sp

    _d, bam, _sg, fasta = inputs
    out = tmp_path / "granted"
    out.mkdir()
    res = sp.run(
        [
            shutil.which("docker"), "run", "--rm", "--cpus=8",
            "-v", "{}:/u:ro".format(REPO / "util"),
            "-v", "{}:/d:ro".format(bam.parent),
            "-v", "{}:/w".format(out),
            DOCKER_IMAGE,
            "python3", "/u/partition_data_by_chromosome.py",
            "--input-bam", "/d/{}".format(bam.name),
            "--chromosomes", *[c[0] for c in CONTIGS],
            "--samtools-threads", "4",
            "--num-workers", "8",
            "--reserved-cpu", "42",   # the caller asked for plenty; the grant is 8
            "--bam-out-dir", "/w/bams",
            "--fasta-out-dir", "/w/fa",
            "--gtf-out-dir", "/w/gtf",
        ],
        capture_output=True,
        text=True,
    )
    if res.returncode != 0 and "Cannot connect to the Docker daemon" in res.stderr:
        pytest.skip("docker daemon unavailable")
    assert res.returncode == 0, res.stderr[-2000:]

    # the GRANT bound it, and the log names which cap did
    assert "cgroup quota(8 core)" in res.stderr, res.stderr[-2000:]
    assert "1 at a time" in res.stderr, res.stderr[-2000:]


def test_the_planned_units_carry_the_thread_count_they_were_given(inputs):
    """The unit's thread field IS the -@ argv, so it must be the adapted value.

    `_extract_one_contig` passes unit[3] straight to `pysam.view("-@", str(threads))`.
    The units are planned before the pool is sized, so a thread count resolved after
    planning would be logged and never applied -- the summary would claim -@ 0 while
    four samtools ran at -@ 4.
    """

    sys.path.insert(0, str(REPO / "util"))
    from partition_data_by_chromosome import _plan_bam_partition

    _d, bam, _sg, _fasta = inputs
    work = _plan_bam_partition(
        str(bam), [c[0] for c in CONTIGS], str(_d / "planned"), "BAM", 3
    )
    assert work, "fixture should plan work"
    assert {unit[3] for unit in work} == {3}


@pytest.mark.skipif(
    not __import__("shutil").which("docker"), reason="needs docker to constrain a cgroup"
)
def test_a_grant_below_the_floor_lowers_the_threads_that_actually_run(inputs, tmp_path):
    """Under an enforced 3-core quota the pool cap alone is not enough.

    One worker at -@ 4 is 5 runnable threads and the two light jobs are 2 more, so
    capping the pool at 1 still runs 7 processes' worth of threads inside a 3-core
    grant. -@ is the only lever with no correctness meaning, so it gives way -- and
    this asserts the PER-CONTIG line, which sits beside the `-@` argv, not just the
    summary, so a value that was logged but not applied would fail here.

    Adapted rather than refused: a grant is a fact about the environment, and aborting
    would fail a run that can still finish.
    """

    import shutil
    import subprocess as sp

    _d, bam, _sg, fasta = inputs
    out = tmp_path / "tight"
    out.mkdir()
    res = sp.run(
        [
            shutil.which("docker"), "run", "--rm", "--cpus=3",
            "-v", "{}:/u:ro".format(REPO / "util"),
            "-v", "{}:/d:ro".format(bam.parent),
            "-v", "{}:/w".format(out),
            DOCKER_IMAGE,
            "python3", "-u", "/u/partition_data_by_chromosome.py",
            "--input-bam", "/d/{}".format(bam.name),
            "--chromosomes", *[c[0] for c in CONTIGS],
            "--samtools-threads", "4",
            "--num-workers", "8",
            "--bam-out-dir", "/w/bams",
            "--fasta-out-dir", "/w/fa",
            "--gtf-out-dir", "/w/gtf",
        ],
        capture_output=True,
        text=True,
    )
    if res.returncode != 0 and "Cannot connect to the Docker daemon" in res.stderr:
        pytest.skip("docker daemon unavailable")
    assert res.returncode == 0, res.stderr[-2000:]

    assert "--samtools-threads reduced from 4 to 0" in res.stderr, res.stderr[-3000:]
    # the line beside the -@ argv: what samtools was ACTUALLY told
    assert "using 0 additional thread(s)" in res.stderr, res.stderr[-3000:]
    # 1 worker + 2 light jobs == 3 processes, exactly the grant
    assert "1 at a time" in res.stderr, res.stderr[-3000:]
