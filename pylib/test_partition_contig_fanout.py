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


def _run(tmp, bam, fasta, workers, out_root, bam_for_sg=None, extra=None):
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
        # named explicitly: the default is a RELATIVE "split_bams_for_sg", so
        # leaving it out writes into the caller's cwd -- the repo root under
        # pytest, and a read-only image directory under `docker run --user`
        "--bam-for-sg-out-dir", str(out_root / "sg_bams"),
        "--fasta-out-dir", str(out_root / "fa"),
        "--gtf-out-dir", str(out_root / "gtf"),
    ]
    if bam_for_sg is not None:
        cmd += ["--bam-for-sg", str(bam_for_sg), "--bam-for-sg-out-dir", str(out_root / "sg")]
    cmd += list(extra or [])
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
    "machine_cores,extra_threads,expected_workers",
    [(4, 0, 2), (8, 0, 6), (16, 0, 14), (30, 0, 28), (16, 4, 2), (4, 1, 1), (4, 4, 1)],
)
def test_the_workers_fill_the_machine_the_task_runs_on(
    machine_cores, extra_threads, expected_workers
):
    """The WDL derives the worker count from the cores of the machine it picked.

    The task owns its VM, so every core it is given is a core it can use. A worker
    needs extra_threads + 1 runnable threads (samtools' -@ is ADDITIONAL threads) and
    the two single-threaded FASTA/GTF jobs run alongside the bam pool, so the pool is
    (cores - 2) // (extra_threads + 1), never below 1. The script is handed the same
    cores as its reservation and holds itself inside them, so a worker count above
    this is capped rather than oversubscribed.
    """

    workers = max(1, (machine_cores - 2) // (extra_threads + 1))
    assert workers == expected_workers
    # what the WDL reserves for exactly this many workers never exceeds the machine
    assert workers * (extra_threads + 1) + 2 <= max(machine_cores, 3 + extra_threads)

    wdl = (REPO / "WDL" / "subwdls" / "Partition_data_by_chromosome.wdl").read_text()
    assert "(c3d_effective_cpu - 2) / (samtools_extra_threads + 1)" in wdl
    # clamped, because the script floors its pool at 1 regardless
    assert "if workers_that_fit_raw < 1 then 1 else workers_that_fit_raw" in wdl
    # unset fills the machine; a given value still wins
    assert "Int? partition_workers" in wdl
    assert "if defined(partition_workers)" in wdl
    assert "--num-workers ~{effective_partition_workers}" in wdl
    # helper threads give way first on a small machine, or the script refuses the
    # reservation (one worker plus the light jobs needs extra_threads + 3 cores)
    assert "extra_threads_that_fit" in wdl


def test_the_single_cell_workflow_sizes_each_partition_to_its_work():
    """Each partition is sized to what it splits, and the shared inputs are split once.

    The initial call splits the whole library plus the genome fasta and the annotation
    and is the wide one. A per-cluster call now splits only that cluster's own bam,
    because the fasta, the annotation and the shared splice-graph bam are split ONCE
    by the parent and handed in by contig name. The failure mode of getting any of
    this wrong is silent (it still runs, just repeating the work), so the wiring is
    pinned rather than inferred from a run.
    """

    sc = (REPO / "WDL" / "LRAA-singlecell.wdl").read_text()
    top = (REPO / "WDL" / "LRAA.wdl").read_text()

    # the single-cell workflow declares the three sizes and forwards each
    assert "Int initial_partition_cpu = 16" in sc
    assert "partition_cpu = initial_partition_cpu" in sc
    assert "Int cluster_partition_cpu = 4" in sc
    assert "cluster_partition_cpu = cluster_partition_cpu" in sc
    assert "Int shared_partition_cpu = 16" in sc
    assert "shared_partition_cpu = shared_partition_cpu" in sc
    # worker counts are overrides now, unset by default, not a tuned default
    assert "Int? initial_partition_workers" in sc
    assert "Int? cluster_partition_workers" in sc

    # LRAA.wdl accepts the cores and the worker override, hands both to the partition,
    # and takes the pre-split shared inputs by name
    assert "Int partition_cpu = 4" in top
    assert "Int? partition_workers" in top
    assert "cpu = partition_cpu" in top
    for name in ("fastas", "gtfs", "sg_bams"):
        assert "Map[String, File]? internal_presplit_{}".format(name) in top, name

    # a shared input that was pre-split is not handed to the per-cluster partition at
    # all: it would be localized and re-split otherwise, which is the cost removed
    assert "if (!defined(internal_presplit_fastas))" in top
    assert "if (!defined(internal_presplit_gtfs) && defined(annot_gtf))" in top
    assert "if (!defined(internal_presplit_sg_bams) && defined(internal_bam_for_sg))" in top
    # and the shard reads its contig's file by NAME from the pre-split map
    assert "select_first([internal_presplit_fastas])[contig_name]" in top

    for name in ("LRAA-cell_cluster_guided.wdl", "LRAA_quant_by_cluster.wdl"):
        text = (REPO / "WDL" / name).read_text()
        assert "Int cluster_partition_cpu = 4" in text, name
        assert "Int? cluster_partition_workers" in text, name
        assert "Int shared_partition_cpu = 16" in text, name
        assert "partition_cpu = cluster_partition_cpu" in text, name
        assert "partition_workers = cluster_partition_workers" in text, name
        assert "cpu = shared_partition_cpu" in text, name
        assert "internal_presplit_fastas =" in text, name
        assert "internal_presplit_gtfs =" in text, name
    # only the final quant has a shared splice-graph bam to split
    qc = (REPO / "WDL" / "LRAA_quant_by_cluster.wdl").read_text()
    assert "internal_presplit_sg_bams = split_shared_inputs.chromosomeBAMsForSGByName" in qc


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
            "--bam-for-sg-out-dir", str(out / "sg_bams"),
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
    # the cores of the machine the tier picked, not the number that was asked for
    assert "--reserved-cpu ~{c3d_effective_cpu}" in wdl


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
            "--bam-for-sg-out-dir", str(out / "sg_bams"),
            "--fasta-out-dir", str(out / "fa"),
            "--gtf-out-dir", str(out / "gtf"),
        ],
        capture_output=True,
        text=True,
    )
    assert res.returncode != 0
    assert "cannot run this task" in res.stderr
    assert "the floor is 7" in res.stderr


def test_the_default_configuration_fits_its_own_machine():
    """The WDL defaults must fit the machine they select, or the script refuses to start.

    One worker plus the two FASTA/GTF jobs needs extra_threads + 3 cores. The default is
    cpu 4 (a C3D tier) with one core per extraction, so 3 <= 4. If either default is
    raised without the other the task would be refused its own reservation.
    """

    wdl = (REPO / "WDL" / "subwdls" / "Partition_data_by_chromosome.wdl").read_text()
    assert "Int cpu = 4" in wdl
    assert "Int samtools_threads = 1" in wdl
    extra_threads, cpu = 0, 4
    assert extra_threads + 3 <= cpu


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
            # so the output is not root-owned: pytest cannot clean tmp_path
            # otherwise, and every run leaks a directory into /tmp
            "--user", "{}:{}".format(os.getuid(), os.getgid()),
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
            "--bam-for-sg-out-dir", "/w/sg_bams",
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
            "--user", "{}:{}".format(os.getuid(), os.getgid()),
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
            "--bam-for-sg-out-dir", "/w/sg_bams",
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


@pytest.mark.parametrize("level", [1, 4, 9])
def test_the_compression_level_reaches_samtools_and_the_output_is_still_readable(
    inputs, tmp_path, level
):
    """--bam-compression-level must be applied, not accepted and ignored.

    Per-contig BAMs are intermediates a shard reads once, so the script lets the caller
    lower the BGZF level. The planned units carry it, and the emitted BAMs must be
    readable with the same records at any level.
    """

    sys.path.insert(0, str(REPO / "util"))
    from partition_data_by_chromosome import _plan_bam_partition

    _d, bam, _sg, _fasta = inputs
    work = _plan_bam_partition(
        str(bam), [c[0] for c in CONTIGS], str(tmp_path / "planned"), "BAM", 1, level
    )
    assert work, "fixture should plan work"
    assert {unit[6] for unit in work} == {level}

    out = tmp_path / "out"
    _run(tmp_path, bam, _fasta, 1, out, extra=["--bam-compression-level", str(level)])
    assert _counts(out / "bams") == {"{}.bam".format(c[0]): c[2] for c in CONTIGS}


def test_the_default_leaves_the_compression_level_to_samtools():
    """Unset means unset: no --output-fmt-option is passed, so existing callers are unchanged."""

    sys.path.insert(0, str(REPO / "util"))
    import inspect

    from partition_data_by_chromosome import _extract_one_contig

    assert inspect.signature(_extract_one_contig).parameters["level"].default is None
