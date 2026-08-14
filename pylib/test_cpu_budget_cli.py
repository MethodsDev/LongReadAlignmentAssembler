#!/usr/bin/env python3

"""The budget as an LRAA run actually reports it.

The arithmetic is asserted directly in test_cpu_budget.py. What is left, and what needs
a real invocation, is that the CLI resolves ONE budget, divides it across the flat queue
it actually built, prints the division on one line, and leaves component forking off
unless asked. Every subprocess here carries a timeout.
"""

import os
import re
import subprocess
import sys

import pysam
import pytest

import CpuBudget


LRAA = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "LRAA")

# Two tiny contigs -> four (contig, strand) units, which is enough to see the flat queue
# divide a budget and to see longest-first order two contigs of unequal read count.
CONTIGS = (("ctgA", 6000), ("ctgB", 6000))
RUN_TIMEOUT_S = 300


def _write_inputs(tmp_path):
    """A genome, an indexed bam and a gtf small enough to quantify in about a second.

    ctgA carries three times ctgB's reads, so the cost proxy has something to order by.
    """
    fasta = str(tmp_path / "genome.fa")
    with open(fasta, "wt") as ofh:
        for contig, length in CONTIGS:
            sequence = ("ACGTTGCA" * (length // 8 + 1))[:length]
            print(">{}".format(contig), file=ofh)
            for i in range(0, length, 60):
                print(sequence[i : i + 60], file=ofh)
    pysam.faidx(fasta)

    gtf = str(tmp_path / "annot.gtf")
    with open(gtf, "wt") as ofh:
        for contig, _length in CONTIGS:
            attrs = 'gene_id "g.{c}"; transcript_id "t.{c}";'.format(c=contig)
            for feature in ("gene", "transcript", "exon"):
                print(
                    "\t".join(
                        (contig, "test", feature, "1001", "2000", ".", "+", ".", attrs)
                    ),
                    file=ofh,
                )

    bam = str(tmp_path / "reads.bam")
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": c, "LN": n} for c, n in CONTIGS],
    }
    read_counts = {"ctgA": 30, "ctgB": 10}
    with pysam.AlignmentFile(bam, "wb", header=header) as ofh:
        for tid, (contig, _length) in enumerate(CONTIGS):
            for i in range(read_counts[contig]):
                aln = pysam.AlignedSegment()
                aln.query_name = "{}_read{}".format(contig, i)
                aln.reference_id = tid
                aln.reference_start = 1000
                aln.mapping_quality = 60
                aln.cigarstring = "1000M"
                aln.query_sequence = "A" * 1000
                aln.query_qualities = pysam.qualitystring_to_array("I" * 1000)
                aln.flag = 0
                ofh.write(aln)
    pysam.index(bam)
    return fasta, bam, gtf


def _run(tmp_path, *extra, budget=None, quant_only=True):
    fasta, bam, gtf = _write_inputs(tmp_path)
    cmd = [
        sys.executable,
        LRAA,
        "--genome", fasta,
        "--bam", bam,
        "--output_prefix", "budget_test",
    ]
    if quant_only:
        cmd += ["--gtf", gtf, "--quant_only"]
    if budget is not None:
        cmd += ["--cpu_budget", str(budget)]
    cmd += list(extra)
    proc = subprocess.run(
        cmd,
        cwd=str(tmp_path),
        capture_output=True,
        text=True,
        timeout=RUN_TIMEOUT_S,
    )
    assert proc.returncode == 0, proc.stderr[-4000:]
    return proc.stderr


def _allocation_line(log):
    lines = [
        ln.strip() for ln in log.splitlines()
        if ln.strip().startswith(CpuBudget.ALLOCATION_LOG_PREFIX)
    ]
    assert len(lines) == 1, "expected exactly one allocation line, got {}".format(lines)
    return lines[0]


def test_startup_log_states_the_allocation_on_one_line(tmp_path):
    """Unexpected wall time must be diagnosable from the log alone.

    Two contigs, both strands, budget 4 is four units on four cores: four workers, one
    native-tool thread each, nothing left over.
    """
    log = _run(tmp_path, budget=4)
    assert _allocation_line(log) == (
        "CPU allocation: budget=4 cores; phase=main; units=4; unit_workers=4; "
        "tool_threads_per_worker=1; component_worker_ceiling=4; unallocated_cores=0"
    )


def test_an_explicit_budget_is_what_gets_divided(tmp_path):
    """Budget 2 over four units: two workers, and the other two units queue behind them.

    This is the whole point of the change -- the flag names a TOTAL that the levels
    divide, not a per-level multiplier.
    """
    line = _allocation_line(_run(tmp_path, budget=2))
    assert "budget=2 cores" in line
    assert "units=4" in line
    assert "unit_workers=2" in line
    assert "tool_threads_per_worker=1" in line


def test_a_budget_above_the_unit_count_lends_the_rest_to_native_tool_steps(tmp_path):
    """Budget 8 over four units: four workers with two tool threads each, and the log
    says outright that the spare cores are not doing assembly work."""
    log = _run(tmp_path, budget=8)
    line = _allocation_line(log)
    assert "unit_workers=4" in line
    assert "tool_threads_per_worker=2" in line
    assert "unallocated_cores=0" in line
    note = [
        ln for ln in log.splitlines()
        if CpuBudget.SHORTFALL_LOG_PREFIX in ln
    ]
    assert note, "a budget the units cannot claim must be reported"
    assert "reachable only by native tool steps" in note[0]


def test_the_default_budget_comes_from_available_cpus(tmp_path):
    """No --cpu_budget means every core this process may run on, so a cpuset is honoured
    without the user restating it."""
    import Util_funcs

    log = _run(tmp_path)
    assert "budget={} cores".format(Util_funcs.available_cpus()) in _allocation_line(log)
    assert "--cpu_budget not set" in log


def test_no_component_reaches_the_spawn_hint_so_nothing_forks(tmp_path):
    """On real ONT chr20 no component reaches the hint of 150 -- the largest measured was
    95 -- so the default outcome is in-process assembly, and the log must say which reason
    applied rather than leaving it to be inferred.

    Discovery mode, because component assembly is the only thing that forks;
    --no_parallelize_contigs so the worker's own log reaches stderr rather than a per-unit
    error file.
    """
    log = _run(tmp_path, "--no_parallelize_contigs", budget=8, quant_only=False)
    assert "No component reaches the spawn hint" in log
    assert "-Component assembly granted" not in log


def test_the_ceiling_can_be_pinned_off_and_then_nothing_forks(tmp_path):
    """--component_workers 0 is the explicit opt-out, and it must survive the whole run."""
    log = _run(
        tmp_path, "--component_workers", "0", "--no_parallelize_contigs",
        budget=8, quant_only=False,
    )
    assert "component_worker_ceiling=0" in _allocation_line(log)
    assert "-Component assembly granted" not in log


def test_a_ceiling_is_reported_as_reserving_nothing(tmp_path):
    """The startup log must not let a ceiling read as cores set aside."""
    log = _run(tmp_path, "--component_workers", "4", budget=8)
    assert "component_worker_ceiling=4" in _allocation_line(log)
    assert "nothing is reserved for it" in log


def test_one_eligible_component_is_refused_a_grant(tmp_path):
    """Forking runs different components concurrently and cannot subdivide one, so a unit
    with a single eligible component must be declined rather than charged for workers.

    The hint is lowered to 1 so the synthetic single-component graph counts as eligible;
    that isolates the backlog rule from the size rule.
    """
    config = tmp_path / "spawn_at_one.json"
    config.write_text('{"min_mpgn_component_size_for_spawn": 1}\n')
    log = _run(
        tmp_path, "--component_workers", "4", "--no_parallelize_contigs",
        "--config_update", str(config), budget=8, quant_only=False,
    )
    assert "only one eligible component outstanding" in log
    assert "-Component assembly granted" not in log


def test_serialising_gives_the_single_worker_the_whole_budget(tmp_path):
    """--no_parallelize_contigs is kept because it is not a multiplier: it caps the queue
    at one worker, and that worker's native tool steps then get the entire budget."""
    line = _allocation_line(_run(tmp_path, "--no_parallelize_contigs", budget=8))
    assert "unit_workers=1" in line
    assert "tool_threads_per_worker=8" in line


def test_launch_order_is_longest_first_on_the_read_count_proxy(tmp_path):
    """ctgA carries three times ctgB's reads, so its units must launch first.

    Construction order would leave the expensive unit for the tail of the run; this is
    the makespan win that applies in single-node mode, where a small reference supplies
    many units of wildly unequal cost.
    """
    log = _run(tmp_path, budget=2)
    assert "longest-first by cost proxy" in log
    launched = re.findall(r"\[(ctg[AB])[+-]\]", log)
    assert launched, "expected per-unit log lines naming the contig"
    assert launched[0] == "ctgA"


def test_the_removed_multiplicative_flags_are_gone(tmp_path):
    """Clean cutover. These are the flags WDL/ passes and must be migrated off."""
    fasta, bam, gtf = _write_inputs(tmp_path)
    for flag, value in (
        ("--num_threads_per_worker", "2"),
        ("--num_parallel_contigs", "4"),
        ("--CPU", "2"),
        ("--no_shuffle_parallel_jobs", None),
    ):
        cmd = [
            sys.executable, LRAA,
            "--genome", fasta, "--bam", bam, "--gtf", gtf, "--quant_only",
            "--output_prefix", "gone", flag,
        ]
        if value is not None:
            cmd.append(value)
        proc = subprocess.run(
            cmd, cwd=str(tmp_path), capture_output=True, text=True,
            timeout=RUN_TIMEOUT_S,
        )
        assert proc.returncode != 0, "{} should no longer be accepted".format(flag)
        assert "unrecognized arguments" in proc.stderr, proc.stderr[-2000:]


def test_a_region_restricted_unit_keys_its_artifacts_on_the_region(tmp_path):
    """Two --region runs sharing an --output_prefix used to resume from each other's
    checkpoint, because the region appeared in no artifact name. Chunk units of one
    contig-strand would collide the same way, on every artifact including the .ok."""
    _run(tmp_path, "--region", "ctgA:1-3000", "--no_cleanup", budget=2)
    tmp_root = tmp_path / "__budget_test.LRAA.quant-only.contigtmp" / "ctgA" / "+"
    names = sorted(p.name for p in tmp_root.iterdir() if p.is_file())
    assert names, "expected per-unit artifacts under {}".format(tmp_root)
    assert all(n.startswith("ctgA.+.chunk0_1-3000.") for n in names), names
    for suffix in (".quant.expr", ".quant.tracking", ".ok"):
        assert "ctgA.+.chunk0_1-3000" + suffix in names
