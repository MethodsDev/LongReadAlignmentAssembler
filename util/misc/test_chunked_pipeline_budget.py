#!/usr/bin/env python3

"""The orchestrator consumes the SAME budget, and cannot multiply against LRAA's.

Before this change the driver took ``--concurrency`` and separately handed LRAA
``--num_parallel_contigs`` and ``--num_threads_per_worker``. Nothing checked the product,
and the product was the only thing that mattered: 4 chunks x 2 contig workers x 2 threads
is 16 cores requested on a host that may have four.

These assertions are on the command the driver BUILDS and on the arithmetic it uses, so
nothing here runs a pipeline.
"""

import argparse
import importlib.util
import os
import sys

import pytest


_HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "..", "pylib"))

import CpuBudget  # noqa: E402


def _load_driver():
    """The pipeline now lives in pylib/ChunkedRun.py; the script here is a front end."""
    import ChunkedRun

    return ChunkedRun


DRIVER = _load_driver()


def _args(**overrides):
    # quant-only, which is the arm whose command shape these tests pin
    base = dict(
        HiFi=False,
        cpu_budget=16,
        discovery=False,
        # lraa_cmd forwards these unconditionally (mapq) or on presence/truthiness
        # (the rest) -- see pylib/test_chunked_quant_token_argv.py's BASE_SETTINGS
        # for the same fixture-drift fix and why cell_barcode_tag/read_umi_tag are
        # deliberately absent rather than defaulted.
        cell_list=None,
        stream_reads=False,
        stream_reads_rescue_unassigned=False,
        stream_reads_rescue_unassigned_to_targets=False,
        min_mapping_quality=0,
        min_mapping_quality_for_final_quant=0,
    )
    base.update(overrides)
    return argparse.Namespace(**base)


def test_the_per_chunk_command_carries_a_budget_and_nothing_multiplicative():
    """Each chunk's LRAA gets its SHARE as --cpu_budget, and none of the removed flags."""
    cmd = DRIVER.lraa_cmd(
        _args(),
        bam_for_quant="chunk.bam",
        bam_for_sg="chunk.norm.bam",
        genome="chunk.fa",
        gtf="chunk.gtf",
        out_prefix="chunk_quant",
        num_total_reads=1000,
        cpu_budget=2,
    )
    assert "--cpu_budget" in cmd
    assert cmd[cmd.index("--cpu_budget") + 1] == "2"
    for gone in ("--num_parallel_contigs", "--num_threads_per_worker", "--CPU"):
        assert gone not in cmd


def test_chunk_concurrency_times_per_chunk_budget_never_exceeds_the_total():
    """The invariant the old --concurrency could not state.

    A chunk is single-contig and strand-pure, so LRAA's internal queue inside it holds one
    unit and its own pool clamps to 1 regardless -- the per-chunk share is not
    double-counted against that inner level.
    """
    for budget in range(1, 33):
        for chunks in (1, 2, 6, 12, 25, 50):
            alloc = CpuBudget.allocate(budget=budget, num_units=chunks)
            assert alloc.unit_workers * alloc.tool_threads <= budget, (budget, chunks)


def test_the_wdl_shape_and_the_chunked_shape_both_fit_sixteen_cores():
    """chr20 at six chunks per strand is twelve units; chr1 at twenty-five is fifty."""
    twelve = CpuBudget.allocate(budget=16, num_units=12)
    assert (twelve.unit_workers, twelve.tool_threads) == (12, 1)
    assert twelve.unit_workers * twelve.tool_threads == 12 <= 16

    fifty = CpuBudget.allocate(budget=16, num_units=50)
    assert (fifty.unit_workers, fifty.tool_threads) == (16, 1)
    assert fifty.unit_workers * fifty.tool_threads == 16


def test_the_driver_no_longer_accepts_the_multiplicative_knobs():
    """Clean cutover: --concurrency is derived now, not supplied."""
    parser_argv = [
        "--bam", "x.bam", "--genome_fa", "x.fa", "--gtf", "x.gtf",
        "--output_dir", "out",
    ]
    for gone in ("--concurrency", "--num_parallel_contigs", "--num_threads_per_worker"):
        with pytest.raises(SystemExit):
            DRIVER.main(parser_argv + [gone, "4"])


def test_a_zero_or_negative_budget_is_refused_rather_than_clamped():
    """A budget the user typed wrong is an error; only derived arithmetic gets clamped."""
    with pytest.raises(SystemExit):
        DRIVER.main([
            "--bam", "x.bam", "--genome_fa", "x.fa", "--gtf", "x.gtf",
            "--output_dir", "out", "--cpu_budget", "0",
        ])


def test_chunk_units_order_longest_first_on_retained_alignments():
    """The extractor already counted retained alignments per chunk, so the proxy is free.

    Span would be the wrong proxy: it does not bound read count, so a short highly
    expressed chunk can outweigh a long quiet one.
    """
    chunks = [
        {"chunk_id": "chr20_plus_00", "chrom": "chr20", "strand": "+", "index": 0,
         "manifest": {"partition_lend": 1, "partition_rend": 10_000_000,
                      "counts": {"alignments_emitted": 5_000}}},
        {"chunk_id": "chr20_plus_01", "chrom": "chr20", "strand": "+", "index": 1,
         "manifest": {"partition_lend": 10_000_001, "partition_rend": 20_000_000,
                      "counts": {"alignments_emitted": 90_000}}},
        {"chunk_id": "chr20_plus_02", "chrom": "chr20", "strand": "+", "index": 2,
         "manifest": {"partition_lend": 20_000_001, "partition_rend": 30_000_000,
                      "counts": {"alignments_emitted": 40_000}}},
    ]
    units = [
        CpuBudget.WorkUnit(
            contig_acc=c["chrom"],
            contig_strand=c["strand"],
            chunk_index=c["index"],
            region=(c["manifest"]["partition_lend"], c["manifest"]["partition_rend"]),
            cost=c["manifest"]["counts"]["alignments_emitted"],
        )
        for c in chunks
    ]
    by_unit = dict(zip(units, chunks))
    order = [by_unit[u]["chunk_id"] for u in CpuBudget.order_longest_first(units)]
    assert order == ["chr20_plus_01", "chr20_plus_02", "chr20_plus_00"]
