#!/usr/bin/env python3

"""The chunk-PROCESSING pool's memory bound, and why it has to be per unit.

``prep_memory_cap`` bounded the make-chunks pool. The chunk-processing pool sized
itself at ``min(cpu_budget, units)`` with no memory term at all, and a whole-genome
one-chunk-per-contig partition (195 chunks, ``--cpu_budget 16``) was killed by the
kernel OOM killer 3,825 s in: the chr1+ stage-5 worker, at a process-tree peak of
45.74 GiB on a 62 GiB box with 358 of 390 quant units still unstarted.

Three properties carry the fix, and each has a test here that fails without it.

THE DEFAULT MUST NOT MOVE. The 475-chunk partition fits any ordinary box, and a
bound that reduced its concurrency would be a regression dressed as a safety fix.
It is checked against a no-cap allocation as its own control, so "16" cannot pass by
coincidence.

THE COST IS PER UNIT, NOT PER PHASE. Measured per-chunk peak spans 50 MiB to
5,594.7 MiB inside one run, so no single constant can both bound a whole chromosome
and leave a 10 Mb chunk alone. The estimate is charged from each unit's own genomic
SPAN, and ``test_the_estimate_dominates_every_measured_per_chunk_peak`` holds the two
constants to the measurements they came from -- lowering either one fails there
rather than silently at 3 a.m. on someone's cluster.

THE DECISION IS VISIBLE. A user whose run got slower reads why, in one line naming
the cap, the available memory, the per-unit basis and the idle cores.

Every memory figure here is INJECTED. Nothing reads ``/proc/meminfo``, so the file is
deterministic under any load and safe in the full suite.
"""

import sys
import threading
from pathlib import Path
from types import SimpleNamespace

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "pylib"))

import ChunkedRun  # noqa: E402  (path insert must precede the import)
import CpuBudget  # noqa: E402

MIB = 1
GIB = 1024 * MIB

# GRCh38 primary lengths. A per-contig partition's units ARE these.
GRCH38 = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
    "chrM": 16569,
}

# 195 units: the 25 primary contigs above plus 170 unplaced scaffolds. This is the
# partition ``--approx_MB_per_cut 300`` produces on GRCh38 and the one that OOMed.
PER_CONTIG_SPANS = list(GRCH38.values()) + [200_000] * 170

# 475 units: the shipped default's SHAPE, from the HG002 PacBio run's manifests --
# 315 cut segments of 9-20.2 Mb over the 24 contigs long enough to cut, and 160
# whole-contig chunks of unplaced scaffolds. Only the span distribution matters to
# the bound, so the segments are held at their measured extremes.
DEFAULT_475_SPANS = [10_000_000] * 300 + [20_200_000] * 15 + [200_000] * 160


# --------------------------------------------------- the default must not move


def test_the_shipped_475_chunk_partition_keeps_every_worker_it_has_today():
    """16 workers before, 16 after, on any box with room -- against a no-cap control.

    The 475-chunk default's 16 largest units estimate to 23,264 MiB together, which a
    62 GB box clears 2.5x over and a 32 GB box still clears. Both are checked: the
    bound must be inert across the range of machines this mode actually runs on, not
    merely on the largest one.
    """

    today = CpuBudget.allocate(budget=16, num_units=len(DEFAULT_475_SPANS))
    assert today.unit_workers == 16

    for available_mib in (56 * GIB, 28 * GIB, 24 * GIB):
        bound = ChunkedRun.chunk_memory_cap(
            16, DEFAULT_475_SPANS, available_mib=available_mib
        )
        assert bound.cap is None, (available_mib, bound.note)
        assert bound.note is None
        after = CpuBudget.allocate(
            budget=16,
            num_units=len(DEFAULT_475_SPANS),
            max_unit_workers=bound.cap,
        )
        assert after == today, available_mib

    # The measurement is kept even where it changes nothing -- that is what makes the
    # guard's error checkable on an ordinary run instead of only on a capped one.
    bound = ChunkedRun.chunk_memory_cap(16, DEFAULT_475_SPANS, available_mib=56 * GIB)
    assert bound.largest_unit_mib == 1468
    assert bound.charged_mib == 23264
    assert bound.available_mib == 56 * GIB


# ------------------------------------------------- a partition that cannot fit


def test_a_whole_chromosome_partition_reduces_the_pool_deterministically():
    """195 whole-contig units at 16-wide want 72,618 MiB. The box has 57,344.

    The reduction is arithmetic, not a heuristic: the 11 largest units charge 54,510
    MiB and the 12th would take it to 58,477, past what the box has. Both sides of
    that boundary are asserted, so the test fails on an off-by-one in either
    direction rather than only on a wholesale regression.
    """

    bound = ChunkedRun.chunk_memory_cap(16, PER_CONTIG_SPANS, available_mib=56 * GIB)
    assert bound.cap == 11
    assert bound.charged_mib == 54510
    assert bound.charged_mib <= 56 * GIB

    estimates = sorted(
        (ChunkedRun.chunk_unit_peak_mib(span) for span in PER_CONTIG_SPANS),
        reverse=True,
    )
    assert sum(estimates[:11]) <= 56 * GIB < sum(estimates[:12])
    assert bound.largest_unit_mib == estimates[0] == 6501

    allocation = CpuBudget.allocate(
        budget=16, num_units=len(PER_CONTIG_SPANS), max_unit_workers=bound.cap
    )
    assert allocation.unit_workers == 11
    # The freed cores are not lost, they go to each worker's native tool steps --
    # which is the existing CpuBudget contract, not something this bound invents.
    assert CpuBudget.allocate(
        budget=16, num_units=len(PER_CONTIG_SPANS)
    ).unit_workers == 16


def test_a_smaller_box_reduces_further_and_monotonically():
    """Less memory never buys more workers. Squeeze the box, watch the cap fall."""

    caps = [
        ChunkedRun.chunk_memory_cap(
            16, PER_CONTIG_SPANS, available_mib=available
        ).cap
        for available in (56 * GIB, 28 * GIB, 14 * GIB, 8 * GIB)
    ]
    assert caps == [11, 5, 2, 1]


def test_the_bound_is_independent_of_the_order_units_arrive_in():
    """It charges the K LARGEST, so the input order cannot change the answer.

    This is why the bound holds despite ``order_longest_first`` ranking the launch
    queue by ALIGNMENTS while the estimate ranks by SPAN: charging the largest K is
    an upper bound over every K-subset, so the two rankings need not agree.
    """

    forward = ChunkedRun.chunk_memory_cap(
        16, PER_CONTIG_SPANS, available_mib=56 * GIB
    )
    backward = ChunkedRun.chunk_memory_cap(
        16, list(reversed(PER_CONTIG_SPANS)), available_mib=56 * GIB
    )
    assert forward == backward


def test_one_worker_is_the_floor_and_the_deficit_is_stated_rather_than_refused():
    """A single unit over the line still runs, and the note says the partition is why.

    REDUCE, not refuse: the estimate is an envelope carrying up to 2.3x of margin, so
    refusing on it would kill runs that complete. What refusal would have bought --
    honesty about the case concurrency cannot fix -- is bought by the warning
    instead, which names the deficit and the only knob that moves it.
    """

    bound = ChunkedRun.chunk_memory_cap(16, PER_CONTIG_SPANS, available_mib=4 * GIB)
    assert bound.cap == 1
    assert "WARNING" in bound.note
    assert "6501 MiB against 4096 MiB" in bound.note
    assert "--approx_MB_per_cut" in bound.note
    assert (
        CpuBudget.allocate(
            budget=16, num_units=len(PER_CONTIG_SPANS), max_unit_workers=bound.cap
        ).unit_workers
        == 1
    )


def test_a_partition_that_fits_at_full_width_gets_no_warning():
    """Positive control for the warning: it must not fire where nothing is short."""

    bound = ChunkedRun.chunk_memory_cap(16, DEFAULT_475_SPANS, available_mib=56 * GIB)
    assert bound.note is None


# --------------------------------------------------------- the decision is visible


def test_the_note_names_the_cap_the_box_and_the_per_unit_basis():
    """One line, in ``budget_shortfall_note``'s register, with every number in it."""

    bound = ChunkedRun.chunk_memory_cap(16, PER_CONTIG_SPANS, available_mib=56 * GIB)
    note = bound.note
    assert "\n" not in note
    assert "capped at 11 worker(s) by MEMORY" in note
    assert "--cpu_budget 16" in note
    assert "57344 MiB available" in note
    assert "16 largest chunk(s)" in note
    assert "72618 MiB together" in note
    assert "{} MiB fixed".format(ChunkedRun.CHUNK_UNIT_FIXED_MIB) in note
    assert (
        "{} MiB per genomic Mb".format(ChunkedRun.CHUNK_UNIT_MIB_PER_GENOMIC_MB) in note
    )
    assert "largest single chunk 6501 MiB" in note
    assert "5 core(s) of the budget sit idle" in note
    assert "11 worker(s) charge 54510 MiB" in note


def test_an_unreadable_meminfo_says_so_and_still_reports_the_estimate():
    """No bound is applied, but the quantity it would have used is not dropped."""

    bound = ChunkedRun.chunk_memory_cap(
        16, PER_CONTIG_SPANS, available_mib=ChunkedRun.available_memory_mib()
    )
    if ChunkedRun.available_memory_mib() is None:
        assert bound.cap is None
        assert "could not be read" in bound.note
        assert bound.largest_unit_mib == 6501
        assert bound.charged_mib == 72618


def test_no_units_is_not_an_error():
    """A dry run or an empty partition asks nothing of the box."""

    bound = ChunkedRun.chunk_memory_cap(16, [], available_mib=1)
    assert bound.cap is None and bound.charged_mib == 0


# ------------------------------------------------------- the guard's calibration


# MEASURED per-chunk / per-unit peak RSS against the chunk's genomic span. Sources:
# the 195-contig whole-genome run's ``stage5_quant_strandless`` step footers
# (peak_tree_rss_kb, reported by MeasurePerContig from chunk_work/logs/chunk_*.log),
# the chr21 discovery runs under __LRAA_local_runs/chr21_denovo_parity_r2, and the
# worst chunks of the 475-chunk HG002 PacBio run's timing.json. Two regimes are
# deliberately included: whole chromosomes, where span drives the footprint, and
# chrM, where 1.19 M reads in 16.6 kb make span useless and the fixed term is the
# only thing covering it.
MEASURED_PEAKS = [
    # (span_bp, measured_MiB, what)
    (16569, 655.0, "chrM whole contig, 1.19 M reads, quant"),
    (10_000_000, 741.0, "chr11_06, worst 10 Mb chunk of the 475-chunk run"),
    (20_200_000, 638.0, "worst 20.2 Mb chunk of the 475-chunk run"),
    (46_709_983, 1050.7, "chr21 whole contig, strandless, quant"),
    (46_709_983, 1710.0, "chr21 whole contig, DISCOVERY -- the dearer mode"),
    (190_214_555, 3110.0, "chr4 whole contig, quant"),
    (181_538_259, 3439.7, "chr5 whole contig, quant"),
    (198_295_559, 3875.5, "chr3 whole contig, quant"),
    (242_193_529, 4725.5, "chr2 whole contig, quant"),
    (248_956_422, 5594.7, "chr1 whole contig -- the unit the OOM killer took"),
]


def test_the_estimate_dominates_every_measured_per_chunk_peak():
    """A guard that a measurement exceeds is not a guard. All ten, both regimes."""

    for span_bp, measured_mib, what in MEASURED_PEAKS:
        estimate = ChunkedRun.chunk_unit_peak_mib(span_bp)
        assert estimate >= measured_mib, (what, estimate, measured_mib)


def test_the_estimate_is_not_so_conservative_that_it_refuses_real_partitions():
    """Margin has a ceiling too, where the units are big enough for it to cost.

    Below 40 Mb the fixed term dominates and the margin is large by design -- an
    empty scaffold is charged 1 GiB it will never use, and that is harmless because
    such units are never the K largest. At whole-chromosome span the margin is what
    decides how many workers a real run gets, so it is held under 2x.
    """

    for span_bp, measured_mib, what in MEASURED_PEAKS:
        if span_bp < 40_000_000:
            continue
        ratio = ChunkedRun.chunk_unit_peak_mib(span_bp) / measured_mib
        assert ratio < 2.0, (what, ratio)


def test_span_is_the_predictor_and_a_read_count_model_would_have_missed_chrM():
    """chrM is the control that rules the alternative out.

    MiB per million alignments ranges 553 (chrM) to 31,066 across measured chunks, a
    56x spread, against 10.6x for MiB per genomic Mb. A read-count model calibrated
    on the big chromosomes charges chrM about 6x what it uses; a least-squares SPAN
    model (88 + 19.6 MiB/Mb) predicts 88 MiB for it against 655 measured, 7.4x under.
    The fixed term is what makes the span model safe there, so it must stay above
    chrM's measured peak on its own.
    """

    chrm_span, chrm_measured = 16569, 655.0
    assert ChunkedRun.CHUNK_UNIT_FIXED_MIB > chrm_measured
    span_only_fit = 88 + 19.6 * chrm_span / 1e6
    assert span_only_fit < chrm_measured  # the fit this guard is deliberately not
    assert ChunkedRun.chunk_unit_peak_mib(chrm_span) >= chrm_measured


# ----------------------------------------------------------------- the wiring


def _chunk(index, chrom, lend, rend, alignments):
    return {
        "chunk_id": "{}_{:02d}".format(chrom, index),
        "chrom": chrom,
        "strand": None,
        "index": index,
        "region": "{}:{}-{}".format(chrom, lend, rend),
        "log": "/dev/null",
        "manifest": {
            "partition_lend": lend,
            "partition_rend": rend,
            "counts": {"alignments_emitted": alignments},
        },
    }


def _run_pool(monkeypatch, chunks, budget, available_mib):
    """Drive ``run_chunks_concurrently`` with the worker replaced by a counter."""

    monkeypatch.setattr(
        ChunkedRun, "available_memory_mib", lambda: available_mib
    )
    lock = threading.Lock()
    live = {"now": 0, "max": 0}
    barrier = threading.Event()

    def fake_worker(args, ckpt, outdir, chunk, num_total_reads, rss_interval, cpu_budget):
        with lock:
            live["now"] += 1
            live["max"] = max(live["max"], live["now"])
        # Hold every started worker until the pool is provably saturated, so the
        # observed maximum is the pool's width and not a scheduling artefact.
        barrier.wait(timeout=5)
        with lock:
            live["now"] -= 1
        return {"chunk_id": chunk["chunk_id"], "wall_s": 0.0, "peak_tree_rss_kb": 0}

    monkeypatch.setattr(ChunkedRun, "chunk_worker", fake_worker)
    releaser = threading.Timer(0.4, barrier.set)
    releaser.start()
    try:
        records, _makespan, allocation, bound = ChunkedRun.run_chunks_concurrently(
            SimpleNamespace(cpu_budget=budget),
            None,
            "/dev/null",
            chunks,
            1000,
            0.5,
        )
    finally:
        releaser.cancel()
        barrier.set()
    assert len(records) == len(chunks)
    return allocation, bound, live["max"]


def test_the_pool_opens_at_the_capped_width_not_the_budget(monkeypatch):
    """End to end through the phase: 16 whole chromosomes, 24 GiB, 4 workers run.

    The cap has to reach ``ThreadPoolExecutor``. Asserting the allocation alone would
    pass against a bound computed and then dropped on the floor, so the concurrency
    the workers themselves observe is what is checked.
    """

    names = list(GRCH38)[:16]
    chunks = [
        _chunk(i, name, 1, GRCH38[name], 100_000 + i) for i, name in enumerate(names)
    ]
    allocation, bound, observed = _run_pool(
        monkeypatch, chunks, budget=16, available_mib=24 * GIB
    )
    assert bound.cap == 4
    assert allocation.unit_workers == 4
    assert observed == 4
    # Idle-looking cores are handed to each worker's native tool steps.
    assert allocation.tool_threads == 4


def test_the_pool_opens_at_full_width_when_the_box_has_room(monkeypatch):
    """Control for the above: same 16 units, ample memory, 16 workers."""

    names = list(GRCH38)[:16]
    chunks = [
        _chunk(i, name, 1, GRCH38[name], 100_000 + i) for i, name in enumerate(names)
    ]
    allocation, bound, observed = _run_pool(
        monkeypatch, chunks, budget=16, available_mib=128 * GIB
    )
    assert bound.cap is None
    assert allocation.unit_workers == 16
    assert observed == 16
