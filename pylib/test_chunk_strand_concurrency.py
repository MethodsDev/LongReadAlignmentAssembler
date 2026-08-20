#!/usr/bin/env python3

"""A strandless chunk's two orientations run concurrently when there is cpu share
to spare, and in series -- byte-for-byte today's behavior -- when there is not.

``chunk_worker`` decides this itself, per chunk, from a SECOND application of
``CpuBudget.allocate`` to its own ``cpu_budget`` share: ``num_units=2`` for a
strandless chunk. Below a share of 2 cores that inner allocation gives
``unit_workers == 1`` and the loop is exactly what it always was. At or above it,
the two orientations run on their own thread each, each getting half the share
for its own LRAA invocation, and each writing to its OWN log file rather than the
shared per-chunk one -- two subprocesses appending the same file at once would
interleave into something unreadable.

This is the scenario named for it: a Terra WDL scatter by chromosome, then
chunked, can leave a shard with far fewer chunks than the box has cores -- one
small chromosome's few chunks against sixteen cores, say -- and CpuBudget.allocate
folds every core the chunk count could not claim into each chunk's OWN share
rather than into more concurrent chunks. Before this, that spare share sat mostly
idle across two sequential single-threaded quant invocations.
"""

import os
import sys
import threading

import pysam
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."))
sys.path.insert(0, os.path.join(REPO, "pylib"))

import ChunkedRun  # noqa: E402  (path insert must precede the import)
import CpuBudget  # noqa: E402


# --------------------------------------------------------------------- fixtures


def write_bam(path, num_forward, num_reverse):
    """A bam of ``num_forward`` forward and ``num_reverse`` reverse records."""

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chunk", "LN": 100000}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as ofh:
        for i in range(num_forward + num_reverse):
            aln = pysam.AlignedSegment(ofh.header)
            aln.query_name = "read{}".format(i)
            aln.reference_id = 0
            aln.reference_start = 10 * i
            aln.mapping_quality = 60
            aln.cigarstring = "50M"
            aln.query_sequence = "A" * 50
            aln.flag = 0 if i < num_forward else 16
            ofh.write(aln)
    return path


def gtf_line(chrom, feature, lend, rend, strand, gene_id, transcript_id=None):
    attrs = 'gene_id "{}";'.format(gene_id)
    if transcript_id:
        attrs += ' transcript_id "{}";'.format(transcript_id)
    return (
        "\t".join([chrom, "test", feature, str(lend), str(rend), ".", strand, ".", attrs])
        + "\n"
    )


def make_chunk(tmp_path, strand="", emitted=(4, 3)):
    """A chunk record as stage 3 builds one, with its mini GTF on disk. Matches
    test_strandless_chunking.make_chunk's shape so it can share chunk_worker's
    wiring; ``emitted`` is (forward, reverse) as the extractor would have
    counted them.
    """

    cdir = tmp_path / "chunks" / "chrT_00"
    cdir.mkdir(parents=True, exist_ok=True)
    prefix = str(cdir / "chunk")
    forward, reverse = emitted

    with open(prefix + ".gtf", "wt") as ofh:
        for orientation, index in (("+", 1), ("-", 2)):
            gene, tid = "G{}".format(index), "T{}".format(index)
            lend, rend = 100 * index, 100 * index + 50
            ofh.write(gtf_line("chunk", "gene", lend, rend, orientation, gene))
            ofh.write(gtf_line("chunk", "transcript", lend, rend, orientation, gene, tid))
            ofh.write(gtf_line("chunk", "exon", lend, rend, orientation, gene, tid))

    manifest = {
        "strand": strand or None,
        "strand_split_required": not strand,
        "offset": 0,
        "window_origin": 0,
        "partition_lend": 1,
        "partition_rend": 100000,
        "dropped_read_names": [],
        "counts": {
            "alignments_emitted": forward + reverse,
            "alignments_emitted_forward": forward,
            "alignments_emitted_reverse": reverse,
            "alignments_dropped_overhang": 0,
            "gtf_transcripts_emitted": 2,
        },
    }
    return {
        "chunk_id": "chrT_00",
        "chrom": "chrT",
        "strand": strand,
        "strandless": not strand,
        "region": "chrT:1-100000",
        "index": 0,
        "order": 0,
        "dir": str(cdir),
        "prefix": prefix,
        "log": str(cdir / "chunk.log"),
        "manifest": manifest,
        "offset": 0,
        "window_origin": 0,
        "upstream_token": "stage3.up_aaaaaaaaaaaa",
        "units": ChunkedRun.chunk_quant_units("chrT_00", str(cdir), prefix, strand, 0, 0),
    }


class Recorder:
    """Stands in for run_step: records (name, cmd, log) and writes what the real
    tool would have, so chunk_worker's own output-existence checks keep biting."""

    def __init__(self, emitted):
        self.calls = []
        self.emitted = emitted
        self.lock = threading.Lock()

    def __call__(self, name, cmd, log, cwd, rss_interval, append=True):
        with self.lock:
            self.calls.append((name, list(cmd), log))
        if name.startswith("stage3b"):
            prefix = cmd[cmd.index("--output_prefix") + 1]
            write_bam(prefix + ".+.bam", self.emitted[0], 0)
            write_bam(prefix + ".-.bam", 0, self.emitted[1])
        elif name.startswith("stage4"):
            out = cmd[cmd.index("--output_bam") + 1]
            write_bam(out, 1, 1)
        elif name.startswith("stage5"):
            out = os.path.join(cwd, cmd[cmd.index("--output_prefix") + 1])
            base = out + "." + ChunkedRun.LRAA_QUANT_ONLY_SUFFIX
            open(base + ".quant.expr", "wt").close()
            open(base + ".quant.tracking.gz", "wt").close()
        return {"step": name, "wall_s": 0.0, "peak_tree_rss_kb": 0}



class RendezvousRecorder(Recorder):
    """Like ``Recorder``, but every ``stage4`` (normalize) call must rendezvous
    with a SECOND one before either is allowed to proceed.

    This is the one property the other assertions in this file cannot rule out
    by themselves: distinct log files and a halved --cpu_budget are equally
    consistent with a secretly-sequential implementation that just happens to
    name its two calls differently. Only two threads genuinely being inside
    stage4 AT THE SAME TIME proves concurrency, and a bounded ``Barrier`` is
    the direct way to observe that -- if the two calls are sequential, the
    first one blocks alone until the wait times out and raises
    ``BrokenBarrierError``, which is recorded rather than left to hang the
    test suite.
    """

    def __init__(self, emitted, timeout=5.0):
        super().__init__(emitted)
        self.stage4_barrier = threading.Barrier(2, timeout=timeout)
        self.stage4_rendezvous_failed = False

    def __call__(self, name, cmd, log, cwd, rss_interval, append=True):
        if name.startswith("stage4"):
            try:
                self.stage4_barrier.wait()
            except threading.BrokenBarrierError:
                self.stage4_rendezvous_failed = True
        return super().__call__(name, cmd, log, cwd, rss_interval, append=append)


def run_worker(tmp_path, monkeypatch, chunk, cpu_budget, recorder_cls=Recorder):
    recorder = recorder_cls(
        (
            chunk["manifest"]["counts"]["alignments_emitted_forward"],
            chunk["manifest"]["counts"]["alignments_emitted_reverse"],
        )
    )
    monkeypatch.setattr(ChunkedRun, "run_step", recorder)
    args = ChunkedRun.default_args(strandless_chunks=not chunk["strand"])
    ChunkedRun.chunk_worker(
        args,
        ChunkedRun.Checkpoints(str(tmp_path / "__ckpt")),
        str(tmp_path),
        chunk,
        1000,
        0.5,
        cpu_budget,
    )
    return recorder


# ------------------------------------------------------------- the gate itself


def test_serial_at_cpu_budget_one_is_unchanged_from_today(tmp_path, monkeypatch):
    """The common case (many chunks, one core each): zero behavior change."""

    chunk = make_chunk(tmp_path)
    recorder = run_worker(tmp_path, monkeypatch, chunk, cpu_budget=1)

    steps = [name.split("_", 2)[0] for name, _, _ in recorder.calls]
    assert steps == ["stage3b", "stage4", "stage5", "stage4", "stage5"]

    # Every step, both orientations, wrote to the ONE shared chunk log -- exactly
    # what a reader following today's behavior would expect to find.
    assert {log for _, _, log in recorder.calls} == {chunk["log"]}

    # Nothing was split off this single core.
    stage5_cmds = [cmd for name, cmd, _ in recorder.calls if name.startswith("stage5")]
    assert len(stage5_cmds) == 2
    for cmd in stage5_cmds:
        assert cmd[cmd.index("--cpu_budget") + 1] == "1"


def test_concurrent_at_cpu_budget_four_splits_the_share_and_the_log(
    tmp_path, monkeypatch
):
    """Enough spare share (few chunks, many cores) runs both orientations at once,
    each on its own half of the budget and its own log file. The concurrency
    itself is proven, not inferred from bookkeeping a sequential implementation
    could equally produce: both stage4 calls must rendezvous at a shared
    barrier, which only two THREADS actually in flight together can satisfy.
    """

    chunk = make_chunk(tmp_path)
    recorder = run_worker(
        tmp_path, monkeypatch, chunk, cpu_budget=4, recorder_cls=RendezvousRecorder
    )

    assert not recorder.stage4_rendezvous_failed, (
        "the two orientations' stage4 (normalize) calls never overlapped -- "
        "the barrier timed out waiting for a second party, which is exactly "
        "what a secretly-sequential implementation would produce"
    )

    # Both orientations completed -- same five steps as the serial case, order
    # unconstrained (concurrent), grouped by stage rather than asserted verbatim.
    names = [name for name, _, _ in recorder.calls]
    assert sorted(names) == sorted(
        [
            "stage3b_strand_split_chrT_00",
            "stage4_normalize_strandless_chrT_00_plus",
            "stage5_quant_strandless_chrT_00_plus",
            "stage4_normalize_strandless_chrT_00_minus",
            "stage5_quant_strandless_chrT_00_minus",
        ]
    )
    # The split itself still runs first, on the shared log, before either
    # orientation's own log exists to write to.
    assert names[0] == "stage3b_strand_split_chrT_00"
    assert recorder.calls[0][2] == chunk["log"]

    # Stage 4/5 calls split into two DISTINCT logs, neither the shared one --
    # this is what keeps two concurrent subprocesses' output from interleaving.
    per_unit_logs = {
        log for name, _, log in recorder.calls if not name.startswith("stage3b")
    }
    assert len(per_unit_logs) == 2
    assert chunk["log"] not in per_unit_logs
    for log in per_unit_logs:
        assert log.startswith(chunk["log"])

    # cpu_budget=4 over 2 units -> CpuBudget.allocate(4, 2) = 2 workers, 2 threads
    # each: HALF the chunk's share per orientation, not the whole share twice.
    stage5_cmds = [cmd for name, cmd, _ in recorder.calls if name.startswith("stage5")]
    assert len(stage5_cmds) == 2
    for cmd in stage5_cmds:
        assert cmd[cmd.index("--cpu_budget") + 1] == "2"


def test_a_strand_first_chunk_never_splits_regardless_of_budget(tmp_path, monkeypatch):
    """A strand-first chunk holds ONE unit; CpuBudget.allocate(budget, 1) can only
    ever give unit_workers == 1, so this path is unreachable for it by construction,
    not merely untested -- a positive control at a wide-open budget."""

    chunk = make_chunk(tmp_path, strand="+", emitted=(4, 0))
    recorder = run_worker(tmp_path, monkeypatch, chunk, cpu_budget=8)

    names = [name for name, _, _ in recorder.calls]
    assert names == ["stage4_normalize_chrT_00", "stage5_quant_chrT_00"]

    # No split step at all for a strand-first chunk, and its single unit gets
    # the WHOLE share -- nothing to divide it with.
    stage5_cmd = recorder.calls[1][1]
    assert stage5_cmd[stage5_cmd.index("--cpu_budget") + 1] == "8"
    assert {log for _, _, log in recorder.calls} == {chunk["log"]}


# -------------------------------------------------- the budget invariant itself


def test_nested_allocation_never_exceeds_the_outer_budget():
    """outer.unit_workers * inner.unit_workers * inner.tool_threads <= budget,
    for every (budget, chunk count) a run could hand chunk_worker -- the product
    CpuBudget's own module docstring promises is never multiplied past the
    original budget, now checked one level deeper than before."""

    for budget in range(1, 33):
        for num_chunks in range(1, 33):
            outer = CpuBudget.allocate(budget=budget, num_units=num_chunks)
            # Every concurrently-running chunk applies the SAME inner allocation
            # to its own share (chunk_worker has no per-chunk state to vary it).
            inner = CpuBudget.allocate(budget=outer.tool_threads, num_units=2)
            total = outer.unit_workers * inner.unit_workers * inner.tool_threads
            assert total <= budget, (budget, num_chunks, outer, inner, total)


# ---------------------------------------- the split's own token reaches both units


def test_both_units_cache_key_off_the_splits_own_token_not_the_stale_one(
    tmp_path, monkeypatch
):
    """A strandless chunk's stage4 cache token has to chain onto stage3b's OWN
    returned token, not ``chunk["upstream_token"]`` -- the value set once when the
    chunk dict was built, before the split ever ran. Those two are the SAME value
    only by coincidence of this fixture's stub; in production a real
    ``split_chunk_by_strand`` computes its own token from the split's actual
    inputs, and if the split step's own parameters ever change, a norm_token still
    keyed on the pre-split value would not move -- silently reusing a normalized
    bam built against a different split. Caught by making the two deliberately
    different here and asserting the cache key follows the RIGHT one.
    """

    chunk = make_chunk(tmp_path)
    assert chunk["upstream_token"] == "stage3.up_aaaaaaaaaaaa"

    real_split = ChunkedRun.split_chunk_by_strand

    def fake_split(args, ckpt, chunk_arg, rss_interval):
        step, _stale_token, counts = real_split(args, ckpt, chunk_arg, rss_interval)
        # The split's OWN token, deliberately distinct from chunk["upstream_token"]
        # -- exactly what a real split computes from its own command/inputs.
        return step, "stage3b.up_bbbbbbbbbbbb", counts

    monkeypatch.setattr(ChunkedRun, "split_chunk_by_strand", fake_split)

    seen_parents = []
    real_chain_token = ChunkedRun.chain_token

    def spy_chain_token(local, parent, *opaque):
        if local.startswith("stage4_norm"):
            seen_parents.append(parent)
        return real_chain_token(local, parent, *opaque)

    monkeypatch.setattr(ChunkedRun, "chain_token", spy_chain_token)

    run_worker(tmp_path, monkeypatch, chunk, cpu_budget=1)

    assert seen_parents, "no stage4_norm cache key was ever computed"
    assert seen_parents == ["stage3b.up_bbbbbbbbbbbb"] * len(seen_parents), (
        "both orientations' cache keys must chain onto the split's OWN returned "
        "token, not the stale chunk['upstream_token'] set before the split ran"
    )
