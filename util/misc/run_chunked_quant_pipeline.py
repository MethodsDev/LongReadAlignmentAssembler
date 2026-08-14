#!/usr/bin/env python3

"""Run the chunked quant-only pipeline end to end, or the whole-contig control.

This is a thin driver over six tools that already exist and are already tested.
It adds no quantification logic of its own; it wires the stages together, runs
the per-chunk work concurrently, and records what each stage cost.

    stage 1  util/separate_bam_by_strand.py         orientation split
    stage 2  util/misc/select_contig_cut_points.py  cut selection
    stage 3  util/misc/extract_contig_region_inputs.py  chunk extraction
    stage 4  util/normalize_bam_by_strand.py        per-chunk normalization
    stage 5  LRAA --quant_only                      per-chunk quantification
    stage 6  merge + coordinate translation         this script

Two arms, selected by ``--arm``, sharing stage 1 and every downstream setting so
that they differ in PARTITIONING ALONE:

  ``chunked``   stages 1-6.
  ``baseline``  stage 1, then the two orientation bams merged back into one
                whole-contig bam, normalized with ``--window_origin 0``, and
                quantified in a single LRAA run. ``--window_origin 0`` is not
                optional: the normalizer's no-option default anchors the
                depth-window grid on the first aligned base seen per contig, and
                no chunk's grid can match that, so a default-path control would
                fail equivalence for a reason that has nothing to do with
                chunking.

RESUMABILITY, stated rather than left accidental: this script IS resumable, via
sentinel files under ``<outdir>/__ckpt/``. A stage whose sentinel exists is
skipped and its recorded timing is carried forward from ``<outdir>/timing.json``
marked ``reused``. Sentinels are named for the stage and, where the stage has
parameters that change its output, for those parameters, so changing a parameter
invalidates the sentinel. To force a rerun, delete the sentinel or the output
directory. Nothing here is deleted on success: every per-chunk log is kept.

ONE CPU BUDGET, not two. ``--cpu_budget`` is the total for the whole pipeline. Chunk
concurrency and the ``--cpu_budget`` handed to each per-chunk LRAA invocation are both
DERIVED from it by ``pylib/CpuBudget.allocate``, so their product cannot exceed it. This
replaces ``--concurrency`` plus a separately-specified ``--num_parallel_contigs`` and
``--num_threads_per_worker``, which multiplied: nothing checked that the chunk pool times
LRAA's own contig pool times its thread count fit the host. A chunk is single-contig and
strand-pure, so LRAA's internal queue inside it is one unit and its own pool clamps to 1
regardless -- the share is not double-counted against that.

MEASUREMENT NOTES. Wall time is measured around each subprocess. Peak RSS is
sampled from ``/proc`` at ``--rss_sample_interval`` over the step's whole
process tree and summed across that tree, so a chunk's figure already includes
whatever workers its LRAA run spawned; the arm-level figure is the same sum over
this driver's entire descendant tree, which is the real concurrent footprint.
Sampling can miss a spike shorter than the interval.
"""

import argparse
import collections
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time

sys.path.insert(
    0,
    os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "..", "..", "pylib"]),
)
import CpuBudget  # noqa: E402  (path insert must precede the import)
import Util_funcs  # noqa: E402
import LRAA_Globals  # noqa: E402

REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..")
)
SEPARATE_BAM = os.path.join(REPO_ROOT, "util", "separate_bam_by_strand.py")
SELECT_CUTS = os.path.join(REPO_ROOT, "util", "misc", "select_contig_cut_points.py")
EXTRACT_CHUNK = os.path.join(
    REPO_ROOT, "util", "misc", "extract_contig_region_inputs.py"
)
NORMALIZE_BAM = os.path.join(REPO_ROOT, "util", "normalize_bam_by_strand.py")
LRAA = os.path.join(REPO_ROOT, "LRAA")

STRAND_TAG = {"+": "plus", "-": "minus"}
PAGE_KB = os.sysconf("SC_PAGE_SIZE") // 1024
# LRAA rewrites --output_prefix to "<prefix>.LRAA.quant-only" in quant-only mode
# (LRAA:_append_lraa_output_mode_suffix), so the files land under that name.
LRAA_QUANT_ONLY_SUFFIX = "LRAA.quant-only"


class PipelineError(RuntimeError):
    pass


# ---------------------------------------------------------------- rss sampling


def _proc_rss_kb(pid):
    try:
        with open("/proc/{}/statm".format(pid), "rt") as fh:
            return int(fh.read().split()[1]) * PAGE_KB
    except (OSError, ValueError, IndexError):
        return 0


def _ppid_map():
    """pid -> ppid for every visible process. One pass, no dependencies."""

    parents = {}
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            with open("/proc/{}/stat".format(entry), "rt") as fh:
                data = fh.read()
        except OSError:
            continue
        # comm may contain spaces and parentheses; ppid is the field after the
        # closing parenthesis of comm.
        close = data.rfind(")")
        if close < 0:
            continue
        fields = data[close + 2 :].split()
        if len(fields) < 2:
            continue
        try:
            parents[int(entry)] = int(fields[1])
        except ValueError:
            continue
    return parents


def _tree_rss_kb(root_pid):
    """Summed RSS of root_pid and every descendant, in KB."""

    parents = _ppid_map()
    children = collections.defaultdict(list)
    for pid, ppid in parents.items():
        children[ppid].append(pid)
    total = 0
    stack = [root_pid]
    seen = set()
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        total += _proc_rss_kb(pid)
        stack.extend(children.get(pid, ()))
    return total


class RssSampler(threading.Thread):
    """Samples the summed RSS of a process tree until stopped."""

    def __init__(self, root_pid, interval):
        super().__init__(daemon=True)
        self.root_pid = root_pid
        self.interval = interval
        self.peak_kb = 0
        self.samples = 0
        # NOT _stop: threading.Thread already owns that name as a method.
        self._halt = threading.Event()

    def run(self):
        while not self._halt.is_set():
            observed = _tree_rss_kb(self.root_pid)
            if observed > self.peak_kb:
                self.peak_kb = observed
            self.samples += 1
            self._halt.wait(self.interval)

    def stop(self):
        self._halt.set()
        self.join(timeout=self.interval * 4)


# ------------------------------------------------------------- process running


def run_step(name, cmd, log_path, cwd, rss_interval, append=True):
    """Run one subprocess, logging to log_path. Returns a timing record.

    Raises PipelineError naming the step and its log on nonzero exit. The log is
    never removed, on success or on failure.
    """

    mode = "at" if append and os.path.exists(log_path) else "wt"
    started = time.time()
    with open(log_path, mode) as log:
        print(
            "\n===== step {} =====\ncwd: {}\ncmd: {}\nstarted: {}\n".format(
                name, cwd, " ".join(cmd), time.strftime("%Y-%m-%d %H:%M:%S")
            ),
            file=log,
            flush=True,
        )
        proc = subprocess.Popen(
            cmd, stdout=log, stderr=subprocess.STDOUT, cwd=cwd, close_fds=True
        )
        sampler = RssSampler(proc.pid, rss_interval)
        sampler.start()
        try:
            returncode = proc.wait()
        finally:
            sampler.stop()
        elapsed = time.time() - started
        print(
            "\n===== step {} exit {} wall {:.1f}s peak_tree_rss {} KB "
            "({} samples) =====".format(
                name, returncode, elapsed, sampler.peak_kb, sampler.samples
            ),
            file=log,
            flush=True,
        )

    record = {
        "step": name,
        "cmd": cmd,
        "log": log_path,
        "wall_s": round(elapsed, 3),
        "peak_tree_rss_kb": sampler.peak_kb,
        "rss_samples": sampler.samples,
        "exit": returncode,
    }
    if returncode != 0:
        raise PipelineError(
            "step {} FAILED with exit {}; log: {}\ncommand: {}".format(
                name, returncode, log_path, " ".join(cmd)
            )
        )
    return record


# ---------------------------------------------------------------- checkpoints


class Checkpoints:
    def __init__(self, root):
        self.root = root
        os.makedirs(root, exist_ok=True)

    def path(self, token):
        return os.path.join(self.root, "{}.ok".format(token))

    def done(self, token):
        return os.path.exists(self.path(token))

    def mark(self, token, note=""):
        with open(self.path(token), "wt") as fh:
            print(note or time.strftime("%Y-%m-%d %H:%M:%S"), file=fh)


# ------------------------------------------------------------ coord translation

_COORD_LIST = re.compile(r"\d+")


def _shift_coord_string(text, offset):
    """Add offset to every integer in the bracketed coordinate list.

    ``exons``/``introns`` in quant.expr look like ``chr20:(+)[[5, 90], [200, 310]]``.
    The contig name ahead of the first ``[`` may itself carry digits, so only the
    bracketed tail is shifted, and a string with no bracketed list is an error
    rather than a silent pass-through.
    """

    if not text:
        return text
    head, sep, tail = text.partition("[")
    if not sep:
        raise PipelineError(
            "coordinate string {!r} has no bracketed coordinate list".format(text)
        )
    return head + _COORD_LIST.sub(
        lambda m: str(int(m.group(0)) + offset), sep + tail
    )


# --------------------------------------------------------------------- stage 1


def stage_strand_split(args, ckpt, outdir, timing, rss_interval):
    split_dir = os.path.join(outdir, "strand_split")
    os.makedirs(split_dir, exist_ok=True)
    prefix = os.path.join(split_dir, "input")
    log = os.path.join(outdir, "logs", "stage1_strand_split.log")
    token = "stage1_strand_split.maxintron_{}".format(args.max_intron_length)
    cmd = [
        sys.executable,
        SEPARATE_BAM,
        "--bam",
        os.path.abspath(args.bam),
        "--output_prefix",
        prefix,
        "--max_intron_length",
        str(args.max_intron_length),
    ]
    if ckpt.done(token):
        timing.setdefault("stages", {})["strand_split"] = {
            "reused": True,
            "cmd": cmd,
            "log": log,
        }
    else:
        record = run_step("stage1_strand_split", cmd, log, outdir, rss_interval)
        timing.setdefault("stages", {})["strand_split"] = record
        ckpt.mark(token)
    strand_bams = {s: "{}.{}.bam".format(prefix, s) for s in ("+", "-")}
    for strand, path in strand_bams.items():
        if not os.path.exists(path):
            raise PipelineError(
                "stage 1 produced no {} bam at {}; log: {}".format(strand, path, log)
            )
    return strand_bams


def count_records(bam):
    out = subprocess.check_output(["samtools", "view", "-c", bam], text=True)
    return int(out.strip())


# --------------------------------------------------------------------- stage 2


def stage_select_cuts(args, ckpt, outdir, timing, strand_bams, rss_interval):
    cut_dir = os.path.join(outdir, "cuts")
    os.makedirs(cut_dir, exist_ok=True)
    selections = {}
    for strand, bam in strand_bams.items():
        tag = STRAND_TAG[strand]
        prefix = os.path.join(cut_dir, tag)
        log = os.path.join(outdir, "logs", "stage2_cuts_{}.log".format(tag))
        # ".sev" is a version marker, not an input: this stage now also emits
        # <prefix>.severed_reads.bam, so a checkpoint written before that would
        # skip the step and leave the file absent while reporting the stage reused.
        # --HiFi raises min_per_id to 97 inside LRAA, and the emitted severed set is
        # filtered on the value the quant step will use, so the thresholds belong in
        # the token: the same cuts with a different min_per_id emit a different bam.
        effective_min_per_id = 97.0 if args.HiFi else LRAA_Globals.config["min_per_id"]
        effective_min_mapq = int(LRAA_Globals.config["min_mapping_quality"])
        token = "stage2_cuts_{}.mb_{}_wig_{}_dw_{}_margin_{}.sev_pid_{}_mq_{}".format(
            tag,
            args.approx_MB_per_cut,
            args.approx_MB_per_cut_wiggle_window,
            args.depth_window,
            args.margin,
            effective_min_per_id,
            effective_min_mapq,
        )
        cmd = [
            sys.executable,
            SELECT_CUTS,
            "--bam",
            bam,
            "--genome_fa",
            os.path.abspath(args.genome_fa),
            "--gtf",
            os.path.abspath(args.gtf),
            "--approx_MB_per_cut",
            str(args.approx_MB_per_cut),
            "--approx_MB_per_cut_wiggle_window",
            str(args.approx_MB_per_cut_wiggle_window),
            "--depth_window",
            str(args.depth_window),
            "--grid_origin",
            "0",
            "--margin",
            str(args.margin),
            "--max_intron_length",
            str(args.max_intron_length),
            # Always, not on request.  A consumer deciding whether a chunk boundary
            # dissolved a read-sharing component needs the severed alignments
            # themselves: names cannot be fetched from a coordinate-indexed bam,
            # and a span cannot answer compatibility, which follows exon blocks.
            # Without it the only alternative is trusting per-chunk components
            # silently, and the set is small by construction -- a cut severing many
            # reads is one the selector rejects.
            "--severed_reads_bam",
            "{}.severed_reads.bam".format(prefix),
            # The values LRAA will apply, not the selector's own defaults: --HiFi
            # changes min_per_id downstream and the selector cannot see that flag.
            "--min_per_id",
            str(effective_min_per_id),
            "--min_mapping_quality",
            str(effective_min_mapq),
            "--output_prefix",
            prefix,
        ]
        # the bam is orientation-pure already, so --strand is omitted (every
        # record counts); the orientation for the region strings comes from which
        # bam this is, recorded below.
        if args.contig:
            cmd += ["--contig", args.contig]
        key = "cuts_{}".format(tag)
        if ckpt.done(token):
            timing.setdefault("stages", {})[key] = {
                "reused": True,
                "cmd": cmd,
                "log": log,
            }
        else:
            timing.setdefault("stages", {})[key] = run_step(
                "stage2_cuts_{}".format(tag), cmd, log, outdir, rss_interval
            )
            ckpt.mark(token)
        with open("{}.cuts.json".format(prefix), "rt") as fh:
            selections[strand] = json.load(fh)
    return selections, cut_dir


# --------------------------------------------------------------------- stage 3


def stage_extract_chunks(args, ckpt, outdir, timing, strand_bams, selections, rss):
    chunk_root = os.path.join(outdir, "chunks")
    os.makedirs(chunk_root, exist_ok=True)
    chunks = []
    for strand in ("+", "-"):
        tag = STRAND_TAG[strand]
        for selection in selections[strand]:
            chrom = selection["chrom"]
            for idx, segment in enumerate(selection["segments"]):
                region = "{}{}:{}-{}".format(
                    chrom, strand, segment["lend"], segment["rend"]
                )
                chunk_id = "{}_{}_{:02d}".format(chrom, tag, idx)
                cdir = os.path.join(chunk_root, chunk_id)
                os.makedirs(cdir, exist_ok=True)
                prefix = os.path.join(cdir, "chunk")
                log = os.path.join(outdir, "logs", "chunk_{}.log".format(chunk_id))
                token = "stage3_extract_{}.margin_{}".format(chunk_id, args.margin)
                cmd = [
                    sys.executable,
                    EXTRACT_CHUNK,
                    "--genome_fa",
                    os.path.abspath(args.genome_fa),
                    "--bam",
                    strand_bams[strand],
                    "--gtf",
                    os.path.abspath(args.gtf),
                    "--region",
                    region,
                    "--output_prefix",
                    prefix,
                    "--margin",
                    str(args.margin),
                    "--secondary_alignments",
                    "exclude",
                ]
                if ckpt.done(token):
                    record = {"reused": True, "cmd": cmd, "log": log}
                else:
                    record = run_step(
                        "stage3_extract_{}".format(chunk_id), cmd, log, outdir, rss
                    )
                    ckpt.mark(token)
                with open("{}.partition.json".format(prefix), "rt") as fh:
                    manifest = json.load(fh)
                chunks.append(
                    {
                        "chunk_id": chunk_id,
                        "chrom": chrom,
                        "strand": strand,
                        "region": region,
                        "index": idx,
                        "dir": cdir,
                        "prefix": prefix,
                        "log": log,
                        "manifest": manifest,
                        "offset": manifest["offset"],
                        "window_origin": manifest["window_origin"],
                        "extract": record,
                    }
                )
    timing.setdefault("stages", {})["extract_chunks"] = [
        c["extract"] for c in chunks
    ]
    return chunks


# ---------------------------------------------------------------- stages 4 + 5


def lraa_cmd(
    args, bam_for_quant, bam_for_sg, genome, gtf, out_prefix, num_total_reads, cpu_budget
):
    """One LRAA invocation. ``cpu_budget`` is this invocation's SHARE of the total.

    A chunk is single-contig and strand-pure, so LRAA's own flat queue inside it has
    exactly one unit: CpuBudget.allocate() gives 1 unit worker and lends the whole share
    to that worker's native tool steps. Nothing is double-counted, because the share is
    already what this pipeline's own chunk concurrency left over.
    """
    cmd = [
        sys.executable,
        LRAA,
        "--genome",
        genome,
        "--bam",
        bam_for_quant,
        "--gtf",
        gtf,
        "--quant_only",
        "--bam_for_sg",
        bam_for_sg,
        "--no_norm",
        "--num_total_reads",
        str(num_total_reads),
        "--cpu_budget",
        str(cpu_budget),
        "--output_prefix",
        out_prefix,
    ]
    if args.HiFi:
        cmd.append("--HiFi")
    return cmd


def chunk_worker(args, ckpt, outdir, chunk, num_total_reads, rss_interval, cpu_budget):
    """Stages 4 and 5 for one chunk. Everything goes to the chunk's own log.

    ``cpu_budget`` is this chunk's share of the total, handed down by the scheduler.
    """

    cid = chunk["chunk_id"]
    log = chunk["log"]
    cdir = chunk["dir"]
    prefix = chunk["prefix"]
    started = time.time()
    steps = []

    norm_bam = os.path.join(cdir, "chunk.norm.bam")
    norm_token = "stage4_norm_{}.maxcov_{}_dw_{}_seed_{}_origin_{}".format(
        cid,
        args.normalize_max_cov_level,
        args.depth_window,
        args.random_seed,
        chunk["window_origin"],
    )
    norm_cmd = [
        sys.executable,
        NORMALIZE_BAM,
        "--input_bam",
        "{}.bam".format(prefix),
        "--output_bam",
        norm_bam,
        "--normalize_max_cov_level",
        str(args.normalize_max_cov_level),
        "--depth_window",
        str(args.depth_window),
        "--random_seed",
        str(args.random_seed),
        "--max_intron_length",
        str(args.max_intron_length),
        "--input_is_single_strand",
        "--window_origin",
        str(chunk["window_origin"]),
    ]
    if ckpt.done(norm_token):
        steps.append({"step": "stage4_normalize", "reused": True, "cmd": norm_cmd})
    else:
        steps.append(
            run_step("stage4_normalize_{}".format(cid), norm_cmd, log, cdir, rss_interval)
        )
        ckpt.mark(norm_token)

    # LRAA derives its scratch roots by string concatenation on --output_prefix
    # ("__{prefix}.contigtmp", "__{prefix}.sgcache"), so an ABSOLUTE prefix would
    # produce nonsense paths like "__/abs/path.contigtmp" rooted at the cwd.
    # Give it a bare name and let cwd place the outputs.
    quant_prefix = os.path.join(cdir, "chunk_quant." + LRAA_QUANT_ONLY_SUFFIX)
    quant_token = "stage5_quant_{}.N_{}_hifi_{}".format(cid, num_total_reads, args.HiFi)
    quant_cmd = lraa_cmd(
        args,
        bam_for_quant="{}.bam".format(prefix),
        bam_for_sg=norm_bam,
        genome="{}.fa".format(prefix),
        gtf="{}.gtf".format(prefix),
        out_prefix="chunk_quant",
        num_total_reads=num_total_reads,
        cpu_budget=cpu_budget,
    )
    if ckpt.done(quant_token):
        steps.append({"step": "stage5_quant", "reused": True, "cmd": quant_cmd})
    else:
        steps.append(
            run_step("stage5_quant_{}".format(cid), quant_cmd, log, cdir, rss_interval)
        )
        ckpt.mark(quant_token)

    for suffix in (".quant.expr", ".quant.tracking"):
        path = quant_prefix + suffix
        if not os.path.exists(path):
            raise PipelineError(
                "chunk {} produced no {}; log: {}".format(cid, path, log)
            )

    chunk["quant_prefix"] = quant_prefix
    chunk["norm_bam"] = norm_bam
    return {
        "chunk_id": cid,
        "region": chunk["region"],
        "log": log,
        "wall_s": round(time.time() - started, 3),
        "peak_tree_rss_kb": max(
            [s.get("peak_tree_rss_kb", 0) for s in steps] or [0]
        ),
        "steps": steps,
    }


def run_chunks_concurrently(args, ckpt, outdir, chunks, num_total_reads, rss_interval):
    """Run every chunk's stages 4-5, scheduling the flat unit queue from ONE budget.

    Concurrency is not a knob here any more. A chunk is one unit of the same flat queue
    LRAA schedules internally, so CpuBudget.allocate() derives both how many chunks run
    at once and the share each one's LRAA invocation is given. Their product cannot
    exceed the budget, which is what --concurrency and LRAA's own contig pool could
    previously do to each other.

    Launch order is longest-first on retained alignments per chunk, which the extractor
    already counted, so no extra pass is needed. Span would be the wrong proxy: it does
    not bound read count, and a short highly expressed chunk can outweigh a long quiet
    one.

    A chunk failure is fatal: the exception names the chunk and its log, and no
    merge is attempted. Chunks already running are allowed to finish so their
    logs are complete, but nothing new is started.
    """

    from concurrent.futures import ThreadPoolExecutor

    units = [
        CpuBudget.WorkUnit(
            contig_acc=chunk["chrom"],
            contig_strand=chunk["strand"],
            chunk_index=chunk["index"],
            region=(chunk["manifest"]["partition_lend"], chunk["manifest"]["partition_rend"]),
            cost=chunk["manifest"]["counts"]["alignments_emitted"],
        )
        for chunk in chunks
    ]
    allocation = CpuBudget.allocate(budget=args.cpu_budget, num_units=len(units))
    print(CpuBudget.format_allocation(allocation, phase="chunks"), flush=True)
    shortfall = CpuBudget.budget_shortfall_note(allocation)
    if shortfall:
        print(shortfall, flush=True)

    by_unit = dict(zip(units, chunks))
    launch_order = [by_unit[u] for u in CpuBudget.order_longest_first(units)]

    results = {}
    failures = []
    lock = threading.Lock()
    abort = threading.Event()

    def task(chunk):
        if abort.is_set():
            return
        try:
            record = chunk_worker(
                args,
                ckpt,
                outdir,
                chunk,
                num_total_reads,
                rss_interval,
                cpu_budget=allocation.tool_threads,
            )
            with lock:
                results[chunk["chunk_id"]] = record
        except Exception as err:  # noqa: BLE001 - reported, then re-raised below
            abort.set()
            with lock:
                failures.append((chunk["chunk_id"], chunk["log"], str(err)))

    started = time.time()
    with ThreadPoolExecutor(max_workers=allocation.unit_workers) as pool:
        list(pool.map(task, launch_order))
    makespan = time.time() - started

    if failures:
        lines = [
            "chunk {} FAILED -- log: {}\n  {}".format(cid, log, msg)
            for cid, log, msg in failures
        ]
        raise PipelineError(
            "{} of {} chunks failed; refusing to merge a partial result.\n{}".format(
                len(failures), len(chunks), "\n".join(lines)
            )
        )
    if len(results) != len(chunks):
        raise PipelineError(
            "only {} of {} chunks reported; refusing to merge a partial "
            "result".format(len(results), len(chunks))
        )
    ordered = [results[c["chunk_id"]] for c in chunks]
    return ordered, makespan, allocation


# --------------------------------------------------------------------- stage 6


def read_tsv(path):
    """Returns (comments, header_fields, rows). Comment lines start with '#'."""

    comments, header, rows = [], None, []
    with open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                comments.append(line)
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                continue
            rows.append(fields)
    if header is None:
        raise PipelineError("{} has no header line".format(path))
    return comments, header, rows


def merge_and_translate(outdir, chunks):
    """Stage 6. Concatenate per-chunk quant output, back in the whole-run frame.

    Three fields carry chunk-local COORDINATES: ``exons`` and ``introns``, and the
    splice hash code, which is a blake2s digest OF the introns string and so has
    to be recomputed rather than shifted.

    One field carries a chunk-local DENOMINATOR: ``TPM`` is
    ``all_reads / (total all_reads emitted by the same job) * 1e6``
    (pylib/Quantify.py:1527), so a chunk's value is relative to that chunk. This
    is not special pleading for the chunked arm: LRAA's own whole-contig merge
    rebases it the same way, recomputing TPM over every merged row from the
    printed ``all_reads`` (LRAA:_merge_quant_expr_files), because its per-job
    outputs have exactly the same problem at contig/strand granularity. Stage 6
    reproduces that arithmetic exactly -- same inputs, same expression, same
    format -- so equality of TPM follows from equality of ``all_reads``. The
    as-emitted chunk-local values are kept alongside, so nothing is hidden.

    ``RPM_total_reads`` needs no rebasing: its denominator is the ``-N`` passed
    identically to every chunk and to the control.

    ``mp_id`` in quant.tracking is a per-process MultiPath counter, not a
    coordinate; it is left as emitted and is not comparable across any two runs.

    Transcript ids come from the supplied GTF, so no id translation is needed.
    """

    merged_dir = os.path.join(outdir, "merged")
    os.makedirs(merged_dir, exist_ok=True)
    expr_out = os.path.join(merged_dir, "chunked.quant.expr")
    track_out = os.path.join(merged_dir, "chunked.quant.tracking")

    expr_header = None
    expr_rows = []
    hash_remap = {}  # (chunk_id, old_hash) -> new_hash
    for chunk in chunks:
        offset = chunk["offset"]
        _, header, rows = read_tsv(chunk["quant_prefix"] + ".quant.expr")
        if expr_header is None:
            expr_header = header
        elif header != expr_header:
            raise PipelineError(
                "chunk {} quant.expr header differs from the first chunk's".format(
                    chunk["chunk_id"]
                )
            )
        col = {name: i for i, name in enumerate(header)}
        for row in rows:
            if len(row) != len(header):
                raise PipelineError(
                    "chunk {} quant.expr row has {} fields, header has {}".format(
                        chunk["chunk_id"], len(row), len(header)
                    )
                )
            row = list(row)
            row[col["exons"]] = _shift_coord_string(row[col["exons"]], offset)
            introns = row[col["introns"]]
            if introns:
                introns = _shift_coord_string(introns, offset)
                row[col["introns"]] = introns
                new_hash = Util_funcs.get_hash_code(introns)
                hash_remap[(chunk["chunk_id"], row[col["splice_hash_code"]])] = new_hash
                row[col["splice_hash_code"]] = new_hash
            expr_rows.append((chunk["chunk_id"], row))

    ecol = {name: i for i, name in enumerate(expr_header)}
    expr_rows.sort(
        key=lambda item: (item[1][ecol["gene_id"]], item[1][ecol["transcript_id"]])
    )

    # TPM rebase, byte-for-byte the arithmetic of LRAA:_merge_quant_expr_files:
    # float() of the printed all_reads, plain sequential sum, /total*1e6, "%.3f".
    chunk_local_tpm = [row[ecol["TPM"]] for _, row in expr_rows]
    total_reported_read_count = 0.0
    for _, row in expr_rows:
        total_reported_read_count += float(row[ecol["all_reads"]] or 0)
    for _, row in expr_rows:
        counts = float(row[ecol["all_reads"]] or 0)
        tpm = (
            counts / total_reported_read_count * 1e6
            if total_reported_read_count > 0
            else 0
        )
        row[ecol["TPM"]] = "{:.3f}".format(tpm)

    with open(expr_out, "wt") as ofh:
        print("\t".join(expr_header), file=ofh)
        for _, row in expr_rows:
            print("\t".join(row), file=ofh)

    # the as-emitted, chunk-relative TPM, kept so the rebase above can be audited
    tpm_audit = expr_out + ".tpm_chunk_local.tsv"
    with open(tpm_audit, "wt") as ofh:
        print(
            "\t".join(
                ["chunk_id", "gene_id", "transcript_id", "all_reads",
                 "TPM_chunk_local", "TPM_merged_scope"]
            ),
            file=ofh,
        )
        for (chunk_id, row), local in zip(expr_rows, chunk_local_tpm):
            print(
                "\t".join(
                    [
                        chunk_id,
                        row[ecol["gene_id"]],
                        row[ecol["transcript_id"]],
                        row[ecol["all_reads"]],
                        local,
                        row[ecol["TPM"]],
                    ]
                ),
                file=ofh,
            )

    track_header = None
    track_rows = []
    for chunk in chunks:
        _, header, rows = read_tsv(chunk["quant_prefix"] + ".quant.tracking")
        if track_header is None:
            track_header = header
        elif header != track_header:
            raise PipelineError(
                "chunk {} quant.tracking header differs from the first "
                "chunk's".format(chunk["chunk_id"])
            )
        col = {name: i for i, name in enumerate(header)}
        for row in rows:
            row = list(row)
            old = row[col["transcript_splice_hash_code"]]
            row[col["transcript_splice_hash_code"]] = hash_remap.get(
                (chunk["chunk_id"], old), old
            )
            track_rows.append(row)

    tcol = {name: i for i, name in enumerate(track_header)}
    track_rows.sort(
        key=lambda r: (
            r[tcol["read_name"]],
            r[tcol["transcript_id"]],
            r[tcol["gene_id"]],
        )
    )
    with open(track_out, "wt") as ofh:
        print("\t".join(track_header), file=ofh)
        for row in track_rows:
            print("\t".join(row), file=ofh)

    return {
        "quant_expr": expr_out,
        "quant_tracking": track_out,
        "tpm_chunk_local_audit": tpm_audit,
        "expr_rows": len(expr_rows),
        "tracking_rows": len(track_rows),
        "splice_hash_codes_recomputed": len(hash_remap),
        "merged_scope_total_all_reads": total_reported_read_count,
    }


# -------------------------------------------------------------- baseline arm


def run_baseline(args, ckpt, outdir, timing, strand_bams, num_total_reads, rss_interval):
    """Whole-contig control. Same substrate, same settings, no partitioning.

    ``--bam`` is the two orientation bams glued back together, which is exactly
    the union of every chunk's mini bam expressed in genome coordinates (no read
    is dropped at a cut on this input), so the arms differ in partitioning alone
    rather than also in which records reach LRAA.
    """

    bdir = os.path.join(outdir, "baseline")
    os.makedirs(bdir, exist_ok=True)
    whole_bam = os.path.join(bdir, "whole.primary.bam")
    log = os.path.join(outdir, "logs", "baseline.log")
    started = time.time()
    steps = []

    merge_token = "baseline_merge.maxintron_{}".format(args.max_intron_length)
    merge_cmd = [
        "samtools",
        "merge",
        "-f",
        "--write-index",
        whole_bam,
        strand_bams["+"],
        strand_bams["-"],
    ]
    if ckpt.done(merge_token):
        steps.append({"step": "baseline_merge", "reused": True, "cmd": merge_cmd})
    else:
        steps.append(run_step("baseline_merge", merge_cmd, log, bdir, rss_interval))
        ckpt.mark(merge_token)

    norm_bam = os.path.join(bdir, "whole.norm.bam")
    norm_token = "baseline_norm.maxcov_{}_dw_{}_seed_{}_origin_0".format(
        args.normalize_max_cov_level, args.depth_window, args.random_seed
    )
    norm_cmd = [
        sys.executable,
        NORMALIZE_BAM,
        "--input_bam",
        whole_bam,
        "--output_bam",
        norm_bam,
        "--normalize_max_cov_level",
        str(args.normalize_max_cov_level),
        "--depth_window",
        str(args.depth_window),
        "--random_seed",
        str(args.random_seed),
        "--max_intron_length",
        str(args.max_intron_length),
        # THE control's grid must be pinned to absolute 0. The default anchors on
        # the first aligned base per contig and no chunk grid can match it.
        "--window_origin",
        "0",
    ]
    if ckpt.done(norm_token):
        steps.append({"step": "baseline_normalize", "reused": True, "cmd": norm_cmd})
    else:
        steps.append(run_step("baseline_normalize", norm_cmd, log, bdir, rss_interval))
        ckpt.mark(norm_token)

    # bare prefix, cwd=bdir: see the note in chunk_worker about LRAA's scratch
    # roots being string concatenations on --output_prefix
    quant_prefix = os.path.join(bdir, "baseline_quant." + LRAA_QUANT_ONLY_SUFFIX)
    quant_token = "baseline_quant.N_{}_hifi_{}".format(num_total_reads, args.HiFi)
    quant_cmd = lraa_cmd(
        args,
        bam_for_quant=whole_bam,
        bam_for_sg=norm_bam,
        genome=os.path.abspath(args.genome_fa),
        gtf=os.path.abspath(args.gtf),
        out_prefix="baseline_quant",
        num_total_reads=num_total_reads,
        # The control is the only thing running, so it gets the whole budget. Its own
        # flat queue is the two strands of one contig, and LRAA divides the budget across
        # those two units itself -- this is exactly the chromosome-shard shape.
        cpu_budget=args.cpu_budget,
    )
    if args.contig:
        quant_cmd += ["--contig", args.contig]
    if ckpt.done(quant_token):
        steps.append({"step": "baseline_quant", "reused": True, "cmd": quant_cmd})
    else:
        steps.append(run_step("baseline_quant", quant_cmd, log, bdir, rss_interval))
        ckpt.mark(quant_token)

    for suffix in (".quant.expr", ".quant.tracking"):
        if not os.path.exists(quant_prefix + suffix):
            raise PipelineError(
                "baseline produced no {}; log: {}".format(quant_prefix + suffix, log)
            )

    timing.setdefault("arms", {})["baseline"] = {
        "wall_s": round(time.time() - started, 3),
        "peak_tree_rss_kb": max([s.get("peak_tree_rss_kb", 0) for s in steps] or [0]),
        "steps": steps,
        "log": log,
    }
    return {
        "quant_expr": quant_prefix + ".quant.expr",
        "quant_tracking": quant_prefix + ".quant.tracking",
    }


# ---------------------------------------------------------------------- driver


def loadavg():
    with open("/proc/loadavg", "rt") as fh:
        return [float(x) for x in fh.read().split()[:3]]


def main(argv=None):

    parser = argparse.ArgumentParser(
        description="run the chunked quant-only pipeline (stages 1-6) or the "
        "whole-contig control, with per-chunk logs, wall time and peak RSS. "
        "Resumable via sentinel files under <output_dir>/__ckpt/.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bam", required=True, help="input bam of aligned reads")
    parser.add_argument("--genome_fa", required=True, help="genome fasta")
    parser.add_argument("--gtf", required=True, help="reference annotation gtf")
    parser.add_argument(
        "--output_dir", required=True, help="output directory; created if absent"
    )
    parser.add_argument(
        "--cpu_budget",
        type=int,
        default=None,
        help="TOTAL cores this pipeline may use. Chunk concurrency and each chunk's own "
        "LRAA budget are DERIVED from it, so the two cannot multiply past it. "
        "Default: all cores this process may run on",
    )
    parser.add_argument(
        "--arm",
        choices=("chunked", "baseline", "both"),
        default="chunked",
        help="chunked runs stages 1-6; baseline runs the whole-contig control "
        "on the same stage-1 substrate; both runs the chunked arm then the "
        "control, sequentially, so their timings do not contend",
    )
    parser.add_argument(
        "--num_total_reads",
        "-N",
        type=int,
        default=None,
        help="TPM denominator, passed identically to every chunk and to the "
        "control. Default: the stage-1 retained record count, which is what "
        "every arm actually sees. A supplied value that disagrees with that "
        "count is an error, not a preference",
    )
    parser.add_argument("--contig", default=None, help="restrict to one contig")
    parser.add_argument("--HiFi", action="store_true", help="pass --HiFi to LRAA")
    parser.add_argument("--approx_MB_per_cut", type=float, default=10)
    parser.add_argument("--approx_MB_per_cut_wiggle_window", type=float, default=1)
    parser.add_argument("--depth_window", type=int, default=100)
    parser.add_argument("--margin", type=int, default=200)
    parser.add_argument("--max_intron_length", type=int, default=200000)
    parser.add_argument("--normalize_max_cov_level", type=int, default=1000)
    parser.add_argument("--random_seed", type=int, default=42)
    parser.add_argument(
        "--rss_sample_interval",
        type=float,
        default=0.5,
        help="seconds between /proc RSS samples; a spike shorter than this can "
        "be missed",
    )

    args = parser.parse_args(argv)

    if args.cpu_budget is None:
        args.cpu_budget = CpuBudget.default_budget()
    elif args.cpu_budget < 1:
        parser.error("--cpu_budget must be >= 1")

    outdir = os.path.abspath(args.output_dir)
    os.makedirs(os.path.join(outdir, "logs"), exist_ok=True)
    ckpt = Checkpoints(os.path.join(outdir, "__ckpt"))

    timing_path = os.path.join(outdir, "timing.json")
    timing = {}
    if os.path.exists(timing_path):
        with open(timing_path, "rt") as fh:
            timing = json.load(fh)
    timing["invocation"] = sys.argv
    timing["resumable"] = True
    timing["checkpoint_dir"] = ckpt.root
    timing["cpu_budget"] = args.cpu_budget
    timing["arm"] = args.arm

    for tool in (SEPARATE_BAM, SELECT_CUTS, EXTRACT_CHUNK, NORMALIZE_BAM, LRAA):
        if not os.path.exists(tool):
            raise PipelineError("missing stage tool: {}".format(tool))
    if shutil.which("samtools") is None:
        raise PipelineError("samtools is not on PATH")

    rss = args.rss_sample_interval

    strand_bams = stage_strand_split(args, ckpt, outdir, timing, rss)
    retained = {s: count_records(b) for s, b in strand_bams.items()}
    retained_total = sum(retained.values())
    timing["stage1_retained_records"] = {
        "plus": retained["+"],
        "minus": retained["-"],
        "total": retained_total,
    }
    if args.num_total_reads is None:
        num_total_reads = retained_total
    else:
        num_total_reads = args.num_total_reads
        if num_total_reads != retained_total:
            raise PipelineError(
                "-N {} disagrees with the stage-1 retained record count {} "
                "({} on + plus {} on -). The denominator has to be the record "
                "set the arms actually see, or TPM is not comparable between "
                "them.".format(
                    num_total_reads, retained_total, retained["+"], retained["-"]
                )
            )
    timing["num_total_reads"] = num_total_reads

    outputs = {"num_total_reads": num_total_reads}

    def flush():
        with open(timing_path, "wt") as fh:
            json.dump(timing, fh, indent=2, sort_keys=True)
            print("", file=fh)

    flush()

    if args.arm in ("chunked", "both"):
        selections, cut_dir = stage_select_cuts(
            args, ckpt, outdir, timing, strand_bams, rss
        )
        flush()
        chunks = stage_extract_chunks(
            args, ckpt, outdir, timing, strand_bams, selections, rss
        )
        flush()
        print(
            "extracted {} chunks: {}".format(
                len(chunks), ", ".join(c["region"] for c in chunks)
            ),
            flush=True,
        )

        load_before = loadavg()
        arm_sampler = RssSampler(os.getpid(), rss)
        arm_sampler.start()
        try:
            chunk_records, makespan, chunk_allocation = run_chunks_concurrently(
                args, ckpt, outdir, chunks, num_total_reads, rss
            )
        finally:
            arm_sampler.stop()
        load_after = loadavg()

        merged = merge_and_translate(outdir, chunks)
        timing.setdefault("arms", {})["chunked"] = {
            "cpu_budget": chunk_allocation.budget,
            "concurrent_chunk_workers": chunk_allocation.unit_workers,
            "cpu_budget_per_chunk": chunk_allocation.tool_threads,
            "unallocated_cores": chunk_allocation.unallocated_cores,
            "makespan_s": round(makespan, 3),
            "summed_wall_s": round(sum(c["wall_s"] for c in chunk_records), 3),
            "chunks": chunk_records,
            "peak_rss_kb_summed_over_chunk_peaks": sum(
                c["peak_tree_rss_kb"] for c in chunk_records
            ),
            "observed_peak_concurrent_tree_rss_kb": arm_sampler.peak_kb,
            "loadavg_before": load_before,
            "loadavg_after": load_after,
        }
        timing["chunk_manifests"] = [
            {
                "chunk_id": c["chunk_id"],
                "region": c["region"],
                "offset": c["offset"],
                "window_origin": c["window_origin"],
                "alignments_emitted": c["manifest"]["counts"]["alignments_emitted"],
                "alignments_dropped_overhang": c["manifest"]["counts"][
                    "alignments_dropped_overhang"
                ],
                "gtf_transcripts_emitted": c["manifest"]["counts"][
                    "gtf_transcripts_emitted"
                ],
                "log": c["log"],
            }
            for c in chunks
        ]
        outputs["chunked"] = merged
        outputs["cut_dir"] = cut_dir
        flush()

    if args.arm in ("baseline", "both"):
        load_before = loadavg()
        outputs["baseline"] = run_baseline(
            args, ckpt, outdir, timing, strand_bams, num_total_reads, rss
        )
        timing["arms"]["baseline"]["loadavg_before"] = load_before
        timing["arms"]["baseline"]["loadavg_after"] = loadavg()
        flush()

    with open(os.path.join(outdir, "outputs.json"), "wt") as fh:
        json.dump(outputs, fh, indent=2, sort_keys=True)
        print("", file=fh)

    print(json.dumps(outputs, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except PipelineError as err:
        print("\nPIPELINE FAILED\n{}".format(err), file=sys.stderr)
        sys.exit(1)
