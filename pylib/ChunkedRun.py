#!/usr/bin/env python3

"""Run the chunked pipeline end to end, or the whole-contig control.

This is a thin driver over six tools that already exist and are already tested.
It adds no quantification or discovery logic of its own; it wires the stages
together, runs the per-chunk work concurrently, and records what each cost.

    stage 1  util/separate_bam_by_strand.py         orientation split
    stage 2  util/misc/select_contig_cut_points.py  cut selection
    stage 3  util/misc/extract_contig_region_inputs.py  chunk extraction
    stage 4  util/normalize_bam_by_strand.py        per-chunk normalization
    stage 5  LRAA, quant-only or discovery          per-chunk work
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

STRANDLESS CHUNKING, ``--strandless_chunks``. An added mode, off by default,
that moves the orientation split out of the whole-bam serial phase and into the
per-chunk parallel phase:

  strand-first  split(whole bam) -> cuts per contig-STRAND -> per chunk: normalize | quant
  strandless    cuts per contig  -> per chunk: extract | SPLIT | (normalize | quant) x2

Stage 1 is then not run at all for the chunked arm, which is the point: a
single-contig run no longer splits every record in the bam to use one contig's
worth. The two orientations of an interval share ONE extraction -- one mini
FASTA, one mini GTF, one pass over the region -- and the split runs on that
mini bam as the chunk's own first step. Stages 4 and 5 are untouched: they
still each receive one orientation-pure bam of one chunk.

ORDER IS LOAD-BEARING HERE. ``separate_bam_by_strand`` REWRITES
``read.is_reverse`` when the orientation it infers disagrees with the aligner,
and the extractor's strand filter reads the RAW flag. So a chunk is extracted
STRANDLESSLY and split AFTERWARDS; extracting ``chr1+:...`` from a bam that has
not been split would assign every flipped read to the wrong chunk and still
produce output that looks entirely normal. Three things stop that reordering:
the extractor refuses a strand-suffixed region over a mixed-orientation bam,
``assert_extracted_strandlessly`` refuses to split a chunk that was extracted
for one orientation, and the post-split record counts are checked against the
extractor's own per-orientation tallies.

DISCOVERY CHUNKING, ``--discovery``. Off by default, and quant-only is unchanged
in every particular while it is off. With it, stage 5 drops ``--quant_only``,
the annotation becomes OPTIONAL -- absent is de novo, present is ref-guided
discovery -- and stage 6 additionally merges the per-chunk GTFs, shifting
coordinates back into the whole-contig frame and namespacing model ids per unit.
The namespacing is not cosmetic: LRAA names a model after its contig and its
component index, every chunk's mini contig carries the SAME contig name, so
``t:chr1:+:comp-1:iso-1`` is emitted by every chunk that has a first component
and an unpatched concatenation fuses unrelated models into one record.

Stage 2 is IDENTICAL in the two modes. An earlier revision had discovery pass
``--require_zero_severed``, refusing any cut that severs a read; that contract
was REJECTED, because as libraries deepen every base ends up covered and the rule
would decline every cut -- silently disabling chunking on exactly the inputs that
need it. Severing is a cost the selector minimises, weighted so a spliced
alignment counts for ``--severed_multiexon_weight`` monoexonic ones, and a severed
read in discovery is treated as in quantification: dropped, counted and named.

Which makes the placement report load-bearing rather than decorative. Severing is
now EXPECTED, so the report is the only place anyone sees what the partition cost:
per contig-strand, targets requested, cuts placed, targets declined for the
annotation (the only remaining decline reason), and per cut how many alignments it
severs split monoexonic against multi-exon. Printed before the expensive phase and
stored in ``timing.json`` and ``outputs.json``. A partition that quietly shrank,
or quietly severed thousands of junctions, is a regression nobody can explain
afterwards.

MAKE-CHUNKS IS A POOL, not two serial loops. Stages 2 and 3 used to be the
serial half of a chunked run: one cut selection per orientation looping every
contig inside one process, then one extraction per chunk in a ``for``. Both scale
with READS rather than with genome length, which is why they dominate --
whole-genome strandless prep models at 1,009 s, and strand-first at 2,010 s,
against a processing phase that has been pooled since it was written.

``run_prep_concurrently`` schedules both as ONE flat queue of (contig,
orientation) selections and (contig, orientation, chunk) extractions, from the
same ``CpuBudget`` and with the same abort-and-refuse mechanics
``run_chunks_concurrently`` uses. Three things about it are not incidental:

  * NO STAGE BARRIER. A contig's extractions are submitted by that contig's own
    selection, so one contig's chunks start while another is still selecting.
    With a barrier the floor is the largest chromosome's whole prep -- chr1 at
    121.6 s, an 8.3x ceiling no number of cores improves. Without one the floor
    is the longest single dependency chain, one selection plus its longest
    extraction, which on an HG002 library is chrM at 73.6 s: 13.7x.
  * WHICH CONTIGS IS NOT A CHOICE. The enumeration is the selector's own
    ``sorted(lengths)`` off the genome fasta, so the cuts a pooled run places are
    the cuts a serial run places. See ``enumerate_prep_contigs``.
  * ORDER IS RESTORED SERIALLY. ``selections[key]`` is reassembled in contig
    order and the ``order`` counter every quant unit carries is assigned in one
    serial pass afterwards, so the merged table's row order cannot depend on
    which selection finished first. See ``planned_chunks``.

A chunk covering its WHOLE contig is not extracted. Offset 0, mini contig named
and sized like the real one, nothing able to overhang: the mini bam would restate
the source. Measured on chrM of a 6.9 GB HG002 bam, writing it costs 73.2 s
against 6.6 s for the same fetch and the same predicates without the write, so
the chunk's bam is the source and stage 3b reads it restricted to the contig.
Strandless only -- a strand-first chunk bam is normalized and quantified
directly, and neither of those can restrict a multi-contig source.

MEASUREMENT NOTES. Wall time is measured around each subprocess. Peak RSS is
sampled from ``/proc`` at ``--rss_sample_interval`` over the step's whole
process tree and summed across that tree, so a chunk's figure already includes
whatever workers its LRAA run spawned; the arm-level figure is the same sum over
this driver's entire descendant tree, which is the real concurrent footprint.
Sampling can miss a spike shorter than the interval.
"""

import argparse
import collections
import csv
import glob
import gzip
import hashlib
import heapq
import json
import operator
import os
import re
import resource
import shutil
import subprocess
import sys
import tempfile
import threading
import time

import pysam

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
import CpuBudget  # noqa: E402  (path insert must precede the import)
import RdnaMask  # noqa: E402
import Util_funcs  # noqa: E402
import LRAA_Globals  # noqa: E402

REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")
)
SEPARATE_BAM = os.path.join(REPO_ROOT, "util", "separate_bam_by_strand.py")
SELECT_CUTS = os.path.join(REPO_ROOT, "util", "misc", "select_contig_cut_points.py")
EXTRACT_CHUNK = os.path.join(
    REPO_ROOT, "util", "misc", "extract_contig_region_inputs.py"
)
NORMALIZE_BAM = os.path.join(REPO_ROOT, "util", "normalize_bam_by_strand.py")
LRAA = os.path.join(REPO_ROOT, "LRAA")


# Set in the environment of everything this pipeline launches. `LRAA --chunk`
# refuses to start when it is present, so a chunk worker cannot re-enter the
# pipeline and fork without bound.
WORKER_ENV = "LRAA_CHUNKED_PIPELINE_WORKER"
STRAND_TAG = {"+": "plus", "-": "minus"}
# The orientation a strandless cut selection or chunk is for is the empty string,
# the same "both" convention parse_region, _strand_matches and find_islands
# already use. STRANDLESS_TAG is what that key is called in a filename.
STRANDLESS_TAG = "strandless"
STRAND_FIRST_MODE = "strand_first"
STRANDLESS_MODE = "strandless"
PIPELINE_MODES = (STRAND_FIRST_MODE, STRANDLESS_MODE)
PAGE_KB = os.sysconf("SC_PAGE_SIZE") // 1024
# LRAA rewrites --output_prefix to "<prefix>.<mode suffix>"
# (LRAA:_get_lraa_output_mode_suffix), so the files land under that name and the
# mode decides which name.
LRAA_QUANT_ONLY_SUFFIX = "LRAA.quant-only"
LRAA_REF_GUIDED_SUFFIX = "LRAA.ref-guided"
LRAA_REF_FREE_SUFFIX = "LRAA.ref-free"


def lraa_output_suffix(discovery, gtf):
    """What LRAA appends to --output_prefix, mirroring its own three-way choice.

    Restated here rather than imported: LRAA is a script, not a module, and
    importing it would run its argument parser. pylib/test_chunked_entry_point.py
    asserts these three against LRAA's own, so a rename there fails a test rather
    than quietly producing paths that nothing writes.
    """

    if not discovery:
        return LRAA_QUANT_ONLY_SUFFIX
    return LRAA_REF_GUIDED_SUFFIX if gtf else LRAA_REF_FREE_SUFFIX

# LRAA gzips quant.tracking since v0.20.0; v0.19.x and earlier wrote it plain.
# Both are accepted on read so a run at either version resolves its own artifact,
# and the .gz form is preferred when a directory somehow carries both.
QUANT_TRACKING_SUFFIXES = (".quant.tracking.gz", ".quant.tracking")


def gtf_index_cache_dir(outdir):
    """Where a GTF's tabix index goes when its own directory is read-only.

    Under the run directory rather than a temp dir, so that a resumed run finds
    the index its earlier stages built instead of rebuilding it.
    """

    return os.path.join(outdir, "gtf_index")


def resolve_tracking(quant_prefix):
    """Path to <quant_prefix>'s quant.tracking, gzipped or not. None if absent."""

    for suffix in QUANT_TRACKING_SUFFIXES:
        path = quant_prefix + suffix
        if os.path.exists(path):
            return path
    return None


def open_text(path):
    """Open a text file that may be gzipped, decided by extension."""

    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


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
            # os.wait4 rather than proc.wait(): rusage is exact and cumulative, and
            # it includes descendants the child reaped. LRAA joins its per-unit
            # workers, so the whole tree rolls up here. The RssSampler below cannot
            # substitute -- at a 0.5 s interval it misses exactly the short-lived
            # per-unit workers this instrumentation exists to account for.
            #
            # wait4 reaps the pid itself, so assign proc.returncode directly and do
            # NOT call proc.wait() afterwards, or Popen.__del__ warns about an
            # unreaped child.
            _, status, ru = os.wait4(proc.pid, 0)
            proc.returncode = returncode = os.waitstatus_to_exitcode(status)
            cpu_s = ru.ru_utime + ru.ru_stime
            max_rss_kb = ru.ru_maxrss
        finally:
            sampler.stop()
        elapsed = time.time() - started
        print(
            "\n===== step {} exit {} wall {:.1f}s cpu {:.1f}s "
            "peak_tree_rss {} KB ({} samples) max_rss {} KB =====".format(
                name,
                returncode,
                elapsed,
                cpu_s,
                sampler.peak_kb,
                sampler.samples,
                max_rss_kb,
            ),
            file=log,
            flush=True,
        )

    record = {
        "step": name,
        "cmd": cmd,
        "log": log_path,
        "wall_s": round(elapsed, 3),
        "cpu_s": round(cpu_s, 3),
        # Two different quantities, deliberately both kept: peak_tree_rss_kb is a
        # SAMPLED sum over the process tree, max_rss_kb is an EXACT peak for the
        # largest single process in it. Neither replaces the other, and the exact
        # one is what distinguishes a completed run from an OOM-killed one.
        "peak_tree_rss_kb": sampler.peak_kb,
        "max_rss_kb": max_rss_kb,
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


def chain_token(local, parent, *opaque):
    """A sentinel name that moves when this step's inputs or any upstream step's do.

    Sentinels are named for the stage and for the parameters that change that
    stage's output, which is only sound while every such parameter is re-listed at
    every stage it reaches.  It was not.  ``baseline_norm`` named the coverage
    target, window, seed and origin but not ``--max_intron_length``, which it hands
    the normalizer, so changing the cap reran the merge into the same
    ``whole.primary.bam`` and then reused a normalized bam built under the old one.
    Neither quant stage named anything at all about the bams it reads, so the same
    change left both quant results untouched as well.

    Re-listing is the part that drifts.  A step now names its own inputs and
    inherits its parent's token rather than repeating what the parent already
    covers, so adding an input upstream invalidates everything downstream without
    an edit to any of those stages.  ``local`` stays in the filename, because a
    sentinel directory that cannot be read at a glance is its own hazard; the
    digest is what carries the chain.  ``opaque`` is for inputs that determine
    contents but cannot go in a filename, such as resolved paths and stat pairs.
    """
    payload = "|".join([parent or "root", local] + [str(f) for f in opaque])
    return "{}.up_{}".format(local, Util_funcs.get_hash_code(payload)[:12])


def argv_digest(argv):
    """A stable digest of the COMMAND a step runs, order-normalised over its flags.

    A sentinel named for a hand-listed set of settings goes stale the moment a
    caller forwards a new flag: the list is maintained by hand, the new flag is
    not in it, and a resumed run serves the old mode's output as the new mode's.
    ``stage5_quant`` named the unit, the read denominator and ``--HiFi`` and so
    named neither quant-only versus discovery nor anything a later forward would
    add.  The argv a step actually issues is the one description of that step
    that cannot fall behind what it does, so it is what the sentinel is keyed on:
    adding, removing or changing any argument moves the digest with no edit here.

    Flags are paired with their values BEFORE sorting, so reordering
    ``--a 1 --b 2`` is inert while swapping to ``--a 2 --b 1`` is not -- sorting
    the bare token stream would lose that distinction.  Leading positionals
    (interpreter, script) keep their relative order, because their position is
    what makes them arguments of each other rather than of a flag.

    Every argument reaches the digest.  A field that does NOT determine the
    step's output is gated by the caller, visibly and one at a time -- see
    ``quant_command_digest`` -- rather than filtered out by a list here: a
    signature naming a non-determining input is as defective as one missing a
    field, because it makes the token unauditable, and that is how missing
    fields survive.
    """
    positionals = []
    flags = []
    i = 0
    n = len(argv)
    while i < n:
        tok = str(argv[i])
        if tok.startswith("--"):
            # A value is any following token that is not itself a long flag; that
            # keeps negative numbers attached to the flag they belong to.
            if i + 1 < n and not str(argv[i + 1]).startswith("--"):
                flags.append("{}\x1f{}".format(tok, argv[i + 1]))
                i += 2
                continue
            flags.append(tok)
        else:
            positionals.append(tok)
        i += 1
    payload = "\x1e".join(["argv0\x1f" + "\x1f".join(positionals)] + sorted(flags))
    return Util_funcs.get_hash_code(payload)


def quant_command_digest(cmd):
    """``argv_digest`` of an LRAA quant command, with ``--cpu_budget`` gated inert.

    ``--cpu_budget`` cannot change what a quant step produces.  It is this
    invocation's SHARE of the total, and a chunk is single-contig and
    strand-pure, so LRAA's internal queue inside it is one unit and its own pool
    clamps to 1 regardless (see this module's docstring, ONE CPU BUDGET).  The
    baseline control divides its share across the two strands of one contig by
    the same allocator.  Neither changes a count.

    Naming it anyway is not the safe direction.  It breaks the recovery path this
    pipeline most needs to keep working: a whole-genome run OOM-killed at
    ``--cpu_budget 16`` and resumed at 8 would discard every finished stage-5
    sentinel and redo units that are already correct on disk.

    The VALUE is gated, not the flag, following the ``onone`` convention the
    normalize cache stem uses for an unset origin: the flag still appears in the
    payload, so dropping ``--cpu_budget`` from ``lraa_cmd`` altogether still moves
    the digest and remains auditable.  Exactly one field is gated, on the stated
    reason above; anything else added here needs its own.  ``--num_total_reads``
    is the TPM denominator and ``--output_prefix`` decides which files exist, so
    neither qualifies, and a test holds every other argument of the real command
    to moving the token.
    """
    gated = []
    i = 0
    n = len(cmd)
    while i < n:
        tok = str(cmd[i])
        gated.append(tok)
        if tok == "--cpu_budget" and i + 1 < n and not str(cmd[i + 1]).startswith("--"):
            gated.append("inert")
            i += 2
            continue
        i += 1
    return argv_digest(gated)


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


def coords_already_whole_contig(units):
    """True when nothing needs shifting, because no unit is a piece of a contig.

    A unit's ``offset`` is where its mini contig begins inside the real one, so
    offset 0 means the mini contig IS the contig: chunk-local coordinates ARE
    whole-contig coordinates, and every shift the merge would apply is ``+ 0``.

    All-zero over the whole unit list is the one-chunk-per-contig whole-genome
    mode, where no contig is cut -- and it is not a corner case even when
    contigs are cut: 171 of the 475 chunks of the shipped HG002 partition
    already spanned their whole contig.

    DERIVED FROM THE UNITS, never configured. There is no flag: a caller can
    neither ask for the shift to be skipped when it would move a coordinate, nor
    ask for it to be applied when it provably would not.
    """

    return all(unit["offset"] == 0 for unit in units)


# --------------------------------------------------------------------- stage 1


def stage_strand_split(args, ckpt, outdir, timing, rss_interval, inputs_token):
    split_dir = os.path.join(outdir, "strand_split")
    os.makedirs(split_dir, exist_ok=True)
    prefix = os.path.join(split_dir, "input")
    log = os.path.join(outdir, "logs", "stage1_strand_split.log")
    token = chain_token(
        "stage1_strand_split.maxintron_{}".format(args.max_intron_length),
        inputs_token,
    )
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
    return strand_bams, token


def count_records(bam):
    out = subprocess.check_output(["samtools", "view", "-c", bam], text=True)
    return int(out.strip())


def ensure_bam_index(bam):
    """A strandless run fetches the RAW bam by region, so it needs its index.

    Strand-first never did.  Its cut selection and extraction read the stage-1
    output, which is written indexed, and nothing else in the pipeline touches the
    input by coordinate -- so an unindexed input has always been fine, and the
    corpora that exercise chunking all ship a committed .bai, which is why nothing
    noticed.  ``LRAA --chunk`` makes it worse: it dispatches before LRAA's own
    "missing index, will try to make it" step, so the one thing that would have
    built the index is bypassed.  Observed as
    ``ValueError: fetch called on bamfile without index`` raised out of pysam
    inside stage 2, naming neither the file nor the reason.

    Built beside the bam, where every reader looks for it.  An input directory that
    cannot be written to is refused with that said, rather than left to surface as
    the same pysam error one stage later.
    """

    for suffix in (".bai", ".csi"):
        if os.path.exists(bam + suffix):
            return
    print("missing index for {}, building it".format(bam), flush=True)
    try:
        pysam.index(bam)
    except Exception as err:
        raise PipelineError(
            "--strandless_chunks needs a coordinate index for {0}, which has "
            "none, and building one beside it failed: {1}. Strandless cut "
            "selection and extraction fetch the raw bam by region. Index it "
            "where it lives (samtools index {0}), or stage the bam somewhere "
            "writable first.".format(bam, err)
        )


# ----------------------------------------------------- make-chunks scheduling


def contig_mapped_counts(bam):
    """Mapped record count per reference, read out of the bam's own index.

    ``samtools idxstats`` touches only the ``.bai``, so this is 0.19 s on a
    6.9 GB bam rather than a pass over it. It is used ONLY as the launch-order
    cost proxy -- ``enumerate_prep_contigs`` says why it must not choose the
    contig set -- so a reference it cannot report just sorts last.
    """

    counts = {}
    try:
        report = pysam.idxstats(bam)
    except Exception as err:  # noqa: BLE001 - a launch-order proxy, not a result
        print(
            "idxstats on {} failed ({}), so make-chunks will launch in "
            "contig-length order rather than longest-first. Scheduling only; no "
            "output depends on it.".format(bam, err),
            flush=True,
        )
        return counts
    for line in report.splitlines():
        fields = line.split("\t")
        if len(fields) < 3 or fields[0] == "*":
            continue
        counts[fields[0]] = int(fields[2])
    return counts


def enumerate_prep_contigs(args):
    """The contigs cut selection runs on, and their lengths. NOT a choice.

    ``select_contig_cut_points.py`` enumerates ``sorted(lengths)`` off the genome
    fasta when it is given one and off the bam header otherwise (``:1335-1351``),
    and this pipeline always gives it the fasta. Splitting that loop across
    processes must not also change which contigs are in it:

      * a SMALLER set drops the reads AND the TPM denominator of everything it
        omits, silently. "Main chromosomes" on a whole-genome HG002 library omits
        85 references holding 504,971 mapped records, 2.78 % of the library.
        "References with mapped reads" is subtler and still wrong: it changes the
        cut manifest against a serial run, which is the artifact the parity gates
        diff.
      * a LARGER set invents chunks.

    So the set is the selector's own, verbatim, and idxstats supplies the launch
    cost only. A reference holding nothing still gets its selection and its one
    empty chunk, exactly as the serial loop gave it one.

    ONE exception, not a contradiction of the rule above: a reference ABSENT from
    the bam header entirely cannot hold anything, by construction -- a bam cannot
    carry a record against a reference its own header does not list. Filtering
    those out is not the heuristic narrowing the rule above forbids (which risks
    dropping a reference that DOES have reads); it can never drop reads, because
    there is nothing there to drop. ``select_contig_cut_points.py`` already
    filters exactly this set for its own no-``--contig`` enumeration
    (``:1761-1791``, following ``00bf5b9``), but its ``--contig X`` branch --
    which is what per-contig seeding below actually uses, for every chunked run
    -- still raised ``SelectionError`` on such an X: by the time that branch
    checks, X has already been filtered out of the very ``lengths`` dict it is
    tested against. GRCh38_no_alt.fa's chrEBV (195 fasta sequences, 194 in a real
    bam) reproduced this exactly: 26 cut selections seeded, 1 failed, the whole
    make-chunks phase refused to extract or merge a partial result. Filtering
    here -- before a subprocess is ever seeded for such a contig -- is symmetric
    with what the no-``--contig`` path already does, and does not touch the
    "exactly one selection per invocation" contract downstream in
    ``run_cut_selection``.

    ``--contig`` reduces the enumeration to a single entry, which leaves today's
    behaviour as a degenerate case of the same pool rather than a second path.
    ``--contigs`` reduces it to a NAMED SET, which is how a workflow's
    "main chromosomes" reaches a chunk partition. Both are applied AFTER the
    header filter above, so a name present in the fasta but absent from the bam
    is reported by that filter rather than mistaken for a typo here.
    """

    with pysam.FastaFile(os.path.abspath(args.genome_fa)) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
    if not lengths:
        raise PipelineError(
            "{} names no contigs, so there is nothing to select cuts on. The "
            "make-chunks phase enumerates the genome fasta because that is what "
            "cut selection enumerates.".format(args.genome_fa)
        )
    with pysam.AlignmentFile(os.path.abspath(args.bam), "rb") as bam:
        bam_contigs = set(bam.references)
    absent = sorted(c for c in lengths if c not in bam_contigs)
    if absent:
        print(
            "NOTE: {} of {} genome reference sequence(s) are absent from the bam "
            "header and are skipped, having nothing to select cuts on: {}{}".format(
                len(absent),
                len(lengths),
                ", ".join(absent[:5]),
                "" if len(absent) <= 5 else ", ...",
            ),
            file=sys.stderr,
        )
        lengths = {c: n for c, n in lengths.items() if c in bam_contigs}
        if not lengths:
            raise PipelineError(
                "every one of {}'s {} contig(s) is absent from {}'s header; "
                "nothing to select cuts on.".format(
                    args.genome_fa, len(absent), args.bam
                )
            )
    if args.contig and getattr(args, "contigs", None):
        raise PipelineError(
            "--contig names one contig and --contigs names a set; pick one"
        )
    if getattr(args, "contigs", None):
        wanted = [c for c in (n.strip() for n in args.contigs.split(",")) if c]
        if not wanted:
            raise PipelineError("--contigs is empty")
        # A name that is in NEITHER the fasta nor the bam is a typo, and a typo
        # must fail rather than silently shrink the run: the partition, the reads
        # and the TPM denominator would all quietly lose whatever it meant to
        # name. A name filtered out by the header check above is reported there
        # and is not an error here -- it holds nothing by construction.
        unknown = [c for c in wanted if c not in lengths and c not in absent]
        if unknown:
            raise PipelineError(
                "--contigs names {} contig(s) absent from both {} and {}'s "
                "header: {}. Known contigs begin: {}".format(
                    len(unknown),
                    args.genome_fa,
                    args.bam,
                    ", ".join(unknown),
                    ", ".join(sorted(lengths)[:10]),
                )
            )
        kept = [c for c in wanted if c in lengths]
        if not kept:
            raise PipelineError(
                "every contig --contigs names is absent from the bam header, so "
                "there is nothing to select cuts on: {}".format(", ".join(wanted))
            )
        lengths = {c: n for c, n in lengths.items() if c in kept}
        return sorted(kept), lengths
    if args.contig:
        if args.contig not in lengths:
            raise PipelineError(
                "--contig {} is absent from {}{}. Known contigs begin: {}".format(
                    args.contig,
                    args.genome_fa,
                    " or from the bam" if args.contig in absent else "",
                    ", ".join(sorted(lengths)[:10]),
                )
            )
        return [args.contig], lengths
    return sorted(lengths), lengths


def _extractor_module():
    """``extract_contig_region_inputs``, imported into this process.

    For ``ensure_gtf_index``, and for reusing ``_attribute`` (the extractor's
    GTF attribute parser) in ``split_chunk_gtf_by_strand`` -- a pure,
    side-effect-free helper, not the extraction work itself. Everything that
    actually EXTRACTS this pipeline still asks of the extractor through a
    subprocess, so that a per-chunk failure lands in that chunk's own log
    rather than killing the driver; that isolation is about extraction's own
    failure modes (a malformed region, a missing index), not a regex search.
    """

    path = os.path.join(REPO_ROOT, "util", "misc")
    if path not in sys.path:
        sys.path.insert(0, path)
    import extract_contig_region_inputs

    return extract_contig_region_inputs


def _selector_module():
    """``select_contig_cut_points``, imported into this process.

    For ``spanning_read_names`` alone -- the selector's own answer to "which
    reads does a cut at this position sever", reused rather than reimplemented so
    that ``verify_severed_accounting`` keeps checking the extractor against a
    SECOND implementation of the geometry. Rewriting the span test here would
    turn that check into the extractor agreeing with itself.

    Selection itself still goes through a subprocess, for the same isolation
    reason ``_extractor_module`` gives.
    """

    path = os.path.join(REPO_ROOT, "util", "misc")
    if path not in sys.path:
        sys.path.insert(0, path)
    import select_contig_cut_points

    return select_contig_cut_points


def warm_gtf_index(args, outdir, timing):
    """Build the shared GTF tabix index ONCE, serially, before the pool opens.

    ``ensure_gtf_index`` has no lock: it checks a stamp, misses, and builds.
    MEASURED, with ``_build_gtf_index`` marked per pid -- six concurrent callers
    against a COLD cache produced SIX independent builds; against a WARM cache,
    zero. Nothing is corrupted, because each build writes pid-scoped temps and
    ``os.replace``s them into place, but on a 1.49 GB GENCODE GTF one build is
    36.2 s and spills ~1.5 GB of temp, so sixteen-wide cold is sixteen times
    both. This is the single largest hazard the parallel phase introduces, and
    one serial call removes it.

    Reported rather than silent, so a cold run says what it paid.
    """

    extractor = _extractor_module()
    started = time.time()
    index_path = extractor.ensure_gtf_index(
        os.path.abspath(args.gtf), cache_dir=gtf_index_cache_dir(outdir)
    )
    elapsed = round(time.time() - started, 3)
    record = {"wall_s": elapsed, "index": index_path}
    if index_path is None:
        # Not fatal: the selector and the extractor both fall back to scanning
        # the GTF per contig. Said out loud, because that fallback is the slow
        # path and silence would read as a successful pre-warm.
        record["note"] = (
            "no tabix index could be built for {}; every unit will scan the "
            "whole GTF instead of seeking into it".format(args.gtf)
        )
    print(
        "gtf tabix index ready in {:.2f}s: {}".format(
            elapsed, index_path or "NONE -- every unit will scan the gtf"
        ),
        flush=True,
    )
    timing.setdefault("stages", {})["gtf_index_warm"] = record
    return index_path


# Worst per-unit peak RSS across the make-chunks stages, with margin. MEASURED by
# exact wait4 ru_maxrss, not by the sampler, over a WHOLE-GENOME run of this pool:
# 195 contigs of GRCh38 against GENCODE v39, every unit's peak recorded by
# ``run_step``.
#
#   cut selection      43.0 MiB (chrM) - 263.8 MiB (chr1), median 44.5.
#                      chr1 264, chr2 220, chr3 191, chr11 179, chr12 166.
#                      The dominant stage, and it scales with ANNOTATION size,
#                      not read count: chrM carries 1.19 M reads, 6.6 % of the
#                      library, and has the SMALLEST footprint of the 195.
#   chunk extraction   41.5 - 80.5 MiB, and flat in concurrency
#                      (80.5 / 80.6 / 80.3 MiB at 1 / 4 / 8 wide).
#
# So K workers cost about K times this. Measured concurrent tree RSS over the same
# run: 0.30 GiB at 1 worker, 0.83 GiB at 4 -- 181 MiB marginal per worker, under
# the worst case because the queue mixes cheap contigs with expensive ones.
#
# 300 MiB is the measured worst (263.8) plus 14 %. It is a GUARD, not an
# estimator: the MiB-per-GTF-MB ratio varies by 2.0x across contigs (1.49 on
# chr1, 0.74 on chr21), which is not enough to model from, so nothing here
# pretends to model it.
PREP_UNIT_PEAK_MIB = 300


def available_memory_mib():
    """``MemAvailable`` in MiB, or None where /proc/meminfo cannot be read."""

    try:
        with open("/proc/meminfo", "rt") as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) // 1024
    except (OSError, ValueError, IndexError):
        return None
    return None


def prep_memory_cap(budget, available_mib=None):
    """Concurrency this box can hold in the make-chunks phase. ``(cap, note)``.

    Bounded by MEMORY as well as by cores, because the two bounds belong to
    different parties: the core bound is the user's ``--cpu_budget``, the memory
    bound is the box's. Both have to hold, and only one of them was being
    checked.

    What the numbers are, at ``PREP_UNIT_PEAK_MIB``:

      | workers | worst-case prep RSS | 62 GB box | 16 GB box |
      |       8 |            2.34 GiB |       26x |      6.8x |
      |      16 |            4.69 GiB |       13x |      3.4x |
      |      24 |            7.03 GiB |      8.8x |      2.3x |
      |      32 |            9.38 GiB |      6.6x | 1.7x, capped |

    So a 62 GB box runs any width the budget asks for and a 16 GB box runs to
    24-wide; the cap starts doing work on a small box at high budget, and on a
    container whose memory has been squeezed. The processing phase's own measured
    peak is 3.575 GB on a whole chr1, and the two phases do not overlap -- prep
    completes before processing starts -- so a wide make-chunks pool becomes the
    run's high-water MAX rather than adding to it.

    ``None`` cap means unbounded by memory. The cap never falls below 1: running
    one unit at a time is what the serial loop being replaced did, and is
    strictly better than refusing to run.
    """

    if available_mib is None:
        available_mib = available_memory_mib()
    if available_mib is None:
        return None, (
            "make-chunks concurrency is bounded by --cpu_budget alone: "
            "/proc/meminfo could not be read, so the {} MiB-per-unit memory "
            "bound could not be applied on this box.".format(PREP_UNIT_PEAK_MIB)
        )
    affordable = max(1, int(available_mib) // PREP_UNIT_PEAK_MIB)
    if affordable >= budget:
        return None, None
    return affordable, (
        "make-chunks concurrency capped at {} worker(s) by MEMORY rather than by "
        "--cpu_budget {}: {} MiB available, and one make-chunks unit peaks at up "
        "to {} MiB (measured worst case: cut selection against a whole-genome "
        "GTF). {} core(s) of the budget sit idle in this phase deliberately, "
        "because {} x {} MiB is {} MiB and the box does not have it.".format(
            affordable,
            budget,
            available_mib,
            PREP_UNIT_PEAK_MIB,
            budget - affordable,
            budget,
            PREP_UNIT_PEAK_MIB,
            budget * PREP_UNIT_PEAK_MIB,
        )
    )


# Per-CHUNK peak RSS in the chunk-processing phase, as a linear function of the
# chunk's own genomic SPAN. Unlike the make-chunks phase, this one cannot be
# guarded by a single constant: measured per-chunk peak runs from 50 MiB on an
# empty scaffold to 5,594.7 MiB on a whole chr1, a 112x range inside one run, so a
# constant is either useless at whole-chromosome span or destroys the 475-chunk
# default's concurrency. The cost has to be charged per unit, from something known
# before the unit starts.
#
# SPAN, not read count, and the two disagree sharply. MEASURED over the 195-contig
# whole-genome run (per-contig quant peaks reported by MeasurePerContig from the
# ``stage5_quant_strandless`` step footers in ``chunk_work/logs/chunk_*.log``) plus
# 1,480 per-chunk records from the 475-chunk HG002 PacBio and ONT runs, the chr1
# 50-chunk sweep, and the chr21 discovery runs:
#
#   MiB per genomic Mb, chunks >= 5 Mb: p50 32.2, p95 51.9, max 74.0 -- 10.6x range
#   MiB per million alignments:         p50 12,552, max 31,066     -- 56x range
#
# A span fit over the 110 contigs carrying reads reaches R2 0.970 against 0.850 for
# a read-count fit. The reason span wins is that memory tracks the splice graph and
# the multipath structures built over the unit's interval, not the record count:
# chrM carries 1.19 M reads in 16.6 kb and peaks at 655 MiB, the SMALLEST of the
# large-contig set, while chr4 carries 366 k reads over 190 Mb and peaks at 3,110.
#
# The two constants are an upper ENVELOPE, not that fit. A least-squares span fit
# (88 + 19.6 MiB/Mb) predicts 88 MiB for chrM and is 7.4x UNDER its measured 655;
# the 1 GiB fixed term is what covers that second regime, which is why no read term
# is needed. Envelope verified against all 1,487 measured points: zero exceed it.
# Tightest margins are the ones that matter -- whole chr1 1.16x (6,501 estimated
# against 5,594.7 measured), whole chr21 in DISCOVERY mode 1.20x (2,051 against
# 1,710). It is loosest where it is harmless: 2.3x on a 20 Mb quant chunk, up to
# 20x on an empty scaffold, where the fixed term dominates.
#
# How wrong it can be, said plainly: it over-charges a quant-only 10 Mb chunk by
# about 1.7x, so the cap starts binding on the 475-chunk default at roughly 22 GiB
# available rather than the ~9 GiB that partition actually needs at 16-wide. That
# is the price of an envelope that also holds at whole-chromosome span and in
# discovery mode. It is calibrated on human/GRCh38 with GENCODE v39 at PacBio and
# ONT depth; a corpus with much deeper coverage per Mb would need it re-measured,
# and the run report carries the estimate beside the observed peak so that
# re-measurement is a subtraction rather than a new experiment.
CHUNK_UNIT_FIXED_MIB = 1024
CHUNK_UNIT_MIB_PER_GENOMIC_MB = 22


ChunkMemoryBound = collections.namedtuple(
    "ChunkMemoryBound",
    "cap note available_mib largest_unit_mib charged_mib",
)


def chunk_unit_peak_mib(span_bp):
    """Guard on one chunk's peak RSS, from its genomic span. See the constants."""

    span_mb = max(0, int(span_bp)) / 1e6
    return CHUNK_UNIT_FIXED_MIB + int(CHUNK_UNIT_MIB_PER_GENOMIC_MB * span_mb)


def chunk_memory_cap(
    budget, unit_spans_bp, available_mib=None, strand_concurrency_multiplier=1
):
    """Concurrency this box can hold in the chunk-processing phase.

    Same mechanism as ``prep_memory_cap`` -- ``MemAvailable`` against a measured
    per-unit guard, handed to ``CpuBudget.allocate`` as ``max_unit_workers`` -- and
    the same reason: the core bound is the user's ``--cpu_budget``, the memory bound
    is the box's, and only the first was being checked here. ``min(cpu_budget,
    units)`` put 16 whole-chromosome quant units resident at once and the kernel OOM
    killer took the chr1+ worker 3,825 s into a 195-chunk whole-genome run, at a
    process-tree peak of 45.74 GiB on a 62 GiB box with 358 of 390 units still to go.

    What differs from prep is that the per-unit cost is not one number: it is
    ``chunk_unit_peak_mib`` of each unit's own span. So the cap is the largest K for
    which the K MOST EXPENSIVE units fit together. Charging the K largest rather
    than the K that happen to launch first is deliberate: it is an upper bound over
    every K-subset of the queue, which makes it independent of the launch order --
    and the launch order here is longest-first on ALIGNMENTS, a different ranking
    from the span-based estimate, so a bound that assumed the two agreed would not
    hold.

    REDUCE, not refuse, and the cap never falls below 1 -- the same choice as
    ``prep_memory_cap:766``. Refusing a partition whose single largest unit does not
    fit was considered and rejected: the estimate is an envelope carrying up to 2.3x
    of deliberate margin, so refusal would kill runs that would have completed, and
    the OOM this fixes came from 16-way RESIDENCY rather than from one unit -- that
    unit's 4.59 GiB fits a 62 GiB box 13 times over. What a refusal would buy is
    honesty about the case concurrency cannot fix, and that is bought instead by
    ``note`` saying so explicitly when even one unit is over the line.

    A cap of ``None`` means unbounded by memory. ``largest_unit_mib`` and
    ``charged_mib`` are returned whether or not the cap binds, so the run report
    carries the estimate that was compared against the box even on the runs where it
    changed nothing.

    ``strand_concurrency_multiplier`` charges each unit's estimate as if that many
    quant subprocesses could be resident inside ONE chunk-worker slot at once, not
    just one. Pass 2 for a strandless run whose per-chunk cpu share is wide enough
    for ``chunk_worker`` to run its two orientations concurrently instead of in
    series (see ``run_chunks_concurrently``): a chunk-worker slot's real peak is
    then bounded by two resident quant units, not one, and a cap sized on one would
    admit twice the concurrent residency the box was checked to hold. Left at 1,
    this is byte-for-byte the cap a run without that optimization would compute.
    """

    budget = max(1, int(budget))
    estimates = sorted(
        (
            chunk_unit_peak_mib(span) * strand_concurrency_multiplier
            for span in unit_spans_bp
        ),
        reverse=True,
    )
    if not estimates:
        return ChunkMemoryBound(None, None, available_mib, 0, 0)
    largest = estimates[0]
    considered = estimates[: min(budget, len(estimates))]
    if available_mib is None:
        available_mib = available_memory_mib()
    if available_mib is None:
        return ChunkMemoryBound(
            None,
            "chunk-processing concurrency is bounded by --cpu_budget alone: "
            "/proc/meminfo could not be read, so the memory bound could not be "
            "applied on this box. The largest chunk of this partition is "
            "estimated at {} MiB and {} of them would be resident at "
            "once.".format(largest, len(considered)),
            None,
            largest,
            sum(considered),
        )
    available_mib = int(available_mib)
    affordable = 0
    charged = 0
    for estimate in considered:
        if charged + estimate > available_mib:
            break
        charged += estimate
        affordable += 1
    if affordable >= len(considered):
        return ChunkMemoryBound(None, None, available_mib, largest, charged)
    affordable = max(1, affordable)
    charged = sum(considered[:affordable])
    note = (
        "chunk-processing concurrency capped at {} worker(s) by MEMORY rather than "
        "by --cpu_budget {}: {} MiB available, and the {} largest chunk(s) of this "
        "partition are estimated at {} MiB together ({} MiB fixed + {} MiB per "
        "genomic Mb of chunk span, measured; largest single chunk {} MiB). {} "
        "core(s) of the budget sit idle in this phase deliberately; {} worker(s) "
        "charge {} MiB, which the box does have.".format(
            affordable,
            budget,
            available_mib,
            len(considered),
            sum(considered),
            CHUNK_UNIT_FIXED_MIB,
            CHUNK_UNIT_MIB_PER_GENOMIC_MB,
            largest,
            budget - affordable,
            affordable,
            charged,
        )
    )
    if largest > available_mib:
        # Concurrency has nothing left to give at one worker, so say what the user
        # would otherwise learn from dmesg: the PARTITION is the problem, not the
        # width. Not a refusal -- see the docstring -- but not silent either.
        note += (
            " WARNING: one chunk alone is estimated at {} MiB against {} MiB "
            "available, so even a single worker may not fit. A smaller "
            "--approx_MB_per_cut is the only thing that reduces this; the "
            "estimate carries up to 2.3x margin, so the run is attempted "
            "rather than refused.".format(largest, available_mib)
        )
    return ChunkMemoryBound(affordable, note, available_mib, largest, charged)


# Extraction refuses some alignments for reasons no cut selection can predict: a
# record whose aligned end lies past the END OF ITS CONTIG overhangs the last chunk
# without any cut being responsible for it. Those reads are genuinely absent from
# the chunk inputs, so the baseline has to subtract them too or the two arms stop
# consuming the same records. ``verify_severed_accounting`` writes them here, and
# the glob below absorbs them for the baseline while the accounting check excludes
# them from the SELECTION side -- otherwise a file this pipeline wrote would be read
# back as a prediction cut selection never made.
EXTRACTION_ONLY_DROPS = "__extraction_only.dropped_reads.txt"


def severed_read_names(cut_dir, exclude=()):
    """Every read the extractor drops at a cut, across both orientations.

    Cut selection writes one of these per strand.  Their union is exactly the
    difference between the two arms' inputs: measured on a 3.4 Mb contig cut into
    seven segments per strand, the whole-contig bam held 2,266 records, the chunk
    bams held 2,261, and these files named those 5 and nothing else.

    The names rather than the severed bam beside them, because that bam is
    deliberately narrower -- it holds only the severed alignments quantification
    would have used, which is the right set for asking what a cut cost and the
    wrong one for making the two arms consume the same records.

    ``exclude`` drops basenames from the union. Callers asking "what did SELECTION
    predict" pass ``EXTRACTION_ONLY_DROPS``; callers asking "what must the baseline
    subtract" pass nothing, because for that question the two sources are equal.
    """
    names = set()
    for path in sorted(glob.glob(os.path.join(cut_dir, "*.dropped_reads.txt"))):
        if os.path.basename(path) in exclude:
            continue
        with open(path, "rt") as fh:
            names.update(line.strip() for line in fh if line.strip())
    return names


def write_bam_excluding(source_bam, names, dest_bam):
    """Copy ``source_bam`` without the named reads, indexed. Returns kept count."""
    kept = 0
    with pysam.AlignmentFile(source_bam, "rb") as reader:
        with pysam.AlignmentFile(dest_bam, "wb", template=reader) as writer:
            for read in reader.fetch(until_eof=True):
                if read.query_name in names:
                    continue
                writer.write(read)
                kept += 1
    pysam.index(str(dest_bam))
    return kept


# --------------------------------------------------------------------- stage 2


# The depth-window frame every cut is aligned to, pinned to ABSOLUTE 0 by this
# pipeline. ``select_contig_cut_points.py:32-42`` admits a 1-based cut ``b`` only
# when ``(b - grid_origin) % depth_window == 0``, and pinning the frame is what
# lets a chunk's own window origin be ``segment.lend - 1``: per-chunk
# normalization then reproduces whole-contig normalization instead of putting a
# depth window across every boundary. Named here rather than spelled "0" inline
# because a shared cut plan RECORDS it, and a plan whose frame is unstated cannot
# be checked against the run applying it.
SELECTOR_GRID_ORIGIN = 0


def cut_blocking_annotation_models(gtf_path, cuts_by_contig, margin):
    """Models in ``gtf_path`` that a plan's cut positions COLLIDE with.

    THE DIRECT FORM of the property a shared plan must satisfy, and what replaced
    a per-contig sha256 of the annotation the plan was selected against.

    A cut at 1-based ``p`` splits a contig into ``[lend, p]`` and
    ``[p + 1, rend]``, and the extractor emits an annotated LOCUS whole or not at
    all (``extract_contig_region_inputs.Annotation.genes_contained``). So a locus
    is SEVERED by a cut ``p`` exactly when ``locus.lend <= p < locus.rend``. One
    that is severed is contained by neither neighbour, BOTH omit it, and the
    model is quantified by nobody while every chunk reports success. That is the
    failure this began as a guard against.

    THE MARGIN IS THE OTHER HALF, and without it this check was WEAKER than the
    extraction it protects. ``admissibility_offenders`` refuses a boundary whose
    clearance from an annotated locus is under ``--margin``, not merely one that
    severs a locus, and it raises ``ExtractionError`` when it does. So a locus at
    ``[p + 50, p + 500]`` against a cut at ``p`` severs nothing, passed the
    straddle test, and then killed every chunk touching ``p`` at margin 200 --
    minutes into a run, which is exactly what validating the plan up front exists
    to prevent. The predicate here is now the extractor's OWN ``blocks_cut``,
    called with THIS run's ``--margin`` -- the value ``extract_cmd`` hands the
    extractor, not the extractor's ``DEFAULT_MARGIN`` -- so the two cannot answer
    differently and cannot drift apart later. Reused rather than reimplemented
    deliberately: ``lend <= cut + margin and rend >= cut - margin + 1`` is
    asymmetric because a cut falls BETWEEN bases ``p`` and ``p + 1``, and a
    second copy of that arithmetic is a second chance to get it wrong.

    WHY THE DIGEST WENT. It was a proxy for this, and it could not tell a
    legitimately evolved annotation from an incompatible one -- it refused both.
    ONE plan is now emitted for an entire run, before phase 1, so:

      * de novo single-cell emits with NO annotation while every consumer has
        one (per-cluster discovery reads the init GTF, final quant the
        consolidated GTF). Presence equality refused that outright;
      * ref-guided emits against the REFERENCE while phase 2 quantifies the
        CONSOLIDATED GTF, which contains phase-1's novel models. A digest can
        never match, yet the geometry is fine: discovery ran INSIDE those chunks,
        so a novel model cannot span a cut by construction. MEASURED on the
        existing smoke runs -- 0 of 994 consolidated ref-guided models and 0 of
        865 de novo models straddle any of the 9 cut positions.

    Both are accepted here, and a plan whose cut really does land inside a model
    this run quantifies is refused by name -- which the digest could not
    distinguish from either of the above.

    THE GENE LOCUS decides, the TRANSCRIPT is named. Containment is tested per
    gene because the gene is the unit the extractor emits: a gene whose two
    transcripts sit either side of a cut straddles it even though neither
    transcript does, and both transcripts are then lost. Naming the transcript is
    what makes the refusal actionable, so a transcript that collides with the cut
    itself is named when there is one -- preferring one the cut severs -- and the
    gene's first transcript otherwise, with the distinction stated.

    ONE streaming pass over the gtf, restricted to the contigs asked for --
    ``Util_funcs.file_identity_token`` could not serve the digest it replaces for
    the same reason it cannot serve here: Cromwell localizes the same annotation
    to a new path with a fresh mtime in every task.

    Returns one record per (locus, cut) collision, in coordinate order. ``kind``
    is ``"straddle"`` when the cut severs the locus and ``"margin"`` when it
    merely comes too close; a ``"margin"`` record also carries ``side`` and
    ``gap``, the bases strictly between the locus and the cut junction, which is
    under ``margin`` exactly when ``blocks_cut`` fires on either side.
    """

    wanted = {
        chrom: sorted({int(position) for position in positions})
        for chrom, positions in cuts_by_contig.items()
        if positions
    }
    if not wanted:
        return []

    extractor = _extractor_module()
    attribute = extractor._attribute
    genes = {chrom: {} for chrom in wanted}
    transcripts = {chrom: {} for chrom in wanted}
    with extractor._open_gtf_text(str(gtf_path)) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 9:
                continue
            chrom = fields[0].strip()
            if chrom not in wanted:
                continue
            try:
                lend = int(fields[3])
                rend = int(fields[4])
            except ValueError:
                continue
            if rend < lend:
                lend, rend = rend, lend
            gene_id = attribute(fields[8], "gene_id")
            if gene_id is None:
                # Not placeable in a locus, so not emitted by gene containment
                # and not a model this partition can lose.
                continue
            span = genes[chrom].setdefault(gene_id, [lend, rend])
            span[0] = min(span[0], lend)
            span[1] = max(span[1], rend)
            transcript_id = attribute(fields[8], "transcript_id")
            if transcript_id is None:
                continue
            tspan = transcripts[chrom].setdefault(
                transcript_id, [gene_id, lend, rend]
            )
            tspan[1] = min(tspan[1], lend)
            tspan[2] = max(tspan[2], rend)

    by_gene = {chrom: {} for chrom in wanted}
    for chrom, table in transcripts.items():
        for transcript_id, (gene_id, lend, rend) in table.items():
            by_gene[chrom].setdefault(gene_id, []).append((lend, rend, transcript_id))
    for table in by_gene.values():
        for members in table.values():
            members.sort()

    blocks_cut = extractor.blocks_cut
    offences = []
    for chrom in sorted(wanted):
        for position in wanted[chrom]:
            for gene_id, (gene_lend, gene_rend) in genes[chrom].items():
                if not blocks_cut(gene_lend, gene_rend, position, margin):
                    continue
                severed = gene_lend <= position < gene_rend
                members = by_gene[chrom].get(gene_id) or []
                straddling = [m for m in members if m[0] <= position < m[1]]
                blocking = [
                    m for m in members if blocks_cut(m[0], m[1], position, margin)
                ]
                chosen = (straddling or blocking or members or [None])[0]
                if chosen is None:
                    model = "gene {} (whose records name no transcript_id)".format(
                        gene_id
                    )
                    lend, rend = gene_lend, gene_rend
                else:
                    lend, rend, transcript_id = chosen
                    model = "transcript {} of gene {}".format(
                        transcript_id, gene_id
                    )
                    if severed and not straddling:
                        model += (
                            " (no single transcript of that gene spans the cut, "
                            "but the gene LOCUS does, and the locus is what is "
                            "emitted whole or not at all)"
                        )
                    elif not severed and not blocking:
                        model += (
                            " (no single transcript of that gene is that close to "
                            "the cut, but the gene LOCUS is, and the locus is what "
                            "the extractor's boundary check consults)"
                        )
                # Bases strictly between the locus and the cut JUNCTION, which
                # falls between bases ``position`` and ``position + 1``. It is
                # under ``margin`` exactly when ``blocks_cut`` fires, on BOTH
                # sides -- which is why the predicate's asymmetric form is not an
                # asymmetric rule, and why one number states the shortfall.
                if severed:
                    side, gap = None, None
                elif gene_rend <= position:
                    side, gap = "before", position - gene_rend
                else:
                    side, gap = "after", gene_lend - position - 1
                offences.append(
                    {
                        "chrom": chrom,
                        "position": position,
                        "gene_id": gene_id,
                        "gene_lend": gene_lend,
                        "gene_rend": gene_rend,
                        "model": model,
                        "lend": lend,
                        "rend": rend,
                        "kind": "straddle" if severed else "margin",
                        "side": side,
                        "gap": gap,
                        "margin": int(margin),
                    }
                )
    offences.sort(key=lambda o: (o["chrom"], o["position"], o["gene_id"]))
    return offences


def cut_sources(args, strand_bams, inputs_token, split_token):
    """The bams cut selection runs over, and what each run's sentinel chains onto.

    Strand-first: one selection per orientation, over that orientation's bam,
    chained onto the split that produced it.

    Strandless: ONE selection over the RAW bam, chained onto the inputs
    directly, because no split has happened yet. A strandless selection blocks
    on the UNION of both orientations' annotated loci and counts every record,
    so each cut it places serves both orientations -- which is what makes one
    extraction and two quant units possible.

    The key is the orientation the selection is FOR: ``"+"``/``"-"``
    strand-first, ``""`` strandless.

    THE STRANDLESS SELECTION SOURCE IS AN INPUT, not always ``--bam``. Cut
    placement is read-dependent -- blocked positions and the spanning-read cost
    are both computed from the bam handed here -- so two callers with different
    reads get different cut coordinates. That is fatal for the shape
    ``--bam_for_sg`` exists for: WDL/LRAA_quant_by_cluster.wdl hands ONE merged
    normalized bam to ~29 per-cluster quant jobs whose ``--bam`` differs, and
    per-cluster geometry would slice that shared evidence at a different set of
    boundaries in every cluster, leaving each with a different splice graph at
    every locus near a cut.

    ``--bam_for_cut_selection`` names the bam the cuts are chosen from, and it
    MUST be an unthinned SUPERSET of every ``--bam`` the plan is applied to --
    in practice the pre-partition input the cluster bams were split out of. A cut
    safe in the superset is safe in every subset, so every caller gets identical
    geometry AND none of them loses a read.

    The thinned sg bam is NOT a legal source, and that is the trap worth naming:
    normalization REMOVES reads, so a position spanned by nothing in the sg bam
    can still be spanned by a raw read in a caller's own bam. Extraction would
    drop that read, no selector would have named it, and the chunked arm would
    lose a read the unchunked run counts -- silently, and against the parity the
    whole chunked design rests on. ``_run_inner`` refuses that composition, and
    refuses ``--bam_for_sg`` with no explicit source at all rather than picking
    one of the two wrong answers.

    Defaulting to ``--bam`` is today's behaviour exactly, which is why an
    existing run's cuts and sentinels are unchanged by any of this.
    """

    if args.strandless_chunks:
        source = getattr(args, "bam_for_cut_selection", None) or args.bam
        return [("", STRANDLESS_TAG, os.path.abspath(source), inputs_token)]
    return [
        (strand, STRAND_TAG[strand], strand_bams[strand], split_token)
        for strand in ("+", "-")
    ]


def resolve_min_per_id(args):
    """The min_per_id this run actually filters on.

    ``--HiFi`` raises it to 97 inside LRAA and nothing else propagates that, so
    every stage that filters alignments has to derive it the same way. This exists
    as one function because it previously did not: stage 2 computed 97 while stage
    4 passed nothing and the normalizer defaulted to 0, even though
    ``normalize_bam_by_strand.py``'s own help says the value "must match the
    consumer's min_per_id". The two arms then disagreed on one TSS by 4 bp on
    chr21 -- reproducible with a control that placed no cuts at all, so it was the
    pipeline path rather than any boundary.
    """
    # The resolved config already encodes the HiFi preset AND any
    # --config_update layered on top, because LRAA's chunked driver resolves
    # before dispatching (_resolve_chunked_driver_config). Re-deriving the
    # preset here would ignore that file and put prep back out of step with
    # the workers -- the divergence that resolution exists to remove.
    return LRAA_Globals.config["min_per_id"]


def apply_standalone_hifi_preset(args):
    """Seed the HiFi identity floor when THIS MODULE is the entry point.

    ``resolve_min_per_id`` above reads the resolved config rather than
    re-deriving the preset, which is right on the driver route: ``LRAA --chunk``
    applies the preset and any ``--config_update`` and only then calls ``run``.
    Invoked directly, nothing has applied anything, so the config still held the
    library default of 80 while ``--HiFi`` was on the command line.

    The result was two floors in one run rather than a coarser one. Prep selected
    cuts, priced severed reads and normalized at 80; every stage-5 worker gets
    ``--HiFi`` forwarded (see ``lraa_cmd``) and filtered at 97. Nothing reported
    it until a by_chunk single-cell final quant -- whose make_chunks task runs
    this module directly -- refused a shared cut plan its LRAA-driven emitter had
    selected at 97, which is ``validate_cut_plan_geometry`` working.

    Called from ``main`` ONLY. The driver route must keep the config it already
    resolved: re-deriving there would put a ``--config_update`` min_per_id back
    under the preset and flip that precedence in chunked mode alone.

    Only ``min_per_id`` is seeded, because it is the only HiFi-owned key this
    module RESOLVES. Every other preset key is consumed by the LRAA workers,
    which derive it from their own forwarded ``--HiFi``, and several of LRAA's
    depend on flags this parser does not have.
    """

    if not getattr(args, "HiFi", False):
        return
    LRAA_Globals.config["HiFi"] = True
    LRAA_Globals.config["min_per_id"] = LRAA_Globals.HIFI_MIN_PER_ID


def resolve_min_mapping_quality(args):
    """The min_mapping_quality this run's cut selection and normalization filter on.

    Mirrors LRAA's own quant-only/discovery resolution
    (``_normalize_bam_for_splice_graph``): a quant-only chunk has every
    downstream reader of the normalized bam inside the final-quant swap, so the
    final-quant floor is the only one that applies. A discovery chunk builds its
    splice graph from this bam OUTSIDE that swap, at the permissive floor, while
    the same chunk's own final quant reads it inside -- so, like LRAA's own
    resolution, take the more permissive of the two rather than withholding
    evidence discovery wants. Previously this read
    ``LRAA_Globals.config["min_mapping_quality_for_final_quant"]`` directly,
    unconditionally of mode, and that global is never populated from the
    caller's args in chunked mode at all: the config write in ``LRAA`` runs
    after the chunked dispatch already exited, so a non-default value from the
    command line never reached it. This reads the chunked pipeline's own args
    instead, which are now populated by ``lraa_cmd``'s caller either way.
    """
    discovery_mapq = int(args.min_mapping_quality)
    final_mapq = int(args.min_mapping_quality_for_final_quant)
    if not args.discovery:
        return final_mapq
    return min(discovery_mapq, final_mapq)


def cut_selection_plan(args, outdir, cut_dir, source, contig):
    """One per-contig cut selection: where it writes, what it runs, its sentinel.

    Lifted out of the serial loop with the COMMAND unchanged. Three things move,
    and all three name the contig: the output prefix, the log and the sentinel.
    Nothing about how a cut is chosen changes, because the selector already loops
    contigs internally (``select_contig_cut_points.py:1342-1351``) and
    ``--contig`` already restricts it -- per-contig invocation is a refactor of
    the caller and zero lines of the selector.

    What per-contig invocation DOES change is aggregation, and the prefix is
    where it bites. ``.cuts.json``, ``.cuts.txt``, ``.dropped_reads.txt``,
    ``.dropped_reads.tsv`` and ``--severed_reads_bam`` are per-INVOCATION
    artifacts, not per-contig ones -- ``write_severed_alignments_bam``'s own
    docstring says so -- so N invocations sharing one prefix would leave only the
    last contig's on disk. One prefix per contig, and the three that are consumed
    downstream are consumed by a glob-and-union (``severed_read_names``) that N
    files already satisfy.
    """

    key, tag, bam, parent_token = source
    prefix = os.path.join(cut_dir, "{}.{}".format(contig, tag))
    log = os.path.join(
        outdir, "logs", "stage2_cuts_{}_{}.log".format(contig, tag)
    )
    # ".sev" is a version marker, not an input: this stage now also emits
    # <prefix>.severed_reads.bam, so a checkpoint written before that would
    # skip the step and leave the file absent while reporting the stage reused.
    # --HiFi raises min_per_id to 97 inside LRAA, and the emitted severed set is
    # filtered on the value the quant step will use, so the thresholds belong in
    # the token: the same cuts with a different min_per_id emit a different bam.
    effective_min_per_id = resolve_min_per_id(args)
    # Stage 5 runs LRAA --quant_only, which swaps
    # min_mapping_quality_for_final_quant into min_mapping_quality before it
    # filters (LRAA:4201-4204). Reading the discovery key unconditionally would be
    # wrong for every run that raises the final-quant threshold; both default to
    # 0, so nothing would have failed until someone set it. Previously read the
    # global directly rather than through this resolver, and unconditionally of
    # discovery mode -- see resolve_min_mapping_quality's docstring for why both
    # were wrong.
    effective_min_mapq = resolve_min_mapping_quality(args)
    # ``tag`` carries the mode into the sentinel filename -- "strandless"
    # rather than "plus"/"minus" -- and the parent token differs as well
    # (the inputs, not a split), so neither mode can read the other's cuts.
    # ``contig`` is in ``local`` and not merely in the prefix, because two
    # contigs sharing a sentinel would let one contig's completion mark the
    # other's work done.
    # An explicit ``--bam_for_cut_selection`` decides WHICH BAM the selection
    # read (see ``cut_sources``) and therefore the cut coordinates themselves.
    # Appended only when set, so an existing output directory's sentinels are
    # byte-identical without it: a resumed run must never serve cuts computed
    # from one caller's own reads to a run that asked for cuts computed from the
    # shared superset, or the reverse.
    cut_source_bam = getattr(args, "bam_for_cut_selection", None)
    token = chain_token(
        "stage2_cuts_{}_{}.mb_{}_wig_{}_dw_{}_margin_{}.sev_pid_{}_mq_{}"
        ".mxw_{}_annot_{}{}".format(
            contig,
            tag,
            args.approx_MB_per_cut,
            args.approx_MB_per_cut_wiggle_window,
            args.depth_window,
            args.margin,
            effective_min_per_id,
            effective_min_mapq,
            # Two more things that decide the cut COORDINATES. The severed
            # multi-exon weight is the objective itself, so changing it moves
            # cuts; and whether an annotation constrains placement at all is
            # now a property of the run rather than a given. Neither was in
            # the token, and a stale hit here reuses one geometry's cuts under
            # another's name. ``discovery`` is deliberately NOT here: stage 2
            # is identical in the two modes, so including it would force a
            # needless re-selection rather than prevent a wrong reuse.
            args.severed_multiexon_weight,
            bool(args.gtf),
            ".cutsrc_" + Util_funcs.file_identity_token(cut_source_bam)
            if cut_source_bam
            else "",
        ),
        parent_token,
    )
    gtf_args = (
        [
            "--gtf",
            os.path.abspath(args.gtf),
            # Both stages name the same fallback, so a read-only reference
            # directory costs one index build for the run rather than one per
            # invocation. It is built once and SERIALLY by warm_gtf_index before
            # this pool opens: concurrent callers all miss the same cold stamp and
            # each build their own copy, measured six for six.
            "--gtf_index_cache_dir",
            gtf_index_cache_dir(outdir),
        ]
        if args.gtf
        # De novo discovery has no annotation, and the selector treats the
        # annotation constraint as optional: with none supplied every
        # grid-aligned position in the window is admissible on that axis.
        else []
    )
    cmd = [
        sys.executable,
        SELECT_CUTS,
        "--bam",
        bam,
        "--genome_fa",
        os.path.abspath(args.genome_fa),
        *gtf_args,
        "--approx_MB_per_cut",
        str(args.approx_MB_per_cut),
        "--approx_MB_per_cut_wiggle_window",
        str(args.approx_MB_per_cut_wiggle_window),
        "--depth_window",
        str(args.depth_window),
        "--grid_origin",
        str(SELECTOR_GRID_ORIGIN),
        "--margin",
        str(args.margin),
        "--max_intron_length",
        str(args.max_intron_length),
        # The severing cost's shape, and therefore the cut coordinates. Passed
        # explicitly rather than left to the selector's config read so the
        # value that chose the cuts is the value this run recorded.
        "--severed_multiexon_weight",
        str(args.severed_multiexon_weight),
        # Always, not on request.  A consumer deciding whether a chunk boundary
        # dissolved a read-sharing component needs the severed alignments
        # themselves: names cannot be fetched from a coordinate-indexed bam,
        # and a span cannot answer compatibility, which follows exon blocks.
        # Without it the only alternative is trusting per-chunk components
        # silently, and the set is small by construction -- a cut severing many
        # reads is one the selector rejects.
        #
        # One per contig, deliberately NOT merged: nothing in the pipeline reads
        # these bams (verify_severed_accounting and run_baseline's pruning both
        # consume the .txt names), and concatenating record streams from
        # concurrent writers is the one aggregation with nothing to gain.
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
        # The whole point. The flag was already wired for --contig runs; now
        # every run passes it, and a --contig invocation is the one-entry case
        # of the same enumeration rather than a separate branch.
        "--contig",
        contig,
    ]
    # Nothing mode-specific here: severing is priced the same way for
    # quantification and for discovery, because a rule refusing to sever
    # anything declines every cut once every base is covered.
    if not key:
        # Declared, not inferred. Omitting --strand ALSO means "this bam is
        # already orientation-pure, count every record", which is the
        # strand-first case; --strandless says the bam holds both, so the
        # selector emits strandless region strings and reports each
        # orientation's retained count separately.
        cmd.append("--strandless")
    return {
        "kind": "cuts",
        "key": key,
        "tag": tag,
        "bam": bam,
        "contig": contig,
        "unit_id": "cuts_{}_{}".format(contig, tag),
        "prefix": prefix,
        "log": log,
        "cmd": cmd,
        "token": token,
    }


def run_cut_selection(plan, ckpt, outdir, rss_interval):
    """Run (or reuse) one contig's cut selection. ``(timing record, selection)``.

    The selection is read back out of the contig's own ``.cuts.json``, which the
    selector writes as a one-element list because it selected on one contig.
    Unwrapped here so the caller reassembles a list over contigs in a FIXED
    order rather than in completion order.
    """

    if ckpt.done(plan["token"]):
        record = {"reused": True, "cmd": plan["cmd"], "log": plan["log"]}
    else:
        record = run_step(
            "stage2_cuts_{}_{}".format(plan["contig"], plan["tag"]),
            plan["cmd"],
            plan["log"],
            outdir,
            rss_interval,
        )
        ckpt.mark(plan["token"])
    record["contig"] = plan["contig"]
    record["tag"] = plan["tag"]
    path = "{}.cuts.json".format(plan["prefix"])
    with open(path, "rt") as fh:
        payload = json.load(fh)
    if len(payload) != 1:
        raise PipelineError(
            "cut selection for {} wrote {} selection(s) to {}, expected exactly "
            "one: the invocation passed --contig {}, so anything else means the "
            "file is from another run. Use a fresh --output_dir.".format(
                plan["contig"], len(payload), path, plan["contig"]
            )
        )
    selection = payload[0]
    if selection["chrom"] != plan["contig"]:
        raise PipelineError(
            "cut selection for {} wrote a selection for {} to {}. The per-contig "
            "prefixes would then be crossed and the reassembled selection list "
            "would be in the wrong order.".format(
                plan["contig"], selection["chrom"], path
            )
        )
    return record, selection


# --------------------------------------------------------------------- stage 3


def chunk_quant_units(
    chunk_id,
    cdir,
    prefix,
    strand,
    offset,
    order,
    # Defaulted to the quant-only shape this function has always produced, so a
    # caller that predates discovery builds the same units it always did. The
    # real caller passes both explicitly, and a caller that forgot would get
    # paths its LRAA invocation does not write -- a missing file, loudly, not a
    # wrong number.
    lraa_suffix=LRAA_QUANT_ONLY_SUFFIX,
    has_gtf=True,
    has_sg_bam=False,
    has_priors_bam=False,
):
    """The LRAA units one extracted chunk feeds, and where each one's files go.

    Strand-first: one unit, at exactly the paths and sentinel names this stage
    has always used, so an existing output directory still resumes.

    Strandless: two, one per orientation, reading the bam AND the annotation
    that the chunk's own stage 3b writes. They share the mini FASTA, which has
    no orientation and was extracted once, and that is the saving. They cannot
    share the GTF: stage 5 quantifies every transcript the GTF names, so a unit
    handed both orientations' models emits a zero row for each of the other
    orientation's -- 1,110 rows where the strand-first arm emits 555, on the
    same 36 rows carrying reads. Stage 5 consumes exactly what it consumes
    today: one orientation's bam against one orientation's models.

    ``has_gtf`` is False for annotation-free discovery, and the unit then has NO
    annotation rather than an empty one -- an empty GTF handed to LRAA is
    ref-guided discovery against nothing, which is a different run.
    ``lraa_suffix`` is what LRAA will append to --output_prefix, which differs by
    mode, so the paths named here are the paths that mode actually writes.
    ``has_sg_bam`` says the chunk carries caller-supplied splice-graph evidence
    and ``has_priors_bam`` that it carries caller-supplied theta priors, each
    sliced to it by the extractor and split by its own stage 3b. Each unit then
    names the orientation of the slice matching its own reads. ``has_sg_bam``
    additionally means stage 4 is skipped for the unit -- that evidence arrived
    pre-normalized, so normalizing again would compose two acceptance rates. The
    priors slice arrives pre-normalized too but decides nothing about stage 4: a
    unit with priors and no sg evidence still normalizes its own reads, because
    it still needs a splice graph and its priors are not one.
    """

    if strand:
        if has_sg_bam or has_priors_bam:
            raise PipelineError(
                "chunk {} is a strand-FIRST chunk carrying caller-supplied "
                "{}, and there is no orientation-pure slice of it: stage 3b, "
                "which splits those bams, runs only for a strandless chunk. "
                "Reachable only by a caller bypassing the CLI, which refuses "
                "--bam_for_sg and --bam_for_priors with --chunk_by_strand. "
                "Raised rather than returning None, because None here means "
                "'normalize your own' and 'estimate from the sg slice', and the "
                "unit would silently be quantified against a graph built from "
                "its own reads or apportioned by a pooled theta.".format(
                    chunk_id,
                    " and ".join(
                        label
                        for label, present in (
                            ("splice-graph evidence", has_sg_bam),
                            ("theta priors", has_priors_bam),
                        )
                        if present
                    ),
                )
            )
        units = [
            (
                strand,
                chunk_id,
                "chunk.norm.bam",
                "chunk_quant",
                "{}.bam".format(prefix),
                "{}.gtf".format(prefix) if has_gtf else None,
                None,
                None,
            )
        ]
    else:
        units = [
            (
                s,
                "{}_{}".format(chunk_id, STRAND_TAG[s]),
                "chunk.{}.norm.bam".format(STRAND_TAG[s]),
                "chunk_quant_{}".format(STRAND_TAG[s]),
                # both written by stage 3b, inside this chunk's own work
                "{}.strand.{}.bam".format(prefix, s),
                "{}.strand.{}.gtf".format(prefix, s) if has_gtf else None,
                # and so are these two, by the same tool at the same
                # --max_intron_length, so a unit's reads, its evidence and its
                # priors are the same orientation of the same chunk.
                "{}.sg.strand.{}.bam".format(prefix, s) if has_sg_bam else None,
                "{}.priors.strand.{}.bam".format(prefix, s)
                if has_priors_bam
                else None,
            )
            for s in ("+", "-")
        ]
    return [
        {
            "unit_id": unit_id,
            # what the sentinel is named for. The mode is spelled out rather
            # than left implicit in the id shape ("chr1_00_plus" against
            # "chr1_plus_00"), because a checkpoint directory that cannot be
            # read at a glance is its own hazard.
            "sentinel_id": unit_id if strand else "strandless_" + unit_id,
            "chunk_id": chunk_id,
            "strand": s,
            "offset": offset,
            "order": order,
            "bam": bam,
            "gtf": gtf,
            "norm_bam": os.path.join(cdir, norm_name),
            # The pre-normalized evidence stage 5 builds its splice graph from,
            # or None when the run supplied none. Set means stage 4 is SKIPPED
            # for this unit; see _process_unit.
            "sg_bam": sg_bam,
            # The pre-normalized reads stage 5 estimates its FIRST-PASS theta
            # from, or None. Set means this unit's priors are its own rather
            # than the sg bam's, which under a shared sg bam is the pool's.
            "priors_bam": priors_bam,
            "quant_name": quant_name,
            "quant_prefix": os.path.join(cdir, quant_name + "." + lraa_suffix),
        }
        for s, unit_id, norm_name, quant_name, bam, gtf, sg_bam, priors_bam in units
    ]


def planned_chunks_for_selection(source, selection):
    """The chunks ONE contig's selection defines, in segment order, WITHOUT order.

    Split out of ``planned_chunks`` because the global ``order`` counter is the
    one thing a per-contig unit cannot know: it runs over every contig's segments
    in sequence, so it is only defined once every selection is in hand. The pool
    submits extractions from here as each contig finishes and the ``order`` is
    assigned afterwards, serially, by ``planned_chunks``. That is the order
    normalisation -- structural rather than a sort bolted on the end.
    """

    key, tag, bam, _parent_token = source
    chrom = selection["chrom"]
    contig_length = selection.get("contig_length")
    for idx, segment in enumerate(selection["segments"]):
        if key:
            chunk_id = "{}_{}_{:02d}".format(chrom, tag, idx)
            region = "{}{}:{}-{}".format(
                chrom, key, segment["lend"], segment["rend"]
            )
        else:
            chunk_id = "{}_{:02d}".format(chrom, idx)
            region = "{}:{}-{}".format(chrom, segment["lend"], segment["rend"])
        yield {
            "key": key,
            "tag": tag,
            "bam": bam,
            "chrom": chrom,
            "index": idx,
            "chunk_id": chunk_id,
            "region": region,
            "lend": segment["lend"],
            "rend": segment["rend"],
            "contig_length": contig_length,
            # One segment covering [1, contig_length] is a chunk that is not a
            # chunk: the contig was never divided. Recorded here so the extractor
            # call and the plan agree about which chunks those are.
            "spans_whole_contig": (
                contig_length is not None
                and segment["lend"] == 1
                and segment["rend"] == contig_length
            ),
        }


def planned_chunks(sources, selections):
    """Every chunk the run will extract, in extraction order.

    One definition of what a chunk is called and which region it covers, shared
    by stage 3 and by ``--dry_run``, so a printed plan cannot describe a
    partition the real run would not build. Strandless ids and regions carry no
    orientation -- ``chr1_00`` over ``chr1:1-9700000`` -- which is also what the
    extractor needs in order to keep both orientations.

    ``order`` is assigned HERE and only here, over ``sources`` x
    ``selections[key]`` in list order, which is what stage 6 concatenates in. So
    the merged table's row order is set by the order ``selections[key]`` was
    reassembled in, never by the order the parallel selections completed.
    """

    order = 0
    for source in sources:
        for selection in selections[source[0]]:
            for planned in planned_chunks_for_selection(source, selection):
                planned["order"] = order
                yield planned
                order += 1


def extraction_plan(args, outdir, chunk_root, planned, parent_token):
    """One chunk extraction: its directory, its command, its sentinel.

    Strand-first: one chunk per contig-STRAND segment, feeding one quant unit.

    Strandless: one chunk per CONTIG segment, holding both orientations, feeding
    two quant units that this chunk's own strand split produces later. Halving
    the number of extractions is the second saving of the mode, after skipping
    stage 1: the mini FASTA and the mini GTF are written once for the pair
    instead of twice, and the region is read once instead of twice.

    THE WHOLE-CONTIG CASE IS NOT EXTRACTED. A strandless chunk covering
    ``[1, contig_length]`` is a chunk of a contig that was never divided: offset
    0, mini contig named and sized like the real one, nothing able to overhang
    either boundary, so the mini bam it would write is the source's own records
    on that contig, unmoved, minus what ``retained_for_extraction`` rejects --
    and the tool that reads it next applies exactly that filter itself.

    That is not a rare shape. ``cut_targets`` keys on LENGTH, so a contig shorter
    than the segment span gets zero cut targets however deep it is: chrM is
    16,569 bp carrying 1,199,182 mapped reads, 6.6 % of an HG002 library, and it
    is one chunk. MEASURED, that one extraction is 67.5 s of chrM's 73.6 s of
    prep, and the same fetch with the same predicates and no write is 5.45 s --
    92 % of the cost is a copy of something that already exists. Every reference
    holding nothing is the same shape, and a whole-genome fasta has 92 of them.

    So ``--reuse_source_bam``: the region is still read and still counted, so the
    manifest states what a full extraction would have stated and states it from
    its own measurement, but no bam is written and ``files.bam`` names the
    source. Stage 3b then reads the source restricted to this contig. Strandless
    only, because there the chunk bam has exactly ONE consumer -- the split --
    whereas a strand-first chunk bam is normalized and quantified directly, and
    the source is not restricted to one contig for either of those.
    """

    key = planned["key"]
    strandless = not key
    chunk_id = planned["chunk_id"]
    region = planned["region"]
    # THE BAM THIS CHUNK'S READS COME FROM, and strandless is not free to take it
    # from ``planned``. ``cut_sources`` puts the SELECTION source in there, which
    # with --bam_for_cut_selection or --chunk_plan is a shared SUPERSET of every
    # caller's reads -- in the single-cell shape the pre-partition library the ~29
    # cluster bams were split out of. Extracting that would hand every cluster the
    # whole library's reads under its own cluster's name: identical quant for all
    # 29, no error anywhere, and the per-input severed accounting still passes
    # because a superset drops a superset of the names derived from --bam. MEASURED
    # on the defect: `extraction_plan` built `--bam /tmp/superset.bam` for both
    # chunks of a run whose --bam was /tmp/cluster.bam.
    #
    # Strand-first keeps ``planned["bam"]``, which there IS this run's reads -- its
    # own stage-1 orientation split -- and is not the whole bam.
    extraction_bam = os.path.abspath(args.bam) if strandless else planned["bam"]
    reuse_source = (
        strandless
        and planned["spans_whole_contig"]
        and not getattr(args, "no_reuse_source_bam", False)
    )
    chunk_plan_path = getattr(args, "chunk_plan", None)
    sg_bam = getattr(args, "bam_for_sg", None)
    priors_bam = getattr(args, "bam_for_priors", None)
    local = "stage3_extract_{}{}.margin_{}_maxintron_{}{}{}{}{}".format(
        "strandless_" if strandless else "",
        chunk_id,
        args.margin,
        args.max_intron_length,
        # In the sentinel because it decides what is on disk: a directory whose
        # chunk was extracted must not be reported as reusable by a run that
        # would have reused the source, and the reverse.
        ".srcbam" if reuse_source else "",
        # Same reason for the sg slice, which is another FILE IN this directory:
        # a run wanting splice-graph evidence must not be served a directory
        # extracted without it, nor one extracted against a different sg bam.
        ".sg_" + Util_funcs.file_identity_token(sg_bam) if sg_bam else "",
        # And for the priors slice, which is a third file in the same directory
        # and decides this chunk's theta. Empty when no priors bam was given, so
        # a directory extracted before this input existed still resumes.
        ".priors_" + Util_funcs.file_identity_token(priors_bam)
        if priors_bam
        else "",
        # And the PLAN, when the geometry came from one. Cut positions are then not
        # a function of anything else in this token: two plans over the same bam at
        # the same --margin place different boundaries, so without this a resumed
        # run would serve chunks extracted on one plan's geometry to a run applying
        # another's, and every offset downstream would be silently wrong. Named by
        # the plan FILE's identity rather than by "a plan was used", for the same
        # reason the sg slice is.
        ".plan_" + Util_funcs.file_identity_token(chunk_plan_path)
        if chunk_plan_path
        else "",
    )
    cdir = os.path.join(chunk_root, chunk_id)
    os.makedirs(cdir, exist_ok=True)
    prefix = os.path.join(cdir, "chunk")
    log = os.path.join(outdir, "logs", "chunk_{}.log".format(chunk_id))
    gtf_args = (
        [
            "--gtf",
            os.path.abspath(args.gtf),
            "--gtf_index_cache_dir",
            gtf_index_cache_dir(outdir),
        ]
        if args.gtf
        # No annotation to partition. The extractor then writes no mini GTF
        # and reports zero transcripts emitted, which is what a de novo chunk
        # is: a mini contig and its reads.
        else []
    )
    cmd = [
        sys.executable,
        EXTRACT_CHUNK,
        "--genome_fa",
        os.path.abspath(args.genome_fa),
        "--bam",
        extraction_bam,
        *gtf_args,
        "--region",
        region,
        "--output_prefix",
        prefix,
        "--margin",
        str(args.margin),
        # The run's value, at all three of selection, extraction and split,
        # so the three see ONE record set. The extractor defaults to the
        # configured cap instead, which is invisible strand-first -- the bam
        # reaching it was already filtered at the run's value -- and a
        # divergence strandless, where the raw bam still holds the records
        # the split would have discarded and extraction is the first place
        # they are seen.
        "--max_intron_length",
        str(args.max_intron_length),
        "--secondary_alignments",
        "exclude",
    ]
    if reuse_source:
        cmd.append("--reuse_source_bam")
    if sg_bam:
        cmd += ["--sg_bam", os.path.abspath(sg_bam)]
    if priors_bam:
        cmd += ["--priors_bam", os.path.abspath(priors_bam)]
    return {
        "kind": "extract",
        "planned": planned,
        "chunk_id": chunk_id,
        "unit_id": "extract_{}".format(chunk_id),
        "dir": cdir,
        "prefix": prefix,
        "log": log,
        "cmd": cmd,
        "token": chain_token(local, parent_token),
        "reuse_source_bam": reuse_source,
    }


def run_extraction(plan, ckpt, outdir, rss_interval):
    """Extract (or reuse) one chunk. ``(timing record, manifest)``.

    The two manifest refusals live here rather than in the assembly pass so that
    a chunk whose directory came from another run names ITSELF in the pool's
    failure list, beside its own log.
    """

    chunk_id = plan["chunk_id"]
    key = plan["planned"]["key"]
    if ckpt.done(plan["token"]):
        record = {"reused": True, "cmd": plan["cmd"], "log": plan["log"]}
    else:
        record = run_step(
            "stage3_extract_{}".format(chunk_id),
            plan["cmd"],
            plan["log"],
            outdir,
            rss_interval,
        )
        ckpt.mark(plan["token"])
    record["chunk_id"] = chunk_id
    with open("{}.partition.json".format(plan["prefix"]), "rt") as fh:
        manifest = json.load(fh)
    # A manifest extracted for the other mode would be silently wrong rather
    # than absent: a strand-pure chunk bam split again is a no-op on one
    # orientation and an empty bam on the other.
    manifest_strand = manifest.get("strand") or ""
    if manifest_strand != key:
        raise PipelineError(
            "chunk {} carries a manifest for strand {!r}, but this run "
            "extracted it for {!r}. The chunk directory is from another "
            "run; use a fresh --output_dir.".format(chunk_id, manifest_strand, key)
        )
    # The other half of the extractor's refusal, and it applies to BOTH
    # modes: a nonzero count means a strand filter discarded records the
    # source bam still held, which is the raw-flag mistake the ordering
    # exists to prevent -- strand-first because the source should have been
    # orientation-pure already, strandless because no strand filter should
    # have run at all.
    excluded = manifest["counts"].get("opposite_orientation_excluded", 0)
    if excluded:
        raise PipelineError(
            "chunk {} excluded {} record(s) of the opposite orientation at "
            "extraction. Its source bam holds both orientations, so the "
            "records that reached the strand filter carried RAW aligner "
            "flags -- any read the split would have flipped is now in the "
            "wrong chunk. Extract strandlessly and split "
            "afterwards.".format(chunk_id, excluded)
        )
    # Where the chunk's reads live has to be agreed by the manifest and this run,
    # because stage 3b reads that file and would otherwise split the wrong bam
    # with no error of its own -- which is exactly the defect this check caught
    # on a coordinate-remapped corpus before it could reach a number.
    #
    # The two directions are NOT symmetric, and the asymmetry is the point.
    #
    # PLANNED A COPY, GOT A REUSE is impossible to reach honestly and stays
    # fatal: nothing downgrades an extraction into a reuse, so the only way to
    # observe it is a chunk directory left by a different plan.
    #
    # PLANNED A REUSE, GOT A COPY is a legitimate outcome the extractor chooses
    # for itself. Reuse is sound only for a chunk that drops nothing, and whether
    # anything is dropped cannot be known until the region has been read, so the
    # extractor withdraws reuse and materialises the mini bam when it finds an
    # alignment reaching past the contig end. See the REUSE block in
    # util/misc/extract_contig_region_inputs.py at the drop site. Accepting the
    # downgrade here is what lets such a corpus run at all; the alternative is
    # refusing a bam for a defect in the bam rather than in the partition.
    #
    # It is accepted, not ignored: files.bam names the real chunk bam either way,
    # so every later stage follows the manifest and nothing keys on the request.
    reused = bool(manifest.get("bam_reused_from_source"))
    if reused and not plan["reuse_source_bam"]:
        raise PipelineError(
            "chunk {}'s manifest reports bam_reused_from_source=True but this run "
            "planned a full extraction. Nothing upgrades a copy into a reuse, so "
            "the chunk directory predates the plan; use a fresh "
            "--output_dir.".format(chunk_id)
        )
    if plan["reuse_source_bam"] and not reused:
        print(
            "  chunk {}: planned to reuse the source bam, extractor materialised "
            "the chunk bam instead (an alignment reaches past the contig end); "
            "proceeding with the manifest's bam.".format(chunk_id),
            flush=True,
        )
    # AND THE STATE THAT MUST NOT EXIST: reused AND dropping. The extractor
    # withdraws reuse precisely to prevent it, so this is unreachable today -- but
    # that is a property of the current extractor rather than of the contract, and
    # the contract is what later stages rely on. A chunk naming the source while
    # reporting a drop is the SILENT failure: verify_severed_accounting writes the
    # dropped names to EXTRACTION_ONLY_DROPS so the baseline subtracts them, while
    # the reused source still feeds them to the chunked arm, and the two arms then
    # quantify different record sets with both reporting success.
    #
    # Checked here because here it is cheap and precise. verify_chunk_split would
    # catch it a stage later -- reuse is strandless-only, so the split always runs
    # -- but it surfaces as an unexplained record-count mismatch after the chunk
    # has been split, rather than as the contradiction it is, before the expensive
    # phase starts.
    overhang = manifest["counts"].get("alignments_dropped_overhang", 0)
    if reused and overhang:
        raise PipelineError(
            "chunk {} reports bam_reused_from_source=True and {} alignment(s) "
            "dropped for overhang. Those cannot both hold: the manifest says the "
            "chunk holds {} record(s) while naming the source bam, which "
            "restricted to this contig still holds the dropped ones too. Every "
            "stage reading a reused chunk would consume them, and the parity "
            "baseline subtracts them, so the two arms would quantify different "
            "record sets and both report success. Reuse is sound only for a chunk "
            "that dropped nothing.".format(
                chunk_id, overhang, manifest["counts"].get("alignments_emitted")
            )
        )
    return record, manifest


def assemble_chunks(args, timing, sources, selections, extractions):
    """Chunk records in PLANNED order, from manifests the pool produced.

    Every order-sensitive decision of the make-chunks phase is made here, in one
    serial pass over ``planned_chunks``: the ``order`` counter each quant unit
    embeds, the ``chunks`` list stage 6 concatenates through, and the
    ``extract_chunks`` timing list. The pool contributes only a mapping from
    chunk id to what that chunk's extraction did, which has no order in it.
    """

    chunks = []
    for planned in planned_chunks(sources, selections):
        chunk_id = planned["chunk_id"]
        key = planned["key"]
        if chunk_id not in extractions:
            raise PipelineError(
                "chunk {} was planned but never extracted; refusing to build a "
                "partial chunk list".format(chunk_id)
            )
        plan, record, manifest = extractions[chunk_id]
        chunks.append(
            {
                "chunk_id": chunk_id,
                "chrom": planned["chrom"],
                "strand": key,
                "strandless": not key,
                "region": planned["region"],
                "index": planned["index"],
                "order": planned["order"],
                "dir": plan["dir"],
                "prefix": plan["prefix"],
                "log": plan["log"],
                "manifest": manifest,
                "offset": manifest["offset"],
                "window_origin": manifest["window_origin"],
                "extract": record,
                # what stages 3b-5 chain their own sentinels onto
                "upstream_token": plan["token"],
                # Whether this chunk has an annotation at all: stage 3b splits a
                # GTF that exists, and de novo discovery has none.
                "has_gtf": bool(args.gtf),
                "units": chunk_quant_units(
                    chunk_id,
                    plan["dir"],
                    plan["prefix"],
                    key,
                    manifest["offset"],
                    planned["order"],
                    lraa_output_suffix(args.discovery, args.gtf),
                    bool(args.gtf),
                    bool(manifest["files"].get("sg_bam")),
                    bool(manifest["files"].get("priors_bam")),
                ),
            }
        )
    timing.setdefault("stages", {})["extract_chunks"] = [
        c["extract"] for c in chunks
    ]
    return chunks


def run_prep_concurrently(
    args,
    ckpt,
    outdir,
    timing,
    sources,
    rss_interval,
    select_only=False,
    chunk_plan=None,
):
    """Stages 2 and 3 as ONE flat queue over (contig, strand) and (contig, chunk).

    The make-chunks phase used to be two serial loops: one cut selection per
    orientation over every contig at once, then one extraction per chunk. Both
    scale with READS rather than with genome length, and on a whole-genome
    strandless run they model at 1,009 s together -- the Amdahl floor of a
    chunked run, paid before the already-pooled processing phase starts.

    ONE POOL, NOT TWO PHASES. When a contig's cut selection finishes, that
    contig's extractions are submitted to the same executor, so contig A's chunks
    start while contig B is still selecting. With a stage barrier the floor is the
    largest chromosome's whole prep (chr1, 121.6 s measured, an 8.3x ceiling no
    number of cores improves); without one the floor is the longest single
    dependency chain, which is one selection plus its longest extraction. On an
    HG002 library that chain belongs to chrM -- 16,569 bp holding 6.6 % of the
    library, zero cut targets because targeting keys on length -- at 73.6 s, i.e.
    13.7x. A 5-contig 33-unit prototype of exactly this shape measured 150.6 s of
    summed wall against a 41.5 s makespan at ``--cpu_budget 4``: 91 % of ideal.

    Mechanics are ``run_chunks_concurrently``'s, deliberately: one
    ``ThreadPoolExecutor``, ``order_longest_first`` for the launch order, a
    ``failures`` list under a lock with a ``threading.Event`` abort, a refusal to go
    on with a partial result, and a pool bounded by MEMORY as well as by the budget
    (``prep_memory_cap`` here, ``chunk_memory_cap`` there). One thing is this
    phase's own: the queue GROWS while it drains, which the wait loop below has to
    account for.

    ``select_only`` runs the cut selections and stops, for ``--dry_run`` and for
    ``--emit_cut_plan``: the selection is real work whose result the plan
    describes, and extraction is exactly the part neither may perform. ``chunks``
    comes back empty.

    ``chunk_plan`` is the opposite half: a loaded plan whose geometry REPLACES the
    selection, so this phase seeds the extractions themselves. That is what makes
    ~29 per-cluster runs comparable -- one selection over the pre-partition bam,
    identical bounds in every caller, each caller extracting its OWN reads and its
    own slice of the shared splice-graph evidence. With both, there is nothing to
    select and nothing to extract, so the pool is empty by design.

    Returns ``(selections, cut_dir, chunks, makespan, allocation)``.
    """

    from concurrent.futures import ThreadPoolExecutor, wait

    cut_dir = os.path.join(outdir, "cuts")
    os.makedirs(cut_dir, exist_ok=True)
    chunk_root = os.path.join(outdir, "chunks")
    os.makedirs(chunk_root, exist_ok=True)

    contigs, lengths = enumerate_prep_contigs(args)
    if args.gtf:
        # De novo discovery has no annotation, so there is no index to warm and
        # nothing in the pool will seek into one. The pre-warm exists only to
        # keep concurrent callers from each building the same index.
        warm_gtf_index(args, outdir, timing)

    source_by_key = {source[0]: source for source in sources}
    # One idxstats per distinct source bam, not per unit: strand-first selects
    # over two bams whose per-contig counts differ, and a cost proxy taken from
    # the wrong one would order the queue by the wrong orientation's depth.
    costs = {}
    for source in sources:
        costs[source[0]] = contig_mapped_counts(source[2])

    plan_selections = None
    if chunk_plan is not None:
        # CUT SELECTION IS SKIPPED, which is the whole saving: one selection over
        # the pre-partition bam serves every caller, and each caller pays only its
        # own extractions. The geometry is validated against this run's arguments,
        # genome and bam header first -- see selections_from_chunk_plan -- because
        # an unchecked plan does not fail, it produces plausible output at bounds
        # nobody selected.
        plan_selections = selections_from_chunk_plan(
            chunk_plan, getattr(args, "chunk_plan", None), args, contigs
        )

    seed_kind = "selection" if plan_selections is None else "extraction"
    selection_plans = []
    selection_units = []
    if plan_selections is None:
        for source in sources:
            for contig in contigs:
                plan = cut_selection_plan(args, outdir, cut_dir, source, contig)
                selection_plans.append(plan)
                selection_units.append(
                    CpuBudget.WorkUnit(
                        contig_acc=contig,
                        contig_strand=source[0] or STRANDLESS_TAG,
                        # Reads, not bases. Prep scales with the record count: chrM is
                        # 604x shorter than a 10 Mb segment span and still the most
                        # expensive unit on the queue. A reference idxstats could not
                        # report falls back to length, which at least separates a
                        # 250 Mb contig from a 5 kb scaffold.
                        cost=costs[source[0]].get(contig, lengths.get(contig, 0)),
                    )
                )
    elif not select_only:
        # Nothing to select, so the EXTRACTIONS are the seed rather than the
        # fan-out. Same pool, same failure handling, same longest-first launch: the
        # queue simply does not grow, because no task creates another.
        for planned in planned_chunks(sources, plan_selections):
            source = source_by_key[planned["key"]]
            selection_plans.append(
                # Chained onto the INPUTS token (source[3] for a strandless
                # source), which is what stage 2 would have chained onto. The
                # plan's own identity is folded in by extraction_plan, so a
                # directory extracted on another plan's geometry is not reused.
                extraction_plan(args, outdir, chunk_root, planned, source[3])
            )
            span = planned["rend"] - planned["lend"] + 1
            contig_length = planned["contig_length"] or span
            selection_units.append(
                CpuBudget.WorkUnit(
                    contig_acc=planned["chrom"],
                    # The chunk, not the orientation: a strandless plan has one
                    # source and several chunks per contig, so the contig alone
                    # would collapse them onto one key.
                    contig_strand=planned["chunk_id"],
                    # The contig's records apportioned by span, since idxstats
                    # counts per REFERENCE and an extraction reads one segment of
                    # it. Still reads rather than bases, for the reason above.
                    cost=int(
                        costs[planned["key"]].get(planned["chrom"], contig_length)
                        * span
                        / max(contig_length, 1)
                    ),
                )
            )
    if not selection_units and not (plan_selections is not None and select_only):
        raise PipelineError(
            "no make-chunks units to run: {} contig(s) x {} source(s), seeding "
            "{}s".format(len(contigs), len(sources), seed_kind)
        )
    # Keyed on the pair that identifies the unit rather than on the WorkUnit
    # itself: WorkUnit is a namedtuple, so two units agreeing on every field --
    # including a cost of 0, which every reference holding no reads has -- would
    # be the same dict key and one plan would be lost.
    plan_by_unit = {
        (u.contig_acc, u.contig_strand): p
        for u, p in zip(selection_units, selection_plans)
    }
    if len(plan_by_unit) != len(selection_plans):
        raise PipelineError(
            "{} make-chunks unit(s) collapsed to {} distinct (contig, "
            "orientation) key(s). Two units sharing a key would share an output "
            "prefix and a sentinel, so one contig's work would overwrite the "
            "other's and the run would look "
            "complete.".format(len(selection_plans), len(plan_by_unit))
        )
    launch_order = [
        plan_by_unit[(u.contig_acc, u.contig_strand)]
        for u in CpuBudget.order_longest_first(selection_units)
    ]

    seeded = len(selection_plans)
    mem_cap, mem_note = prep_memory_cap(args.cpu_budget)
    # The queue GROWS, so the seed length is a LOWER bound on it and sizing the
    # pool from the seed alone would run a single-contig invocation's extractions
    # one at a time -- which is exactly the serial stage-3 loop this phase exists
    # to remove, and the only part of it a WDL per-contig shard can use (a chr1
    # shard is 25 extractions and 103.9 s of them, measured).
    allocation = CpuBudget.allocate(
        budget=args.cpu_budget,
        num_units=max(seeded, args.cpu_budget),
        max_unit_workers=mem_cap,
    )
    print(CpuBudget.format_allocation(allocation, phase="make_chunks"), flush=True)
    shortfall = CpuBudget.budget_shortfall_note(allocation)
    if shortfall:
        print(shortfall, flush=True)
    if mem_note:
        print(mem_note, flush=True)
    if seed_kind == "selection":
        print(
            "make_chunks queue: {} cut selection(s) seeded over {} contig(s) x {} "
            "source bam(s), longest-first on mapped reads; each contig's chunk "
            "extractions join the SAME pool when its own selection finishes, so the "
            "queue grows past the seed and there is no stage barrier.".format(
                seeded, len(contigs), len(sources)
            ),
            flush=True,
        )
    else:
        print(
            "make_chunks queue: {} extraction(s) seeded over {} contig(s) on the "
            "geometry in {}; cut selection is SKIPPED, so every caller sharing "
            "this plan extracts at identical bounds.".format(
                seeded, len(contigs), getattr(args, "chunk_plan", None)
            ),
            flush=True,
        )

    selections_by_contig = {}
    cuts_records = {}
    extractions = {}
    failures = []
    lock = threading.Lock()
    abort = threading.Event()
    futures = []

    sampler = RssSampler(os.getpid(), rss_interval)
    sampler.start()
    started = time.time()
    try:
        with ThreadPoolExecutor(max_workers=allocation.unit_workers) as pool:

            def submit(fn, payload):
                # The lock is held ACROSS the submit so the drain loop cannot
                # observe a future that exists inside the pool but not yet in
                # ``futures`` -- which would let it exit while an extraction was
                # still queued.
                with lock:
                    futures.append(pool.submit(fn, payload))

            def extraction_task(plan):
                if abort.is_set():
                    return
                try:
                    record, manifest = run_extraction(
                        plan, ckpt, outdir, rss_interval
                    )
                    with lock:
                        extractions[plan["chunk_id"]] = (plan, record, manifest)
                except Exception as err:  # noqa: BLE001 - collected, raised below
                    abort.set()
                    with lock:
                        failures.append(
                            (plan["unit_id"], plan["log"], str(err))
                        )

            def selection_task(plan):
                if abort.is_set():
                    return
                try:
                    record, selection = run_cut_selection(
                        plan, ckpt, outdir, rss_interval
                    )
                    with lock:
                        cuts_records[(plan["key"], plan["contig"])] = record
                        selections_by_contig[(plan["key"], plan["contig"])] = (
                            selection
                        )
                except Exception as err:  # noqa: BLE001 - collected, raised below
                    abort.set()
                    with lock:
                        failures.append(
                            (plan["unit_id"], plan["log"], str(err))
                        )
                    return
                if select_only or abort.is_set():
                    return
                # Guarded too, and not out of caution: a future's exception is
                # never retrieved here (``wait`` does not raise), so an error
                # thrown while BUILDING the fan-out -- an unwritable chunk
                # directory, say -- would vanish and surface later as
                # "planned but never extracted" with nothing named.
                try:
                    source = source_by_key[plan["key"]]
                    for planned in planned_chunks_for_selection(source, selection):
                        submit(
                            extraction_task,
                            extraction_plan(
                                args, outdir, chunk_root, planned, plan["token"]
                            ),
                        )
                except Exception as err:  # noqa: BLE001 - collected, raised below
                    abort.set()
                    with lock:
                        failures.append(
                            (
                                plan["unit_id"] + " (fanning out its chunks)",
                                plan["log"],
                                str(err),
                            )
                        )

            for plan in launch_order:
                submit(
                    selection_task if seed_kind == "selection" else extraction_task,
                    plan,
                )

            # Drain: every round waits for the futures known at its start, and
            # only a task in that set can create new ones. Exiting the ``with``
            # block before the queue is empty would call shutdown() and make a
            # still-running selection's submit raise instead of extracting.
            waited = set()
            while True:
                with lock:
                    fresh = [f for f in futures if f not in waited]
                if not fresh:
                    break
                waited.update(fresh)
                wait(fresh)
    finally:
        sampler.stop()
    makespan = time.time() - started

    if failures:
        lines = [
            "make_chunks unit {} FAILED -- log: {}\n  {}".format(uid, log, msg)
            for uid, log, msg in failures
        ]
        raise PipelineError(
            "{} make-chunks unit(s) failed out of {} cut selection(s) and {} "
            "extraction(s) that reported; refusing to extract or merge a partial "
            "result.\n{}".format(
                len(failures),
                len(cuts_records),
                len(extractions),
                "\n".join(lines),
            )
        )

    if plan_selections is None:
        missing = [
            (key, contig)
            for key in source_by_key
            for contig in contigs
            if (key, contig) not in selections_by_contig
        ]
        if missing:
            raise PipelineError(
                "{} of {} cut selection(s) never reported: {}. Refusing to "
                "assemble a selection list with holes in it -- the missing "
                "contigs would simply be absent from the run.".format(
                    len(missing),
                    seeded,
                    ", ".join("{}{}".format(c, k or "") for k, c in missing[:5]),
                )
            )

    # THE ORDER NORMALISATION. Per-contig selections complete in whatever order
    # the pool finishes them in; ``selections[key]`` is rebuilt in ``contigs``
    # order, which is the selector's own ``sorted(lengths)``. That list feeds
    # planned_chunks -> the ``order`` counter -> chunk_quant_units ->
    # ordered_units -> merge_and_translate, so this one sort is what keeps the
    # merged table's row order identical to a serial run's.
    if plan_selections is not None:
        # Already in ``contigs`` order, and the same order in every caller sharing
        # the plan, which is what makes their ``order`` counters -- and therefore
        # their merged row order -- agree as well as their bounds.
        selections = plan_selections
        timing.setdefault("stages", {})["cuts_" + STRANDLESS_TAG] = {
            "skipped": "geometry supplied by --chunk_plan {}; cut selection "
            "would have chosen per-caller positions from this run's own "
            "bam".format(getattr(args, "chunk_plan", None))
        }
    else:
        selections = {}
        for source in sources:
            key, tag = source[0], source[1]
            selections[key] = [selections_by_contig[(key, c)] for c in contigs]
            timing.setdefault("stages", {})["cuts_" + tag] = [
                cuts_records[(key, c)] for c in contigs
            ]

    if select_only:
        chunks = []
        reused_source = 0
    else:
        chunks = assemble_chunks(args, timing, sources, selections, extractions)
        reused_source = sum(
            1 for c in chunks if c["manifest"].get("bam_reused_from_source")
        )
    timing.setdefault("arms", {}).setdefault("make_chunks", {}).update(
        {
            "cpu_budget": allocation.budget,
            "concurrent_unit_workers": allocation.unit_workers,
            "memory_cap_unit_workers": mem_cap,
            "memory_note": mem_note,
            "unit_peak_rss_budget_mib": PREP_UNIT_PEAK_MIB,
            "contigs": len(contigs),
            # ``seeded`` is what the pool was handed, which is selections in the
            # selecting shape and extractions in the plan-driven one; keeping the
            # two names apart stops a plan-driven run from reporting 300 "cut
            # selections" it never ran.
            "cut_selections": seeded if plan_selections is None else 0,
            "geometry_from_chunk_plan": (
                None if plan_selections is None else getattr(args, "chunk_plan", None)
            ),
            "select_only": bool(select_only),
            "extractions": len(chunks),
            "extractions_reusing_source_bam": reused_source,
            "makespan_s": round(makespan, 3),
            "summed_wall_s": round(
                sum(
                    r.get("wall_s", 0.0)
                    for r in list(cuts_records.values())
                    + [rec for _p, rec, _m in extractions.values()]
                ),
                3,
            ),
            "observed_peak_concurrent_tree_rss_kb": sampler.peak_kb,
        }
    )
    print(
        "make_chunks: {} cut selection(s) + {} extraction(s) in {:.1f}s at {} "
        "worker(s); {} chunk(s) reused the source bam instead of copying "
        "it".format(
            seeded if plan_selections is None else 0,
            len(chunks),
            makespan,
            allocation.unit_workers,
            reused_source,
        ),
        flush=True,
    )
    return selections, cut_dir, chunks, makespan, allocation


def ordered_units(chunks):
    """Every quant unit of the run, in the order stage 6 concatenates them.

    Grouped by orientation and then by extraction order, which for the
    strand-first path is the order the chunks were built in and so leaves that
    path's merged output byte-identical. A strandless run interleaves the two
    orientations chunk by chunk, and this puts them back into the same shape, so
    the two modes' merged tables are comparable row for row rather than only as
    sets.
    """

    rank = {"+": 0, "-": 1}
    return sorted(
        (unit for chunk in chunks for unit in chunk["units"]),
        key=lambda u: (rank[u["strand"]], u["order"]),
    )


# Bumped from 1 when the plan grew its ``geometry`` block. ONE artifact, one
# writer, one version key: a shared cut plan and a scattered leaf's plan describe
# the same partition from different ends, so forking them would let a run be
# handed a file that looks like a plan and answers a different question.
CHUNK_PLAN_VERSION = 2


def cut_geometry_params(args):
    """Every argument that MOVES A CUT, and nothing that does not.

    This is the equality a shared plan is checked on. A plan applied by a run
    configured differently is not a coarser or finer answer, it is a DIFFERENT
    PARTITION: the chunk bounds no longer sit where the reads were selected
    against, so boundary overhang drops differ per caller and the shared
    splice-graph evidence gets sliced at bounds nobody chose. Refused by name and
    value rather than resolved, because nothing downstream can detect it.

    Membership is the same list ``cut_selection_plan``'s sentinel token carries,
    for the same reason -- these decide the coordinates -- plus one the token
    cannot express because it only ever compared a directory to itself:

      ``strandless``           strand-first selects over the stage-1 split bams,
                               which are per-run files; its geometry is not
                               shareable and the mode has to be part of the match.

    THE ANNOTATION IS DELIBERATELY NOT A MEMBER, in any form. It was, as
    ``annotation_present`` here plus a per-contig digest in
    ``selections_from_chunk_plan``, and both were PROXIES for the property that
    actually matters -- that no cut falls inside a model the consuming run
    quantifies -- and both refused correct runs in order to enforce it. ONE plan
    is emitted for a whole run, before phase 1, so in the de novo single-cell
    shape it carries no annotation at all while every consumer has one:
    per-cluster discovery reads the init GTF and final quant the consolidated
    one. Ref-guided is the same story one step later, since phase 2's
    consolidated GTF never equals phase 1's reference. Strict presence refused
    the first and the digest refused the second, both for geometry that is
    perfectly fine. The property is now checked DIRECTLY, against whatever gtf
    each consumer actually holds; see ``cut_blocking_annotation_models``.

    ``discovery`` is deliberately not a member and is still enforced: it is what
    ``resolve_min_mapping_quality`` resolves the effective floor through, so a mode
    disagreement surfaces as a min_mapping_quality mismatch, named with both
    values. Adding it directly would refuse plans that are genuinely identical --
    stage 2 is the same in both modes at the same effective floor.
    """

    return {
        "strandless": bool(args.strandless_chunks),
        "approx_MB_per_cut": float(args.approx_MB_per_cut),
        "approx_MB_per_cut_wiggle_window": float(
            args.approx_MB_per_cut_wiggle_window
        ),
        "depth_window": int(args.depth_window),
        # The window-origin frame. Every chunk's own ``window_origin`` is
        # ``segment.lend - 1``, carried by the segments below rather than repeated
        # here, and is only meaningful against this frame.
        "grid_origin": SELECTOR_GRID_ORIGIN,
        "margin": int(args.margin),
        "max_intron_length": int(args.max_intron_length),
        "severed_multiexon_weight": float(args.severed_multiexon_weight),
        "min_per_id": resolve_min_per_id(args),
        "min_mapping_quality": resolve_min_mapping_quality(args),
    }


def write_chunk_plan(
    path, args, sources, selections, num_total_reads, chunks_extracted
):
    """The partition, written once and applied by other runs. TWO readers.

    ``--only_chunk`` reads ``chunks`` and ``lraa_suffix``: the minimum a
    separately-scheduled chunk task needs that is NOT in its own mini inputs.
    Everything derivable is DERIVED at load time rather than copied here --
    offsets, window origins and counts come from each chunk's own
    ``chunk.partition.json`` and the unit list from ``chunk_quant_units`` -- so
    this file cannot disagree with the chunk directory it describes.

    ``--chunk_plan`` reads ``geometry``: the cut positions themselves, the segment
    spans they produce, and the parameters that chose them. That is what makes ONE
    partition applicable to ~29 per-cluster bams, which is the point -- cut
    placement READS a bam (blocked positions, spanning-read cost), so per-caller
    selection gives per-caller bounds and the shared splice-graph evidence is then
    sliced differently in every cluster.

    ``chunks`` is built from ``planned_chunks`` -- the SAME generator stage 3
    iterates and ``assemble_chunks`` builds its records from -- rather than from
    the assembled chunk list, so one writer serves both callers and the plan
    cannot describe a partition the run would not build. ``chunks_extracted`` says
    whether the directories those entries name exist yet; ``--emit_cut_plan``
    writes geometry with nothing extracted, and ``--only_chunk`` refuses such a
    plan rather than failing later on an absent manifest.

    GEOMETRY ONLY, on purpose. ``num_total_reads`` is the TPM denominator and
    ``discovery`` is the mode, and both belong to the run doing the quantifying:
    they are recorded for the leaf path and for provenance, and
    ``selections_from_chunk_plan`` never reads either.
    """

    entries = []
    for planned in planned_chunks(sources, selections):
        entries.append(
            {
                "chunk_id": planned["chunk_id"],
                "chrom": planned["chrom"],
                "strand": planned["key"],
                "strandless": not planned["key"],
                "region": planned["region"],
                "index": planned["index"],
                "order": planned["order"],
                "has_gtf": bool(args.gtf),
            }
        )
    plan = {
        "version": CHUNK_PLAN_VERSION,
        "chunks_extracted": bool(chunks_extracted),
        # None on an --emit_cut_plan run, which never reaches quantification and
        # so has no denominator to state. Written as null rather than 0, because a
        # zero denominator is a number and would divide.
        "num_total_reads": None if num_total_reads is None else int(num_total_reads),
        "discovery": bool(args.discovery),
        "lraa_suffix": lraa_output_suffix(args.discovery, args.gtf),
        "geometry": {
            "mode": STRANDLESS_MODE if args.strandless_chunks else STRAND_FIRST_MODE,
            "params": cut_geometry_params(args),
            # One entry per selection SOURCE, in ``sources`` order, each holding
            # that source's per-contig selections in ``contigs`` order. Strandless
            # has exactly one, keyed "" ; strand-first has two, and its geometry is
            # recorded for the record rather than for sharing (``--chunk_plan``
            # refuses it: its sources are per-run stage-1 files).
            "by_source": [
                {
                    "key": source[0],
                    "tag": source[1],
                    # PROVENANCE, not something a consumer reproduces: this is the
                    # unthinned superset the cuts were chosen on, and the runs
                    # applying the plan hold subsets of it by construction. The
                    # stat-pair identity is enough to say which file it was on the
                    # machine that selected.
                    "bam": source[2],
                    "bam_identity": Util_funcs.file_identity_token(source[2]),
                    # The selector's OWN output, verbatim: cut positions with their
                    # targets and severing costs, and the segments they produce.
                    # Stored whole rather than distilled so a consuming run rebuilds
                    # the identical `selections` structure stage 3 already consumes,
                    # and so `cut_placement_report` can describe a plan-driven run
                    # exactly as it describes a selecting one.
                    "selections": selections[source[0]],
                }
                for source in sources
            ],
        },
        "chunks": entries,
    }
    with open(path, "wt") as fh:
        json.dump(plan, fh, indent=2, sort_keys=True)
    return plan


def load_chunk_plan(path):
    if not os.path.exists(path):
        raise PipelineError("chunk plan {} does not exist".format(path))
    with open(path, "rt") as fh:
        plan = json.load(fh)
    if plan.get("version") != CHUNK_PLAN_VERSION:
        raise PipelineError(
            "chunk plan {} is version {!r}, this build writes {}".format(
                path, plan.get("version"), CHUNK_PLAN_VERSION
            )
        )
    return plan


def _format_geometry_value(value):
    return repr(value)


def validate_cut_plan_geometry(plan, path, args):
    """Refuse a plan this run's own arguments would not have produced.

    Every check names both values. A mismatched plan does not fail the run
    downstream -- it produces plausible output at bounds nobody selected -- so the
    only place it can be caught is here, against the arguments, before a single
    region is read.
    """

    geometry = plan.get("geometry")
    if not geometry:
        raise PipelineError(
            "chunk plan {} carries no geometry block, so it states no cut "
            "positions. Only a plan written by --emit_cut_plan or "
            "--stop_after_make_chunks can be applied with "
            "--chunk_plan.".format(path)
        )
    mode = STRANDLESS_MODE if args.strandless_chunks else STRAND_FIRST_MODE
    if geometry.get("mode") != mode:
        raise PipelineError(
            "chunk plan {} was selected in {} mode and this run is {}. The two "
            "cut different substrates -- strand-first selects over the stage-1 "
            "orientation split, strandless over the raw bam -- so the positions "
            "are not interchangeable.".format(path, geometry.get("mode"), mode)
        )
    want = cut_geometry_params(args)
    have = geometry.get("params") or {}
    for name in sorted(want):
        if name not in have:
            raise PipelineError(
                "chunk plan {} states no {}, so there is nothing to check this "
                "run's {!r} against. A plan that does not record every parameter "
                "deciding a cut cannot be applied.".format(
                    path, name, want[name]
                )
            )
        if have[name] != want[name]:
            raise PipelineError(
                "chunk plan {} was selected with {} = {} and this run has {}."
                " A plan applied at different geometry is not a smaller answer, "
                "it is a different partition: the chunk bounds no longer sit "
                "where the reads were selected against.".format(
                    path,
                    name,
                    _format_geometry_value(have[name]),
                    _format_geometry_value(want[name]),
                )
            )
    by_source = geometry.get("by_source") or []
    if len(by_source) != 1 or by_source[0].get("key") != "":
        raise PipelineError(
            "chunk plan {} carries {} selection source(s) {}, and a shared plan "
            "has exactly one strandless source. Strand-first geometry is recorded "
            "for the record, not for sharing.".format(
                path,
                len(by_source),
                [entry.get("key") for entry in by_source],
            )
        )
    return geometry


def selections_from_chunk_plan(plan, path, args, contigs):
    """The ``selections`` mapping stage 3 consumes, taken from a shared plan.

    Cut selection is SKIPPED: this returns exactly what
    ``run_prep_concurrently`` would have produced by running it, so
    ``planned_chunks`` assigns the same ids, regions and ``order`` counters, and
    the extractor is handed the same bounds -- against THIS run's own ``--bam``
    and ``--bam_for_sg``.

    THE CONTIGS THIS RUN WILL PROCESS are checked, not every contig the plan
    names. A plan is genome-wide and a caller need not be: ``LRAA_runner.wdl``
    passes ``--chunk`` unconditionally, so a by_chromosome shard is itself a
    chunked run holding one chromosome's bam and one chromosome's fasta, and
    refusing it for the 300 contigs it does not touch would break the mode this
    plan exists to repair. A plan contig this run never reaches is not an error.

    Each contig it DOES reach must be in the plan, in this run's bam header, and
    in its genome fasta at the SAME length. The length is the silent one: the
    final segment of every contig runs to ``contig_length``, so a plan built on a
    different assembly places the last boundary where this run's reads are not.
    And a contig this run selects that the plan does not partition is refused
    rather than selected for locally -- inventing geometry for it would give this
    caller chunk bounds its siblings do not have, which is the whole thing sharing
    a plan prevents.

    THE ANNOTATION IS CHECKED DIRECTLY rather than by identity: no locus in THIS
    run's ``--gtf`` may span a cut the plan places on a contig this run
    processes, nor sit within ``--margin`` of one -- the pair the extractor
    itself refuses. See ``cut_blocking_annotation_models`` for what that catches
    which a digest of the annotation could not, and which correct runs the digest
    refused.
    """

    geometry = validate_cut_plan_geometry(plan, path, args)
    selections = geometry["by_source"][0]["selections"]
    by_chrom = {}
    for selection in selections:
        chrom = selection["chrom"]
        if chrom in by_chrom:
            raise PipelineError(
                "chunk plan {} names contig {} twice; one contig has one "
                "partition.".format(path, chrom)
            )
        by_chrom[chrom] = selection

    with pysam.FastaFile(os.path.abspath(args.genome_fa)) as fasta:
        fasta_lengths = dict(zip(fasta.references, fasta.lengths))
    with pysam.AlignmentFile(os.path.abspath(args.bam), "rb") as bam:
        bam_lengths = dict(zip(bam.references, bam.lengths))

    for chrom in contigs:
        if chrom not in by_chrom:
            raise PipelineError(
                "this run processes contig {}, which the plan {} does not "
                "partition; it names {} contig(s), beginning {}. Selecting cuts "
                "for it locally would give this caller chunk bounds its siblings "
                "do not have, which is what sharing a plan exists to "
                "prevent.".format(
                    chrom, path, len(by_chrom), ", ".join(sorted(by_chrom)[:5])
                )
            )
        planned_length = by_chrom[chrom].get("contig_length")
        if chrom not in fasta_lengths:
            raise PipelineError(
                "this run processes contig {} on the plan {}'s geometry, and it "
                "is absent from the genome fasta {}. The plan and this run "
                "describe different assemblies.".format(chrom, path, args.genome_fa)
            )
        if chrom not in bam_lengths:
            raise PipelineError(
                "this run processes contig {} on the plan {}'s geometry, and it "
                "is absent from the bam header {}. A bam cannot hold a record "
                "against a reference its header does not list, so there is "
                "nothing here to extract the plan's chunks from.".format(
                    chrom, path, args.bam
                )
            )
        # LENGTH is compared against the FASTA only, though presence is required in
        # both. The fasta is what decides geometry -- the last segment of a contig
        # runs to the length the extractor will fetch sequence for, and that fetch
        # is a fasta fetch. A bam header may legitimately state a different length
        # for the same name: aligning against a full assembly and then analysing
        # against a subsetted fasta is a supported configuration, and the
        # single-cell fixture is exactly it -- reads aligned to GRCh38 chr19
        # (58,617,616 bp) analysed against a 2,000,000 bp slice. Comparing the plan
        # to the header refused every shard of that run for a disagreement that
        # predates the plan and that the unchunked path tolerates.
        if planned_length != fasta_lengths[chrom]:
            raise PipelineError(
                "chunk plan {} places contig {}'s cuts on a {} bp contig and "
                "{} says it is {} bp. The last segment of every contig runs "
                "to its end, so applying this plan would put the final "
                "boundary where this run's reads are not.".format(
                    path, chrom, planned_length, args.genome_fa, fasta_lengths[chrom]
                )
            )

    if args.gtf:
        # THE DIRECT CHECK. Cut positions come from the PLAN, models from THIS
        # run's --gtf, scoped to the contigs this run processes -- for a
        # by_chromosome shard, the single contig its slice holds. A plan contig
        # this run never reaches is not checked, for the same reason its length
        # is not: the run cannot sever a model it does not quantify.
        offences = cut_blocking_annotation_models(
            args.gtf,
            {
                chrom: [
                    int(cut["position"])
                    for cut in (by_chrom[chrom].get("cuts") or [])
                ]
                for chrom in contigs
            },
            # THIS RUN'S margin, which is also the plan's -- a disagreement was
            # already refused above by ``validate_cut_plan_geometry`` -- and the
            # value ``extract_cmd`` passes the extractor. Checking against the
            # extractor's DEFAULT_MARGIN instead would validate a run configured
            # at another margin against 200 and then extract it at its own.
            args.margin,
        )
        if offences:
            severing = [o for o in offences if o["kind"] == "straddle"]
            # A severed locus is named first when there is one: it is the graver
            # diagnosis and its remedy is the narrower. Both kinds are refused,
            # because the extractor refuses both.
            first = (severing or offences)[0]
            if first["kind"] == "straddle":
                raise PipelineError(
                    "chunk plan {} places a cut at {}:{} and {} in this run's --gtf "
                    "{} spans it: the model lies at {}:{}-{} and its gene locus at "
                    "{}:{}-{}, so the locus is contained by neither the chunk ending "
                    "at {} nor the one starting at {}. The extractor emits an "
                    "annotated locus whole or not at all, so BOTH neighbours omit it "
                    "and the model is quantified by nobody while every chunk reports "
                    "success. {} model(s) of this annotation straddle a cut and {} "
                    "block one at this run's --margin of {} bp; select the plan "
                    "against an annotation that contains them, or quantify with the "
                    "one it was selected against.".format(
                        path,
                        first["chrom"],
                        first["position"],
                        first["model"],
                        args.gtf,
                        first["chrom"],
                        first["lend"],
                        first["rend"],
                        first["chrom"],
                        first["gene_lend"],
                        first["gene_rend"],
                        first["position"],
                        first["position"] + 1,
                        len(severing),
                        len(offences),
                        args.margin,
                    )
                )
            if first["side"] == "before":
                holder = "ending at {}".format(first["position"])
                neighbour = "starting at {}".format(first["position"] + 1)
            else:
                holder = "starting at {}".format(first["position"] + 1)
                neighbour = "ending at {}".format(first["position"])
            raise PipelineError(
                "chunk plan {} places a cut at {}:{} that {} in this run's --gtf {} "
                "clears by only {} bp: the model lies at {}:{}-{} and its gene locus "
                "at {}:{}-{}, so the chunk {} holds the locus WHOLE -- nothing is "
                "severed -- but the boundary between that chunk and the one {} sits "
                "closer to the locus than this run's --margin of {} bp. The "
                "extractor refuses such a boundary outright "
                "(extract_contig_region_inputs.admissibility_offenders, over the "
                "same blocks_cut predicate checked here), so both chunks either side "
                "of {}:{} would die on it mid-extraction instead. {} locus/loci of "
                "this annotation block a cut at this margin, {} of them by straddling "
                "one; re-emit the plan against an annotation containing this locus -- "
                "cut selection never places a cut within --margin of one, so a plan "
                "selected with this --gtf clears it -- or quantify with the "
                "annotation the plan was selected against.".format(
                    path,
                    first["chrom"],
                    first["position"],
                    first["model"],
                    args.gtf,
                    first["gap"],
                    first["chrom"],
                    first["lend"],
                    first["rend"],
                    first["chrom"],
                    first["gene_lend"],
                    first["gene_rend"],
                    holder,
                    neighbour,
                    first["margin"],
                    first["chrom"],
                    first["position"],
                    len(offences),
                    len(severing),
                )
            )

    rebuilt = []
    for chrom in contigs:
        selection = dict(by_chrom[chrom])
        # THIS RUN'S length, not the plan's, even though the two were just proven
        # equal. ``planned_chunks_for_selection`` derives ``spans_whole_contig``
        # from it, and that is what decides whether the chunk reuses the source
        # bam instead of extracting one -- a decision about THIS run's files, which
        # must never be inherited from a value another run wrote down.
        selection["contig_length"] = fasta_lengths[chrom]
        rebuilt.append(selection)
    return {"": rebuilt}


def rebuild_chunk_record(plan, chunk_id, outdir):
    """One chunk record for ``chunk_worker``, rebuilt from the plan entry plus the
    chunk's OWN ``chunk.partition.json``.

    Paths are re-derived from this machine's ``outdir``; the manifest's own
    ``files`` entries are absolute paths from the machine that extracted and are
    deliberately overwritten. ``upstream_token`` is a fresh constant because
    checkpoint sentinels key on resolved path plus size plus mtime
    (``Util_funcs.file_identity_token``) and so never reproduce across machines --
    a scattered chunk always runs its stages rather than resuming them.
    """

    if plan.get("chunks_extracted") is False:
        # A geometry-only plan from --emit_cut_plan names chunks that do not exist
        # yet: the point of that flag is to pay cut selection ONCE and extraction
        # per caller. Said here, where the plan is identified, rather than left to
        # surface as an absent chunk.partition.json for chunk 00 of 300.
        #
        # Keyed on the field being explicitly False rather than on it being
        # truthy, because the field is a POSITIVE CLAIM and its absence is not the
        # opposite claim: a plan that does not state it still fails, on its own
        # missing manifest, naming the chunk and the path.
        raise PipelineError(
            "the chunk plan names {} chunk(s) but reports chunks_extracted=false, "
            "so it carries cut GEOMETRY only and no chunk directory exists to "
            "process. Extract first: --chunk_plan <that file> "
            "--stop_after_make_chunks.".format(len(plan.get("chunks") or []))
        )
    if plan.get("num_total_reads") is None:
        # The TPM denominator, forwarded to stage 5 as --num_total_reads. A leaf
        # cannot default it: it holds one chunk and has no way to count a library.
        raise PipelineError(
            "the chunk plan states no num_total_reads, so this leaf has no TPM "
            "denominator. Only a plan written by a run given -N carries one."
        )
    entries = [c for c in plan["chunks"] if c["chunk_id"] == chunk_id]
    if not entries:
        raise PipelineError(
            "chunk {} is not in the plan; it names: {}".format(
                chunk_id, ", ".join(c["chunk_id"] for c in plan["chunks"])
            )
        )
    entry = entries[0]

    cdir = os.path.join(outdir, "chunks", chunk_id)
    prefix = os.path.join(cdir, "chunk")
    manifest_path = "{}.partition.json".format(prefix)
    if not os.path.exists(manifest_path):
        raise PipelineError("chunk {} has no {}".format(chunk_id, manifest_path))
    with open(manifest_path, "rt") as fh:
        manifest = json.load(fh)
    if manifest.get("bam_reused_from_source"):
        raise PipelineError(
            "chunk {} was extracted with source-bam reuse, so its manifest names "
            "the whole input bam rather than a mini bam. Re-run make-chunks with "
            "--no_reuse_source_bam".format(chunk_id)
        )
    if manifest.get("sg_bam_reused_from_source"):
        raise PipelineError(
            "chunk {} was extracted with sg-source reuse, so its manifest names "
            "the whole splice-graph bam rather than this chunk's slice of it. "
            "Re-run make-chunks with --no_reuse_source_bam".format(chunk_id)
        )
    if manifest.get("priors_bam_reused_from_source"):
        raise PipelineError(
            "chunk {} was extracted with priors-source reuse, so its manifest "
            "names the whole priors bam rather than this chunk's slice of it. "
            "Re-run make-chunks with --no_reuse_source_bam".format(chunk_id)
        )
    # Read BEFORE files is rebuilt below, which is the point: whether this chunk
    # carries caller-supplied evidence or caller-supplied priors is a fact about
    # what was extracted, while the rebuilt dict is only about where those files
    # live on THIS machine. A rebuild that omitted either key would not fail --
    # the leaf would build the graph from its own reads, or estimate theta from
    # the pool, which is the very defect these inputs exist to close, visible
    # nowhere on the driver's machine.
    had_sg_bam = bool(manifest["files"].get("sg_bam"))
    had_priors_bam = bool(manifest["files"].get("priors_bam"))
    manifest["files"] = {
        "fasta": "{}.fa".format(prefix),
        "bam": "{}.bam".format(prefix),
        "gtf": "{}.gtf".format(prefix) if entry["has_gtf"] else None,
        "dropped_reads": "{}.dropped_reads.txt".format(prefix),
        "sg_bam": "{}.sg.bam".format(prefix) if had_sg_bam else None,
        "priors_bam": "{}.priors.bam".format(prefix) if had_priors_bam else None,
    }

    units = chunk_quant_units(
        chunk_id,
        cdir,
        prefix,
        entry["strand"],
        manifest["offset"],
        entry["order"],
        lraa_suffix=plan["lraa_suffix"],
        has_gtf=entry["has_gtf"],
        has_sg_bam=had_sg_bam,
        has_priors_bam=had_priors_bam,
    )
    return {
        "chunk_id": chunk_id,
        "chrom": entry["chrom"],
        "strand": entry["strand"],
        "strandless": entry["strandless"],
        "region": entry["region"],
        "index": entry["index"],
        "order": entry["order"],
        "dir": cdir,
        "prefix": prefix,
        "log": os.path.join(outdir, "logs", "chunk_{}.log".format(chunk_id)),
        "manifest": manifest,
        "extract": {"reused": True, "cmd": [], "log": None, "chunk_id": chunk_id},
        "offset": manifest["offset"],
        "window_origin": manifest["window_origin"],
        "upstream_token": "scattered_{}".format(chunk_id),
        "has_gtf": entry["has_gtf"],
        "units": units,
    }


# ---------------------------------------------------------------- stages 4 + 5


def resolve_oversimplify_for_chunk(args, chunk):
    """The ``--oversimplify`` value THIS chunk's stage-5 run should get, or None.

    ``--oversimplify`` names contigs as the genome fasta spells them, but each
    chunk is extracted onto a mini contig the extractor RENAMED, so the name a
    caller wrote never appears in the fasta stage 5 runs against. Forwarding the
    caller's list verbatim would match nothing and leave the mode silently off --
    accepted, plumbed, inert, which is the failure the chunked tests exist to
    catch.

    So the decision is made per chunk against the chunk's OWN contig, and the
    value returned is the mini name that contig became. Whitespace around a name
    is tolerated because a workflow builds this list by joining a user's field.
    """

    wanted = getattr(args, "oversimplify", None)
    if not wanted:
        return None
    names = {n.strip() for n in str(wanted).split(",") if n.strip()}
    if chunk.get("chrom") not in names:
        return None
    # Falls back to the original name only when the manifest does not carry a
    # mini name, which is the case where the two are the same thing.
    return (chunk.get("manifest") or {}).get("mini_contig_name") or chunk["chrom"]


def lraa_cmd(
    args,
    bam_for_quant,
    bam_for_sg,
    genome,
    gtf,
    out_prefix,
    num_total_reads,
    cpu_budget,
    bam_for_priors=None,
    oversimplify=None,
):
    """One LRAA invocation. ``cpu_budget`` is this invocation's SHARE of the total.

    A chunk is single-contig and strand-pure, so LRAA's own flat queue inside it has
    exactly one unit: CpuBudget.allocate() gives 1 unit worker and lends the whole share
    to that worker's native tool steps. Nothing is double-counted, because the share is
    already what this pipeline's own chunk concurrency left over.

    Quant-only and discovery differ by exactly two tokens -- whether ``--quant_only``
    is passed, and whether there is a ``--gtf`` to pass -- and in nothing else. The
    splice graph is still built from the per-chunk NORMALIZED bam while the reads are
    counted against the full one, which is what ``--bam_for_sg`` with ``--no_norm``
    says. That is the same composition LRAA performs internally when it normalizes
    for itself, moved out to stage 4 so it runs per chunk and in parallel.

    ``bam_for_priors`` is the third bam and the third ROLE: the slice this unit
    estimates its first-pass theta from. Forwarded only when the caller supplied
    one, so a run without it produces the identical argv it always did. With a
    SHARED ``--bam_for_sg`` and no priors slice, that estimate is taken over
    pooled evidence and every caller's ambiguous-read apportionment depends on
    every other caller's expression -- see LRAA's ``_first_pass_assignment_bam``.

    ``--min_mapping_quality``/``--min_mapping_quality_for_final_quant`` are forwarded
    RAW, not through ``resolve_min_mapping_quality`` -- that resolver is for THIS
    pipeline's own shared stage-2/stage-4 preprocessing, which has one bam to filter
    for two possible consumers. The worker is a full LRAA invocation with its own
    quant-only/discovery swap already built in (LRAA:``_normalize_bam_for_splice_graph``,
    run_quant_only), so it needs both floors distinctly, exactly as the unchunked path
    does, not one pre-collapsed value.

    Single-cell tags follow LRAA's own SUPPRESS convention: forwarded only when the
    caller set one, so an unforwarded default cannot be mistaken for a user's choice.
    ``--cell_list`` is forwarded as the absolute path ``default_args``/``parse_args``
    already resolved it to, since a chunk worker's cwd is its own chunk directory.
    """
    cmd = [
        sys.executable,
        LRAA,
        "--genome",
        genome,
        "--bam",
        bam_for_quant,
        # Explicit, not omitted: --chunk now defaults to True on this SAME LRAA
        # script. A worker's job is to process exactly one already-extracted
        # chunk plain, so it must never re-chunk -- and if it inherited a bare
        # default of True, WORKER_ENV would immediately refuse it (see
        # _run_chunked_mode's re-entry guard), killing every chunked run at its
        # first worker. Omission stopped being safe the moment the default did.
        "--no_chunk",
    ]
    if gtf:
        cmd += ["--gtf", gtf]
    if not args.discovery:
        cmd.append("--quant_only")
    cmd += [
        "--bam_for_sg",
        bam_for_sg,
        "--no_norm",
        "--num_total_reads",
        str(num_total_reads),
        "--cpu_budget",
        str(cpu_budget),
        "--output_prefix",
        out_prefix,
        "--min_mapping_quality",
        str(args.min_mapping_quality),
        "--min_mapping_quality_for_final_quant",
        str(args.min_mapping_quality_for_final_quant),
    ]
    # ONLY when the caller supplied one: absent, LRAA's own pass-1 resolution is
    # unchanged and this argv is byte-identical to the one every chunked run has
    # always produced.
    if bam_for_priors:
        cmd += ["--bam_for_priors", bam_for_priors]
    if args.HiFi:
        cmd.append("--HiFi")
    if getattr(args, "no_rdna_mask", False):
        cmd.append("--no_rdna_mask")
    elif getattr(args, "rdna_mask_fasta", None):
        cmd += ["--rdna_mask_fasta", args.rdna_mask_fasta]
    if hasattr(args, "cell_barcode_tag"):
        cmd += ["--cell_barcode_tag", args.cell_barcode_tag]
    if hasattr(args, "read_umi_tag"):
        cmd += ["--read_umi_tag", args.read_umi_tag]
    if args.cell_list:
        cmd += ["--cell_list", args.cell_list]
    # Explicit either way, not merely omitted when False: --stream_reads and
    # --stream_reads_rescue_unassigned now both default to True on LRAA's own
    # parser (streaming on by default; rescue-inside-streaming on whenever
    # transcriptome rescue itself is on). A chunk worker is a fresh LRAA
    # invocation that resolves its OWN defaults when a flag is simply absent --
    # omitting the flag here would let that worker re-derive True even after
    # this run explicitly asked for False via ChunkedRun's own --no_stream_reads
    # / --no_stream_reads_rescue_unassigned. args.stream_reads_rescue_unassigned
    # is always a resolved bool by here: LRAA's main() resolves the sentinel
    # before dispatching into this pipeline, and this parser's own
    # parse_args()/default_args() resolve it too (see build_parser).
    if args.stream_reads:
        cmd.append("--stream_reads")
        cmd.append(
            "--stream_reads_rescue_unassigned"
            if args.stream_reads_rescue_unassigned
            else "--no_stream_reads_rescue_unassigned"
        )
        if args.stream_reads_rescue_unassigned_to_targets:
            cmd.append("--stream_reads_rescue_unassigned_to_targets")
    else:
        cmd.append("--no_stream_reads")
    # The MASTER rescue opt-out, forwarded on the same principle as the streaming
    # one above and for the same reason: a worker that is not handed it re-derives
    # rescue as True. Placed AFTER the stream if/else, not inside it -- rescue is
    # independent of streaming, and it applies to a --no_stream_reads worker too.
    # getattr, because a caller reaching this through an older namespace may not
    # carry the dest; True keeps that caller's behaviour unchanged.
    if not getattr(args, "rescue_unassigned_reads_via_transcriptome_alignment", True):
        cmd.append("--no_rescue_unassigned_reads_via_transcriptome_alignment")
    # The threshold settings. Forwarded as one --config_update file rather than as
    # per-flag arguments because the config keys a user can set outnumber the flags
    # that name them (min_TSS_iso_fraction and the containment/PolyA thresholds
    # have no flag at all), and because this list is exactly the allowlist that
    # silently omitted every threshold before: a new tunable added upstream now
    # arrives here without anyone remembering to extend this function.
    worker_config_json = getattr(args, "worker_config_json", None)
    if worker_config_json:
        cmd += ["--config_update", worker_config_json]

    # Already resolved to THIS chunk's mini contig name by
    # resolve_oversimplify_for_chunk, or None when the caller's list does not name
    # this chunk's contig.
    if oversimplify:
        cmd += ["--oversimplify", oversimplify]

    return cmd


def assert_extracted_strandlessly(chunk):
    """Refuse to split a chunk that was not extracted with both orientations.

    ``separate_bam_by_strand`` REWRITES ``read.is_reverse`` when the orientation
    it infers disagrees with the aligner, and the extractor's strand filter
    reads the RAW flag. Extract-then-split is therefore the only correct order:
    extracting ``chr1+:...`` from a bam that has not been split puts every
    flipped read in the wrong chunk and produces output that looks entirely
    normal.

    Checked here rather than at the call site because this is the one place the
    order cannot be got wrong by accident -- the chunk in hand either was
    extracted strandlessly or it was not, and no rearrangement of the stages
    changes which. The extractor refuses the same mistake from the other side.
    """

    manifest = chunk["manifest"]
    if manifest.get("strand"):
        raise PipelineError(
            "REFUSING to strand-split chunk {}: it was extracted for strand {!r}, "
            "so its bam was already filtered on the RAW aligner flag. Any read "
            "whose inferred orientation disagrees with its flag is in the wrong "
            "chunk, and splitting now would not move it back. Extract the chunk "
            "strandlessly and split afterwards.".format(
                chunk["chunk_id"], manifest["strand"]
            )
        )
    # The extractor states this outright for a strandless chunk. A manifest
    # written before the key existed does not carry it, and the strand check
    # above is the same statement, so absence is not itself a failure.
    if manifest.get("strand_split_required", True) is not True:
        raise PipelineError(
            "REFUSING to strand-split chunk {}: its manifest reports the bam has "
            "already been strand-separated. Splitting it again would empty one "
            "orientation.".format(chunk["chunk_id"])
        )


def split_chunk_by_strand(args, ckpt, chunk, rss_interval):
    """Stage 3b. Separate ONE chunk's mini bam and mini GTF into two orientations.

    The step the strandless mode exists for: the tool stage 1 runs once over the
    whole bam, run instead over a chunk's own reads, concurrently with every
    other chunk. Same command and therefore the same filters -- a per-chunk
    difference in filtering would be a difference in the record set rather than
    a difference in scheduling, and the arms would no longer be comparable.

    The annotation is partitioned here too, in process, because it is a column-7
    filter on a file the extractor has already rebased -- no coordinates move
    and nothing is re-derived. It has to happen: stage 5 quantifies every model
    the GTF names, so a unit given both orientations reports the other
    orientation's transcripts as zero rows.

    A chunk whose extraction was SKIPPED -- one covering its whole contig, whose
    manifest names the source bam rather than a mini copy -- is split straight out
    of that source, restricted to the contig. Sound, not merely cheaper:
    ``separate_bam_by_strand`` retains exactly what ``retained_for_extraction``
    retains (primary, non-duplicate, non-qcfail, intron within the same cap), so
    the records it writes are the records the skipped copy would have held, and
    ``verify_chunk_split`` still checks its output against the extractor's own
    count of them. What differs is the accounting the splitter PRINTS: it now
    reports the secondary, supplementary and long-intron records it discarded,
    where a mini bam had none left to discard. Same record set, said one stage
    earlier.

    The restriction is a CONTIG NAME rather than a region on purpose. The reuse
    only ever applies to a whole contig, and a name cannot carry the
    1-based/0-based boundary error a region can.

    THE AUXILIARY SLICES ARE SPLIT HERE TOO -- the splice-graph evidence and the
    priors -- by the same tool at the same cap and under the same reuse rule, in
    ONE loop rather than a copy per role. A unit's reads, its evidence and its
    priors have to be the same orientation of the same chunk; any divergence in
    how the files are split is a splice graph built from the other strand's
    reads, or a theta estimated from them.

    Returns ``(step record, token stages 4-5 chain onto, counts)``.
    """

    assert_extracted_strandlessly(chunk)

    cid = chunk["chunk_id"]
    prefix = chunk["prefix"]
    manifest = chunk["manifest"]
    split_prefix = "{}.strand".format(prefix)
    reused = bool(manifest.get("bam_reused_from_source"))
    source_bam = manifest["files"]["bam"] if reused else "{}.bam".format(prefix)
    # ``.get`` on a possibly-absent ``files`` for the same reason line 3340 only
    # indexes it under ``reused``: a manifest that never named files is a legal
    # shape here (hand-built fixtures, and any chunk extracted before the
    # auxiliary slices existed), and an unconditional index turned that into a
    # KeyError three stages after the manifest was read.
    files = manifest.get("files") or {}
    aux = [
        {
            "role": role,
            "source": files.get(role + "_bam"),
            "reused": bool(manifest.get(role + "_bam_reused_from_source")),
            "split_prefix": "{}.{}.strand".format(prefix, role),
            "what": what,
        }
        for role, what in (("sg", "sg slice"), ("priors", "priors slice"))
    ]
    aux = [entry for entry in aux if entry["source"]]
    token = chain_token(
        "stage3b_split_{}.maxintron_{}{}{}".format(
            cid,
            args.max_intron_length,
            # In the sentinel because it decides WHICH FILE was read: a resumed
            # run must not report a split of a mini bam as a split of the source.
            ".srcbam_" + chunk["chrom"] if reused else "",
            # Each auxiliary slice is split by this same stage, so its identity
            # decides this stage's output too. Named here rather than merely
            # inherited from the extraction token, because on the ``--only_chunk``
            # path the upstream token is a constant naming the chunk and nothing
            # else. An sg-only chunk still gets exactly ".sg_<token>", so a
            # directory split before the priors slice existed still resumes.
            "".join(
                ".{}_{}".format(
                    entry["role"], Util_funcs.file_identity_token(entry["source"])
                )
                for entry in aux
            ),
        ),
        chunk["upstream_token"],
    )
    cmd = [
        sys.executable,
        SEPARATE_BAM,
        "--bam",
        source_bam,
        "--output_prefix",
        split_prefix,
        "--max_intron_length",
        str(args.max_intron_length),
    ]
    if reused:
        cmd += ["--contig", chunk["chrom"]]
    for entry in aux:
        entry["cmd"] = [
            sys.executable,
            SEPARATE_BAM,
            "--bam",
            entry["source"],
            "--output_prefix",
            entry["split_prefix"],
            "--max_intron_length",
            str(args.max_intron_length),
        ]
        if entry["reused"]:
            entry["cmd"] += ["--contig", chunk["chrom"]]
    if ckpt.done(token):
        step = {"step": "stage3b_strand_split", "reused": True, "cmd": cmd}
        for entry in aux:
            step[entry["role"] + "_split"] = {"reused": True, "cmd": entry["cmd"]}
    else:
        step = run_step(
            "stage3b_strand_split_{}".format(cid),
            cmd,
            chunk["log"],
            chunk["dir"],
            rss_interval,
        )
        for entry in aux:
            step[entry["role"] + "_split"] = run_step(
                "stage3b_strand_split_{}_{}".format(entry["role"], cid),
                entry["cmd"],
                chunk["log"],
                chunk["dir"],
                rss_interval,
            )
        ckpt.mark(token)

    counts = verify_chunk_split(chunk, split_prefix)
    for entry in aux:
        # Against the extractor's OWN counts for that role, never the read
        # counts: an auxiliary bam is normalized and the reads are not, so the
        # two accountings are different numbers by design and a check that mixed
        # them would fail on correct output. A partition that does not add up
        # here has lost records, and lost records are a quietly sparser splice
        # graph or a quietly different theta -- no error either way.
        counts.update(
            verify_chunk_split(
                chunk,
                entry["split_prefix"],
                counts_prefix=entry["role"] + "_",
                what=entry["what"],
            )
        )
    # Rewritten on every pass, sentinel or not: it is a few hundred kB of text
    # and a resumed run must find the files, not the sentinel that says they
    # were once written.
    # Defaulted True: a chunk record built before annotation-free discovery
    # existed always had a mini GTF, because the pipeline required one.
    if chunk.get("has_gtf", True):
        counts.update(split_chunk_gtf_by_strand(chunk, split_prefix))
    step["counts"] = counts
    return step, token, counts


def split_chunk_gtf_by_strand(chunk, split_prefix):
    """Partition the chunk's mini GTF into one file per orientation.

    Column 7 and nothing else decides which file a line goes to: every line the
    extractor emitted for a gene carries that gene's orientation, and the
    coordinates were rebased when it was written, so this selects lines and
    moves nothing. The transcript counts are checked against the extractor's
    own tally, which is the same discipline the bam split gets -- a partition
    that does not add up is a partition that lost a model, and a lost model is
    a row missing from the merged table rather than an error anyone would see.

    Counted by DISTINCT transcript_id, not by a literal ``transcript`` feature
    row: a transcript is identifiable from its exon lines' shared
    ``transcript_id`` alone, a GTF is not required to also carry a summary
    ``transcript`` row for it, and the extractor's own
    ``gtf_transcripts_emitted`` (util/misc/extract_contig_region_inputs.py's
    ``_GtfIngest``) counts registered transcript ids, not row types. Counting
    row types here disagreed with it on any GTF containing an exon-only
    transcript -- 5 of 654 on ``testing/single_cells/data/chr19.gtf`` alone --
    and failed the chunk outright. Reuses ``_attribute``, the extractor's own
    parser, rather than a second implementation free to disagree with it again.
    """

    source = "{}.gtf".format(chunk["prefix"])
    seen = {"+": set(), "-": set()}
    handles = {}
    attribute = _extractor_module()._attribute
    try:
        for strand in ("+", "-"):
            handles[strand] = open("{}.{}.gtf".format(split_prefix, strand), "wt")
        with open(source, "rt") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.split("\t")
                strand = fields[6] if len(fields) > 7 else None
                if strand not in handles:
                    raise PipelineError(
                        "chunk {} annotation line has orientation {!r}, which is "
                        "neither + nor -: {}".format(
                            chunk["chunk_id"], strand, line.rstrip()
                        )
                    )
                handles[strand].write(line)
                if len(fields) > 8:
                    transcript_id = attribute(fields[8], "transcript_id")
                    if transcript_id is not None:
                        seen[strand].add(transcript_id)
    finally:
        for handle in handles.values():
            handle.close()

    written = {strand: len(ids) for strand, ids in seen.items()}
    emitted = chunk["manifest"]["counts"]["gtf_transcripts_emitted"]
    total = written["+"] + written["-"]
    if total != emitted:
        raise PipelineError(
            "chunk {} annotation split accounting: {} + {} = {} distinct "
            "transcript(s) across the two orientations, but extraction emitted "
            "{}. Every model the chunk holds has to reach exactly one of the "
            "two quant units.".format(
                chunk["chunk_id"], written["+"], written["-"], total, emitted
            )
        )

    return {
        "gtf_transcripts_plus": written["+"],
        "gtf_transcripts_minus": written["-"],
        "gtf_transcripts_emitted": emitted,
    }


def verify_chunk_split(chunk, split_prefix, counts_prefix="", what="mini bam"):
    """Prove the in-chunk split partitioned the chunk and discarded nothing.

    Two identities, both exact, both cheap, and checked on every chunk of every
    strandless run rather than sampled:

    1. ``+`` records plus ``-`` records equals ``alignments_emitted``. The split
       CANNOT discard here, because extraction already applied the same record
       filter -- primary, non-duplicate, non-qcfail, intron within the cap, at
       the same ``--max_intron_length`` -- so every record in the chunk bam is
       one the split retains. A nonzero loss is a filter divergence between the
       two tools, not an expected cost, and it would silently shrink the chunked
       arm's input relative to the control.

    2. Each orientation's record count equals the extractor's own tally for that
       orientation. The extractor counts by raw flag at the moment it writes;
       the split partitions by raw flag, because this pipeline never passes
       ``--infer_read_orient``. Two tools, one number. If inference were ever
       turned on the counts would differ by exactly
       ``num_records_strand_flipped`` -- by design, not by fault -- and this
       check would have to become the sum alone.

    ``counts_prefix`` selects WHICH of the extractor's accountings to check
    against: ``""`` for the read bam, ``"sg_"`` for the splice-graph slice,
    ``"priors_"`` for the theta-priors slice. They are separate keys in the
    manifest precisely so that this check can be run once per role without any
    set of numbers being able to excuse another.
    """

    counts = chunk["manifest"]["counts"]
    emitted = counts[counts_prefix + "alignments_emitted"]
    observed = {}
    for strand in ("+", "-"):
        bam = "{}.{}.bam".format(split_prefix, strand)
        if not os.path.exists(bam):
            raise PipelineError(
                "chunk {} {} strand split produced no {}".format(
                    chunk["chunk_id"], what, bam
                )
            )
        observed[strand] = count_records(bam)

    total = observed["+"] + observed["-"]
    if total != emitted:
        raise PipelineError(
            "chunk {} {} strand split accounting: {} + {} = {} records across "
            "the two orientations, but extraction emitted {}. The split lost {} "
            "record(s) the {} held. Extraction and the split must apply the same "
            "record filter at the same --max_intron_length; they did not.".format(
                chunk["chunk_id"],
                what,
                observed["+"],
                observed["-"],
                total,
                emitted,
                emitted - total,
                what,
            )
        )

    expected = (
        counts.get(counts_prefix + "alignments_emitted_forward"),
        counts.get(counts_prefix + "alignments_emitted_reverse"),
    )
    if None not in expected and (observed["+"], observed["-"]) != expected:
        raise PipelineError(
            "chunk {} {} strand split disagrees with extraction on which records "
            "are which: the split wrote {}/{} (+/-), extraction counted {}/{} by "
            "raw flag. The totals match, so nothing was lost -- records moved "
            "between orientations, which only happens under "
            "--infer_read_orient, and this pipeline does not pass it.".format(
                chunk["chunk_id"],
                what,
                observed["+"],
                observed["-"],
                expected[0],
                expected[1],
            )
        )

    return {
        counts_prefix + "alignments_emitted": emitted,
        counts_prefix + "records_plus": observed["+"],
        counts_prefix + "records_minus": observed["-"],
        counts_prefix + "records_total": total,
        counts_prefix + "records_lost": emitted - total,
    }


def _process_unit(
    args,
    ckpt,
    chunk,
    unit,
    upstream_token,
    log,
    cdir,
    num_total_reads,
    rss_interval,
    cpu_budget,
):
    """Stages 4 and 5 for ONE quant unit (one orientation of one chunk).

    Split out of ``chunk_worker`` so a strandless chunk's two orientations can
    run either in series on the chunk's own log (today's default, when the
    chunk's cpu share is too thin to split) or concurrently on their own logs
    and half that share (see ``chunk_worker`` and
    ``run_chunks_concurrently``'s memory-aware two-pass allocation). Returns
    this unit's own steps list; raises ``PipelineError`` naming THIS unit's log,
    whichever one the caller gave it.
    """
    cid = chunk["chunk_id"]
    prefix = chunk["prefix"]
    steps = []

    uid = unit["sentinel_id"]
    sg_bam = unit.get("sg_bam")
    priors_bam = unit.get("priors_bam")
    if sg_bam:
        # STAGE 4 IS SKIPPED, only here and deliberately, and what it skips is
        # normalizing THIS UNIT'S OWN READS to build the splice graph from --
        # which is precisely what ``--bam_for_sg`` says not to do, because the
        # caller already supplied that evidence normalized. Thinning it a second
        # time COMPOSES two acceptance rates: the splice graph divides each
        # read's acceptance probability back out through XW (LRAA:165-192), so a
        # second pass leaves the graph's weights disagreeing with the reads that
        # built it. LRAA refuses to re-normalize a --bam_for_sg for exactly this
        # reason (LRAA:2329-2340); this is that same decision taken one stage
        # earlier, where the subprocess would otherwise have run.
        #
        # The priors slice does not enter this decision and never did: it too
        # arrives pre-normalized, so there is nothing for stage 4 to do for it
        # either, and a unit carrying priors but no sg evidence still normalizes
        # its own reads below -- it needs a splice graph, and its priors are not
        # one.
        norm_bam = None
        # Stage 5 chains onto whatever this branch produces, so the skip still
        # has to name the evidence's identity. Without it a resumed unit could be
        # served a quant built against DIFFERENT evidence purely because no
        # normalization ran for the sentinel to key on.
        quant_upstream_token = chain_token(
            "stage4_skipped_presupplied_sg_{}.sg_{}".format(
                uid, Util_funcs.file_identity_token(sg_bam)
            ),
            upstream_token,
        )
        steps.append(
            {
                "step": "stage4_normalize",
                "unit": unit["unit_id"],
                "skipped": "splice-graph evidence was supplied pre-normalized as "
                "--bam_for_sg; normalizing it again would compose two acceptance "
                "rates and leave the graph's XW weights disagreeing with the "
                "reads it was built from",
                "sg_bam": sg_bam,
            }
        )
    else:
        # Same genome LRAA itself will be handed below (genome="{}.fa".format(prefix)
        # at the lraa_cmd call further down) -- built here too because stage 4's
        # normalize_bam_by_strand.py subprocess needs a BED, not the in-process mask
        # the per-chunk LRAA invocation builds for itself. RdnaMask.build_rdna_mask_bed
        # caches on (genome, cassette, pad) identity, so calling it again for the
        # chunk's OTHER orientation, or again inside that LRAA invocation's own
        # main(), is a cache hit rather than a second minimap2 run. Built inside
        # this branch because stage 4 is its only consumer: the sg path above runs
        # no normalizer, and LRAA builds its own mask in process.
        rdna_mask_bed = None
        if not getattr(args, "no_rdna_mask", False):
            rdna_fasta = RdnaMask.resolve_rdna_fasta(
                getattr(args, "rdna_mask_fasta", None)
            )
            rdna_mask_bed = RdnaMask.build_rdna_mask_bed(
                "{}.fa".format(prefix),
                rdna_fasta,
                cache_dir=os.path.join(cdir, "__rdna_mask_cache"),
            )

        norm_bam = unit["norm_bam"]
        quant_upstream_token = chain_token(
            # maxintron is passed to the normalizer below, so it decides these
            # contents and has to name them. So does min_per_id: the normalizer
            # discards alignments below it before measuring depth, so the same
            # input at a different threshold yields a different normalized bam,
            # and a cache written before it was forwarded must not be reused.
            # Same for min_mapping_quality, added alongside it below --
            # previously absent from both this token and norm_cmd, so a
            # non-default value from the command line silently normalized as
            # though it were 0.
            "stage4_norm_{}.maxcov_{}_dw_{}_seed_{}_origin_{}_maxintron_{}"
            ".minperid_{}_minmapq_{}_rdnamask_{}".format(
                uid,
                args.normalize_max_cov_level,
                args.depth_window,
                args.random_seed,
                chunk["window_origin"],
                args.max_intron_length,
                resolve_min_per_id(args),
                resolve_min_mapping_quality(args),
                Util_funcs.file_identity_token(rdna_mask_bed)
                if rdna_mask_bed
                else "none",
            ),
            upstream_token,
        )
        norm_cmd = [
            sys.executable,
            NORMALIZE_BAM,
            "--input_bam",
            unit["bam"],
            "--output_bam",
            norm_bam,
            "--normalize_max_cov_level",
            str(args.normalize_max_cov_level),
            "--depth_window",
            str(args.depth_window),
            "--random_seed",
            str(args.random_seed),
            # The normalizer discards alignments below this before measuring
            # depth, so it has to be the same value the rest of the run filters
            # on. Its own help says "must match the consumer's min_per_id"; it
            # did not, and the chunked and unchunked arms disagreed on a TSS by
            # 4 bp as a result.
            "--min_per_id",
            str(resolve_min_per_id(args)),
            "--min_mapping_quality",
            str(resolve_min_mapping_quality(args)),
            "--max_intron_length",
            str(args.max_intron_length),
            # true of both modes' units: strand-first because the whole bam was
            # split before extraction, strandless because stage 3b split this
            # chunk above.
            "--input_is_single_strand",
            "--window_origin",
            str(chunk["window_origin"]),
        ]
        if rdna_mask_bed:
            norm_cmd += ["--rdna_mask_bed", rdna_mask_bed]
        if ckpt.done(quant_upstream_token):
            steps.append(
                {
                    "step": "stage4_normalize",
                    "unit": unit["unit_id"],
                    "reused": True,
                    "cmd": norm_cmd,
                }
            )
        else:
            steps.append(
                run_step(
                    "stage4_normalize_{}".format(uid),
                    norm_cmd,
                    log,
                    cdir,
                    rss_interval,
                )
            )
            ckpt.mark(quant_upstream_token)

    if priors_bam:
        # The priors slice decides stage 5's THETA, so a resumed unit must not be
        # served a quant estimated from different priors. The stage-5 sentinel
        # keys on argv, and argv cannot say so: the slice path is stable across
        # runs while its contents are not. Chained here rather than inside either
        # branch above, because a unit can carry priors with or without
        # caller-supplied evidence.
        quant_upstream_token = chain_token(
            "stage5_priors_{}.priors_{}".format(
                uid, Util_funcs.file_identity_token(priors_bam)
            ),
            quant_upstream_token,
        )

    # LRAA derives its scratch roots by string concatenation on
    # --output_prefix ("__{prefix}.contigtmp", "__{prefix}.sgcache"), so an
    # ABSOLUTE prefix would produce nonsense paths like
    # "__/abs/path.contigtmp" rooted at the cwd. Give it a bare name and let
    # cwd place the outputs.
    quant_prefix = unit["quant_prefix"]
    quant_cmd = lraa_cmd(
        args,
        bam_for_quant=unit["bam"],
        bam_for_sg=sg_bam or norm_bam,
        # None unless the caller supplied one. Set, it is the ONLY thing pass 1
        # reads, and the unit's theta is then estimated from its own normalized
        # reads rather than from the evidence it shares with every sibling.
        bam_for_priors=priors_bam,
        # ONE mini contig for the pair -- sequence has no orientation and
        # the extraction wrote it once -- but each unit's OWN annotation,
        # because stage 5 quantifies every model its GTF names.
        genome="{}.fa".format(prefix),
        gtf=unit["gtf"],
        out_prefix=unit["quant_name"],
        num_total_reads=num_total_reads,
        cpu_budget=cpu_budget,
        oversimplify=resolve_oversimplify_for_chunk(args, chunk),
    )
    # Built BEFORE the sentinel because the sentinel is keyed on it. The old
    # key named the unit, the read denominator and --HiFi, so quant-only
    # versus discovery was absent from it and so was every flag a caller
    # might later forward -- a resumed run could serve one mode's output for
    # the other, which no downstream check can detect. args.discovery needs
    # no field of its own: it IS the presence of --quant_only in this argv.
    quant_token = chain_token(
        "stage5_quant_{}.argv_{}".format(uid, quant_command_digest(quant_cmd)[:12]),
        quant_upstream_token,
    )
    if ckpt.done(quant_token):
        steps.append(
            {
                "step": "stage5_quant",
                "unit": unit["unit_id"],
                "reused": True,
                "cmd": quant_cmd,
            }
        )
    else:
        steps.append(
            run_step(
                "stage5_quant_{}".format(uid), quant_cmd, log, cdir, rss_interval
            )
        )
        ckpt.mark(quant_token)

    expr_path = quant_prefix + ".quant.expr"
    if not os.path.exists(expr_path):
        raise PipelineError(
            "chunk {} unit {} produced no {}; log: {}".format(
                cid, unit["unit_id"], expr_path, log
            )
        )
    if resolve_tracking(quant_prefix) is None:
        raise PipelineError(
            "chunk {} unit {} produced no {}{{{}}}; log: {}".format(
                cid,
                unit["unit_id"],
                quant_prefix,
                ",".join(QUANT_TRACKING_SUFFIXES),
                log,
            )
        )
    if args.discovery:
        gtf_path = quant_prefix + ".gtf"
        if not os.path.exists(gtf_path):
            raise PipelineError(
                "chunk {} unit {} produced no {}; discovery is what this run "
                "is for, so a missing model set is a failure, not an empty "
                "result. log: {}".format(cid, unit["unit_id"], gtf_path, log)
            )

    return steps


def chunk_worker(args, ckpt, outdir, chunk, num_total_reads, rss_interval, cpu_budget):
    """Stages 3b, 4 and 5 for one chunk. Everything goes to the chunk's own log.

    ``cpu_budget`` is this chunk's share of the total, handed down by the scheduler.

    A strand-first chunk is one quant unit -- ``CpuBudget.allocate`` below gives
    it ``unit_workers == 1`` regardless of ``cpu_budget``, so it always takes
    the serial branch. A strandless chunk is TWO, and the orientation split
    that produces them runs HERE -- inside the chunk, against the chunk's own
    reads, concurrently with every other chunk -- instead of once over the
    whole bam before any of this starts. The two units then run concurrently
    on their own half of this chunk's share when that share is 2+ cores, or
    one after the other on the whole share when it is 1 -- see
    ``run_chunks_concurrently``'s docstring for why the memory bound has
    already accounted for this before handing down a share wide enough to
    trigger it. Either way they share the mini FASTA and mini GTF the single
    extraction wrote.
    """

    from concurrent.futures import ThreadPoolExecutor

    cid = chunk["chunk_id"]
    log = chunk["log"]
    cdir = chunk["dir"]
    prefix = chunk["prefix"]
    started = time.time()
    steps = []

    upstream_token = chunk["upstream_token"]
    if chunk["strandless"]:
        split_step, upstream_token, split_counts = split_chunk_by_strand(
            args, ckpt, chunk, rss_interval
        )
        steps.append(split_step)
        chunk["split_counts"] = split_counts

    units = chunk["units"]
    inner = CpuBudget.allocate(budget=cpu_budget, num_units=len(units))
    unit_logs = [log]
    if inner.unit_workers > 1:
        # Enough of this chunk's own share is spare (fewer chunks than the
        # box has cores -- see run_chunks_concurrently) to run both
        # orientations at once instead of one after the other. Each gets its
        # OWN log while it runs: run_step appends to a single file, and two
        # subprocesses writing the same one at once would interleave into
        # something unreadable. The two are merged back into the chunk's ONE
        # log below, in declared unit order, so "everything goes to the
        # chunk's own log" -- this function's own opening line -- still holds
        # for any reader or test that looks at chunk["log"] alone; only the
        # WRITE-time interleaving hazard is what the split log avoids.
        unit_logs = ["{}.{}".format(log, unit["sentinel_id"]) for unit in units]
        with ThreadPoolExecutor(max_workers=inner.unit_workers) as pool:
            futures = [
                pool.submit(
                    _process_unit,
                    args,
                    ckpt,
                    chunk,
                    unit,
                    upstream_token,
                    ulog,
                    cdir,
                    num_total_reads,
                    rss_interval,
                    inner.tool_threads,
                )
                for unit, ulog in zip(units, unit_logs)
            ]
            # Collected before merging or raising: a failure in one unit must
            # not orphan the OTHER unit's log content in a per-unit file
            # nobody reads by default -- the merge below runs regardless.
            unit_results = []
            first_error = None
            for future in futures:
                try:
                    unit_results.append(future.result())
                except Exception as err:  # noqa: BLE001 - merged in below, re-raised after
                    unit_results.append(None)
                    if first_error is None:
                        first_error = err
        # The per-unit logs are RETAINED, not removed: run_step's own promise
        # ("the log is never removed, on success or on failure") applies here
        # too. The merge above gives chunk["log"] the complete picture for a
        # reader who only looks at one file; these are the same content
        # attributed to its own orientation, kept as an independent artifact
        # in case the merge itself is what needs debugging.
        with open(log, "at") as ofh:
            for ulog in unit_logs:
                if os.path.exists(ulog):
                    with open(ulog, "rt") as ifh:
                        ofh.write(ifh.read())
        if first_error is not None:
            raise first_error
        for result in unit_results:
            steps.extend(result)
    else:
        for unit in units:
            steps.extend(
                _process_unit(
                    args,
                    ckpt,
                    chunk,
                    unit,
                    upstream_token,
                    log,
                    cdir,
                    num_total_reads,
                    rss_interval,
                    cpu_budget,
                )
            )

    return {
        "chunk_id": cid,
        "region": chunk["region"],
        "log": log,
        "unit_logs": unit_logs,
        "wall_s": round(time.time() - started, 3),
        "peak_tree_rss_kb": max(
            [s.get("peak_tree_rss_kb", 0) for s in steps] or [0]
        ),
        "steps": steps,
    }


def run_chunks_concurrently(args, ckpt, outdir, chunks, num_total_reads, rss_interval):
    """Run every chunk's stages 3b-5, scheduling the flat unit queue from ONE budget.

    Concurrency is not a knob here any more. A chunk is one unit of the same flat queue
    LRAA schedules internally, so CpuBudget.allocate() derives both how many chunks run
    at once and the share each one's LRAA invocation is given. Their product cannot
    exceed the budget, which is what --concurrency and LRAA's own contig pool could
    previously do to each other.

    The scheduled unit is the CHUNK, not the quant job. A strandless chunk holds
    two quant jobs. When this chunk's own cpu share (``tool_threads``) is only
    1 core, they still run in series on it, exactly as before -- splitting a
    single core between two invocations buys nothing. When the share is 2+
    cores -- the common shape once chunk count drops below core count, e.g. one
    Terra-scattered chromosome producing fewer chunks than the shard's cores --
    ``chunk_worker`` runs both orientations concurrently on their own half of
    that share (see ``CpuBudget.allocate`` applied a second, nested time
    there). The per-chunk product still cannot exceed this chunk's own budget,
    which is what keeps the whole-run product bounded by ``args.cpu_budget``.

    Launch order is longest-first on retained alignments per chunk, which the extractor
    already counted, so no extra pass is needed. Span would be the wrong proxy: it does
    not bound read count, and a short highly expressed chunk can outweigh a long quiet
    one. For a strandless chunk that cost is both orientations together, which is
    what the chunk's worker will actually do.

    The pool is bounded by MEMORY as well as by the budget (``chunk_memory_cap``),
    for the same reason the make-chunks pool is: this phase's per-unit footprint is
    multi-GiB at whole-chromosome span, and ``min(cpu_budget, units)`` alone let 16
    of them go resident at once.

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
    # Both bounds, not just the budget: a partition whose units are whole
    # chromosomes puts every worker's multi-GiB quant footprint resident at once,
    # which is what the kernel OOM killer ended a 195-chunk whole-genome run over.
    #
    # Two passes when strandless: the first pass's own allocation tells us whether
    # a chunk-worker slot will run one orientation at a time or two concurrently
    # (chunk_worker makes that same call from its OWN cpu share, so this has to
    # predict it before knowing the final answer). Below tool_threads 2 nothing
    # changes -- one resident quant unit per slot, same charge as always. At or
    # above it, a slot's real peak is bounded by TWO resident units, so the memory
    # bound is recomputed charging double before the concurrency it protects
    # against is allowed to run. This can only ever REDUCE unit_workers relative to
    # the first pass, which by construction cannot lower tool_threads back below 2
    # for a fixed budget -- so the second pass is self-consistent with no further
    # iteration needed.
    spans = [unit.region[1] - unit.region[0] + 1 for unit in units]
    bound = chunk_memory_cap(args.cpu_budget, spans)
    allocation = CpuBudget.allocate(
        budget=args.cpu_budget,
        num_units=len(units),
        max_unit_workers=bound.cap,
    )
    if args.strandless_chunks and allocation.tool_threads >= 2:
        bound = chunk_memory_cap(
            args.cpu_budget, spans, strand_concurrency_multiplier=2
        )
        allocation = CpuBudget.allocate(
            budget=args.cpu_budget,
            num_units=len(units),
            max_unit_workers=bound.cap,
        )
    print(CpuBudget.format_allocation(allocation, phase="chunks"), flush=True)
    shortfall = CpuBudget.budget_shortfall_note(allocation)
    if shortfall:
        print(shortfall, flush=True)
    if bound.note:
        print(bound.note, flush=True)

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
    return ordered, makespan, allocation, bound


# --------------------------------------------------------------------- stage 6


def read_tsv(path):
    """Returns (comments, header_fields, rows). Comment lines start with '#'."""

    comments, header, rows = [], None, []
    with open_text(path) as fh:
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


def iter_tsv(path):
    """``read_tsv`` without the list: returns (header_fields, row_iterator).

    Same file grammar as ``read_tsv`` -- first non-comment line is the header,
    ``#`` lines are skipped wherever they appear, no header is a
    ``PipelineError`` -- but rows arrive one at a time and are never
    accumulated. Comments are dropped rather than collected, because the one
    caller that streams (the stage-6 tracking merge) discards them anyway.

    The handle is closed when the iterator is exhausted or closed, so the caller
    MUST either drain it or let it be garbage collected; holding a half-read
    iterator holds the file open.
    """

    fh = open_text(path)
    header = None
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith("#"):
            continue
        header = line.split("\t")
        break
    if header is None:
        fh.close()
        raise PipelineError("{} has no header line".format(path))

    def rows():
        try:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    continue
                yield line.split("\t")
        finally:
            fh.close()

    return header, rows()


# gene_id / transcript_id in a GTF attribute column, quoted.
_GTF_ID_ATTR = re.compile(r'(gene_id|transcript_id)\s+"([^"]*)"')


# Separates a chunk's unit id from the chunk-local model id it qualifies. NOT "|":
# gffcompare's .tracking file delimits the subfields of its query column with an
# unescaped pipe -- `qJ:gene_id|transcript_id|num_exons|FPKM|TPM|cov|len` -- so a
# model id containing one shatters that record positionally. MEASURED on the
# single-cell fixture: 182 tracking rows parsed to keys like "chr19_00_plus" (the
# unit id alone) and "g:chr19:+:comp-44" (the model id stripped of its namespace),
# neither of which is an id anything downstream holds, so every gene symbol went
# unassigned and incorporate_gene_symbols_in_sc_features.py refused the run. The
# pipe was invisible until then because chunked discovery's own merge round-trips
# it fine; only a THIRD-party parser exposed it.
#
# "@" because it is unclaimed everywhere these ids travel: GTF attribute values
# are quoted, TSV needs only tab-freedom, gffcompare treats it as ordinary, and it
# is not regex-special the way "^" is for the R/Seurat feature-name patterns these
# same ids reach. Deliberately not "." "-" or "_", all of which already occur
# inside both unit ids and model ids and so could not be split back out.
NAMESPACE_SEP = "@"


def _namespace_id(unit_id, value):
    """A chunk-local model id made unique across the run."""

    return "{}{}{}".format(unit_id, NAMESPACE_SEP, value)


VERSION_COMMENT_PREFIX = "# LRAA version "


def merged_provenance_header(units, suffix):
    """The provenance lines a MERGED artifact must carry, read from its inputs.

    A chunked run's merged output is the file that SHIPS, and it used to carry no
    version at all: the GTF merge wrote only its own unit-count line, the
    quant.expr merge wrote nothing, while every per-chunk input carried both
    ``# LRAA version`` and ``# LRAA CMD:``. That left the work tree as the sole
    provenance for a chunked deliverable -- so reclaiming that tree, which a
    benchmark harness does per row, made the shipped file unattributable to a
    build. Unchunked outputs never had this problem, which is how three different
    behaviours grew for what is one contract.

    Read from the UNITS rather than from this process, deliberately: the answer has
    to describe what actually produced the rows, not what is doing the merging.
    That also makes a mixed-version merge visible, and it is REFUSED rather than
    resolved, because a merged file naming one version while half its rows came
    from another is worse than a file naming none.

    A unit with no version line is refused for the same reason: the alternative is
    asserting a version this merge could not verify. Every current writer emits it
    (LRAA prints it on the first line of both artifacts), so its absence means a
    truncated or foreign file rather than an old one.

    Returned so that ``head -1`` answers "which build made this" the same way for
    chunked and unchunked outputs. Callers that also write their own provenance
    line put it AFTER these.
    """

    if not units:
        raise PipelineError("no units to merge; cannot establish provenance")

    seen = {}
    for unit in units:
        path = unit["quant_prefix"] + suffix
        version = None
        with open_text(path) as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                if line.startswith(VERSION_COMMENT_PREFIX):
                    version = line.rstrip("\n")
                    break
        if version is None:
            raise PipelineError(
                "unit {} carries no '{}' line in {}; refusing to write a merged "
                "header asserting a version this merge could not read".format(
                    unit["unit_id"], VERSION_COMMENT_PREFIX.strip(), path
                )
            )
        seen.setdefault(version, []).append(unit["unit_id"])

    if len(seen) > 1:
        detail = "; ".join(
            "{} from {} unit(s) e.g. {}".format(v, len(u), ", ".join(u[:3]))
            for v, u in sorted(seen.items())
        )
        raise PipelineError(
            "units disagree on LRAA version, refusing to merge them into one "
            "artifact: {}".format(detail)
        )

    # VERSION ONLY, no `# LRAA CMD:` line. A CMD line would have to come from this
    # process's argv, and that differs between the CLI entry point and a direct
    # in-process call -- which would break the byte-identity those two are tested
    # against (TestCliMatchesDirectCall). The merge's own shape is already recorded
    # by the merge-provenance line each caller writes after these, and the version
    # is the part that makes a shipped artifact attributable to a build.
    return [next(iter(seen))]


def merge_discovery_gtf(merged_dir, units):
    """Concatenate the per-chunk MODEL gtfs, back in the whole-run frame.

    Two rewrites, both mandatory.

    COORDINATES. The extractor rebases each chunk onto a mini contig that starts
    at 1 and is NAMED after the real contig, so a chunk's models carry chunk-local
    coordinates under a contig name that looks absolute. Columns 4 and 5 take the
    chunk offset; nothing else in a GTF line is a coordinate. Uncut contigs are
    exempt: ``coords_already_whole_contig`` means both columns already hold
    whole-contig coordinates and they are written through as read.

    IDS. LRAA names a model after its contig, its strand and a per-run component
    index -- ``t:chr1:+:comp-1:iso-1`` -- and every chunk of a contig shares the
    contig name and restarts the component counter at 1. Concatenating unpatched
    therefore FUSES unrelated models from different chunks into single records
    spanning the chromosome. That is not hypothetical: it happened, and produced
    37 spurious chromosome-crossing "models" that had to be diagnosed before any
    conclusion could be drawn from a chunked de novo comparison
    (FINDINGS.chr21_denovo_parity.md, section 2). Prefixing every gene_id and
    transcript_id with the unit id is the smallest change that cannot collide,
    and quant.expr and quant.tracking get the same prefix so the three artifacts
    name the same models.
    """

    gtf_out = os.path.join(merged_dir, "chunked.gtf")
    lines_written = 0
    transcripts = 0
    translate = not coords_already_whole_contig(units)
    # Version and CMD first, then this merge's own line: a consumer reading line 1
    # gets the same answer here as from an unchunked output.
    provenance = merged_provenance_header(units, ".gtf")

    with open(gtf_out, "wt") as ofh:
        for line in provenance:
            print(line, file=ofh)
        print(
            "# LRAA chunked discovery merge: {} unit(s); {}, model ids "
            "namespaced per unit".format(
                len(units),
                "coordinates translated to the whole-contig frame"
                if translate
                else "every unit spans its whole contig, so coordinates were "
                "already in the whole-contig frame and none were translated",
            ),
            file=ofh,
        )
        for unit in units:
            offset = unit["offset"]
            path = unit["quant_prefix"] + ".gtf"
            with open(path, "rt") as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 9:
                        raise PipelineError(
                            "unit {} gtf line carries {} field(s), not 9: "
                            "{}".format(
                                unit["unit_id"], len(fields), line.rstrip()
                            )
                        )
                    if translate:
                        fields[3] = str(int(fields[3]) + offset)
                        fields[4] = str(int(fields[4]) + offset)
                    fields[8] = _GTF_ID_ATTR.sub(
                        lambda m: '{} "{}"'.format(
                            m.group(1), _namespace_id(unit["unit_id"], m.group(2))
                        ),
                        fields[8],
                    )
                    print("\t".join(fields), file=ofh)
                    lines_written += 1
                    if fields[2] == "transcript":
                        transcripts += 1

    return {
        "gtf": gtf_out,
        "gtf_lines": lines_written,
        "gtf_transcripts": transcripts,
    }


# Rows held in one in-memory sorted run of the stage-6 tracking sort.
#
# WHY A CAP AT ALL: the merged quant.tracking is the largest table the pipeline
# produces -- one row per (read, assigned transcript) -- and materialising it to
# sort it was MEASURED at 588 B of Python object per 140 B of data, 4.2x, which
# put the serial merge at 41.6 GiB on a 63.8M-row PBMC arm and made stage 6 the
# peak of the whole run in all four measured arms. Capping the only container
# that grows with row count replaces that with a peak set by this constant.
#
# WHY 500k AND NOT SMALLER: peak resident rows is max(cap, rows/cap), so the cap
# IS the peak for any table below cap**2 -- 2.5e11 rows, against the 63.8M of
# the largest arm ever measured. Below about 20k the trade reverses: each
# spilled run is an open GzipFile during the merge, and MEASURED on 240k
# synthetic rows the allocation peak rose from 7.7 MiB at cap 20k to 9.5 MiB at
# cap 2k, the run buffer shrinking by 10x and the readers more than eating it.
#
# At 500k the same measurement gives 383 B per buffered row against the
# materialising implementation's 568 B, so ~0.19 GiB here against ~34 GiB there
# at PBMC arm B's row count. It is also FASTER: 29.4 s against 36.7 s on 2.4M
# rows, because five 500k sorts under `operator.itemgetter` beat one 2.4M sort
# under a key lambda by more than the spill round trip costs.
_TRACKING_SORT_RUN_ROWS = 500_000

# (sort_key, output_line) pairs are sorted on the key alone. NOT a bare
# `list.sort()`: that would fall through to comparing the line on equal keys and
# so reorder rows the materialising implementation left in input order.
_by_sort_key = operator.itemgetter(0)


def _read_sorted_run(path, key_columns, live):
    """Stream a spilled run back as (sort_key, line), counting what is resident.

    ``live`` is a two-slot list [current, peak]. ``heapq.merge`` holds exactly
    one pulled item per source at a time, so incrementing on the way out of the
    yield and decrementing when the source is resumed makes ``live[1]`` the
    MEASURED high-water mark of rows the merge held, rather than an assertion
    about heapq's internals. It is what the returned
    ``tracking_merge_peak_resident_rows`` reports.
    """

    read_i, tx_i, gene_i = key_columns
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            fields = line.split("\t")
            live[0] += 1
            if live[0] > live[1]:
                live[1] = live[0]
            yield (fields[read_i], fields[tx_i], fields[gene_i]), line
            live[0] -= 1


def _merge_tracking_streaming(units, hash_remap, discovery, track_out):
    """Stage 6's quant.tracking merge, as an external sort rather than a list.

    Byte-for-byte the previous implementation's output, reached differently:
    read every unit's rows applying the per-row transforms, buffer at most
    ``_TRACKING_SORT_RUN_ROWS`` of them, sort and spill each full buffer, then
    ``heapq.merge`` the spilled runs straight into the gzip writer.

    THE PREREQUISITE DOES NOT HOLD, which is why the sort is here rather than a
    plain k-way merge over the unit files. A per-unit quant.tracking is emitted
    transcript-major -- ``Quantify.report_quant_results`` walks transcripts by
    descending read count and only the read names WITHIN one multipath are
    sorted -- and LRAA's own shard merge concatenates those files without
    reordering (``LRAA:_append_without_first_line``). Checked on the 14
    strandless-parity gate chunks: 11 of the 14 are not in
    ``(read_name, transcript_id, gene_id)`` order, and the three that are have
    0, 1 and 2 rows. So the rows have to be sorted; only the container they are
    sorted in is negotiable.

    TRANSFORMS RUN ON THE WAY IN, not in a second pass, because the sort key is
    taken from the namespaced ids: discovery prefixes gene_id and transcript_id,
    two thirds of the key, so a merge over untransformed rows would order the
    table on strings that never get written.

    STABILITY, which byte-identity depends on. ``list.sort`` is stable, so the
    old implementation left equal-keyed rows in (unit order, within-unit file
    order). Reproduced by two facts: each run's buffer is sorted stably on the
    key alone, and ``heapq.merge`` breaks ties by source position. Runs are
    appended in read order, so source position IS (unit, offset in unit).

    Returns (track_header, n_rows, peak_resident_rows).
    """

    track_header = None
    tcol = None
    key_columns = None
    gene_i = tx_i = hash_i = read_i = None

    buf = []
    run_paths = []
    tmp_dir = None
    n_rows = 0
    buf_peak = 0

    def spill():
        nonlocal buf, tmp_dir
        if tmp_dir is None:
            # Beside the output rather than in /tmp: a node's scratch is often
            # far smaller than the space already provisioned for the merged
            # table, and these runs are the same order of magnitude as it.
            tmp_dir = tempfile.mkdtemp(
                dir=os.path.dirname(track_out), prefix="__tracking_sort."
            )
        buf.sort(key=_by_sort_key)
        path = os.path.join(tmp_dir, "run.{:05d}.tsv.gz".format(len(run_paths)))
        # compresslevel=1: written once, read once, deleted. The cheapest
        # setting that keeps the spill from costing more disk than the output.
        with gzip.open(path, "wt", compresslevel=1) as ofh:
            for _, line in buf:
                ofh.write(line)
                ofh.write("\n")
        run_paths.append(path)
        buf = []

    try:
        for unit in units:
            header, rows = iter_tsv(resolve_tracking(unit["quant_prefix"]))
            if track_header is None:
                track_header = header
                tcol = {name: i for i, name in enumerate(header)}
                gene_i = tcol["gene_id"]
                tx_i = tcol["transcript_id"]
                hash_i = tcol["transcript_splice_hash_code"]
                read_i = tcol["read_name"]
                key_columns = (read_i, tx_i, gene_i)
            elif header != track_header:
                raise PipelineError(
                    "unit {} quant.tracking header differs from the first "
                    "unit's".format(unit["unit_id"])
                )
            unit_id = unit["unit_id"]
            for row in rows:
                if discovery:
                    row[gene_i] = _namespace_id(unit_id, row[gene_i])
                    row[tx_i] = _namespace_id(unit_id, row[tx_i])
                old = row[hash_i]
                row[hash_i] = hash_remap.get((unit_id, old), old)
                buf.append(
                    ((row[read_i], row[tx_i], row[gene_i]), "\t".join(row))
                )
                n_rows += 1
                if len(buf) >= _TRACKING_SORT_RUN_ROWS:
                    buf_peak = _TRACKING_SORT_RUN_ROWS
                    spill()

        # The buffer only grows between spills, so its high-water mark is the
        # cap if it ever filled and the leftover otherwise. No per-row check.
        buf_peak = max(buf_peak, len(buf))

        live = [0, 0]
        if not run_paths:
            # Nothing spilled: the whole table fit inside one run, which is
            # already the bounded case. Sort it and write it, no temp files and
            # no merge -- the path every run smaller than the cap takes.
            buf.sort(key=_by_sort_key)
            merged = iter(buf)
        else:
            if buf:
                spill()
            merged = heapq.merge(
                *[_read_sorted_run(p, key_columns, live) for p in run_paths],
                key=_by_sort_key,
            )

        with gzip.open(track_out, "wt") as ofh:
            print("\t".join(track_header), file=ofh)
            for _, line in merged:
                print(line, file=ofh)
    finally:
        if tmp_dir is not None:
            shutil.rmtree(tmp_dir, ignore_errors=True)

    return track_header, n_rows, max(buf_peak, live[1])


def rebase_tpm(rows, header, total_override=None):
    """Rewrite TPM in place over ``rows``, byte-for-byte the arithmetic of
    ``LRAA:_merge_quant_expr_files``: float() the printed ``all_reads``, plain
    sequential sum, ``/total*1e6``, formatted ``"%.3f"``.

    ``rows`` are mutated in place (each a list of column strings matching
    ``header``) so a caller already holding references to them, as
    ``merge_and_translate`` does via ``expr_rows``, sees the update without a
    second pass.

    ``total_override``, given, is used INSTEAD of summing ``all_reads`` over
    ``rows``. This is what makes the function reusable for a merge that only
    ever sees part of the run: ``combine_grouped_expr`` below sums each GROUP's
    already-reported total rather than re-reading every group's rows to get
    one, and passes that sum in here. Returns the total actually used, so a
    caller that let this function compute it (the whole-run case) can record
    it exactly as before.
    """

    i_ar = header.index("all_reads")
    i_tpm = header.index("TPM")
    if total_override is not None:
        total = total_override
    else:
        total = 0.0
        for row in rows:
            total += float(row[i_ar] or 0)
    for row in rows:
        counts = float(row[i_ar] or 0)
        tpm = counts / total * 1e6 if total > 0 else 0
        row[i_tpm] = "{:.3f}".format(tpm)
    return total


def combine_grouped_expr(result_paths, out_path):
    """Recombine N grouped ``merge_and_translate`` results into ONE correct quant.expr.

    Grouping ``quant.tracking`` for a scattered merge (see
    ``util/misc/merge_chunk_outputs.py``'s ``--group``) is free: each group's
    output is already sorted on the global key, so a k-way merge of the N
    outputs reproduces the single global merge byte-for-byte with no re-sort.

    ``quant.expr`` does NOT get that for free, because its TPM is rebased over
    ``total_reported_read_count`` -- the SCOPE of whatever was merged in one
    call. A grouped call's scope is that group, not the run, so a grouped
    ``quant.expr``'s TPM column is correct only within its own group and wrong
    if read as a run-wide value.

    The fix needs no new arithmetic. `merge_and_translate` already reports
    ``merged_scope_total_all_reads`` for exactly this reason -- "so nothing is
    hidden" -- and already preserves each row's ``all_reads`` untouched (only
    ``TPM`` is overwritten). So: sum the N groups' already-reported totals to
    get the run's true total (proven exact -- summing partitions of a sum is
    associative), concatenate their quant.expr rows, and call ``rebase_tpm``
    ONE more time over the union with that sum as ``total_override``. VERIFIED
    byte-identical to a single global merge on the strandless-parity gate
    corpus, 555 of 555 rows, by ``test_combine_grouped_expr_matches_global_merge``.

    ``result_paths`` are the ``--result`` JSON files ``merge_chunk_outputs.py``
    writes per group; each carries ``quant_expr`` (a path) and
    ``merged_scope_total_all_reads`` (a float), which is everything this needs.
    """

    if not result_paths:
        raise PipelineError("no group results to combine; cannot establish provenance")

    header = None
    rows = []
    global_total = 0.0
    versions = {}
    for rp in result_paths:
        with open(rp, "rt") as fh:
            result = json.load(fh)
        global_total += float(result["merged_scope_total_all_reads"])
        # read_tsv rather than lines[0]/lines[1:]: a merged quant.expr opens with
        # the provenance comments this function has to carry forward, so the column
        # row is no longer guaranteed to be the first line.
        comments, this_header, these_rows = read_tsv(result["quant_expr"])
        # EVERY group must carry one, checked here per group rather than in
        # aggregate: a global "did any group have a version" test would let a group
        # with none inherit a sibling's, stamping the combined table with a version
        # that does not cover its rows -- the exact defect this header exists to
        # remove. Same contract as merged_provenance_header.
        version = None
        for comment in comments:
            if comment.startswith(VERSION_COMMENT_PREFIX):
                version = comment
                break
        if version is None:
            raise PipelineError(
                "group {} carries no '{}' line in {}; refusing to write a combined "
                "table asserting a version this combine could not read".format(
                    rp, VERSION_COMMENT_PREFIX.strip(), result["quant_expr"]
                )
            )
        versions.setdefault(version, []).append(rp)
        if header is None:
            header = this_header
        elif this_header != header:
            raise PipelineError(
                "{} quant.expr header differs from the first group's".format(rp)
            )
        rows.extend(these_rows)

    # Same contract as merge_discovery_gtf and the whole-run quant merge: the
    # combined table is the artifact that ships, so it names the build that made
    # its rows, and a disagreement is refused rather than resolved.
    if len(versions) > 1:
        raise PipelineError(
            "groups disagree on LRAA version, refusing to combine them into one "
            "artifact: {}".format(", ".join(sorted(versions)))
        )

    rebase_tpm(rows, header, total_override=global_total)

    with open(out_path, "wt") as ofh:
        print(next(iter(versions)), file=ofh)
        print(
            "# LRAA chunked quant combine: {} group(s); TPM rebased to the "
            "whole-run scope".format(len(result_paths)),
            file=ofh,
        )
        print("\t".join(header), file=ofh)
        for row in rows:
            print("\t".join(row), file=ofh)

    return {
        "quant_expr": out_path,
        "expr_rows": len(rows),
        "combined_from_groups": len(result_paths),
        "merged_scope_total_all_reads": global_total,
    }


# The counters a merged TOTAL row sums. Names match what LRAA's own
# read_assignment.summary.tsv writer emits.
_SUMMARY_TOTAL_KEYS = (
    "reads_total",
    "reads_kept_genome",
    "reads_selected_tx_total",
    "reads_selected_tx_missing_genome",
    "reads_selected_tx_failed_genome",
    "reads_rescue_requested",
    "reads_rescue_rescued",
    "reads_rescue_unrescued",
    "reads_rescue_requested_failed_genome",
    "reads_rescue_requested_unassigned_quant",
    "reads_rescue_declined_locality",
    "reads_rescue_displaced_locality",
    "alignments_rescue_rejected_locality",
)


def merge_read_assignment_summaries(merged_dir, units):
    """Every unit's read-assignment summary, as ONE table for the whole run.

    Nothing did this before. Chunking is unconditional, so the task-level
    ``<output_prefix>.read_assignment.summary.tsv`` that ``LRAA_runner.wdl``
    declares was never written by a chunked run: the per-unit files sat in the
    chunk directories and the declared output resolved to nothing. MEASURED on
    the single-cell fixtures: 931 per-unit summaries on disk, zero delocalized,
    and the workflow's own merge -- handed an empty list by ``select_all`` --
    emitted a correctly-schema'd table whose TOTAL row read 0 reads against a
    per-unit file reporting 11,768. Present, plausible and wrong, which is worse
    than absent.

    Input rows of ``row_type == "TOTAL"`` are SKIPPED when aggregating. Each
    unit's file is itself the output of a complete LRAA run and therefore already
    carries both a ``worker`` row and its own ``TOTAL``; summing every row would
    count each unit twice. The TOTAL here is recomputed from the ``worker`` rows --
    ENFORCED, not assumed: a row that is neither ``worker`` nor ``TOTAL`` is refused.
    Summing every non-TOTAL row would have let a third row type into the grand total
    and dropping every non-``worker`` row would have let one vanish from it, and this
    table IS the run's read accounting, so neither may happen quietly.

    ``rescue_alignment_rejections`` is a comma-separated ``key=count`` list, so
    it is summed BY KEY rather than concatenated -- concatenation would repeat
    the same reason once per unit and make the field unreadable.

    ``chunk_id`` is FILLED IN here rather than by the writer. LRAA emits the column
    empty (``LRAA._write_read_assignment_summary``) because a work unit knows only
    its contig and strand -- on a chunked single-contig run every worker row then
    reads identically (``worker chr19 +``) and the per-chunk detail, though present,
    is unusable. The value comes from the unit's ``chunk_id`` and NOT from
    ``unit_id``, which appends the orientation (``chr19_00_plus``) and so would name
    a strand rather than a cut. A unit record without one -- a manifest-driven stage
    6 whose manifest predates the field -- leaves the column empty rather than
    guessing at it. The TOTAL row's ``chunk_id`` is empty by construction: it spans
    every chunk, so no single chunk id is true of it.

    Every unit must contribute. A missing per-unit summary is refused rather than
    skipped -- see the comment on ``absent`` below for what that replaced -- so the
    merged ``reads_total`` is the whole run's read population and not an unstated
    subset of it.
    """

    summaries = [
        (u, "{}.read_assignment.summary.tsv".format(u["quant_prefix"])) for u in units
    ]
    present = [(u, p) for u, p in summaries if os.path.exists(p)]
    absent = [(u["unit_id"], p) for u, p in summaries if not os.path.exists(p)]

    # EVERY unit has to have one. LRAA writes a summary on every path now: the three
    # ``run_quant_only`` early returns that quantify nothing write their own, and a
    # backstop where a work unit FINISHES covers the rest -- discovery's nine early
    # returns and the oversimplify paths, which took no summary argument at all and
    # so left chrM out of the accounting of every run. Each counts the reads it SAW
    # and reports them with zero assigned
    # (LRAA._write_unquantified_read_assignment_summary). A missing file therefore
    # no longer means "this unit quantified nothing"; it means the unit did not
    # finish, and merging around it publishes a total short by its reads with
    # nothing in the table saying so.
    #
    # This replaces an evidence-based tolerance: accept a missing summary when the
    # unit's quant.expr held no data rows, refuse it otherwise. That rule existed
    # only to cover those early returns, and it covered them badly -- a shard taking
    # the third one can hold a bam full of reads while quantifying nothing, so the
    # tolerated case was exactly the undercounting one. With the returns writing,
    # the evidence is never needed, and keeping it would let a genuinely broken unit
    # through whenever its expr happened to be empty.
    if absent:
        raise PipelineError(
            "{} of {} quant unit(s) have no read-assignment summary, so a merged "
            "total would undercount by their reads: {}{}. Every unit writes one, "
            "including one that quantified nothing, so a missing file means that "
            "unit did not complete. Expected "
            "<quant_prefix>.read_assignment.summary.tsv".format(
                len(absent),
                len(summaries),
                "; ".join("{} ({})".format(uid, path) for uid, path in absent[:5]),
                "" if len(absent) <= 5 else "; ...",
            )
        )

    if not present:
        # No units at all: the standalone merge over an empty manifest. There is
        # nothing to describe, so nothing is published. Every other shape of
        # "nothing to merge" is refused above.
        return None

    rows = []
    in_fieldnames = None
    totals = {k: 0 for k in _SUMMARY_TOTAL_KEYS}
    rejections = collections.OrderedDict()
    for unit, path in present:
        unit_id = unit["unit_id"]
        # The CUT this unit's reads came from. Optional on the unit record so that a
        # manifest-driven stage 6 (util/misc/merge_chunk_outputs.py) still merges when
        # its manifest carries no chunk_id.
        #
        # ABSENT and FALSY are different questions and ``or ""`` answered the wrong
        # one: a real chunk id of 0 rendered as an empty column, indistinguishable
        # from the deliberate "this unit record carries no chunk id" case above and
        # from the TOTAL row, whose blank means "spans every chunk". Ids are strings
        # like "chr19_00" today so nothing reaches it; a sentinel confused with a
        # value in a provenance column is worth closing before something does.
        chunk_id_value = unit.get("chunk_id")
        chunk_id = "" if chunk_id_value is None else str(chunk_id_value)
        with open(path, "rt", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if in_fieldnames is None:
                in_fieldnames = list(reader.fieldnames or [])
            elif list(reader.fieldnames or []) != in_fieldnames:
                # Columns are summed by NAME, so a differing schema would drop or
                # misplace counters rather than fail. Refuse instead.
                raise PipelineError(
                    "read-assignment summaries disagree on their columns: {} has "
                    "{} field(s) against {} in the first unit merged".format(
                        unit_id, len(reader.fieldnames or []), len(in_fieldnames)
                    )
                )
            for row in reader:
                row_type = (row.get("row_type") or "").strip()
                if row_type == "TOTAL":
                    continue
                if row_type != "worker":
                    # LRAA writes exactly two row types -- "worker" (LRAA:5694) and
                    # "TOTAL" (LRAA:5948) -- so no in-tree writer reaches this. If a
                    # third ever appears, summing it into the run total and dropping
                    # it from the run total are both wrong answers that no column in
                    # the output would contradict, which is why neither is guessed.
                    raise PipelineError(
                        "read-assignment summary {} carries a row of row_type {!r}; "
                        "only 'worker' (summed into the merged TOTAL) and 'TOTAL' "
                        "(skipped) are known. Guessing would either double-count the "
                        "new row into the run total or silently drop it".format(
                            path, row_type
                        )
                    )
                row["chunk_id"] = chunk_id
                rows.append(row)
                for key in _SUMMARY_TOTAL_KEYS:
                    totals[key] += int(row.get(key) or 0)
                for item in (row.get("rescue_alignment_rejections") or "").split(","):
                    item = item.strip()
                    if not item or "=" not in item:
                        continue
                    reason, _, count = item.partition("=")
                    try:
                        rejections[reason] = rejections.get(reason, 0) + int(count)
                    except ValueError:
                        continue

    if not in_fieldnames:
        raise PipelineError(
            "every read-assignment summary stage 6 found was empty of rows: {}".format(
                ", ".join(p for _, p in present[:5])
            )
        )

    # chunk_id sits immediately after row_type. Inserted rather than assumed, so a
    # per-unit summary written by a build that predates the column still merges into
    # a table carrying it.
    fieldnames = list(in_fieldnames)
    if "chunk_id" not in fieldnames:
        fieldnames.insert(
            fieldnames.index("row_type") + 1 if "row_type" in fieldnames else 0,
            "chunk_id",
        )

    total_reads = totals["reads_total"]

    def frac(key):
        if total_reads <= 0:
            return "0.000000"
        return "{:.6f}".format(float(totals[key]) / float(total_reads))

    total_row = {f: "" for f in fieldnames}
    total_row.update(
        {"row_type": "TOTAL", "contig_acc": "TOTAL", "contig_strand": "."}
    )
    for key in _SUMMARY_TOTAL_KEYS:
        if key in total_row:
            total_row[key] = str(totals[key])
    # Fractions are recomputed against the merged denominator, never carried
    # over: a chunk-local fraction means nothing once the rows are pooled.
    for field in fieldnames:
        if not field.startswith("frac_"):
            continue
        counter = "reads_" + field[len("frac_") :]
        if counter in totals:
            total_row[field] = frac(counter)
        elif field == "frac_rescue_rescued_of_requested":
            req = totals["reads_rescue_requested"]
            total_row[field] = (
                "{:.6f}".format(float(totals["reads_rescue_rescued"]) / float(req))
                if req > 0
                else "0.000000"
            )
        elif field == "frac_rescue_unrescued_of_requested":
            req = totals["reads_rescue_requested"]
            total_row[field] = (
                "{:.6f}".format(float(totals["reads_rescue_unrescued"]) / float(req))
                if req > 0
                else "0.000000"
            )
    if "rescue_alignment_rejections" in total_row:
        total_row["rescue_alignment_rejections"] = ",".join(
            "{}={}".format(k, v) for k, v in sorted(rejections.items())
        )

    out_path = os.path.join(merged_dir, "chunked.read_assignment.summary.tsv")
    with open(out_path, "wt", newline="") as ofh:
        writer = csv.DictWriter(ofh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({f: row.get(f, "") for f in fieldnames})
        writer.writerow(total_row)
    # Path AND coverage. ``units_absent`` is 0 on every run that gets here, because
    # a unit without a summary is refused above -- the keys are kept because a
    # consumer reading this table still has to be able to see its own denominator's
    # coverage rather than assume it, and because that is what makes the refusal
    # visible as a refusal: the numbers say "every unit" instead of leaving the
    # question unasked.
    return {
        "path": out_path,
        "units_merged": len(present),
        "units_absent": len(absent),
        "units_total": len(summaries),
        "complete": not absent,
    }


def merge_and_translate(outdir, units, discovery=False):
    """Stage 6. Concatenate per-chunk quant output, back in the whole-run frame.

    Three fields carry chunk-local COORDINATES: ``exons`` and ``introns``, and the
    splice hash code, which is a blake2s digest OF the introns string and so has
    to be recomputed rather than shifted.

    NONE OF THAT APPLIES WHEN NO CONTIG WAS CUT. ``coords_already_whole_contig``
    asks the units, and all-zero offsets mean every unit's mini contig IS its
    contig: each shift would be ``+ 0``, and the digest would be a digest of an
    unchanged string. Both are skipped and the three fields pass through exactly
    as emitted, which is what makes the skip an identity rather than a shortcut.
    ``coordinates_translated`` in the returned dict says which path ran, and
    ``splice_hash_codes_recomputed`` is 0 in the skipping one because nothing was
    recomputed.

    One field carries a chunk-local DENOMINATOR: ``TPM`` is
    ``all_reads / (total all_reads emitted by the same job) * 1e6``
    (pylib/Quantify.py:2133), so a chunk's value is relative to that chunk. This
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

    THE TWO TABLES ARE MERGED DIFFERENTLY, deliberately. quant.expr is one row
    per model and fits anywhere, so it is sorted as a list here. quant.tracking
    is one row per (read, assigned transcript) and does not: it was the peak of
    the whole run in all four measured arms. It goes through
    ``_merge_tracking_streaming``, an external sort, and never exists as a list.

    Quantifying, transcript ids come from the supplied GTF and are already
    global, so no id translation is needed. DISCOVERING they are chunk-local and
    collide, so ``discovery`` additionally namespaces every gene_id and
    transcript_id by unit and merges the per-chunk model GTFs; see
    ``merge_discovery_gtf``. A quant-only run takes neither branch and emits
    exactly the bytes it emitted before.

    A quant UNIT rather than a chunk is what is concatenated here: strand-first
    they are the same thing, and strandless each chunk contributes two, one per
    orientation, from the one extraction. ``ordered_units`` puts them back into
    the strand-first order so the two modes' merged tables line up row for row.
    """

    merged_dir = os.path.join(outdir, "merged")
    os.makedirs(merged_dir, exist_ok=True)
    expr_out = os.path.join(merged_dir, "chunked.quant.expr")
    # Gzipped to match what LRAA itself emits per chunk, so a consumer opens the
    # merged table and a per-arm table the same way.
    track_out = os.path.join(merged_dir, "chunked.quant.tracking.gz")

    # PREFLIGHT, before a single merged byte is written: resolve provenance for
    # every artifact this call will produce. The gtf merge runs last, so leaving its
    # check to `merge_discovery_gtf` would refuse only AFTER quant.expr, tracking
    # and the read-assignment summary were already on disk -- and a partial merged
    # dir is exactly what a harness stages and ships. Discarding the gtf result is
    # deliberate: `merge_discovery_gtf` re-resolves it, staying self-sufficient for
    # its own callers, and the reread costs one leading comment line per unit.
    expr_provenance = merged_provenance_header(units, ".quant.expr")
    if discovery:
        merged_provenance_header(units, ".gtf")

    expr_header = None
    expr_rows = []
    hash_remap = {}  # (unit_id, old_hash) -> new_hash
    translate = not coords_already_whole_contig(units)
    for unit in units:
        offset = unit["offset"]
        _, header, rows = read_tsv(unit["quant_prefix"] + ".quant.expr")
        if expr_header is None:
            expr_header = header
        elif header != expr_header:
            raise PipelineError(
                "unit {} quant.expr header differs from the first unit's".format(
                    unit["unit_id"]
                )
            )
        col = {name: i for i, name in enumerate(header)}
        for row in rows:
            if len(row) != len(header):
                raise PipelineError(
                    "unit {} quant.expr row has {} fields, header has {}".format(
                        unit["unit_id"], len(row), len(header)
                    )
                )
            row = list(row)
            if discovery:
                # The model ids are this chunk's, and every chunk names a
                # comp-1. The merged GTF applies the identical prefix, so the
                # table and the models keep naming the same things.
                row[col["gene_id"]] = _namespace_id(
                    unit["unit_id"], row[col["gene_id"]]
                )
                row[col["transcript_id"]] = _namespace_id(
                    unit["unit_id"], row[col["transcript_id"]]
                )
            if translate:
                row[col["exons"]] = _shift_coord_string(row[col["exons"]], offset)
                introns = row[col["introns"]]
                if introns:
                    introns = _shift_coord_string(introns, offset)
                    row[col["introns"]] = introns
                    new_hash = Util_funcs.get_hash_code(introns)
                    hash_remap[(unit["unit_id"], row[col["splice_hash_code"]])] = (
                        new_hash
                    )
                    row[col["splice_hash_code"]] = new_hash
            expr_rows.append((unit["unit_id"], row))

    ecol = {name: i for i, name in enumerate(expr_header)}
    expr_rows.sort(
        key=lambda item: (item[1][ecol["gene_id"]], item[1][ecol["transcript_id"]])
    )

    # TPM rebase, byte-for-byte the arithmetic of LRAA:_merge_quant_expr_files.
    # Shared with `combine_grouped_expr` below, which is the SAME arithmetic
    # applied a second time over the union of grouped merges -- see that
    # function's docstring for why one recompute suffices for both cases.
    chunk_local_tpm = [row[ecol["TPM"]] for _, row in expr_rows]
    total_reported_read_count = rebase_tpm(
        [row for _, row in expr_rows], expr_header
    )

    with open(expr_out, "wt") as ofh:
        # This table used to start at the column row, carrying no provenance at all
        # while every per-chunk input carried both lines. Version first, so `head -1`
        # answers the same question here as on an unchunked quant.expr, then the
        # merge's own line for what the unchunked form has no equivalent of.
        for line in expr_provenance:
            print(line, file=ofh)
        print(
            "# LRAA chunked quant merge: {} unit(s); TPM rebased to the merged "
            "scope".format(len(units)),
            file=ofh,
        )
        print("\t".join(expr_header), file=ofh)
        for _, row in expr_rows:
            print("\t".join(row), file=ofh)

    # the as-emitted, chunk-relative TPM, kept so the rebase above can be audited
    tpm_audit = expr_out + ".tpm_chunk_local.tsv"
    with open(tpm_audit, "wt") as ofh:
        print(
            "\t".join(
                ["quant_unit", "gene_id", "transcript_id", "all_reads",
                 "TPM_chunk_local", "TPM_merged_scope"]
            ),
            file=ofh,
        )
        for (unit_id, row), local in zip(expr_rows, chunk_local_tpm):
            print(
                "\t".join(
                    [
                        unit_id,
                        row[ecol["gene_id"]],
                        row[ecol["transcript_id"]],
                        row[ecol["all_reads"]],
                        local,
                        row[ecol["TPM"]],
                    ]
                ),
                file=ofh,
            )

    _track_header, track_row_count, track_peak_rows = _merge_tracking_streaming(
        units, hash_remap, discovery, track_out
    )

    merged_summary = merge_read_assignment_summaries(merged_dir, units)

    merged = {
        "quant_expr": expr_out,
        "quant_tracking": track_out,
        "tpm_chunk_local_audit": tpm_audit,
        "expr_rows": len(expr_rows),
        "tracking_rows": track_row_count,
        "splice_hash_codes_recomputed": len(hash_remap),
        "merged_scope_total_all_reads": total_reported_read_count,
        "coordinates_translated": translate,
        # The high-water mark of tracking rows the merge held at once, MEASURED
        # rather than derived (see ``_read_sorted_run``). Reported because this
        # table is what set the run's memory ceiling, and a number that can be
        # read off a run is the only way to notice if it starts climbing again.
        "tracking_merge_peak_resident_rows": track_peak_rows,
        # Stays a PATH (or None), because that is what LRAA's publication step and
        # the workflow merge both consume. Coverage rides alongside rather than
        # inside it.
        "read_assignment_summary": (
            merged_summary["path"] if merged_summary else None
        ),
        # Whether that table covers every quant unit. It always does on a run that
        # completes -- a unit without a summary is refused now that every LRAA path
        # writes one -- so these are the numbers a consumer checks rather than a
        # tolerance it has to interpret. False here means no units at all.
        "read_assignment_summary_units_merged": (
            merged_summary["units_merged"] if merged_summary else 0
        ),
        "read_assignment_summary_units_absent": (
            merged_summary["units_absent"] if merged_summary else 0
        ),
        "read_assignment_summary_complete": (
            bool(merged_summary["complete"]) if merged_summary else False
        ),
    }
    if discovery:
        merged.update(merge_discovery_gtf(merged_dir, units))
    return merged


def shared_cut_positions(chunks):
    """The interior boundaries this partition actually has, per contig.

    Read off the chunk MANIFESTS rather than the selector's output, so what gets
    checked is the geometry the chunk directories were built with rather than the
    geometry something intended. A cut is the LAST base of a segment, which is
    what ``spanning_alignments`` means by ``position``: a read spans it when
    ``start <= position < end``.

    Strandless only. Strand-first has one segment series per contig-STRAND, so
    the boundaries are not a single per-contig sequence and a re-derivation there
    would have to name the orientation too -- unnecessary, because a distinct
    selection source is refused in that mode.
    """

    spans = {}
    for chunk in chunks:
        if not chunk["strandless"]:
            raise PipelineError(
                "per-input severed accounting is defined for strandless chunks "
                "only, and chunk {} carries orientation {!r}".format(
                    chunk["chunk_id"], chunk["strand"]
                )
            )
        manifest = chunk["manifest"]
        spans.setdefault(chunk["chrom"], set()).add(
            (manifest["partition_lend"], manifest["partition_rend"])
        )
    return {
        chrom: [rend for _lend, rend in sorted(ordered)[:-1]]
        for chrom, ordered in spans.items()
    }


def derived_severed_names(bam_path, cuts, max_intron_length):
    """Which reads of THIS bam the shared cuts sever, by the selector's own test.

    ``spanning_read_names`` is the selector's implementation, reused rather than
    copied, and that is the whole value of the check it feeds: the extractor's
    containment rule (``start < lend or end > rend``, applied per chunk) and the
    selector's span rule (``start <= position < end``, applied per cut) are two
    different statements of one geometry. Reimplementing either one here would
    reduce the check to the extractor agreeing with itself.

    Orientation-blind (``strand=""``), matching both how a strandless chunk is
    extracted and how a strandless selection is run.
    """

    selector = _selector_module()
    names = set()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom, positions in sorted(cuts.items()):
            for position in positions:
                names.update(
                    selector.spanning_read_names(
                        bam, chrom, "", position, max_intron_length
                    )
                )
    return names


def verify_severed_accounting(cut_dir, chunks, inputs=None, max_intron_length=None):
    """The set the control subtracts must BE the set the chunks dropped.

    WHAT IS COUNTED, exactly. ``run_baseline`` prunes the whole-contig bam by the
    read NAMES cut selection wrote to ``<cut_dir>/*.dropped_reads.txt``, so the
    two arms consume the same records only if those names are precisely the reads
    extraction then refused to place. Three quantities, deliberately kept apart:

      selected            names cut selection predicted it would sever. One file
                          per contig-STRAND strand-first; ONE FILE PER CONTIG
                          strandless, holding the reads a shared cut severs on
                          EITHER orientation, because a strandless cut is one
                          position serving both.
      dropped             names extraction actually refused, unioned over chunks.
                          A SUPERSET of ``selected``, not an equal: extraction also
                          refuses alignments reaching past the end of their contig,
                          which no cut is responsible for.
      mentions            the same drops counted per chunk rather than per read.
                          A read spanning a cut overlaps BOTH neighbouring
                          chunks and is listed by both, so this is about twice
                          ``dropped`` and is not a read count. It is reported so
                          that nobody reaches for the per-chunk sum -- which is
                          what ``alignments_dropped_overhang`` totals to -- as if
                          it were the severed set.

    The opening statement is still the contract, but it is enforced by making the
    control subtract more rather than by refusing the run: reads dropped without
    being named are written to ``EXTRACTION_ONLY_DROPS`` in the same directory, which
    the baseline's own glob then picks up. Only the opposite direction -- named and
    NOT dropped -- remains fatal, because that one needs the two tools to disagree
    about identical geometry and no bookkeeping can decide which is right.

    Strand-first never had to check this: one selection and one chunk series per
    orientation meant the same coordinates decided both sets twice over.
    Strandless breaks that symmetry -- the selection is per contig and the drops
    are per orientation within a shared chunk -- so the identity is verified on
    every run instead of assumed, in both modes, because a check that only runs
    in the new mode proves nothing about the old one.

    ``inputs`` IS THE SHARED-CUTS CASE and changes which equality is asserted.
    With ``--bam_for_cut_selection`` the cuts are chosen on a SUPERSET of the bam
    being extracted, so a name selection wrote can simply be absent from this
    caller's reads -- legitimate, and the one-source equality above would fail on
    correct geometry. The equality is re-established PER INPUT instead: each
    input's own severed set, re-derived at these same cuts, against that input's
    own drops. Both directions are exact again, and a genuine geometry error --
    a chunk boundary that does not sit at the cut its reads were selected against
    -- still fails, because no bam can drop fewer reads than the cuts sever.
    Each entry is ``{"label", "bam", "names_key"}``, ``names_key`` naming the
    manifest field holding that input's drops (``dropped_read_names`` for the
    reads, ``sg_dropped_read_names`` for the splice-graph slice).
    """

    selected = severed_read_names(cut_dir, exclude=(EXTRACTION_ONLY_DROPS,))
    dropped = set()
    mentions = 0
    for chunk in chunks:
        names = chunk["manifest"]["dropped_read_names"]
        mentions += len(names)
        dropped.update(names)

    per_input = {}
    if inputs:
        cuts = shared_cut_positions(chunks)
        for spec in inputs:
            derived = derived_severed_names(spec["bam"], cuts, max_intron_length)
            actual = set()
            for chunk in chunks:
                actual.update(chunk["manifest"].get(spec["names_key"]) or [])
            missed = sorted(derived - actual)
            if missed:
                raise PipelineError(
                    "severed-read accounting is inexact for the {} input ({}): at "
                    "the {} shared cut position(s) this partition has, {} "
                    "alignment(s) span a cut and extraction placed {} of them "
                    "whole (e.g. {}). The cut positions and the chunk boundaries "
                    "disagree, so a chunk holds a record crossing its own "
                    "boundary.".format(
                        spec["label"],
                        spec["bam"],
                        sum(len(v) for v in cuts.values()),
                        len(derived),
                        len(missed),
                        ", ".join(missed[:5]),
                    )
                )
            per_input[spec["label"]] = {
                "bam": spec["bam"],
                "severed_at_shared_cuts": len(derived),
                "dropped_by_extraction": len(actual),
                # Drops no cut is responsible for: an alignment ending past its
                # contig. Reported per input rather than folded into the run-level
                # number below, because for the sg slice there is no baseline arm
                # to subtract them from and the count is the only record of them.
                "dropped_beyond_cuts": len(actual - derived),
            }

    # Deliberately NOT a place where a nonzero drop fails the run, in either mode.
    # An earlier revision raised here whenever discovery dropped anything, because
    # selection had promised zero by refusing every severing position. That promise
    # is gone: severing is priced, not forbidden, so a nonzero drop is the expected
    # outcome on any dense input.
    #
    # The two directions of a mismatch are NOT the same defect, and only one of them
    # is about the tools disagreeing.
    #
    # NAMED BUT NOT DROPPED is fatal. Selection promised a cut would sever a read and
    # extraction placed it whole, so the two disagree about the same geometry. The
    # baseline would then lose records the chunks did see, and no bookkeeping here can
    # reconstruct which side is right.
    unrealized = sorted(selected - dropped)
    if unrealized and not inputs:
        raise PipelineError(
            "severed-read accounting is inexact: cut selection named {} read(s) as "
            "severed and extraction placed {} of them whole (e.g. {}). The two tools "
            "disagree about the same cut geometry, so the parity comparison would "
            "lose records the chunks did see.".format(
                len(selected), len(unrealized), ", ".join(unrealized[:5])
            )
        )

    # DROPPED BUT NOT NAMED is tolerated, because it does not require the tools to
    # disagree about anything. An alignment whose end lies past the end of its contig
    # overhangs the final chunk with no cut responsible for it, and coordinate-remapped
    # bams carry these routinely -- 453 of 2,507 records in testing/sep_contigs do,
    # one ending at 199,427 on a 77,313 bp contig. Refusing the run over a defect in
    # the INPUT would make such a bam unchunkable for a reason that has nothing to do
    # with chunking.
    #
    # Tolerated is not ignored. These reads are absent from the chunk inputs, so the
    # baseline must subtract them or the arms stop consuming the same records and the
    # parity comparison silently compares different sets. Persisting them where
    # ``severed_read_names`` looks is what keeps that true across separate
    # invocations of the two arms, which is the whole reason the baseline reads this
    # directory from disk rather than being handed a set in memory.
    extraction_only = sorted(dropped - selected)
    if extraction_only:
        path = os.path.join(cut_dir, EXTRACTION_ONLY_DROPS)
        with open(path, "wt") as fh:
            for name in extraction_only:
                fh.write(name + "\n")
        print(
            "NOTE: extraction dropped {} read(s) no cut selection named, which is "
            "what an alignment reaching past the end of its contig looks like. They "
            "are named in {} so the baseline subtracts them too; the run continues. "
            "e.g. {}".format(
                len(extraction_only), path, ", ".join(extraction_only[:3])
            ),
            flush=True,
        )

    return {
        "severed_reads": len(selected),
        "named_by_cut_selection": len(selected),
        "dropped_by_extraction": len(dropped),
        "per_chunk_drop_mentions": mentions,
        "dropped_not_named": len(extraction_only),
        "sets_identical": not extraction_only,
        # Nonzero only with a shared selection source, where it is expected: the
        # selector saw reads of every caller and this one holds a subset. Fatal
        # without one, and raised above rather than reported.
        "named_absent_from_extraction_input": len(unrealized),
        "shared_cut_selection": bool(inputs),
        "per_input": per_input,
    }


def roll_up_split_accounting(chunks):
    """Per-chunk split counts, summed, for a strandless run's timing record.

    ``verify_chunk_split`` already refused anything inconsistent, chunk by chunk,
    inside the parallel phase. This is the run-level statement of what survived:
    the number of alignment records the chunked arm quantifies, and the number
    the in-chunk splits lost on the way, which is zero or the run failed.
    """

    total = collections.Counter()
    for chunk in chunks:
        for key, value in chunk.get("split_counts", {}).items():
            total[key] += value
    rolled = {
        "intervals_split": sum(1 for c in chunks if "split_counts" in c),
        "alignments_emitted": total["alignments_emitted"],
        "records_plus": total["records_plus"],
        "records_minus": total["records_minus"],
        "records_quantified": total["records_total"],
        "records_lost_in_split": total["records_lost"],
    }
    # Only when the run had splice-graph evidence to split. Reported at all
    # because evidence lost in the split is not visible in any quant number: the
    # graph is simply sparser, so the run succeeds with fewer junctions.
    # ``in`` rather than indexing, because a Counter would answer 0 and a zero
    # here would read as "no evidence records" instead of "no evidence".
    if "sg_alignments_emitted" in total:
        rolled.update(
            {
                "sg_alignments_emitted": total["sg_alignments_emitted"],
                "sg_records_plus": total["sg_records_plus"],
                "sg_records_minus": total["sg_records_minus"],
                "sg_records_quantified": total["sg_records_total"],
                "sg_records_lost_in_split": total["sg_records_lost"],
            }
        )
    return rolled


# -------------------------------------------------------------- baseline arm


def run_baseline(
    args,
    ckpt,
    outdir,
    timing,
    strand_bams,
    num_total_reads,
    rss_interval,
    split_token,
    severed_names=None,
):
    """Whole-contig control. Same substrate, same settings, no partitioning.

    ``--bam`` is the two orientation bams glued back together, MINUS the reads a
    cut severs, which is exactly the union of every chunk's mini bam expressed in
    genome coordinates -- so the arms differ in partitioning alone rather than
    also in which records reach LRAA.

    The subtraction is the point.  A severed read is dropped from both
    neighbouring chunks and is an accepted loss, but the whole-contig arm has no
    cut to lose it at, so leaving it in makes the control read a record set the
    chunked arm never saw and any difference between the two is confounded by
    exactly those reads.  This claimed the arms already matched while doing
    nothing to make them: on a 3.4 Mb contig cut seven ways it quantified 2,266
    records against the chunked arm's 2,261.

    ``severed_names`` absent means no cut selection has run in this output
    directory, so there is no chunked arm to be comparable with and the whole
    input is the right input.
    """

    bdir = os.path.join(outdir, "baseline")
    os.makedirs(bdir, exist_ok=True)
    whole_bam = os.path.join(bdir, "whole.primary.bam")
    log = os.path.join(outdir, "logs", "baseline.log")
    started = time.time()
    steps = []

    merge_token = chain_token(
        "baseline_merge.maxintron_{}".format(args.max_intron_length), split_token
    )
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

    quant_bam = whole_bam
    input_token = merge_token
    if severed_names:
        quant_bam = os.path.join(bdir, "whole.parity.bam")
        # The names decide the contents, so they name the artifact.  A digest of
        # the sorted set rather than the set itself: it is unbounded, and the
        # count alone would collide between two runs severing different reads.
        input_token = chain_token(
            "baseline_parity.severed_{}".format(len(severed_names)),
            merge_token,
            Util_funcs.get_hash_code("|".join(sorted(severed_names))),
        )
        prune_step = {
            "step": "baseline_parity_prune",
            "excluded_read_names": len(severed_names),
            "output": quant_bam,
        }
        if ckpt.done(input_token):
            prune_step["reused"] = True
        else:
            # This step runs IN-PROCESS -- write_bam_excluding is a call, not a
            # subprocess -- so os.wait4 has no child to report and it is the one
            # step run_step cannot instrument. Measure the driver's own CPU across
            # the span instead, and label the source, because a self-rusage delta
            # and a child rusage total are not interchangeable: this one includes
            # any other work this process does concurrently, which here is none.
            _ru0 = resource.getrusage(resource.RUSAGE_SELF)
            _t0 = time.time()
            prune_step["kept_records"] = write_bam_excluding(
                whole_bam, severed_names, quant_bam
            )
            _ru1 = resource.getrusage(resource.RUSAGE_SELF)
            prune_step["wall_s"] = round(time.time() - _t0, 3)
            prune_step["cpu_s"] = round(
                (_ru1.ru_utime - _ru0.ru_utime) + (_ru1.ru_stime - _ru0.ru_stime), 3
            )
            prune_step["cpu_source"] = "rusage_self_delta"
            prune_step["max_rss_kb"] = _ru1.ru_maxrss
            ckpt.mark(input_token)
        steps.append(prune_step)

    norm_bam = os.path.join(bdir, "whole.norm.bam")
    # Built against the WHOLE genome, not per chunk -- this arm has no chunks --
    # but from the same cassette and padding, so a locus masked in the chunked
    # arm is masked here too. Anything else would make the two arms disagree for
    # a reason that has nothing to do with partitioning, which is the one axis
    # this control exists to isolate.
    baseline_rdna_mask_bed = None
    if not getattr(args, "no_rdna_mask", False):
        baseline_rdna_fasta = RdnaMask.resolve_rdna_fasta(
            getattr(args, "rdna_mask_fasta", None)
        )
        baseline_rdna_mask_bed = RdnaMask.build_rdna_mask_bed(
            os.path.abspath(args.genome_fa),
            baseline_rdna_fasta,
            cache_dir=os.path.join(bdir, "__rdna_mask_cache"),
        )
    norm_token = chain_token(
        # maxintron is passed to the normalizer below; it was missing here, which
        # is what let a new cap rerun the merge and then reuse this bam.
        # min_mapping_quality was missing the same way: the chunk arm's own
        # stage-4 token gained it for the identical reason (see chunk_worker),
        # and the baseline control has to name the same inputs or a non-default
        # value would make the two arms disagree without either being wrong.
        "baseline_norm.maxcov_{}_dw_{}_seed_{}_origin_0_maxintron_{}_minmapq_{}"
        "_rdnamask_{}".format(
            args.normalize_max_cov_level,
            args.depth_window,
            args.random_seed,
            args.max_intron_length,
            resolve_min_mapping_quality(args),
            Util_funcs.file_identity_token(baseline_rdna_mask_bed)
            if baseline_rdna_mask_bed
            else "none",
        ),
        input_token,
    )
    norm_cmd = [
        sys.executable,
        NORMALIZE_BAM,
        "--input_bam",
        quant_bam,
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
        "--min_mapping_quality",
        str(resolve_min_mapping_quality(args)),
        # THE control's grid must be pinned to absolute 0. The default anchors on
        # the first aligned base per contig and no chunk grid can match it.
        "--window_origin",
        "0",
    ]
    if baseline_rdna_mask_bed:
        norm_cmd += ["--rdna_mask_bed", baseline_rdna_mask_bed]
    if ckpt.done(norm_token):
        steps.append({"step": "baseline_normalize", "reused": True, "cmd": norm_cmd})
    else:
        steps.append(run_step("baseline_normalize", norm_cmd, log, bdir, rss_interval))
        ckpt.mark(norm_token)

    # bare prefix, cwd=bdir: see the note in chunk_worker about LRAA's scratch
    # roots being string concatenations on --output_prefix
    quant_prefix = os.path.join(bdir, "baseline_quant." + LRAA_QUANT_ONLY_SUFFIX)
    quant_cmd = lraa_cmd(
        args,
        bam_for_quant=quant_bam,
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
    # AFTER the --contig append, and after the command is complete: the sentinel
    # is keyed on the argv this arm actually runs, so the restriction is part of
    # the key rather than a setting the old hand-listed name happened to omit.
    # Same defect and same fix as the chunk arm's stage-5 token; leaving either
    # one enumerated leaves the whole class reachable.
    quant_token = chain_token(
        "baseline_quant.argv_{}".format(quant_command_digest(quant_cmd)[:12]), norm_token
    )
    if ckpt.done(quant_token):
        steps.append({"step": "baseline_quant", "reused": True, "cmd": quant_cmd})
    else:
        steps.append(run_step("baseline_quant", quant_cmd, log, bdir, rss_interval))
        ckpt.mark(quant_token)

    expr_path = quant_prefix + ".quant.expr"
    if not os.path.exists(expr_path):
        raise PipelineError(
            "baseline produced no {}; log: {}".format(expr_path, log)
        )
    baseline_tracking = resolve_tracking(quant_prefix)
    if baseline_tracking is None:
        raise PipelineError(
            "baseline produced no {}{{{}}}; log: {}".format(
                quant_prefix, ",".join(QUANT_TRACKING_SUFFIXES), log
            )
        )

    timing.setdefault("arms", {})["baseline"] = {
        "wall_s": round(time.time() - started, 3),
        "peak_tree_rss_kb": max([s.get("peak_tree_rss_kb", 0) for s in steps] or [0]),
        "steps": steps,
        "log": log,
    }
    return {
        "quant_expr": quant_prefix + ".quant.expr",
        "quant_tracking": baseline_tracking,
    }


# ---------------------------------------------------------------------- driver


def loadavg():
    with open("/proc/loadavg", "rt") as fh:
        return [float(x) for x in fh.read().split()[:3]]


def build_parser():
    """The chunked pipeline's own CLI.

    Kept as a parser rather than folded into LRAA's so that the parity harness
    and ``LRAA --chunk`` cannot drift apart on a default: both end up calling
    ``run`` with a namespace this produced.
    """

    parser = argparse.ArgumentParser(
        description="run the chunked quant-only pipeline (stages 1-6) or the "
        "whole-contig control, with per-chunk logs, wall time and peak RSS. "
        "Resumable via sentinel files under <output_dir>/__ckpt/.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # Required on the command line, but declared optional so that a fully
    # defaulted namespace can be built for `LRAA --chunk` to fill in.
    # `parse_args` enforces the requirement for CLI callers, `run` for both.
    parser.add_argument("--bam", default=None, help="input bam of aligned reads")
    parser.add_argument(
        "--bam_for_sg",
        default=None,
        help="OPTIONAL pre-normalized splice-graph evidence, supplied by the "
        "caller and partitioned per chunk so every unit sees the same evidence "
        "the unchunked run would have seen. --bam remains the READS: they are "
        "what is counted and what sets the TPM denominator. This exists for the "
        "single-cell shape where ~29 per-cluster quant jobs must share ONE splice "
        "graph (WDL/LRAA_quant_by_cluster.wdl); chunking used to manufacture its "
        "own sg bam per unit, which silently destroyed that. Stage 4 is skipped "
        "for a unit that has one and stage 5 gets --bam_for_sg <slice> --no_norm. "
        "REQUIRES --bam_for_cut_selection: without it the cuts come from --bam, "
        "so every caller slices the shared evidence at its own boundaries",
    )
    parser.add_argument(
        "--bam_for_priors",
        default=None,
        help="OPTIONAL pre-normalized bam the per-chunk LRAA estimates its "
        "FIRST-PASS theta from, supplied by the caller and partitioned per chunk "
        "on this run's own geometry exactly as --bam_for_sg is. Three roles, "
        "three files: --bam_for_sg is ONE splice graph shared across cell "
        "clusters, this is THIS cluster's own normalized reads, and --bam is the "
        "full cluster bam the second pass assigns. Without it a cluster's theta "
        "is estimated over the SHARED sg slice, so every cluster's ambiguous-read "
        "apportionment depends on every other cluster's expression -- pooled "
        "priors reported as cluster-local quant. Stage 4 is skipped for a unit "
        "that has one (it arrived normalized) and stage 5 gets "
        "--bam_for_priors <slice>. Strandless chunking only, because stage 3b, "
        "which splits it, runs only for a strandless chunk",
    )
    parser.add_argument(
        "--bam_for_cut_selection",
        default=None,
        help="the bam CUT POSITIONS are chosen from. Default --bam, which is "
        "today's behaviour exactly. It must be an unthinned SUPERSET of every "
        "--bam this plan is applied to -- in practice the pre-partition input the "
        "per-caller bams were split out of. A cut safe in the superset is safe in "
        "every subset, so every caller gets identical geometry AND none of them "
        "loses a read. Strandless chunking only",
    )
    parser.add_argument("--genome_fa", default=None, help="genome fasta")
    parser.add_argument(
        "--gtf",
        default=None,
        help="reference annotation gtf. Required for quant-only; OPTIONAL with "
        "--discovery, where absent means de novo and present means ref-guided "
        "discovery",
    )
    parser.add_argument(
        "--discovery",
        action="store_true",
        help="run isoform DISCOVERY per chunk instead of quantification against "
        "a supplied annotation: stage 5 drops --quant_only and stage 6 also "
        "merges the per-chunk model gtfs. Cut selection is IDENTICAL to "
        "quant-only's -- severing is a cost the selector minimises, weighted by "
        "--severed_multiexon_weight, and a severed read is dropped, counted and "
        "named in both modes",
    )
    parser.add_argument(
        "--output_dir", default=None, help="output directory; created if absent"
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
        "count is an error, not a preference. REQUIRED with "
        "--strandless_chunks --arm chunked, which does not run stage 1 and so "
        "has no count to default to",
    )
    parser.add_argument("--contig", default=None, help="restrict to one contig")
    parser.add_argument(
        "--oversimplify",
        # MUST match LRAA's own --oversimplify default ("chrM,M"). A chunk worker
        # that receives no flag falls back to that default and applies it to its
        # MINI contig, whose name is not chrM -- so leaving this None would turn
        # the default off in chunked mode only, for the same command. Carrying the
        # same value here lets resolve_oversimplify_for_chunk translate it into the
        # mini name the worker will actually match.
        default="chrM,M",
        metavar="CHR1,CHR2,...",
        help="run these contigs in LRAA's oversimplify mode, named as the genome "
        "fasta spells them. Resolved PER CHUNK by "
        "resolve_oversimplify_for_chunk: a chunk whose contig is named here passes "
        "--oversimplify to its stage-5 LRAA run, spelled as the mini contig the "
        "extractor renamed it to, because the name written here never appears in a "
        "chunk's own fasta",
    )
    parser.add_argument(
        "--contigs",
        default=None,
        metavar="CHR1,CHR2,...",
        help="restrict the run to these contigs, comma-separated. The partition "
        "is otherwise every reference the genome fasta and the bam header agree "
        "on; this filters that set, which is how a workflow's main-chromosomes "
        "list reaches a chunk partition. A name absent from both the fasta and "
        "the bam header is a typo and is refused rather than silently dropped. "
        "Mutually exclusive with --contig, which restricts to exactly one",
    )
    parser.add_argument("--HiFi", action="store_true", help="pass --HiFi to LRAA")
    parser.add_argument(
        "--rdna_mask_fasta",
        type=str,
        default=LRAA_Globals.config["rdna_mask_fasta"],
        help="passed identically to LRAA's own --rdna_mask_fasta for every chunk",
    )
    parser.add_argument(
        "--no_rdna_mask",
        action="store_true",
        default=not LRAA_Globals.config["rdna_mask_enabled"],
        help="passed identically to LRAA's own --no_rdna_mask for every chunk",
    )
    # DEFAULT ON, and it must agree with LRAA's own default -- default_args() is
    # documented as the single place both routes get their defaults from, so a
    # disagreement here would make `LRAA --chunk` and this driver run different
    # pipelines from the same arguments. LRAA_Globals.config["strandless_chunks"]
    # is the value; it is not read directly here because this parser is usable
    # without that module loaded.
    parser.add_argument(
        "--strandless_chunks",
        dest="strandless_chunks",
        action="store_true",
        default=True,
        help="cut and extract STRANDLESS chunks and run the orientation split "
        "inside each chunk, concurrently with every other chunk, instead of "
        "splitting the whole bam up front. THE DEFAULT; the flag is accepted so "
        "a command can say so explicitly. The chunked arm then skips stage 1 "
        "entirely and extracts once per interval rather than once per "
        "contig-strand. The control still needs the split -- it IS the "
        "strand-split whole bam -- so --arm baseline/both still runs stage 1. "
        "An output directory serves one mode or the other, never both",
    )
    parser.add_argument(
        "--chunk_by_strand",
        dest="strandless_chunks",
        action="store_false",
        help="split the whole bam by orientation FIRST, then cut and process each "
        "contig-STRAND. The pre-default ordering, kept so a regression can be "
        "bisected against it. Must stay spelled the same as LRAA's flag: the two "
        "parsers are one interface as far as anyone reading a command is concerned",
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="print the plan and stop. Cut selection IS the plan, so it runs "
        "(and is checkpointed for the real run); nothing is extracted, "
        "normalized or quantified. Strand-first also needs stage 1 first, "
        "because its cut selection reads the split bams -- which is itself the "
        "difference the strandless mode exists to remove",
    )
    parser.add_argument(
        "--no_reuse_source_bam",
        action="store_true",
        help="always extract a mini bam, even for a strandless chunk spanning its "
        "whole contig. The reuse optimization names the SOURCE bam in the chunk "
        "manifest, which a chunk processed on another machine cannot open; this "
        "makes every chunk directory self-contained at the cost of the copy "
        "reuse exists to avoid",
    )
    parser.add_argument(
        "--stop_after_make_chunks",
        action="store_true",
        help="run cut selection and chunk extraction, write <output_dir>/chunk_plan.json, "
        "and stop before stages 3b-6. For a scattered workflow that processes each "
        "chunk as its own task; see --only_chunk",
    )
    parser.add_argument(
        "--only_chunk",
        default=None,
        metavar="CHUNK_ID",
        help="process ONE already-extracted chunk (stages 3b-5) and stop, writing "
        "<output_dir>/chunks/<CHUNK_ID>/units.json. Requires --output_dir holding "
        "chunk_plan.json and chunks/<CHUNK_ID>/ from a --stop_after_make_chunks run. "
        "Skips stages 1-3 entirely, so --bam, --genome_fa and --gtf are not needed",
    )
    parser.add_argument(
        "--emit_cut_plan",
        default=None,
        metavar="PATH",
        help="run stage 2 cut selection ONLY, write the partition to PATH and "
        "stop. Nothing is extracted -- which is the point, and the whole "
        "difference from --stop_after_make_chunks: that flag extracts every chunk "
        "of the bam it is given (~300 of them on a whole genome), and here the "
        "extraction belongs to the callers. Emit ONE plan from the pre-partition "
        "bam, then hand it to every per-cluster run with --chunk_plan: selection "
        "is paid once, extraction once per caller, and all of them cut at "
        "identical positions. GEOMETRY ONLY -- no -N is needed or recorded as "
        "authoritative, because the TPM denominator belongs to whoever quantifies",
    )
    parser.add_argument(
        "--chunk_plan",
        default=None,
        metavar="PATH",
        help="apply the geometry in a plan from --emit_cut_plan instead of "
        "selecting cuts. Extraction runs against THIS run's --bam and "
        "--bam_for_sg, and this run's own -N is the TPM denominator; the plan "
        "supplies nothing but where the boundaries are. Refused if any parameter "
        "that decides a cut disagrees with the plan's, or if a contig this run "
        "processes is missing from the plan or differs in length -- a plan "
        "applied at other geometry is a different partition, not a smaller "
        "answer. Satisfies --bam_for_sg's shared-geometry requirement, so "
        "--bam_for_cut_selection is neither needed nor accepted with it",
    )
    # Defaults come from the canonical definition rather than a second copy. These
    # were hardcoded here and happened to agree with LRAA_Globals; a change to
    # either would have diverged silently, and worse than silently, because the
    # values are baked into the stage-2 cache token above -- a run would have kept
    # the old cut geometry while the token asserted the new one.
    parser.add_argument(
        "--approx_MB_per_cut",
        type=float,
        default=LRAA_Globals.config["approx_MB_per_cut"],
    )
    parser.add_argument(
        "--approx_MB_per_cut_wiggle_window",
        type=float,
        # The MAXIMUM search window, and an ABSOLUTE default: not derived from
        # --approx_MB_per_cut. The config comment on the key carries the
        # measurement that rules out a proportional form.
        default=LRAA_Globals.config["approx_MB_per_cut_wiggle_window"],
    )
    parser.add_argument(
        "--severed_multiexon_weight",
        type=int,
        default=LRAA_Globals.config["chunk_severed_multiexon_weight"],
        help="what a severed MULTI-EXON alignment costs cut selection, against 1 "
        "for a monoexonic one. 1 makes the cost a plain alignment count",
    )
    parser.add_argument(
        "--depth_window",
        type=int,
        default=LRAA_Globals.config["chunk_depth_window"],
    )
    parser.add_argument(
        "--margin", type=int, default=LRAA_Globals.config["chunk_margin"]
    )
    parser.add_argument(
        "--max_intron_length",
        type=int,
        default=LRAA_Globals.config["max_intron_length"],
    )
    parser.add_argument(
        "--normalize_max_cov_level",
        type=int,
        default=LRAA_Globals.config["normalize_max_cov_level"],
    )
    parser.add_argument(
        "--random_seed",
        type=int,
        default=LRAA_Globals.config["chunk_random_seed"],
    )
    # Mirrors LRAA's own flags exactly -- same names, same SUPPRESS-vs-real-default
    # split -- because `lraa_cmd` forwards them to a per-chunk LRAA invocation
    # unchanged, and `default_args`'s override validation and the presence checks
    # below both depend on that. See LRAA_Globals.py and LRAA's "Single cell
    # settings" / "--stream_reads" argument groups for the values these mirror.
    parser.add_argument(
        "--cell_barcode_tag",
        type=str,
        default=argparse.SUPPRESS,
        help="bam tag for cell barcode (default: {})".format(
            LRAA_Globals.config["cell_barcode_tag"]
        ),
    )
    parser.add_argument(
        "--read_umi_tag",
        type=str,
        default=argparse.SUPPRESS,
        help="bam tag for read umi (default: {})".format(
            LRAA_Globals.config["read_umi_tag"]
        ),
    )
    parser.add_argument(
        "--cell_list",
        type=str,
        default=None,
        help="file of cell barcodes considered real cells, one per line. "
        "Resolved to an absolute path immediately, because a chunk worker's "
        "cwd is its own chunk directory, not the caller's",
    )
    parser.add_argument(
        "--stream_reads",
        dest="stream_reads",
        action="store_true",
        default=LRAA_Globals.config["stream_reads"],
        help="pass --stream_reads to every chunk worker, quant-only or discovery "
        "alike. Each worker is a full LRAA invocation and validates the combination "
        "itself (a thinner first-pass bam, rescue settings, --tag_bam), exactly as "
        "the unchunked path would. ON BY DEFAULT, matching LRAA's own default; see "
        "--no_stream_reads",
    )
    parser.add_argument(
        "--no_stream_reads",
        dest="stream_reads",
        action="store_false",
        help="pass --no_stream_reads to every chunk worker instead",
    )
    parser.add_argument(
        "--rescue_unassigned_reads_via_transcriptome_alignment",
        dest="rescue_unassigned_reads_via_transcriptome_alignment",
        action="store_true",
        default=LRAA_Globals.config[
            "rescue_unassigned_reads_via_transcriptome_alignment"
        ],
        help="pass transcriptome rescue on to every chunk worker. ON BY DEFAULT, "
        "matching LRAA's own default; see "
        "--no_rescue_unassigned_reads_via_transcriptome_alignment",
    )
    parser.add_argument(
        # Declared here so it exists on BOTH routes into the pipeline, which is the
        # whole point of this parser holding every default. Without it, turning
        # rescue off did not disable rescue -- it killed the run: the master flag
        # reached the shard and stopped there, each worker re-derived rescue as True
        # from its own default, and LRAA's guard then refused the resulting
        # rescue-on + --stream_reads pair, naming the very flag the worker was never
        # handed. The config-override route cannot substitute: that guard reads
        # `args` and fires ~140 lines before --config_update is applied.
        "--no_rescue_unassigned_reads_via_transcriptome_alignment",
        dest="rescue_unassigned_reads_via_transcriptome_alignment",
        action="store_false",
        help="pass --no_rescue_unassigned_reads_via_transcriptome_alignment to "
        "every chunk worker instead",
    )
    parser.add_argument(
        "--stream_reads_rescue_unassigned",
        dest="stream_reads_rescue_unassigned",
        action="store_true",
        default=None,
        help="unset (the default) tracks --stream_reads AND transcriptome rescue: "
        "on only when both are on, which is how LRAA itself resolves it -- see "
        "--no_stream_reads_rescue_unassigned",
    )
    parser.add_argument(
        "--no_stream_reads_rescue_unassigned",
        dest="stream_reads_rescue_unassigned",
        action="store_false",
        help="pass --no_stream_reads_rescue_unassigned to every chunk worker instead",
    )
    parser.add_argument(
        "--stream_reads_rescue_unassigned_to_targets",
        action="store_true",
        default=LRAA_Globals.config["stream_reads_rescue_unassigned_to_targets"],
    )
    parser.add_argument(
        "--min_mapping_quality",
        type=int,
        default=LRAA_Globals.config["min_mapping_quality"],
        help="forwarded to every chunk worker AND to stage 2 cut selection and "
        "stage 4 normalization, resolved the same way LRAA itself resolves it "
        "(see resolve_min_mapping_quality)",
    )
    parser.add_argument(
        "--min_mapping_quality_for_final_quant",
        type=int,
        default=LRAA_Globals.config["min_mapping_quality_for_final_quant"],
    )
    parser.add_argument(
        "--rss_sample_interval",
        type=float,
        default=0.5,
        help="seconds between /proc RSS samples; a spike shorter than this can "
        "be missed",
    )
    parser.add_argument(
        "--worker_config_json",
        default=None,
        help="path to a JSON file of LRAA_Globals.config overrides to forward to "
        "EVERY chunk worker as --config_update. This is how a threshold the user "
        "set on `LRAA --chunk` reaches the workers that actually apply it: the "
        "outer driver exits before its own config resolution runs, and lraa_cmd "
        "forwards an explicit allowlist that no threshold flag was ever on, so "
        "without this a chunked run silently used library defaults",
    )

    return parser


def parse_args(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    _resolve_stream_reads_rescue_unassigned(args)

    required = ["output_dir"] if args.only_chunk else ["bam", "genome_fa", "output_dir"]
    if not args.discovery and not args.only_chunk:
        # Quant-only has nothing to quantify without it. Discovery does: no gtf
        # is de novo, which is the mode this option exists to reach. An
        # --only_chunk leaf needs neither: its inputs are the chunk's own mini
        # files and its mode comes from chunk_plan.json.
        required.append("gtf")
    for name in required:
        if not getattr(args, name, None):
            parser.error("--{} is required".format(name))
    if args.cpu_budget is not None and args.cpu_budget < 1:
        parser.error("--cpu_budget must be >= 1")
    # A chunk worker's cwd is its own chunk directory (run_step passes cwd=cdir),
    # not the caller's -- resolved here, once, rather than trusting every caller
    # of this parser to remember to.
    if args.cell_list:
        args.cell_list = os.path.abspath(args.cell_list)

    return args


def _resolve_stream_reads_rescue_unassigned(args):
    """Resolve the ``None`` sentinel left by an unset ``--stream_reads_rescue_unassigned``.

    Byte-for-byte the rule at ``LRAA:1855`` -- ``stream_reads and rescue`` -- now
    that this parser carries the rescue toggle too. It previously tracked
    ``--stream_reads`` alone, because rescue was never forwarded to a chunk worker
    at all; a run with rescue off therefore resolved streaming rescue back ON here,
    the opposite of what the caller asked for. The two resolutions must agree or a
    chunked run and an unchunked one rescue differently from the same command line.
    """

    if args.stream_reads_rescue_unassigned is None:
        args.stream_reads_rescue_unassigned = bool(
            args.stream_reads
            and args.rescue_unassigned_reads_via_transcriptome_alignment
        )


def default_args(**overrides):
    """A fully defaulted namespace, with ``overrides`` applied.

    This is how ``LRAA --chunk`` reaches the pipeline: every default lives in
    one parser, so a flag added here is present on both routes rather than only
    on the one whose argparse got edited.
    """

    parser = build_parser()
    args = parser.parse_args([])
    # Checked against the parser's declared destinations, not `hasattr(args,
    # key)` on the freshly parsed namespace: a SUPPRESS-defaulted flag (the
    # single-cell tags) has no attribute at all until something sets one, so a
    # presence check here would reject the exact override that flag exists for.
    valid_dests = {action.dest for action in parser._actions}
    for key, value in overrides.items():
        if key not in valid_dests:
            raise PipelineError(
                "no such chunked-pipeline setting: {}".format(key)
            )
        setattr(args, key, value)
    # AFTER overrides, not before: an override of `stream_reads` with
    # `stream_reads_rescue_unassigned` left unstated must resolve against the
    # OVERRIDDEN stream_reads, not the pre-override default. An override that
    # named `stream_reads_rescue_unassigned` explicitly (True or False) already
    # left it non-None, so this is a no-op for that caller.
    _resolve_stream_reads_rescue_unassigned(args)
    if getattr(args, "cell_list", None):
        args.cell_list = os.path.abspath(args.cell_list)
    return args


def main(argv=None):

    args = parse_args(argv)
    # This is the ONE path the driver never takes -- LRAA --chunk calls run()
    # with a config it already resolved -- so it is where a directly-invoked
    # run's own --HiFi has to reach the config.
    apply_standalone_hifi_preset(args)
    return run(args)


def claim_pipeline_mode(ckpt, mode):
    """One output directory serves one chunking mode, and records which one.

    ``merged/``, ``timing.json`` and ``outputs.json`` sit at fixed paths and
    stage 6 is not checkpointed, so a second run in the other mode would
    overwrite the first's results rather than add to them. Every per-stage
    sentinel already differs between the modes -- the two chains have different
    roots -- and that is exactly what would make the collision quiet: each stage
    would find itself not done, rerun, and replace a merged table that looked
    finished. So the directory carries an explicit claim, and the second mode is
    refused instead of served.
    """

    for other in PIPELINE_MODES:
        if other != mode and ckpt.done("mode_" + other):
            raise PipelineError(
                "this output directory has already served the {} pipeline "
                "(sentinel {}). The two modes share merged/, timing.json and "
                "outputs.json, so running {} here would overwrite those results "
                "rather than add to them, while every per-stage sentinel "
                "reported the work as new. Use a different "
                "--output_dir.".format(other, ckpt.path("mode_" + other), mode)
            )
    ckpt.mark("mode_" + mode, mode)


def cut_placement_report(selections, discovery):
    """What the partition was asked for, what it got, and what it cost.

    Two things a run must not do quietly. Produce fewer chunks than its geometry
    implies -- six months later that is a performance regression with no visible
    cause. And sever reads -- which is now the EXPECTED outcome rather than a
    refused one, so this report is the only place anyone sees which reads a cut
    cost, and whether they carried junctions.

    Hence per contig-strand: requested, placed, declined for the ANNOTATION (the
    only remaining decline reason -- a locus straddling a boundary is emitted by
    neither chunk, and unlike a read it cannot be dropped and counted), otherwise
    unplaced, tail-merged, chunks. Then per placed cut, the alignments it severs
    split monoexonic against multi-exon, because a cut that dropped three reads of
    depth and one that dropped three junctions are not the same event.

    Printed in BOTH modes, and identical in both: severing is priced the same way
    for quantification and for discovery.

    Returns ``(text, summary)``: the text is printed, the summary is stored in
    timing.json and outputs.json.
    """

    totals = collections.Counter()
    per_selection = []
    lines = [
        "",
        "CUT PLACEMENT ({} mode): severing is a COST, never a veto -- a severed "
        "read is dropped, counted and named, and a multi-exon one costs more than "
        "a monoexonic one. Only the annotation can decline a target.".format(
            "discovery" if discovery else "quant-only"
        ),
        "",
    ]
    for key in sorted(selections):
        for selection in selections[key]:
            counts = selection["counts"]
            label = "{}{}".format(selection["chrom"], key or "")
            unplaced = selection["unplaced_targets"]
            declined = [i for i in unplaced if i.get("declined_annotation")]
            other = [i for i in unplaced if not i.get("declined_annotation")]
            severed = counts.get("alignments_dropped_at_cuts") or 0
            monoexonic = counts.get("alignments_dropped_monoexonic") or 0
            multiexon = counts.get("alignments_dropped_multiexon") or 0
            totals["targets"] += counts["targets"]
            totals["placed"] += counts["cuts_placed"]
            totals["declined"] += len(declined)
            totals["unplaced_other"] += len(other)
            totals["tail_merged"] += counts["targets_tail_merged"]
            totals["chunks"] += counts["segments"]
            totals["severed"] += severed
            totals["severed_monoexonic"] += monoexonic
            totals["severed_multiexon"] += multiexon
            lines.append(
                "  {:<14} {} requested, {} placed, {} declined for annotation, {} "
                "otherwise unplaced, {} tail-merged -> {} chunk(s); {} alignment(s) "
                "severed ({} monoexonic, {} multi-exon)".format(
                    label,
                    counts["targets"],
                    counts["cuts_placed"],
                    len(declined),
                    len(other),
                    counts["targets_tail_merged"],
                    counts["segments"],
                    severed,
                    monoexonic,
                    multiexon,
                )
            )
            for cut in selection["cuts"]:
                if not cut.get("spanning_alignments_dropped"):
                    continue
                lines.append(
                    "      cut at {} (target {}, {:+d}) severs {} alignment(s): {} "
                    "monoexonic, {} multi-exon; searched {} bp".format(
                        cut["position"],
                        cut["target"],
                        cut["offset_from_target"],
                        cut["spanning_alignments_dropped"],
                        cut.get("severed_monoexonic", 0),
                        cut.get("severed_multiexon", 0),
                        cut.get("search_radius"),
                    )
                )
            for item in declined:
                lines.append(
                    "      DECLINED target {}: {}".format(
                        item["target"], item["reason"]
                    )
                )
            for item in other:
                lines.append(
                    "      UNPLACED target {}: {}".format(
                        item["target"], item["reason"]
                    )
                )
            per_selection.append(
                {
                    "label": label,
                    "chrom": selection["chrom"],
                    "strand": selection["strand"],
                    "targets": counts["targets"],
                    "cuts_placed": counts["cuts_placed"],
                    "cuts_declined_annotation": len(declined),
                    "targets_unplaced_other": len(other),
                    "targets_tail_merged": counts["targets_tail_merged"],
                    "chunks": counts["segments"],
                    "alignments_severed": severed,
                    "alignments_severed_monoexonic": monoexonic,
                    "alignments_severed_multiexon": multiexon,
                    "cuts": [
                        {
                            "target": cut["target"],
                            "position": cut["position"],
                            "severed": cut["spanning_alignments_dropped"],
                            "severed_monoexonic": cut.get("severed_monoexonic"),
                            "severed_multiexon": cut.get("severed_multiexon"),
                            "search_radius": cut.get("search_radius"),
                        }
                        for cut in selection["cuts"]
                    ],
                    "declined": [
                        {
                            "target": item["target"],
                            "best_spanning_in_window": item.get(
                                "best_spanning_in_window"
                            ),
                            "reason": item["reason"],
                        }
                        for item in declined
                    ],
                }
            )
    lines.append("")
    lines.append(
        "  TOTAL {} cut(s) requested, {} placed, {} declined for annotation, {} "
        "otherwise unplaced, {} tail-merged -> {} chunk(s); {} alignment(s) severed "
        "({} monoexonic, {} multi-exon)".format(
            totals["targets"],
            totals["placed"],
            totals["declined"],
            totals["unplaced_other"],
            totals["tail_merged"],
            totals["chunks"],
            totals["severed"],
            totals["severed_monoexonic"],
            totals["severed_multiexon"],
        )
    )
    if totals["declined"]:
        lines.append(
            "  {} cut(s) were DECLINED because the annotation left their windows "
            "with no admissible position. The chunks they would have separated stay "
            "joined, so this run is slower than its geometry suggests: a read can "
            "be dropped and counted, a locus cannot.".format(totals["declined"])
        )
    if totals["severed"]:
        lines.append(
            "  {} alignment(s) were severed and are dropped, counted and named "
            "(<cuts>.dropped_reads.tsv, <cuts>.severed_reads.bam). {} of them "
            "carried junctions, which is the part that can cost a model.".format(
                totals["severed"], totals["severed_multiexon"]
            )
        )
    lines.append("")

    summary = {
        "mode": "discovery" if discovery else "quant_only",
        "targets": totals["targets"],
        "cuts_placed": totals["placed"],
        "cuts_declined_annotation": totals["declined"],
        "targets_unplaced_other": totals["unplaced_other"],
        "targets_tail_merged": totals["tail_merged"],
        "chunks": totals["chunks"],
        "alignments_severed": totals["severed"],
        "alignments_severed_monoexonic": totals["severed_monoexonic"],
        "alignments_severed_multiexon": totals["severed_multiexon"],
        "per_selection": per_selection,
    }
    return "\n".join(lines), summary


def format_plan(args, mode, sources, selections):
    """The ``--dry_run`` report: what would run, how much of it, and over what.

    Built from ``planned_chunks``, the same generator stage 3 iterates, so the
    plan cannot describe a partition the real run would not build.
    """

    planned = list(planned_chunks(sources, selections))
    strandless = args.strandless_chunks
    per_chunk = 2 if strandless else 1
    quant_units = len(planned) * per_chunk
    noun = "interval" if strandless else "contig-strand chunk"

    def row(stage, what, detail):
        return "  {:<9} {:<24} {}".format(stage, what, detail)

    lines = [
        "",
        "PLAN (--dry_run): {} chunking, {}, arm={}".format(
            mode, "discovery" if args.discovery else "quant-only", args.arm
        ),
        "",
    ]
    if strandless:
        lines.append(
            row(
                "stage 1",
                "whole-bam strand split",
                "SKIPPED for the chunked arm"
                + (
                    "; run for the control, which IS the strand-split bam"
                    if args.arm in ("baseline", "both")
                    else ""
                ),
            )
        )
    else:
        lines.append(
            row("stage 1", "whole-bam strand split", "1 run over the whole bam")
        )
    lines += [
        row(
            "stage 2",
            "cut selection",
            "{} run(s) over {}".format(
                len(sources), "the RAW bam" if strandless else "each orientation bam"
            ),
        ),
        row(
            "stage 3",
            "chunk extraction",
            "{} extraction(s), one per {}".format(len(planned), noun),
        ),
        row(
            "stage 3b",
            "strand split IN CHUNK",
            "{} run(s), after extraction, inside the parallel phase".format(
                len(planned)
            )
            if strandless
            else "not run -- the whole bam was split at stage 1",
        ),
        row(
            "stage 4",
            "normalize",
            "{}{}".format(
                quant_units,
                " = {} x 2 orientations".format(len(planned)) if strandless else "",
            ),
        ),
        row(
            "stage 5",
            "discovery" if args.discovery else "quant",
            str(quant_units),
        ),
        row(
            "stage 6",
            "merge",
            "{} unit(s){}".format(
                quant_units, ", incl. model gtfs" if args.discovery else ""
            ),
        ),
        "",
        "  {:<20} {:<28} {:>9}  {}".format(
            "chunk", "region", "span Mb", "quant units"
        ),
    ]
    for chunk in planned:
        span = (chunk["rend"] - chunk["lend"] + 1) / 1e6
        if strandless:
            units = "{0}_plus, {0}_minus".format(chunk["chunk_id"])
        else:
            units = chunk["chunk_id"]
        lines.append(
            "  {:<20} {:<28} {:>9.2f}  {}".format(
                chunk["chunk_id"], chunk["region"], span, units
            )
        )
    lines += [
        "",
        "  {} {}s, {} extraction(s), {} quant unit(s)".format(
            len(planned), noun, len(planned), quant_units
        ),
    ]
    if strandless:
        lines.append(
            "  strand-first over the same substrate extracts once per quant "
            "unit, so {} extractions rather than {}".format(
                quant_units, len(planned)
            )
        )
    lines.append("")
    return "\n".join(lines)


def run(args):
    """Execute the pipeline described by ``args``. Returns the outputs mapping.

    Every subprocess below inherits ``WORKER_ENV``, which is how a per-chunk
    LRAA refuses ``--chunk``: the chunk commands this module builds never pass
    it, but a recursive pipeline would fork until the box died, so the guard is
    a marker in the environment rather than a promise about a command line.

    The marker is RESTORED on the way out. In a real run this makes no difference
    -- the process spawns its children and exits -- but this function is also
    called in-process, and leaving the variable set marked the whole caller as a
    chunk worker for everything that followed. That turned test ORDER into a
    correctness input: a later in-process run, or any subprocess inheriting the
    caller's environment, would be refused by the recursion guard with nothing to
    say why. Set-and-restore keeps the guard's meaning ("my parent is a chunked
    pipeline") true for the children that need it and false for everyone else.
    """

    had_worker_env = WORKER_ENV in os.environ
    previous_worker_env = os.environ.get(WORKER_ENV)
    os.environ[WORKER_ENV] = "1"
    try:
        return _run_inner(args)
    finally:
        if had_worker_env:
            os.environ[WORKER_ENV] = previous_worker_env
        else:
            os.environ.pop(WORKER_ENV, None)


def _run_inner(args):
    """The body of :func:`run`, with the worker marker already in place."""

    if args.cpu_budget is None:
        args.cpu_budget = CpuBudget.default_budget()
    elif args.cpu_budget < 1:
        raise PipelineError("cpu_budget must be >= 1")

    # No wiggle resolution here: the MAXIMUM search window is an ABSOLUTE default,
    # deliberately not derived from the spacing. See the config comment on
    # approx_MB_per_cut_wiggle_window for the measurement behind that.
    if args.severed_multiexon_weight < 1:
        raise PipelineError(
            "severed_multiexon_weight must be >= 1: 0 would make severing a "
            "spliced alignment free"
        )
    if args.only_chunk and args.stop_after_make_chunks:
        raise PipelineError(
            "--only_chunk processes an already-extracted chunk and "
            "--stop_after_make_chunks produces them; pick one"
        )
    sg_bam = getattr(args, "bam_for_sg", None)
    priors_bam = getattr(args, "bam_for_priors", None)
    cut_bam = getattr(args, "bam_for_cut_selection", None)
    emit_plan = getattr(args, "emit_cut_plan", None)
    chunk_plan_path = getattr(args, "chunk_plan", None)
    if emit_plan and chunk_plan_path:
        raise PipelineError(
            "--emit_cut_plan produces a partition and --chunk_plan applies one; "
            "pick one. Re-emitting a plan from a plan would select nothing and "
            "write a copy under a name that claims this run chose it."
        )
    if emit_plan and args.stop_after_make_chunks:
        raise PipelineError(
            "--emit_cut_plan stops BEFORE extraction and --stop_after_make_chunks "
            "stops after it; pick one. Extracting here is the cost the split "
            "exists to avoid: on a whole-genome bam that is ~300 extractions of "
            "the pre-partition library, and the callers applying the plan each "
            "extract their own reads anyway."
        )
    for name, value in (
        ("--emit_cut_plan", emit_plan),
        ("--chunk_plan", chunk_plan_path),
    ):
        if not value:
            continue
        if args.only_chunk:
            raise PipelineError(
                "{} decides the PARTITION and --only_chunk processes one chunk of "
                "an existing one; a leaf reads its geometry from the chunk "
                "directory it was handed.".format(name)
            )
        if not args.strandless_chunks:
            raise PipelineError(
                "{} is strandless-chunking only. Strand-first selects cuts over "
                "the stage-1 orientation split, which is a per-run file, so its "
                "geometry is not shareable.".format(name)
            )
    if chunk_plan_path and cut_bam:
        raise PipelineError(
            "--chunk_plan supplies the cut positions and --bam_for_cut_selection "
            "says which bam to choose them from; with a plan there is no "
            "selection to run, so passing both asks for one partition and names "
            "the source of another."
        )
    if emit_plan and args.arm != "chunked":
        raise PipelineError(
            "--emit_cut_plan has no baseline arm: it stops at cut selection, so "
            "the whole-contig control would be skipped silently. Run "
            "--arm chunked."
        )
    if sg_bam and not cut_bam and not chunk_plan_path:
        # ACCEPTED with --chunk_plan, and that is not a loosening: a shared plan
        # IS the shared geometry this refusal was demanding, and a stronger form
        # of it -- the positions themselves rather than a bam to re-derive them
        # from, checked against every parameter that could have moved them.
        #
        # REFUSED rather than resolved, because both available answers are wrong
        # and neither is detectable downstream.
        #
        # Selecting on --bam gives every caller its OWN cut coordinates: cut
        # placement reads the bam (blocked positions, spanning-read cost), so the
        # ~29 per-cluster jobs sharing one sg bam would each slice that shared
        # evidence at different boundaries and each build a different splice graph
        # near every cut. That is the cross-caller incomparability --bam_for_sg
        # exists to remove.
        #
        # Selecting on the sg bam is worse. Normalization REMOVES reads, so a
        # position spanned by nothing in the sg bam can still be spanned by a raw
        # read in a caller's bam; extraction drops that read, no selector named
        # it, and the chunked arm loses a read the unchunked run counts. Silent,
        # and against the parity the whole chunked design rests on.
        raise PipelineError(
            "--bam_for_sg needs --bam_for_cut_selection naming an unthinned "
            "SUPERSET of every --bam this plan is applied to (in practice the "
            "pre-partition input the per-caller bams were split out of). "
            "Selecting cuts on --bam instead would give each caller its own cut "
            "geometry, so the shared splice-graph evidence would be sliced at "
            "different boundaries per caller -- the guarantee --bam_for_sg "
            "exists to restore. Selecting them on the sg bam is not the "
            "alternative: it is coverage-normalized, so a cut can look "
            "unspanned there while a raw read in --bam spans it, and that read "
            "is then dropped by extraction and named by nobody. "
            "--chunk_plan PATH satisfies this too, and is the cheaper shape: one "
            "--emit_cut_plan run over that same superset, then every caller "
            "applies the plan instead of re-selecting on it."
        )
    if (sg_bam or priors_bam) and args.arm in ("baseline", "both"):
        raise PipelineError(
            "--bam_for_sg and --bam_for_priors have no baseline arm: the control "
            "normalizes the whole contig for itself, so it would build its splice "
            "graph from evidence the chunked arm never saw and estimate theta "
            "from reads the chunked arm never used, and the two arms would not be "
            "comparable. Run --arm chunked."
        )
    for name, value in (
        ("--bam_for_sg", sg_bam),
        ("--bam_for_priors", priors_bam),
        ("--bam_for_cut_selection", cut_bam),
    ):
        if not value:
            continue
        if not args.strandless_chunks:
            raise PipelineError(
                "{} is strandless-chunking only. Strand-first splits the whole "
                "bam by orientation before cutting, so there is no such split of "
                "a caller-supplied bam: the shared evidence has none to cut "
                "against, and the priors slice has no orientation-pure form for a "
                "unit to read. Stage 3b, which produces both, runs only for a "
                "strandless chunk.".format(name)
            )
        if not os.path.exists(value):
            raise PipelineError("{} {} does not exist".format(name, value))
        if args.bam and os.path.realpath(value) == os.path.realpath(args.bam):
            if name == "--bam_for_sg":
                raise PipelineError(
                    "--bam_for_sg and --bam resolve to the same file, so the "
                    "splice graph would be built from the reads being counted "
                    "with no normalization anywhere. LRAA refuses the same "
                    "composition (LRAA:1945-1957). Drop the flag and let each "
                    "chunk normalize its own reads."
                )
            if name == "--bam_for_priors":
                raise PipelineError(
                    "--bam_for_priors and --bam resolve to the same file, so pass "
                    "1 would estimate theta over the very population pass 2 "
                    "assigns -- the no-priors configuration, spelled with an "
                    "extra slice extracted per chunk. LRAA refuses it downstream "
                    "too: --bam must carry no XW weight and --bam_for_priors "
                    "must. Drop the flag, or name this caller's own NORMALIZED "
                    "reads."
                )
    if sg_bam and cut_bam and os.path.realpath(sg_bam) == os.path.realpath(cut_bam):
        raise PipelineError(
            "--bam_for_cut_selection and --bam_for_sg resolve to the same file. "
            "The sg bam is coverage-normalized, so it is a SUBSET of the reads "
            "being extracted: a cut it shows as unspanned can still be spanned "
            "by a raw read, which extraction then drops with no selector having "
            "named it -- read loss the chunked arm reports as success. The "
            "selection source must be an unthinned superset."
        )
    if sg_bam and priors_bam and Util_funcs.paths_name_one_file(sg_bam, priors_bam):
        raise PipelineError(
            "--bam_for_sg and --bam_for_priors are ONE file. That is the POOLED "
            "configuration --bam_for_priors exists to discourage: the sg bam is ONE "
            "splice graph shared by every cell cluster, so a theta estimated over it "
            "apportions this caller's ambiguous reads by every other caller's "
            "expression -- and it looks like it worked, because each caller still "
            "reports its own read totals. Compared by device and inode, so a second "
            "path, a symlink or a hard link is caught; a byte-identical COPY is not, "
            "and cannot be without hashing the bams, so passing this does not rule "
            "the pooled shape out. Name this caller's OWN normalized reads, or pass "
            "no --bam_for_priors and get today's behaviour."
        )

    required = ["output_dir"] if args.only_chunk else ["bam", "genome_fa", "output_dir"]
    if not args.discovery and not args.only_chunk:
        required.append("gtf")
    for name in required:
        if not getattr(args, name, None):
            raise PipelineError(
                "the chunked pipeline needs {}".format(name)
            )
    if args.discovery and args.arm in ("baseline", "both"):
        # The control arm is a whole-contig QUANT run against a supplied
        # annotation; it has no discovery counterpart here. Refused rather than
        # served, because silently running a quant baseline beside a discovery
        # chunked arm would produce two artifacts that look comparable and are
        # not.
        raise PipelineError(
            "--discovery has no baseline arm: the control is a whole-contig "
            "quantification against a supplied annotation, which is not the "
            "comparison a discovery run wants. Run --arm chunked, and produce "
            "the unchunked discovery control with a plain LRAA invocation."
        )

    mode = STRANDLESS_MODE if args.strandless_chunks else STRAND_FIRST_MODE
    outdir = os.path.abspath(args.output_dir)
    os.makedirs(os.path.join(outdir, "logs"), exist_ok=True)
    ckpt = Checkpoints(os.path.join(outdir, "__ckpt"))
    claim_pipeline_mode(ckpt, mode)

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
    timing["chunk_order"] = mode
    timing["dry_run"] = bool(args.dry_run)
    timing["discovery"] = bool(args.discovery)
    timing["severed_multiexon_weight"] = args.severed_multiexon_weight
    timing["emit_cut_plan"] = emit_plan
    timing["chunk_plan"] = chunk_plan_path

    shared_plan = None
    if chunk_plan_path:
        if os.path.realpath(chunk_plan_path) == os.path.realpath(
            os.path.join(outdir, "chunk_plan.json")
        ):
            # --stop_after_make_chunks writes the leaf plan to exactly that path,
            # so this composition would overwrite the shared plan every sibling
            # run is reading -- with a file describing one caller's extracted
            # chunks. Refused rather than renamed, because the two are different
            # artifacts and the caller has to say which it means.
            raise PipelineError(
                "--chunk_plan {} is the path --stop_after_make_chunks writes its "
                "own leaf plan to. A shared plan is read by every caller, so a "
                "run that rewrote it in place would change its siblings' "
                "geometry mid-flight. Localize it under another name.".format(
                    chunk_plan_path
                )
            )
        shared_plan = load_chunk_plan(chunk_plan_path)
        # Checked HERE, before the bam is indexed and before a single region is
        # read: every refusal in it is about the arguments, so it costs nothing
        # and a mismatched plan should not fail minutes into a run.
        validate_cut_plan_geometry(shared_plan, chunk_plan_path, args)

    for tool in (SEPARATE_BAM, SELECT_CUTS, EXTRACT_CHUNK, NORMALIZE_BAM, LRAA):
        if not os.path.exists(tool):
            raise PipelineError("missing stage tool: {}".format(tool))
    if shutil.which("samtools") is None:
        raise PipelineError("samtools is not on PATH")

    rss = args.rss_sample_interval

    if args.only_chunk:
        plan = load_chunk_plan(os.path.join(outdir, "chunk_plan.json"))
        args.discovery = bool(plan["discovery"])
        chunk = rebuild_chunk_record(plan, args.only_chunk, outdir)
        record = chunk_worker(
            args,
            ckpt,
            outdir,
            chunk,
            plan["num_total_reads"],
            rss,
            cpu_budget=args.cpu_budget,
        )
        units_path = os.path.join(chunk["dir"], "units.json")
        with open(units_path, "wt") as fh:
            json.dump(
                {
                    "chunk_id": chunk["chunk_id"],
                    "units": [
                        {
                            "unit_id": u["unit_id"],
                            "strand": u["strand"],
                            "offset": u["offset"],
                            "order": u["order"],
                            "quant_prefix": u["quant_prefix"],
                        }
                        for u in chunk["units"]
                    ],
                },
                fh,
                indent=2,
                sort_keys=True,
            )
        print("chunk {} complete: {}".format(chunk["chunk_id"], units_path), flush=True)
        return {"only_chunk": chunk["chunk_id"], "units": units_path, "chunk": record}

    # The chain root. Every artifact below derives from these three files and the
    # contig restriction, and none of the per-stage sentinels named any of them:
    # pointing the same outdir at a different bam, genome or gtf reused the whole
    # pipeline. Identity is the resolved path with size and mtime, not a content
    # hash, because these run to gigabytes.
    inputs_token = chain_token(
        "inputs",
        None,
        Util_funcs.file_identity_token(args.bam),
        Util_funcs.file_identity_token(args.genome_fa),
        # De novo discovery has no annotation, and a null identity is not the
        # identity of some other file: it is the absence of an input.
        Util_funcs.file_identity_token(args.gtf) if args.gtf else "no_gtf",
        args.contig or "",
    )

    # Skipping stage 1 IS the strandless mode; it is not made cheaper, it is not
    # run. The CONTROL still needs it, because the whole-contig baseline is
    # precisely the strand-split bam glued back together, so a parity run pays
    # for it once and a pure chunked run does not pay at all.
    needs_split = (not args.strandless_chunks) or args.arm in ("baseline", "both")
    strand_bams = split_token = None
    if needs_split:
        strand_bams, split_token = stage_strand_split(
            args, ckpt, outdir, timing, rss, inputs_token
        )
        retained = {s: count_records(b) for s, b in strand_bams.items()}
        retained_total = sum(retained.values())
        timing["stage1_retained_records"] = {
            "plus": retained["+"],
            "minus": retained["-"],
            "total": retained_total,
        }
        if args.num_total_reads is None:
            num_total_reads = retained_total
            timing["num_total_reads_source"] = "stage 1 retained record count"
        else:
            # AUTHORITATIVE, and only warned about when it differs. It used to be a
            # hard refusal on the grounds that "the denominator has to be the record
            # set the arms actually see, or TPM is not comparable between them" --
            # but that argument is about not inventing a DIFFERENT default per arm.
            # One stated value goes to both arms, so they stay comparable with each
            # other, and it is the only way to stay comparable with an UNCHUNKED run:
            # unchunked counts -F 0x904, which retains MORE than stage 1 does (the
            # splitter also drops duplicate, qcfail and long-intron records), so the
            # two can legitimately differ on a perfectly good bam. Refusing made the
            # whole-library denominator unreachable on this arm -- MEASURED on a SIRV
            # bam with 2,095 of 104,766 primary records flagged duplicate, `-N 104766`
            # exited 1.
            #
            # Safe to trust because nothing DERIVES this value any more: LRAA leaves
            # it unset for this arm precisely so stage 1 owns the default, and the
            # only routes that set it are a caller's own -N.
            num_total_reads = args.num_total_reads
            timing["num_total_reads_source"] = "supplied via -N"
            if num_total_reads != retained_total:
                # Both channels, because this module has no logger: stderr so the
                # caller sees it, and timing.json so it survives into the run record
                # a benchmark reads back.
                warning = (
                    "WARNING: -N {} differs from the stage-1 retained record count "
                    "{} ({} on + plus {} on -); USING THE SUPPLIED VALUE. Unset -N "
                    "to use the retained count instead. The two differ when the bam "
                    "holds records the splitter drops but -F 0x904 keeps: "
                    "duplicate, qcfail, or long-intron.".format(
                        num_total_reads, retained_total, retained["+"], retained["-"]
                    )
                )
                print(warning, file=sys.stderr)
                timing["num_total_reads_warning"] = warning
    else:
        timing.setdefault("stages", {})["strand_split"] = {
            "skipped": "--strandless_chunks with --arm chunked: the orientation "
            "split runs per chunk, at stage 3b"
        }
        timing["stage1_retained_records"] = None
        # There is no stage-1 count to default to, and inventing one would
        # silently move the RPM_total_reads column relative to a strand-first
        # run of the same substrate -- the one column stage 6 does not rebase.
        # A dry run never reaches quantification, so it does not need the number,
        # and neither does --emit_cut_plan: cut selection has no use for a TPM
        # denominator, the plan records none as authoritative, and requiring one
        # would put a whole-bam counting pass in front of the one phase that
        # exists to be cheap and shared.
        if args.num_total_reads is None:
            if not args.dry_run and not emit_plan:
                raise PipelineError(
                    "-N is required with --strandless_chunks --arm chunked. That "
                    "arm does not run stage 1, so there is no retained record "
                    "count to default to, and RPM_total_reads is computed "
                    "against whatever is passed -- a guess here would look like "
                    "a quantification difference. Use "
                    "stage1_retained_records.total from a strand-first run over "
                    "the same bam, or run --arm both once, which does run "
                    "stage 1."
                )
            num_total_reads = None
        else:
            num_total_reads = args.num_total_reads
        timing["num_total_reads_source"] = (
            "supplied via -N (stage 1 skipped, nothing to cross-check against)"
        )
    timing["num_total_reads"] = num_total_reads

    outputs = {"num_total_reads": num_total_reads, "chunk_order": mode}

    def flush():
        with open(timing_path, "wt") as fh:
            json.dump(timing, fh, indent=2, sort_keys=True)
            print("", file=fh)

    flush()

    if args.arm in ("chunked", "both"):
        # Only the strandless arm reads the input by coordinate, so only it pays for
        # an index build -- imposing one on strand-first would cost a whole-bam pass
        # for an index that arm never opens. After argument validation, not before:
        # a missing -N is a refusal that should arrive immediately rather than
        # behind minutes of indexing.
        if args.strandless_chunks:
            ensure_bam_index(args.bam)
            # Both are fetched by coordinate as well: the selection source by
            # cut_selection_plan and by derived_severed_names, the sg bam by the
            # per-chunk extractor. Indexed here for the same reason --bam is,
            # rather than discovered missing inside a pool worker.
            for extra in (
                getattr(args, "bam_for_cut_selection", None),
                getattr(args, "bam_for_sg", None),
                getattr(args, "bam_for_priors", None),
            ):
                if extra:
                    ensure_bam_index(extra)
        sources = cut_sources(args, strand_bams, inputs_token, split_token)
        # Stages 2 and 3 are ONE pool now rather than two serial loops, so they
        # are entered together. A dry run selects and stops: the selection is real
        # work whose result the printed plan describes, and extraction is the part
        # a dry run must not do. --emit_cut_plan stops in the same place for the
        # opposite reason -- the selection is the DELIVERABLE and the extraction
        # belongs to the callers -- and --chunk_plan replaces the selection
        # entirely, seeding those extractions from geometry someone else chose.
        selections, cut_dir, chunks, prep_makespan, prep_allocation = (
            run_prep_concurrently(
                args,
                ckpt,
                outdir,
                timing,
                sources,
                rss,
                select_only=bool(args.dry_run or emit_plan),
                chunk_plan=shared_plan,
            )
        )
        flush()
        placement_text, placement = cut_placement_report(selections, args.discovery)
        print(placement_text, flush=True)
        timing["cut_placement"] = placement
        outputs["cut_placement"] = placement
        flush()

        if args.dry_run:
            plan = format_plan(args, mode, sources, selections)
            print(plan, flush=True)
            timing["plan"] = plan.splitlines()
            flush()
            # outputs.json is NOT written: a dry run must not overwrite the
            # record of a real run that already finished in this directory.
            outputs["dry_run"] = True
            outputs["cut_dir"] = cut_dir
            return outputs

        if emit_plan:
            write_chunk_plan(
                emit_plan,
                args,
                sources,
                selections,
                # Recorded as null. The denominator belongs to whoever quantifies,
                # and this run does not: a plan that carried one caller's count
                # would hand every sibling the wrong TPM scale.
                None,
                chunks_extracted=False,
            )
            print(
                "cut plan written: {} chunk(s) of geometry at {}; nothing "
                "extracted. Apply it with --chunk_plan {}".format(
                    len(list(planned_chunks(sources, selections))),
                    emit_plan,
                    emit_plan,
                ),
                flush=True,
            )
            timing["emitted_cut_plan"] = emit_plan
            flush()
            # outputs.json is NOT written, for the reason a dry run does not write
            # it: the plan is this run's artifact, and the output directory may
            # already hold the record of a real run that finished in it.
            outputs["cut_plan"] = emit_plan
            outputs["cut_dir"] = cut_dir
            return outputs

        # Before the expensive phase, not after it: this is the check that keeps
        # the control's pruned bam and the chunked arm's inputs the same record
        # set, and it needs nothing but the manifests and the cut selection.
        #
        # With a shared selection source the equality is asserted PER INPUT --
        # each bam's own severed set at these same cuts against its own drops --
        # because the selector saw a superset of any one caller's reads. Without
        # one this is exactly the single-source check it has always been.
        #
        # A --chunk_plan run is the same situation and needs the same treatment,
        # which is what keeps its promise honest: the cuts were chosen on a bam
        # this run does not hold, so a name the selection wrote may be absent here
        # (fine), while every read of THIS bam that spans a shared cut must appear
        # in some chunk's drops. That direction is re-derived and stays fatal, so
        # shared geometry cannot lose a read without saying so. The cut directory
        # is empty on this path -- no selection ran in it -- so the drops become
        # EXTRACTION_ONLY_DROPS and are counted and named there.
        severed_inputs = None
        if cut_bam or chunk_plan_path:
            severed_inputs = [
                {
                    "label": "reads",
                    "bam": os.path.abspath(args.bam),
                    "names_key": "dropped_read_names",
                }
            ]
            if getattr(args, "bam_for_sg", None):
                severed_inputs.append(
                    {
                        "label": "splice_graph_evidence",
                        "bam": os.path.abspath(args.bam_for_sg),
                        "names_key": "sg_dropped_read_names",
                    }
                )
            if getattr(args, "bam_for_priors", None):
                severed_inputs.append(
                    {
                        "label": "theta_priors",
                        "bam": os.path.abspath(args.bam_for_priors),
                        "names_key": "priors_dropped_read_names",
                    }
                )
        timing["severed_read_accounting"] = verify_severed_accounting(
            cut_dir,
            chunks,
            inputs=severed_inputs,
            max_intron_length=args.max_intron_length,
        )
        flush()
        print(
            "extracted {} chunk(s): {}".format(
                len(chunks), ", ".join(c["region"] for c in chunks)
            ),
            flush=True,
        )

        if args.stop_after_make_chunks:
            plan_path = os.path.join(outdir, "chunk_plan.json")
            write_chunk_plan(
                plan_path, args, sources, selections, num_total_reads,
                chunks_extracted=True,
            )
            print(
                "make-chunks complete: {} chunk(s), plan at {}".format(
                    len(chunks), plan_path
                ),
                flush=True,
            )
            timing["stopped_after_make_chunks"] = True
            outputs["stopped_after_make_chunks"] = True
            outputs["chunk_plan"] = plan_path
            outputs["cut_dir"] = cut_dir
            flush()
            with open(os.path.join(outdir, "outputs.json"), "wt") as fh:
                json.dump(outputs, fh, indent=2, sort_keys=True)
                print("", file=fh)
            return outputs

        load_before = loadavg()
        # Spans chunk processing AND the stage-6 merge that follows, not just the
        # former. It used to stop the instant run_chunks_concurrently returned, so
        # "observed_peak_concurrent_tree_rss_kb" undercounted a whole-genome run's
        # real peak by up to 2.67x -- the peak lives in the merge (an external
        # merge sort over every unit's tracking rows), which is SERIAL and ran
        # entirely outside the old window. Renamed without "concurrent": the merge
        # is not concurrent with anything, and the field now covers both phases.
        arm_sampler = RssSampler(os.getpid(), rss)
        arm_sampler.start()
        try:
            chunk_records, makespan, chunk_allocation, chunk_mem = (
                run_chunks_concurrently(
                    args, ckpt, outdir, chunks, num_total_reads, rss
                )
            )
            load_after = loadavg()

            units = ordered_units(chunks)
            # Emit the same unit list stage 6 is about to consume, in the same order,
            # as a manifest util/misc/merge_chunk_outputs.py can be handed directly.
            # Written from `units` rather than rebuilt from `chunks` so the standalone
            # merger and this driver cannot describe different work, and written BEFORE
            # the merge so it survives a merge that dies. Order is load-bearing here:
            # see the comment on `ordered_units`.
            merge_manifest = os.path.join(outdir, "merge_manifest.json")
            try:
                # Imported HERE, not at module scope, because merge_chunk_outputs
                # imports this module for `merge_and_translate` -- there is exactly one
                # merge implementation and it lives here. Deferring to call time breaks
                # the cycle and keeps the manifest schema documented and implemented in
                # the one file that also reads it.
                sys.path.insert(
                    0,
                    os.path.sep.join(
                        [os.path.dirname(os.path.realpath(__file__)), "..", "util", "misc"]
                    ),
                )
                import merge_chunk_outputs as _merge_cli

                _merge_cli.write_manifest(merge_manifest, units)
                print(
                    "stage 6 unit manifest written to {}".format(merge_manifest),
                    flush=True,
                )
            except Exception as _e:
                # A manifest is a convenience for scattering the merge elsewhere; it is
                # not an input to the merge that follows, so failing to write it must
                # not lose a run that has already done all of the work.
                print(
                    "NOTE: could not write {}: {}".format(merge_manifest, _e),
                    flush=True,
                )
            merged = merge_and_translate(outdir, units, discovery=args.discovery)
        finally:
            arm_sampler.stop()
        timing.setdefault("arms", {})["chunked"] = {
            "cpu_budget": chunk_allocation.budget,
            "concurrent_chunk_workers": chunk_allocation.unit_workers,
            "cpu_budget_per_chunk": chunk_allocation.tool_threads,
            "unallocated_cores": chunk_allocation.unallocated_cores,
            "makespan_s": round(makespan, 3),
            "summed_wall_s": round(sum(c["wall_s"] for c in chunk_records), 3),
            "chunks_extracted": len(chunks),
            # The make-chunks phase's own numbers, beside the processing phase's,
            # because the two together are the run and the first used to be the
            # serial half of it.
            "make_chunks_makespan_s": round(prep_makespan, 3),
            "make_chunks_unit_workers": prep_allocation.unit_workers,
            "quant_units": len(units),
            "chunks": chunk_records,
            # peak_rss_kb_summed_over_chunk_peaks (removed): summing each chunk's
            # OWN peak treated concurrent, independent peaks as if they stacked.
            # MEASURED 12.0x over the real concurrent peak on HG002, 17.5x on ONT,
            # and its own control showed zero discriminating power -- chunked and
            # unchunked summed peaks agreed within 4.5% while their real peaks
            # differed 5.2x. Per-chunk peaks remain available above, in `chunks`.
            "observed_peak_tree_rss_kb": arm_sampler.peak_kb,
            # Reported whether or not the bound reduced anything: the estimate beside
            # the observed peak is what makes the guard's error measurable on every
            # run instead of only on the ones it changed.
            "memory_cap_chunk_workers": chunk_mem.cap,
            "memory_note": chunk_mem.note,
            "memory_available_mib_at_dispatch": chunk_mem.available_mib,
            "largest_chunk_estimate_mib": chunk_mem.largest_unit_mib,
            "concurrent_chunk_estimate_mib": chunk_mem.charged_mib,
            "loadavg_before": load_before,
            "loadavg_after": load_after,
        }
        if args.strandless_chunks:
            timing["strandless_split_accounting"] = roll_up_split_accounting(chunks)
        timing["chunk_manifests"] = [
            {
                "chunk_id": c["chunk_id"],
                "region": c["region"],
                "offset": c["offset"],
                "window_origin": c["window_origin"],
                "quant_units": [u["unit_id"] for u in c["units"]],
                "alignments_emitted": c["manifest"]["counts"]["alignments_emitted"],
                "alignments_dropped_overhang": c["manifest"]["counts"][
                    "alignments_dropped_overhang"
                ],
                "gtf_transcripts_emitted": c["manifest"]["counts"][
                    "gtf_transcripts_emitted"
                ],
                "split_counts": c.get("split_counts"),
                "log": c["log"],
            }
            for c in chunks
        ]
        outputs["chunked"] = merged
        outputs["cut_dir"] = cut_dir
        flush()
    elif args.dry_run:
        # --arm baseline has no partition to plan: the control is one merge, one
        # normalize and one quant over the whole contig. Said, not run, because
        # a dry run that quantified would be a run.
        print(
            "\nPLAN (--dry_run): {} chunking, arm=baseline\n\n"
            "  baseline   merge + normalize + quant   1 whole-contig run, minus "
            "the reads any cut selection in this directory named as "
            "severed\n".format(mode),
            flush=True,
        )
        return outputs

    if args.arm in ("baseline", "both"):
        # Read from disk rather than carried down from the chunked arm, so a
        # baseline run into an output directory whose cuts already exist is
        # comparable with them even when the two arms ran as separate invocations.
        severed = severed_read_names(os.path.join(outdir, "cuts"))
        timing["baseline_excluded_severed_reads"] = len(severed)
        load_before = loadavg()
        outputs["baseline"] = run_baseline(
            args,
            ckpt,
            outdir,
            timing,
            strand_bams,
            num_total_reads,
            rss,
            split_token,
            severed_names=severed,
        )
        timing["arms"]["baseline"]["loadavg_before"] = load_before
        timing["arms"]["baseline"]["loadavg_after"] = loadavg()
        flush()

    with open(os.path.join(outdir, "outputs.json"), "wt") as fh:
        json.dump(outputs, fh, indent=2, sort_keys=True)
        print("", file=fh)

    print(json.dumps(outputs, indent=2, sort_keys=True))
    return outputs


if __name__ == "__main__":
    try:
        main()
        sys.exit(0)
    except PipelineError as err:
        print("\nPIPELINE FAILED\n{}".format(err), file=sys.stderr)
        sys.exit(1)
