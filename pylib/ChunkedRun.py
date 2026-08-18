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
import glob
import gzip
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

    A changed ``--cpu_budget`` or interpreter path moves the digest too, and that
    is deliberate rather than tolerated: those are a MISS, costing a rerun, and
    the alternative is an exclusion list that has to be right forever.
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

    ``--contig`` reduces the enumeration to a single entry, which leaves today's
    behaviour as a degenerate case of the same pool rather than a second path.
    """

    with pysam.FastaFile(os.path.abspath(args.genome_fa)) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
    if not lengths:
        raise PipelineError(
            "{} names no contigs, so there is nothing to select cuts on. The "
            "make-chunks phase enumerates the genome fasta because that is what "
            "cut selection enumerates.".format(args.genome_fa)
        )
    if args.contig:
        if args.contig not in lengths:
            raise PipelineError(
                "--contig {} is absent from {}. Known contigs begin: {}".format(
                    args.contig, args.genome_fa, ", ".join(sorted(lengths)[:10])
                )
            )
        return [args.contig], lengths
    return sorted(lengths), lengths


def _extractor_module():
    """``extract_contig_region_inputs``, imported into this process.

    For ``ensure_gtf_index`` alone. Everything else this pipeline asks of the
    extractor it asks through a subprocess, so that a per-chunk failure lands in
    that chunk's own log rather than killing the driver.
    """

    path = os.path.join(REPO_ROOT, "util", "misc")
    if path not in sys.path:
        sys.path.insert(0, path)
    import extract_contig_region_inputs

    return extract_contig_region_inputs


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


def chunk_memory_cap(budget, unit_spans_bp, available_mib=None):
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
    """

    budget = max(1, int(budget))
    estimates = sorted((chunk_unit_peak_mib(span) for span in unit_spans_bp), reverse=True)
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
    """

    if args.strandless_chunks:
        return [("", STRANDLESS_TAG, os.path.abspath(args.bam), inputs_token)]
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
    return 97.0 if args.HiFi else LRAA_Globals.config["min_per_id"]


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
    # filters (LRAA:4201-4204).  Reading the discovery key here would be wrong
    # for every run that raises the final-quant threshold; both default to 0,
    # so nothing would have failed until someone set it.
    effective_min_mapq = int(
        LRAA_Globals.config["min_mapping_quality_for_final_quant"]
    )
    # ``tag`` carries the mode into the sentinel filename -- "strandless"
    # rather than "plus"/"minus" -- and the parent token differs as well
    # (the inputs, not a split), so neither mode can read the other's cuts.
    # ``contig`` is in ``local`` and not merely in the prefix, because two
    # contigs sharing a sentinel would let one contig's completion mark the
    # other's work done.
    token = chain_token(
        "stage2_cuts_{}_{}.mb_{}_wig_{}_dw_{}_margin_{}.sev_pid_{}_mq_{}"
        ".mxw_{}_annot_{}".format(
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
        "0",
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
    """

    if strand:
        units = [
            (
                strand,
                chunk_id,
                "chunk.norm.bam",
                "chunk_quant",
                "{}.bam".format(prefix),
                "{}.gtf".format(prefix) if has_gtf else None,
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
            "quant_name": quant_name,
            "quant_prefix": os.path.join(cdir, quant_name + "." + lraa_suffix),
        }
        for s, unit_id, norm_name, quant_name, bam, gtf in units
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
    reuse_source = strandless and planned["spans_whole_contig"]
    local = "stage3_extract_{}{}.margin_{}_maxintron_{}{}".format(
        "strandless_" if strandless else "",
        chunk_id,
        args.margin,
        args.max_intron_length,
        # In the sentinel because it decides what is on disk: a directory whose
        # chunk was extracted must not be reported as reusable by a run that
        # would have reused the source, and the reverse.
        ".srcbam" if reuse_source else "",
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
        planned["bam"],
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
                ),
            }
        )
    timing.setdefault("stages", {})["extract_chunks"] = [
        c["extract"] for c in chunks
    ]
    return chunks


def run_prep_concurrently(
    args, ckpt, outdir, timing, sources, rss_interval, select_only=False
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

    ``select_only`` runs the cut selections and stops, for ``--dry_run``: the
    selection is real work whose result the plan describes, and extraction is
    exactly the part a dry run must not perform. ``chunks`` comes back empty.

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

    selection_plans = []
    selection_units = []
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
    if not selection_units:
        raise PipelineError(
            "no cut selection units to run: {} contig(s) x {} source(s)".format(
                len(contigs), len(sources)
            )
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
            "{} cut selection unit(s) collapsed to {} distinct (contig, "
            "orientation) key(s). Two units sharing a key would share an output "
            "prefix and a sentinel, so one contig's cuts would overwrite the "
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
    print(
        "make_chunks queue: {} cut selection(s) seeded over {} contig(s) x {} "
        "source bam(s), longest-first on mapped reads; each contig's chunk "
        "extractions join the SAME pool when its own selection finishes, so the "
        "queue grows past the seed and there is no stage barrier.".format(
            seeded, len(contigs), len(sources)
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
                submit(selection_task, plan)

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

    missing = [
        (key, contig)
        for key in source_by_key
        for contig in contigs
        if (key, contig) not in selections_by_contig
    ]
    if missing:
        raise PipelineError(
            "{} of {} cut selection(s) never reported: {}. Refusing to assemble a "
            "selection list with holes in it -- the missing contigs would simply "
            "be absent from the run.".format(
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
            "cut_selections": seeded,
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
            seeded,
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


# ---------------------------------------------------------------- stages 4 + 5


def lraa_cmd(
    args, bam_for_quant, bam_for_sg, genome, gtf, out_prefix, num_total_reads, cpu_budget
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
    """
    cmd = [
        sys.executable,
        LRAA,
        "--genome",
        genome,
        "--bam",
        bam_for_quant,
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
    ]
    if args.HiFi:
        cmd.append("--HiFi")
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

    Returns ``(step record, token stages 4-5 chain onto, counts)``.
    """

    assert_extracted_strandlessly(chunk)

    cid = chunk["chunk_id"]
    prefix = chunk["prefix"]
    manifest = chunk["manifest"]
    split_prefix = "{}.strand".format(prefix)
    reused = bool(manifest.get("bam_reused_from_source"))
    source_bam = manifest["files"]["bam"] if reused else "{}.bam".format(prefix)
    token = chain_token(
        "stage3b_split_{}.maxintron_{}{}".format(
            cid,
            args.max_intron_length,
            # In the sentinel because it decides WHICH FILE was read: a resumed
            # run must not report a split of a mini bam as a split of the source.
            ".srcbam_" + chunk["chrom"] if reused else "",
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
    if ckpt.done(token):
        step = {"step": "stage3b_strand_split", "reused": True, "cmd": cmd}
    else:
        step = run_step(
            "stage3b_strand_split_{}".format(cid),
            cmd,
            chunk["log"],
            chunk["dir"],
            rss_interval,
        )
        ckpt.mark(token)

    counts = verify_chunk_split(chunk, split_prefix)
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

    Column 7 and nothing else: every line the extractor emitted for a gene
    carries that gene's orientation, and the coordinates were rebased when it
    was written, so this selects lines and moves nothing. The transcript counts
    are checked against the extractor's own tally, which is the same discipline
    the bam split gets -- a partition that does not add up is a partition that
    lost a model, and a lost model is a row missing from the merged table rather
    than an error anyone would see.
    """

    source = "{}.gtf".format(chunk["prefix"])
    written = {"+": 0, "-": 0}
    handles = {}
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
                if fields[2] == "transcript":
                    written[strand] += 1
    finally:
        for handle in handles.values():
            handle.close()

    emitted = chunk["manifest"]["counts"]["gtf_transcripts_emitted"]
    total = written["+"] + written["-"]
    if total != emitted:
        raise PipelineError(
            "chunk {} annotation split accounting: {} + {} = {} transcript "
            "line(s) across the two orientations, but extraction emitted {}. "
            "Every model the chunk holds has to reach exactly one of the two "
            "quant units.".format(
                chunk["chunk_id"], written["+"], written["-"], total, emitted
            )
        )

    return {
        "gtf_transcripts_plus": written["+"],
        "gtf_transcripts_minus": written["-"],
        "gtf_transcripts_emitted": emitted,
    }


def verify_chunk_split(chunk, split_prefix):
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
    """

    counts = chunk["manifest"]["counts"]
    emitted = counts["alignments_emitted"]
    observed = {}
    for strand in ("+", "-"):
        bam = "{}.{}.bam".format(split_prefix, strand)
        if not os.path.exists(bam):
            raise PipelineError(
                "chunk {} strand split produced no {}".format(chunk["chunk_id"], bam)
            )
        observed[strand] = count_records(bam)

    total = observed["+"] + observed["-"]
    if total != emitted:
        raise PipelineError(
            "chunk {} strand split accounting: {} + {} = {} records across the "
            "two orientations, but extraction emitted {}. The split lost {} "
            "record(s) the chunk bam held. Extraction and the split must apply "
            "the same record filter at the same --max_intron_length; they did "
            "not.".format(
                chunk["chunk_id"],
                observed["+"],
                observed["-"],
                total,
                emitted,
                emitted - total,
            )
        )

    expected = (
        counts.get("alignments_emitted_forward"),
        counts.get("alignments_emitted_reverse"),
    )
    if None not in expected and (observed["+"], observed["-"]) != expected:
        raise PipelineError(
            "chunk {} strand split disagrees with extraction on which records "
            "are which: the split wrote {}/{} (+/-), extraction counted {}/{} by "
            "raw flag. The totals match, so nothing was lost -- records moved "
            "between orientations, which only happens under "
            "--infer_read_orient, and this pipeline does not pass it.".format(
                chunk["chunk_id"],
                observed["+"],
                observed["-"],
                expected[0],
                expected[1],
            )
        )

    return {
        "alignments_emitted": emitted,
        "records_plus": observed["+"],
        "records_minus": observed["-"],
        "records_total": total,
        "records_lost": emitted - total,
    }


def chunk_worker(args, ckpt, outdir, chunk, num_total_reads, rss_interval, cpu_budget):
    """Stages 3b, 4 and 5 for one chunk. Everything goes to the chunk's own log.

    ``cpu_budget`` is this chunk's share of the total, handed down by the scheduler.

    A strand-first chunk is one quant unit. A strandless chunk is TWO, and the
    orientation split that produces them runs HERE -- inside the chunk, against
    the chunk's own reads, concurrently with every other chunk -- instead of
    once over the whole bam before any of this starts. The two units then run one
    after the other on this chunk's share of the budget rather than competing
    for it, and they share the mini FASTA and mini GTF the single extraction
    wrote.
    """

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

    for unit in chunk["units"]:
        uid = unit["sentinel_id"]
        norm_bam = unit["norm_bam"]
        norm_token = chain_token(
            # maxintron is passed to the normalizer below, so it decides these
            # contents and has to name them. So does min_per_id: the normalizer
            # discards alignments below it before measuring depth, so the same input
            # at a different threshold yields a different normalized bam, and a
            # cache written before it was forwarded must not be reused.
            "stage4_norm_{}.maxcov_{}_dw_{}_seed_{}_origin_{}_maxintron_{}_minperid_{}".format(
                uid,
                args.normalize_max_cov_level,
                args.depth_window,
                args.random_seed,
                chunk["window_origin"],
                args.max_intron_length,
                resolve_min_per_id(args),
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
            # The normalizer discards alignments below this before measuring depth,
            # so it has to be the same value the rest of the run filters on. Its own
            # help says "must match the consumer's min_per_id"; it did not, and the
            # chunked and unchunked arms disagreed on a TSS by 4 bp as a result.
            "--min_per_id",
            str(resolve_min_per_id(args)),
            "--max_intron_length",
            str(args.max_intron_length),
            # true of both modes' units: strand-first because the whole bam was
            # split before extraction, strandless because stage 3b split this
            # chunk above.
            "--input_is_single_strand",
            "--window_origin",
            str(chunk["window_origin"]),
        ]
        if ckpt.done(norm_token):
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
            ckpt.mark(norm_token)

        # LRAA derives its scratch roots by string concatenation on
        # --output_prefix ("__{prefix}.contigtmp", "__{prefix}.sgcache"), so an
        # ABSOLUTE prefix would produce nonsense paths like
        # "__/abs/path.contigtmp" rooted at the cwd. Give it a bare name and let
        # cwd place the outputs.
        quant_prefix = unit["quant_prefix"]
        quant_cmd = lraa_cmd(
            args,
            bam_for_quant=unit["bam"],
            bam_for_sg=norm_bam,
            # ONE mini contig for the pair -- sequence has no orientation and
            # the extraction wrote it once -- but each unit's OWN annotation,
            # because stage 5 quantifies every model its GTF names.
            genome="{}.fa".format(prefix),
            gtf=unit["gtf"],
            out_prefix=unit["quant_name"],
            num_total_reads=num_total_reads,
            cpu_budget=cpu_budget,
        )
        # Built BEFORE the sentinel because the sentinel is keyed on it. The old
        # key named the unit, the read denominator and --HiFi, so quant-only
        # versus discovery was absent from it and so was every flag a caller
        # might later forward -- a resumed run could serve one mode's output for
        # the other, which no downstream check can detect. args.discovery needs
        # no field of its own: it IS the presence of --quant_only in this argv.
        quant_token = chain_token(
            "stage5_quant_{}.argv_{}".format(uid, argv_digest(quant_cmd)[:12]),
            norm_token,
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

    return {
        "chunk_id": cid,
        "region": chunk["region"],
        "units": [u["unit_id"] for u in chunk["units"]],
        "log": log,
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
    two quant jobs and runs them in series on its own share, so the queue is 25
    intervals rather than 50 contig-strand chunks. That coarsens the granularity
    the makespan heuristic has to work with, and it is the one place this mode
    could lose what it gains: with a budget near the interval count, a long tail
    interval now carries both orientations instead of one. It buys the shared
    extraction and a split that is per-chunk rather than whole-bam. Splitting the
    queue back into 50 quant units would need a barrier after every chunk's
    split, or units that appear mid-run; both were rejected as more machinery
    than the granularity is worth at 8-wide over 25 intervals.

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
    bound = chunk_memory_cap(
        args.cpu_budget, [unit.region[1] - unit.region[0] + 1 for unit in units]
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


def _namespace_id(unit_id, value):
    """A chunk-local model id made unique across the run."""

    return "{}|{}".format(unit_id, value)


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
    with open(gtf_out, "wt") as ofh:
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
    }
    if discovery:
        merged.update(merge_discovery_gtf(merged_dir, units))
    return merged


def verify_severed_accounting(cut_dir, chunks):
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
    """

    selected = severed_read_names(cut_dir, exclude=(EXTRACTION_ONLY_DROPS,))
    dropped = set()
    mentions = 0
    for chunk in chunks:
        names = chunk["manifest"]["dropped_read_names"]
        mentions += len(names)
        dropped.update(names)

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
    if unrealized:
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
    return {
        "intervals_split": sum(1 for c in chunks if "split_counts" in c),
        "alignments_emitted": total["alignments_emitted"],
        "records_plus": total["records_plus"],
        "records_minus": total["records_minus"],
        "records_quantified": total["records_total"],
        "records_lost_in_split": total["records_lost"],
    }


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
    norm_token = chain_token(
        # maxintron is passed to the normalizer below; it was missing here, which
        # is what let a new cap rerun the merge and then reuse this bam
        "baseline_norm.maxcov_{}_dw_{}_seed_{}_origin_0_maxintron_{}".format(
            args.normalize_max_cov_level,
            args.depth_window,
            args.random_seed,
            args.max_intron_length,
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
        "baseline_quant.argv_{}".format(argv_digest(quant_cmd)[:12]), norm_token
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
    parser.add_argument("--HiFi", action="store_true", help="pass --HiFi to LRAA")
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
    parser.add_argument(
        "--rss_sample_interval",
        type=float,
        default=0.5,
        help="seconds between /proc RSS samples; a spike shorter than this can "
        "be missed",
    )

    return parser


def parse_args(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    required = ["bam", "genome_fa", "output_dir"]
    if not args.discovery:
        # Quant-only has nothing to quantify without it. Discovery does: no gtf
        # is de novo, which is the mode this option exists to reach.
        required.append("gtf")
    for name in required:
        if not getattr(args, name, None):
            parser.error("--{} is required".format(name))
    if args.cpu_budget is not None and args.cpu_budget < 1:
        parser.error("--cpu_budget must be >= 1")

    return args


def default_args(**overrides):
    """A fully defaulted namespace, with ``overrides`` applied.

    This is how ``LRAA --chunk`` reaches the pipeline: every default lives in
    one parser, so a flag added here is present on both routes rather than only
    on the one whose argparse got edited.
    """

    args = build_parser().parse_args([])
    for key, value in overrides.items():
        if not hasattr(args, key):
            raise PipelineError(
                "no such chunked-pipeline setting: {}".format(key)
            )
        setattr(args, key, value)
    return args


def main(argv=None):

    return run(parse_args(argv))


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

    required = ["bam", "genome_fa", "output_dir"]
    if not args.discovery:
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

    for tool in (SEPARATE_BAM, SELECT_CUTS, EXTRACT_CHUNK, NORMALIZE_BAM, LRAA):
        if not os.path.exists(tool):
            raise PipelineError("missing stage tool: {}".format(tool))
    if shutil.which("samtools") is None:
        raise PipelineError("samtools is not on PATH")

    rss = args.rss_sample_interval

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
            num_total_reads = args.num_total_reads
            timing["num_total_reads_source"] = "supplied via -N"
            if num_total_reads != retained_total:
                raise PipelineError(
                    "-N {} disagrees with the stage-1 retained record count {} "
                    "({} on + plus {} on -). The denominator has to be the record "
                    "set the arms actually see, or TPM is not comparable between "
                    "them.".format(
                        num_total_reads, retained_total, retained["+"], retained["-"]
                    )
                )
    else:
        timing.setdefault("stages", {})["strand_split"] = {
            "skipped": "--strandless_chunks with --arm chunked: the orientation "
            "split runs per chunk, at stage 3b"
        }
        timing["stage1_retained_records"] = None
        # There is no stage-1 count to default to, and inventing one would
        # silently move the RPM_total_reads column relative to a strand-first
        # run of the same substrate -- the one column stage 6 does not rebase.
        # A dry run never reaches quantification, so it does not need the number.
        if args.num_total_reads is None:
            if not args.dry_run:
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
        sources = cut_sources(args, strand_bams, inputs_token, split_token)
        # Stages 2 and 3 are ONE pool now rather than two serial loops, so they
        # are entered together. A dry run selects and stops: the selection is real
        # work whose result the printed plan describes, and extraction is the part
        # a dry run must not do.
        selections, cut_dir, chunks, prep_makespan, prep_allocation = (
            run_prep_concurrently(
                args,
                ckpt,
                outdir,
                timing,
                sources,
                rss,
                select_only=bool(args.dry_run),
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

        # Before the expensive phase, not after it: this is the check that keeps
        # the control's pruned bam and the chunked arm's inputs the same record
        # set, and it needs nothing but the manifests and the cut selection.
        timing["severed_read_accounting"] = verify_severed_accounting(cut_dir, chunks)
        flush()
        print(
            "extracted {} chunk(s): {}".format(
                len(chunks), ", ".join(c["region"] for c in chunks)
            ),
            flush=True,
        )

        load_before = loadavg()
        arm_sampler = RssSampler(os.getpid(), rss)
        arm_sampler.start()
        try:
            chunk_records, makespan, chunk_allocation, chunk_mem = (
                run_chunks_concurrently(
                    args, ckpt, outdir, chunks, num_total_reads, rss
                )
            )
        finally:
            arm_sampler.stop()
        load_after = loadavg()

        units = ordered_units(chunks)
        merged = merge_and_translate(outdir, units, discovery=args.discovery)
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
            "peak_rss_kb_summed_over_chunk_peaks": sum(
                c["peak_tree_rss_kb"] for c in chunk_records
            ),
            "observed_peak_concurrent_tree_rss_kb": arm_sampler.peak_kb,
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
