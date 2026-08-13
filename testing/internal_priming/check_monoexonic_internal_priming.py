#!/usr/bin/env python3

"""End-to-end exercise of the monoexonic internal-priming branch.

chr20 emits zero flagged monoexonic models, so the monoexonic path was never exercised
by any real-data arm. This builds the locus deliberately: a synthetic contig with two
monoexonic pileups that differ ONLY in what the genome does immediately 3' of where the
reads end.

  locus A, reads ending at 1500, followed by a 20 nt genomic A-run
          -> the PolyA candidate is internal priming and must be rejected at
             identification, so no PolyA vertex may exist there
  locus B, reads ending at 4500, followed by ordinary sequence
          -> the candidate must survive, which is the control that proves the rejection
             is about genomic context and not about being monoexonic

Two arms are run over the identical locus, and their outcomes must differ:

  default                          -> a monoexonic model terminating in A-rich context
                                      is DELETED by the transcript stage
  --no_filter_internal_priming     -> the same model is EMITTED, carrying
                                      InternalPriming "True"

The second arm is what makes the first meaningful: it proves the model is reconstructed
at all, so its absence under the default is a deletion rather than a locus that simply
never assembled.

  check_monoexonic_internal_priming.py [--keep] [--lraa /path/to/LRAA]
"""

import argparse
import os
import random
import subprocess
import sys
import tempfile

import pysam

CONTIG = "synth1"
CONTIG_LEN = 6000
A_LOCUS_END = 1500      # reads end here, genomic A-run follows
CLEAN_LOCUS_END = 4500  # reads end here, ordinary sequence follows
READ_LEN = 500
N_READS = 40


def build_genome(path):
    random.seed(20260813)
    bases = []
    for i in range(1, CONTIG_LEN + 1):
        if A_LOCUS_END < i <= A_LOCUS_END + 20:
            bases.append("A")            # the internal-priming template
        elif CLEAN_LOCUS_END < i <= CLEAN_LOCUS_END + 20:
            bases.append("CGTCGTCGTG"[(i - CLEAN_LOCUS_END - 1) % 10])
        else:
            bases.append(random.choice("CGT"))   # no incidental A-runs anywhere else
    seq = "".join(bases)
    assert seq[A_LOCUS_END : A_LOCUS_END + 20] == "A" * 20
    with open(path, "w") as fh:
        fh.write(f">{CONTIG}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i : i + 60] + "\n")
    pysam.faidx(path)
    return seq


def build_bam(path, seq):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}]}
    reads = []
    for locus_end in (A_LOCUS_END, CLEAN_LOCUS_END):
        start = locus_end - READ_LEN + 1          # 1-based inclusive
        for i in range(N_READS):
            a = pysam.AlignedSegment()
            a.query_name = f"r{locus_end}_{i}"
            a.reference_id = 0
            a.reference_start = start - 1          # 0-based
            a.mapping_quality = 60
            a.cigartuples = [(0, READ_LEN)]       # a single M block: monoexonic, no clips
            a.query_sequence = seq[start - 1 : locus_end]
            a.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
            a.flag = 0
            reads.append(a)
    unsorted = path + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for a in reads:
            out.write(a)
    pysam.sort("-o", path, unsorted)
    pysam.index(path)
    os.unlink(unsorted)


def polyA_vertices(workdir):
    """PolyA vertex coordinates from __PolyAsite_info.bed.

    Asserting on the graph artifact rather than on a log line: Splice_graph logs from
    worker processes whose output never reaches the driver's stdout, so the gate's own
    INFO line is not observable from here. The bed file is appended across graph builds
    and strands, so this is the union -- which is what "does a vertex exist here" wants.
    """
    path = os.path.join(workdir, "__PolyAsite_info.bed")
    if not os.path.exists(path):
        return None
    coords = set()
    for line in open(path):
        f = line.split("\t")
        if len(f) >= 3:
            coords.add(int(f[2]))          # bed end == the site coordinate
    return sorted(coords)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lraa", default=os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "LRAA"))
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args()

    workdir = tempfile.mkdtemp(prefix="ip_mono_", dir=os.getcwd())
    print(f"workdir {workdir}")
    genome = os.path.join(workdir, "synth.fa")
    bam = os.path.join(workdir, "synth.bam")
    seq = build_genome(genome)
    build_bam(bam, seq)

    def run(prefix, extra=()):
        # each arm gets its own directory: with --debug, LRAA writes graph artifacts
        # whose names are not qualified by --output_prefix (__MPGN_components_described
        # .*.bed), and it refuses to overwrite them
        rundir = os.path.join(workdir, prefix)
        os.makedirs(rundir, exist_ok=True)
        cmd = [sys.executable, args.lraa,
               "--genome", genome, "--bam", bam,
               "--output_prefix", prefix,
               "--num_threads_per_worker", "1",
               "--min_mapping_quality", "0",
               "--min_mapping_quality_for_final_quant", "0",
               "--HiFi",
               # --debug so the graph writes __PolyAsite_info.bed; the gate's own log
               # line is emitted from a worker process and never reaches this stdout
               "--debug", *extra]
        print("running: " + " ".join(cmd))
        proc = subprocess.run(cmd, cwd=rundir, capture_output=True, text=True)
        log = proc.stdout + proc.stderr
        with open(os.path.join(rundir, "run.log"), "w") as fh:
            fh.write(log)
        if proc.returncode != 0:
            print(log[-4000:])
            sys.exit(f"LRAA exited {proc.returncode}; log in {rundir}/run.log")

    run("ip_mono")

    failures = []

    # --- 1. the gate: no PolyA VERTEX at the A-rich locus, one at the clean locus
    verts = polyA_vertices(os.path.join(workdir, "ip_mono"))
    if verts is None:
        failures.append("no __PolyAsite_info.bed written; cannot inspect PolyA vertices")
    else:
        near_A = [v for v in verts if abs(v - A_LOCUS_END) <= 25]
        near_clean = [v for v in verts if abs(v - CLEAN_LOCUS_END) <= 25]
        print(f"PolyA vertices: {verts}")
        print(f"  near the A-rich locus ({A_LOCUS_END}): {near_A}")
        print(f"  near the clean locus  ({CLEAN_LOCUS_END}): {near_clean}")
        if near_A:
            failures.append(
                f"A-rich locus still has PolyA vertex/vertices {near_A}: the "
                f"identification-stage gate did not reject the internally primed "
                f"candidate")
        if not near_clean:
            failures.append(
                f"CONTROL FAILED: no PolyA vertex near the clean locus either, so the "
                f"absence at {A_LOCUS_END} does not demonstrate the gate -- it may just "
                f"be that no PolyA site was called anywhere")

    # --- 2. the deletion, and the flag that disables it
    def emitted(prefix):
        gtf = os.path.join(workdir, prefix, prefix + ".LRAA.ref-free.gtf")
        if not os.path.exists(gtf):
            failures.append(f"no GTF emitted at {gtf}")
            return None, None
        ends, flags = [], {}
        for line in open(gtf):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "transcript":
                continue
            tid = [kv for kv in f[8].split(";") if "transcript_id" in kv]
            tid = tid[0].split('"')[1] if tid else "?"
            three = int(f[4]) if f[6] == "+" else int(f[3])
            ends.append((tid, f[6], three))
            flags[tid] = 'InternalPriming "True"' in f[8]
        return ends, flags

    # 25 nt is max_dist_between_alt_polyA_sites/2, the capture window of a site
    def near(ends, locus):
        return [e for e in ends if e[1] == "+" and abs(e[2] - locus) <= 25]

    ends, flags = emitted("ip_mono")
    if ends is not None:
        print(f"default arm, emitted models: {ends}")
        print(f"default arm, InternalPriming flags: {flags}")
        if near(ends, CLEAN_LOCUS_END) == []:
            failures.append(
                "control failed: no model terminates near the CLEAN locus, so the test "
                "cannot distinguish deletion from the locus never being assembled")
        surviving = near(ends, A_LOCUS_END)
        if surviving:
            failures.append(
                f"monoexonic model(s) {[e[0] for e in surviving]} terminate in A-rich "
                f"context and were NOT deleted; the transcript stage should remove them")

    # the same locus with deletion disabled: the model must reappear, annotated
    run("ip_mono_nofilt", extra=("--no_filter_internal_priming",))
    ends_nf, flags_nf = emitted("ip_mono_nofilt")
    if ends_nf is not None:
        print(f"--no_filter_internal_priming arm, emitted models: {ends_nf}")
        print(f"--no_filter_internal_priming arm, flags: {flags_nf}")
        kept = near(ends_nf, A_LOCUS_END)
        if not kept:
            failures.append(
                "with --no_filter_internal_priming no model terminates near the A-rich "
                "locus either, so the default arm's absence does not demonstrate the "
                "deletion -- the model may never be reconstructed at all")
        for tid, strand, three in kept:
            if not flags_nf.get(tid):
                failures.append(
                    f"model {tid} ends at {three}, in A-rich context, but is not "
                    f"annotated InternalPriming True -- the annotation was lost")

    if not args.keep:
        import shutil
        shutil.rmtree(workdir, ignore_errors=True)
    else:
        print(f"kept {workdir}")

    if failures:
        for f in failures:
            print("FAIL: " + f)
        sys.exit(1)
    print("PASS: monoexonic internal-priming branch exercised end to end")


if __name__ == "__main__":
    main()
