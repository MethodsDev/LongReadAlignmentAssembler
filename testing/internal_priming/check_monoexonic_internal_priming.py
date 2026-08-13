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

Also asserts that an emitted model whose own terminus sits in A-rich context still
carries InternalPriming "True": the annotation is preserved, only the deletion moved.

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

    cmd = [sys.executable, args.lraa,
           "--genome", genome, "--bam", bam,
           "--output_prefix", "ip_mono",
           "--num_threads_per_worker", "1",
           "--min_mapping_quality", "0",
           "--min_mapping_quality_for_final_quant", "0",
           "--HiFi",
           # --debug so the graph writes __PolyAsite_info.bed; the gate's own log line
           # is emitted from a worker process and never reaches this stdout
           "--debug"]
    print("running: " + " ".join(cmd))
    proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    log = proc.stdout + proc.stderr
    with open(os.path.join(workdir, "run.log"), "w") as fh:
        fh.write(log)
    if proc.returncode != 0:
        print(log[-4000:])
        sys.exit(f"LRAA exited {proc.returncode}; log in {workdir}/run.log")

    failures = []

    # --- 1. the gate: no PolyA VERTEX at the A-rich locus, one at the clean locus
    verts = polyA_vertices(workdir)
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

    # --- 2. no PolyA vertex survives at the A-rich locus, and one does at the clean locus
    gtf = os.path.join(workdir, "ip_mono.LRAA.ref-free.gtf")
    if not os.path.exists(gtf):
        failures.append(f"no GTF emitted at {gtf}")
    else:
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
        print(f"emitted models: {ends}")
        print(f"InternalPriming flags: {flags}")

        # a model terminating within the capture window of the A-rich site would mean the
        # vertex survived; 25 nt is max_dist_between_alt_polyA_sites/2
        snapped_to_A = [e for e in ends if e[1] == "+" and abs(e[2] - A_LOCUS_END) <= 25]
        snapped_to_clean = [
            e for e in ends if e[1] == "+" and abs(e[2] - CLEAN_LOCUS_END) <= 25]
        if snapped_to_clean == []:
            failures.append(
                "control failed: no model terminates near the CLEAN locus, so the test "
                "cannot distinguish rejection from the locus simply not being assembled")
        # annotation must survive for anything that does end in A-rich context
        for tid, strand, three in snapped_to_A:
            if not flags.get(tid):
                failures.append(
                    f"model {tid} ends at {three}, in A-rich context, but is not "
                    f"annotated InternalPriming True -- the annotation was lost")
        if snapped_to_A:
            print(f"note: {len(snapped_to_A)} model(s) still terminate near the A-rich "
                  f"site; they are annotated, which is the intended behaviour now that "
                  f"the transcript stage no longer deletes")

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
