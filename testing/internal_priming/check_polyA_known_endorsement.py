#!/usr/bin/env python3

"""End-to-end exercise of the internal-priming veto's ENDORSEMENT SOURCE (v0.41.0).

The veto rejects a read-derived PolyA candidate whose downstream genome is A-rich, unless
a trusted 3' end endorses cleavage there (spare_polyA_veto_at_known_3prime). Before
v0.41.0 the endorsement source was every --gtf transcript terminus, which let a de-novo
guide (e.g. the cluster-guided init GTF) endorse its OWN A-rich internal-priming termini
and re-bless them as PolyA vertices in the guided build -- exactly the artifacts a plain
de-novo run rejects. v0.41.0: (a) a --gtf terminus flagged InternalPriming no longer
endorses, and (b) a --polyA_known BED/GTF can supply the endorsement set explicitly.

One synthetic locus: monoexonic reads ending where a 20 nt genomic A-run begins, so the
candidate there is internal priming. Arms differ ONLY in the endorsement source:

  ref-free                    -> vetoed (no PolyA vertex at the A-rich locus)   [baseline]
  --gtf guide (clean 3' end)  -> SPARED  (vertex survives)      [mechanism still fires]
  --gtf guide (InternalPriming "True" on the same 3' end)
                              -> vetoed  (vertex gone)          [THE v0.41.0 FIX]
  --polyA_known bed @ locus   -> SPARED  (vertex survives)      [the new knob]

A clean control locus (ordinary downstream sequence) must carry a PolyA vertex in every
arm, so "no vertex at the A-rich locus" means "vetoed" and not "nothing assembled".

  check_polyA_known_endorsement.py [--keep] [--lraa /path/to/LRAA]
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
A_LOCUS_END = 1500      # reads end here, genomic A-run follows -> internal priming
CLEAN_LOCUS_END = 4500  # reads end here, ordinary sequence follows -> real cleavage
N_READS = 40

# Reads are MULTI-EXONIC (one intron), so the candidate goes through the spliced (ME)
# graph where the internal-priming veto actually REJECTS -- under the default
# reject_internally_primed_polyA_sites="spliced_only" the monoexonic (SE) graph would
# instead keep it for deferred filtering, and the graph-level reprieve under test here
# never engages. Each locus: exon1, a 100 nt intron, then a terminal exon ending AT the
# locus.
A_EXONS = [(1000, 1300), (1401, A_LOCUS_END)]           # intron 1301-1400
CLEAN_EXONS = [(4000, 4300), (4401, CLEAN_LOCUS_END)]   # intron 4301-4400


def build_genome(path):
    random.seed(20260924)
    bases = []
    for i in range(1, CONTIG_LEN + 1):
        if A_LOCUS_END < i <= A_LOCUS_END + 20:
            bases.append("A")                                  # internal-priming template
        elif CLEAN_LOCUS_END < i <= CLEAN_LOCUS_END + 20:
            bases.append("CGTCGTCGTG"[(i - CLEAN_LOCUS_END - 1) % 10])
        else:
            bases.append(random.choice("CGT"))                 # no incidental A-runs
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
    for tag, exons in (("A", A_EXONS), ("clean", CLEAN_EXONS)):
        (e1s, e1e), (e2s, e2e) = exons
        e1len, e2len, intron = e1e - e1s + 1, e2e - e2s + 1, e2s - e1e - 1
        matched = seq[e1s - 1 : e1e] + seq[e2s - 1 : e2e]
        for i in range(N_READS):
            a = pysam.AlignedSegment()
            a.query_name = f"r{tag}_{i}"
            a.reference_id = 0
            a.reference_start = e1s - 1                       # 0-based
            a.mapping_quality = 60
            a.cigartuples = [(0, e1len), (3, intron), (0, e2len)]  # M N M
            a.query_sequence = matched
            a.query_qualities = pysam.qualitystring_to_array("I" * len(matched))
            a.flag = 0
            reads.append(a)
    unsorted = path + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for a in reads:
            out.write(a)
    pysam.sort("-o", path, unsorted)
    pysam.index(path)
    os.unlink(unsorted)


def _write_guide(path, internal_primed):
    # Multi-exonic, matching the reads' A-locus structure, so the guide is routed to the
    # spliced (ME) graph round -- the same round the multi-exon-read PolyA candidate is
    # built in. A monoexonic guide would go to the SE round and never endorse it.
    attrs = 'gene_id "refg"; transcript_id "reft";'
    if internal_primed:
        attrs += ' InternalPriming "True";'
    (e1s, e1e), (e2s, e2e) = A_EXONS
    with open(path, "w") as fh:
        print("\t".join([CONTIG, "ref", "transcript", str(e1s), str(e2e), ".",
                         "+", ".", attrs]), file=fh)
        for (s, e) in A_EXONS:
            print("\t".join([CONTIG, "ref", "exon", str(s), str(e), ".",
                             "+", ".", attrs]), file=fh)


def polyA_vertices(workdir):
    """PolyA vertex coordinates (bed end) from the --debug __PolyAsite_info.bed."""
    path = os.path.join(workdir, "__PolyAsite_info.bed")
    if not os.path.exists(path):
        return None
    coords = set()
    for line in open(path):
        f = line.split("\t")
        if len(f) >= 3:
            coords.add(int(f[2]))
    return sorted(coords)


def main():
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap = argparse.ArgumentParser()
    ap.add_argument("--lraa", default=os.path.join(here, "LRAA"))
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args()

    workdir = tempfile.mkdtemp(prefix="polyA_known_", dir=os.getcwd())
    print(f"workdir {workdir}")
    genome = os.path.join(workdir, "synth.fa")
    bam = os.path.join(workdir, "synth.bam")
    seq = build_genome(genome)
    build_bam(bam, seq)

    guide_clean = os.path.join(workdir, "guide_clean.gtf")
    guide_ip = os.path.join(workdir, "guide_ip.gtf")
    _write_guide(guide_clean, internal_primed=False)
    _write_guide(guide_ip, internal_primed=True)

    known_bed = os.path.join(workdir, "known.bed")
    with open(known_bed, "w") as fh:
        # a trusted cleavage site AT the A-rich locus, on '+'
        print("\t".join([CONTIG, str(A_LOCUS_END - 1), str(A_LOCUS_END),
                          f"k:{A_LOCUS_END}:+", ".", "+"]), file=fh)

    def run(prefix, extra=()):
        rundir = os.path.join(workdir, prefix)
        os.makedirs(rundir, exist_ok=True)
        cmd = [sys.executable, args.lraa,
               "--genome", genome, "--bam", bam,
               "--output_prefix", prefix,
               "--cpu_budget", "1",
               "--min_mapping_quality", "0",
               "--min_mapping_quality_for_final_quant", "0",
               "--HiFi", "--no_chunk", "--no_stream_reads", "--debug", *extra]
        proc = subprocess.run(cmd, cwd=rundir, capture_output=True, text=True)
        with open(os.path.join(rundir, "run.log"), "w") as fh:
            fh.write(proc.stdout + proc.stderr)
        if proc.returncode != 0:
            print((proc.stdout + proc.stderr)[-4000:])
            sys.exit(f"LRAA exited {proc.returncode}; log in {rundir}/run.log")
        return rundir

    arms = {
        "reffree":     run("reffree"),
        "guide_clean": run("guide_clean", extra=("--gtf", guide_clean)),
        "guide_ip":    run("guide_ip", extra=("--gtf", guide_ip)),
        "polyA_known": run("polyA_known", extra=("--polyA_known", known_bed)),
    }

    # SPARED means a PolyA vertex survives at the A-rich locus; vetoed means none does.
    expected_spared = {
        "reffree": False,
        "guide_clean": True,    # a clean guide terminus endorses -> mechanism fires
        "guide_ip": False,      # v0.41.0: an InternalPriming-flagged guide must NOT endorse
        "polyA_known": True,    # the explicit trusted list endorses
    }

    failures = []
    for arm, rundir in arms.items():
        verts = polyA_vertices(rundir)
        if verts is None:
            failures.append(f"[{arm}] no __PolyAsite_info.bed written")
            continue
        near_A = [v for v in verts if abs(v - A_LOCUS_END) <= 25]
        near_clean = [v for v in verts if abs(v - CLEAN_LOCUS_END) <= 25]
        print(f"[{arm}] vertices={verts} nearA={near_A} nearClean={near_clean}")
        if not near_clean:
            failures.append(
                f"[{arm}] CONTROL FAILED: no PolyA vertex near the clean locus, so an "
                f"absence at the A-rich locus would not prove a veto")
        spared = bool(near_A)
        if spared != expected_spared[arm]:
            failures.append(
                f"[{arm}] A-rich PolyA vertex spared={spared}, expected "
                f"{expected_spared[arm]} (near_A={near_A})")

    if failures:
        print("\nFAIL:")
        for f in failures:
            print("  - " + f)
        sys.exit(1)

    print("\nPASS: endorsement source behaves as specified "
          "(clean guide & --polyA_known spare; InternalPriming-flagged guide does not).")
    if not args.keep:
        import shutil
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
