#!/usr/bin/env python3
# encoding: utf-8

"""Execution regression test for the whole-genome alignment-mismapping filter
(v0.40.0). See README.md.

Each minigenome contains a real, highly-expressed source gene AND the
alignment/strand-mismapping artifacts derived from it (wrong-strand near-mirrors
for the `s` cases, run-on/near-identical copies for the `x` cases). We run LRAA
de novo TWICE on the frozen, de-identified reads -- once with
--no_filter_mismappings (artifacts present) and once with the filter on -- and
assert that the filter REMOVES the opposite-strand artifacts while RETAINING the
real source-strand models. Truth-annotation-independent: the assertions are the
filter-off vs filter-on delta, so an example needs no reference GTF.
"""

import os
import re
import sys
import glob
import shutil
import argparse
import subprocess
from collections import defaultdict

HERE = os.path.dirname(os.path.realpath(__file__))
LRAA = os.path.normpath(os.path.join(HERE, "..", "..", "LRAA"))

# name, source_strand, min_removed, min_source_models
#   artifacts are on the strand OPPOSITE source_strand.
EXAMPLES = [
    ("s1_mirror", "-", 2, 2),
    ("s2_mirror", "-", 2, 2),
    ("x1_runon", "+", 3, 2),
    ("x2_runon", "+", 5, 3),
]


def run_lraa(name, cond, workdir):
    d = os.path.join(workdir, f"{name}.{cond}")
    if os.path.exists(d):
        shutil.rmtree(d)
    os.makedirs(d)
    cmd = [
        LRAA, "--genome", os.path.join(HERE, f"{name}.fa"),
        "--bam", os.path.join(HERE, f"{name}.bam"),
        "--no_chunk", "--output_prefix", "m",
    ]
    if cond == "off":
        cmd.append("--no_filter_mismappings")
    with open(os.path.join(d, "run.log"), "w") as log:
        rc = subprocess.call(cmd, cwd=d, stdout=log, stderr=subprocess.STDOUT)
    if rc != 0:
        raise RuntimeError(f"{name} filter-{cond}: LRAA exited {rc}; see {d}/run.log")
    return d


def model_strands(gtf):
    counts = defaultdict(int)
    if not os.path.exists(gtf):
        return counts
    for line in open(gtf):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) >= 9 and f[2] == "transcript":
            counts[f[6]] += 1
    return counts


def removed_count(log):
    if not os.path.exists(log):
        return 0
    n = 0
    for line in open(log):
        if line.startswith("#") or line.startswith("transcript_id"):
            continue
        if line.strip():
            n += 1
    return n


def tpm_sum(quant):
    s = 0.0
    for line in open(quant):
        if line.startswith("#") or line.startswith("gene_id"):
            continue
        parts = line.split("\t")
        if len(parts) > 6:
            s += float(parts[6])
    return s


def check_example(name, source_strand, min_removed, min_source, workdir):
    art_strand = "-" if source_strand == "+" else "+"
    fails = []

    d_off = run_lraa(name, "off", workdir)
    d_on = run_lraa(name, "on", workdir)
    off = model_strands(os.path.join(d_off, "m.LRAA.ref-free.gtf"))
    on = model_strands(os.path.join(d_on, "m.LRAA.ref-free.gtf"))
    log = os.path.join(d_on, "m.LRAA.ref-free.gtf.mismapping_filter.log")
    n_removed = removed_count(log)
    quant = os.path.join(d_on, "m.LRAA.ref-free.quant.expr")
    tsum = tpm_sum(quant) if os.path.exists(quant) else -1

    # 1. the artifact reproduces with the filter off
    if off[art_strand] < 1:
        fails.append(
            f"expected >=1 {art_strand}-strand artifact model with filter OFF, "
            f"got {off[art_strand]}"
        )
    # 2. the filter reduces the artifact-strand models
    if on[art_strand] >= off[art_strand]:
        fails.append(
            f"filter did not reduce {art_strand}-strand models: "
            f"off={off[art_strand]} on={on[art_strand]}"
        )
    # 3. it removed at least the expected number of artifacts
    if n_removed < min_removed:
        fails.append(f"mismapping log removed {n_removed} < expected {min_removed}")
    # 4. the real source-strand models are retained (unchanged, and >= expected)
    if on[source_strand] != off[source_strand]:
        fails.append(
            f"source-strand ({source_strand}) model count changed: "
            f"off={off[source_strand]} on={on[source_strand]} (filter must not touch real models)"
        )
    if on[source_strand] < min_source:
        fails.append(
            f"expected >={min_source} retained {source_strand}-strand models, got {on[source_strand]}"
        )
    # 5. quant renormalized to 1e6
    if not (999999.0 <= tsum <= 1000001.0):
        fails.append(f"filtered quant.expr TPM sum {tsum:.1f} != 1e6")

    status = "PASS" if not fails else "FAIL"
    print(
        f"[{status}] {name}: OFF {dict(off)} -> ON {dict(on)}; "
        f"removed={n_removed}; TPMsum={tsum:.0f}"
    )
    for f in fails:
        print(f"    - {f}")
    return not fails


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=os.path.join(HERE, "__test_runs"))
    ap.add_argument("--only", default=None, help="run one example by name")
    args = ap.parse_args()
    os.makedirs(args.workdir, exist_ok=True)

    if not os.path.exists(LRAA):
        sys.exit(f"LRAA driver not found at {LRAA}")

    ok = True
    for name, strand, minrem, minsrc in EXAMPLES:
        if args.only and name != args.only:
            continue
        ok = check_example(name, strand, minrem, minsrc, args.workdir) and ok
    if not ok:
        sys.exit("MISMAPPING FILTER REGRESSION TEST FAILED")
    print("All mismapping-filter minigenome checks passed.")


if __name__ == "__main__":
    main()
