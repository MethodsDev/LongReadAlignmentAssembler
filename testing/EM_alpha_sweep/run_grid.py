#!/usr/bin/env python3
"""Run the EM_alpha x 3'-weighting quant-only grid locally.

Measurement scaffolding.  Nothing here changes LRAA behaviour: every arm is a
plain quant-only invocation differing only in --EM_alpha and the presence of
--no_weight_reads_by_3prime_agreement.

Each finished run leaves exactly one artifact behind -- the quant.expr, hard
linked into the per-sample scoring directory.  The tracking file is the only
output that scales with library size and nothing downstream reads it, so it is
deleted as soon as the run succeeds; 180 runs of it would be tens of GB.

Usage:
  run_grid.py --samples CL_HG002_E2_sirv --alphas 0,0.001,0.01,0.1,1.0 --threep on
  run_grid.py                       # the whole 180-run grid
  run_grid.py --dry_run
"""

import argparse
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep_config as C  # noqa: E402


def run_dir(sample, arm):
    return os.path.join(C.WORK, sample, arm)


def expected_expr(sample, arm):
    return os.path.join(run_dir(sample, arm), f"{sample}.LRAA.{arm}.LRAA.quant-only.quant.expr")


def score_link(sample, arm):
    return os.path.join(
        C.SCORE, sample, "raw_prog_results", f"{sample}.LRAA.{arm}.LRAA.quant-only.quant.expr"
    )


def build_cmd(smp, armd, cpu_budget, bam=None, prefix=None, extra=()):
    cmd = [
        C.LRAA,
        "--genome", smp["genome"],
        "--bam", bam or smp["bam"],
        "--gtf", smp["gtf"],
        "--quant_only",
        "--min_mapping_quality", "0",
        "--min_mapping_quality_for_final_quant", "0",
        "--EM_alpha", repr(armd["alpha"]),
        "--cpu_budget", str(cpu_budget),
        "--output_prefix", prefix or f'{smp["sample"]}.LRAA.{armd["arm"]}',
    ]
    if smp["hifi"]:
        cmd.append("--HiFi")
    if not armd["threep"]:
        cmd.append("--no_weight_reads_by_3prime_agreement")
    # Execution-only flags (e.g. --no_parallelize_contigs), applied uniformly to
    # every arm of a batch so they cannot become a hidden difference between
    # arms that get compared.  The full command line is written to run.log, so
    # what an arm actually ran with is recoverable from the artifact.
    cmd.extend(extra)
    return cmd


def one_run(smp, armd, cpu_budget, keep_tracking=False, extra=()):
    sample, arm = smp["sample"], armd["arm"]
    d = run_dir(sample, arm)
    expr = expected_expr(sample, arm)
    link = score_link(sample, arm)
    if os.path.exists(expr):
        _relink(expr, link)
        return sample, arm, 0.0, "cached"

    os.makedirs(d, exist_ok=True)
    os.makedirs(os.path.dirname(link), exist_ok=True)
    cmd = build_cmd(smp, armd, cpu_budget, extra=extra)
    t0 = time.time()
    with open(os.path.join(d, "run.log"), "w") as fh:
        fh.write("# " + " ".join(cmd) + "\n")
        fh.flush()
        rc = subprocess.call(cmd, cwd=d, stdout=fh, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    if rc != 0 or not os.path.exists(expr):
        return sample, arm, dt, f"FAIL rc={rc}"
    if not keep_tracking:
        for f in os.listdir(d):
            if ".quant.tracking" in f or f.startswith("__"):
                p = os.path.join(d, f)
                shutil.rmtree(p, ignore_errors=True) if os.path.isdir(p) else os.remove(p)
    _relink(expr, link)
    return sample, arm, dt, "ok"


def _relink(src, dst):
    if os.path.exists(dst) or os.path.islink(dst):
        os.remove(dst)
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", help="comma list of sample names (default: all 10)")
    ap.add_argument("--alphas", help="comma list of alpha values (default: the 9-point grid)")
    ap.add_argument("--threep", default="both", choices=["on", "off", "both"])
    ap.add_argument("--cpu_budget", type=int, default=3)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--keep_tracking", action="store_true")
    ap.add_argument("--dry_run", action="store_true")
    ap.add_argument("--extra_flags", default="",
                    help="space-separated LRAA flags appended to EVERY arm in this batch, "
                         "e.g. '--no_parallelize_contigs'")
    args = ap.parse_args()

    samples = (
        [C.SAMPLES_BY_NAME[s] for s in args.samples.split(",")]
        if args.samples
        else list(C.SAMPLES)
    )
    alphas = [float(a) for a in args.alphas.split(",")] if args.alphas else list(C.ALPHAS)
    threeps = {"on": [True], "off": [False], "both": [True, False]}[args.threep]

    jobs = []
    for smp in samples:
        for w in threeps:
            for a in alphas:
                jobs.append((smp, dict(arm=C.arm_name(a, w), alpha=a, threep=w)))

    # Longest first: the two MORF samples dominate the makespan, so starting
    # them last would leave three cores idle at the tail.
    jobs.sort(key=lambda j: 0 if j[0]["corpus"] == "MORF" else 1)

    print(f"{len(jobs)} runs, {args.jobs}-way, cpu_budget={args.cpu_budget}", flush=True)
    if args.dry_run:
        for smp, armd in jobs:
            print(" ".join(build_cmd(smp, armd, args.cpu_budget, extra=args.extra_flags.split())))
        return 0

    t0 = time.time()
    nfail = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [
            ex.submit(one_run, smp, armd, args.cpu_budget, args.keep_tracking,
                      tuple(args.extra_flags.split()))
            for smp, armd in jobs
        ]
        for i, f in enumerate(as_completed(futs), 1):
            sample, arm, dt, status = f.result()
            if status.startswith("FAIL"):
                nfail += 1
            print(
                f"[{i}/{len(jobs)}] {sample} {arm} {dt:7.1f}s {status} "
                f"(elapsed {time.time() - t0:.0f}s)",
                flush=True,
            )
    print(f"DONE in {time.time() - t0:.0f}s, {nfail} failures")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
