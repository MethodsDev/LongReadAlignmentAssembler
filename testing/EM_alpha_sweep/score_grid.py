#!/usr/bin/env python3
"""Score whatever arms are currently staged in each per-sample scoring dir.

Runs bmark_nb_runner.py once per sample.  The scorer needs at least two
registry entries present in a directory (its scatterplot cell indexes the
subplot array, which matplotlib collapses to a bare Axes at n=1), so a sample
with fewer than two staged arms is skipped rather than crashed.

Emits nothing itself beyond the scorer's own outputs; collect_results.py turns
those into the comparable table.
"""

import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep_config as C  # noqa: E402

RUNNER = (
    "/home/unix/bhaas/projects/LRAA_PAPER_Analyses/iso-reconstruct-benchmark2"
    "/Isosceles_bmarking_routines/bmark_nb_runner.py"
)


def score_one(smp, registry):
    d = os.path.join(C.SCORE, smp["sample"])
    raw = os.path.join(d, "raw_prog_results")
    n = len(os.listdir(raw)) if os.path.isdir(raw) else 0
    if n < 2:
        return smp["sample"], 0.0, f"skip ({n} arms staged, scorer needs >=2)"
    cmd = [
        RUNNER,
        "--analysisType", "quant_only",
        "--registry", registry,
        "--truth_gtf", smp["truth_gtf"],
        "--truth_quant", smp["truth_quant"],
    ]
    t0 = time.time()
    with open(os.path.join(d, "score.log"), "w") as fh:
        fh.write("# " + " ".join(cmd) + "\n")
        fh.flush()
        rc = subprocess.call(cmd, cwd=d, stdout=fh, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    ok = os.path.exists(os.path.join(d, "oarfish_style_metrics.tsv"))
    return smp["sample"], dt, ("ok %d arms" % n) if (rc == 0 and ok) else f"FAIL rc={rc}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", help="comma list (default: all staged)")
    ap.add_argument("--registry", default=os.path.join(C.BASE, "registry.yaml"))
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()

    samples = (
        [C.SAMPLES_BY_NAME[s] for s in args.samples.split(",")]
        if args.samples
        else [s for s in C.SAMPLES if os.path.isdir(os.path.join(C.SCORE, s["sample"]))]
    )
    nfail = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(score_one, s, args.registry) for s in samples]
        for f in as_completed(futs):
            sample, dt, status = f.result()
            if status.startswith("FAIL"):
                nfail += 1
            print(f"{sample:24s} {dt:6.1f}s  {status}", flush=True)
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main())
