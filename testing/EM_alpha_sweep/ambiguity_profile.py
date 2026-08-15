#!/usr/bin/env python3
"""What prior does EM_alpha pull toward, and does it resemble the truth?

The M-step is theta_t proportional to (counts_t + alpha * ambiguous_t). As
alpha -> infinity the counts term becomes negligible and theta approaches
ambiguous_t / sum(ambiguous), so:

    alpha interpolates theta between the likelihood solution at alpha=0 and the
    NORMALIZED AMBIGUITY PROFILE at alpha -> infinity.

alpha is therefore a mixing weight between the data's own answer and one
specific named prior, and it should help exactly to the degree that prior
resembles the truth. That is a falsifiable prediction and it needs no sweep to
evaluate: the profile is computable from read compatibility alone.

The profile is built from LRAA's own quant.tracking, not re-derived. For every
read, the set of transcripts it was assigned to; reads whose set has more than
one member contribute 1 to EVERY member. That is exactly what EM.py:284-287
accumulates into ambiguous_read_counts, and it is independent of alpha, since
alpha does not change which transcripts a read is compatible with.

Similarity to truth is measured with the SAME symmetric relative difference the
benchmark uses for MARD, so profile-truth distance is on the footing of the
objective rather than an unrelated divergence, plus Spearman for the ordering.

Two modes:
  --run       produce a tracking file per sample (one arm, alpha irrelevant)
  --analyze   build profiles, compare to truth, write the table
"""

import argparse
import csv
import gzip
import os
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep_config as C  # noqa: E402

PROBE = os.path.join(C.BASE, "probe", "profile")
EPS = 1e-9


def run_one(smp, cpu_budget):
    d = os.path.join(PROBE, smp["sample"])
    trk = os.path.join(d, "prof.LRAA.quant-only.quant.tracking.gz")
    if os.path.exists(trk):
        return smp["sample"], "cached"
    os.makedirs(d, exist_ok=True)
    cmd = [
        C.LRAA, "--genome", smp["genome"], "--bam", smp["bam"], "--gtf", smp["gtf"],
        "--quant_only", "--min_mapping_quality", "0",
        "--min_mapping_quality_for_final_quant", "0",
        "--EM_alpha", "0.01", "--cpu_budget", str(cpu_budget), "--output_prefix", "prof",
    ]
    if smp["hifi"]:
        cmd.append("--HiFi")
    with open(os.path.join(d, "run.log"), "w") as fh:
        fh.write("# " + " ".join(cmd) + "\n")
        fh.flush()
        rc = subprocess.call(cmd, cwd=d, stdout=fh, stderr=subprocess.STDOUT)
    return smp["sample"], ("ok" if rc == 0 and os.path.exists(trk) else f"FAIL rc={rc}")


def profile_from_tracking(path):
    """ambiguous_read_counts per transcript, exactly as EM accumulates it.

    Grouped by read rather than by mp_id: mp ids are only unique within a
    component, and a read's candidate set is what defines ambiguity anyway.
    """
    read_to_tx = defaultdict(set)
    with gzip.open(path, "rt") as fh:
        rows = (ln for ln in fh if not ln.startswith("#"))
        for r in csv.DictReader(rows, delimiter="\t"):
            read_to_tx[r["read_name"]].add(r["transcript_id"])
    amb = defaultdict(float)
    n_amb_reads = n_reads = 0
    for _read, txs in read_to_tx.items():
        n_reads += 1
        if len(txs) > 1:
            n_amb_reads += 1
            for t in txs:
                amb[t] += 1.0
    return amb, n_reads, n_amb_reads


def load_truth(path):
    d = pd.read_csv(path, sep="\t")
    idc, vc = d.columns[0], d.columns[1]
    d = d[d[vc] > 0]
    return dict(zip(d[idc].astype(str), d[vc].astype(float)))


def symmetric_mard(pred, true):
    """mean(|p-t| / (p+t+eps)) over the union, the benchmark's own MARD form."""
    keys = sorted(set(pred) | set(true))
    p = np.array([pred.get(k, 0.0) for k in keys], dtype=float)
    t = np.array([true.get(k, 0.0) for k in keys], dtype=float)
    p = p / p.sum() * 1e6 if p.sum() else p
    t = t / t.sum() * 1e6 if t.sum() else t
    return float(np.mean(np.abs(p - t) / (p + t + EPS))), p, t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--analyze", action="store_true")
    ap.add_argument("--samples")
    ap.add_argument("--cpu_budget", type=int, default=3)
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--out", default=os.path.join(C.BASE, "results", "ambiguity_profile.tsv"))
    args = ap.parse_args()

    samples = (
        [C.SAMPLES_BY_NAME[s] for s in args.samples.split(",")]
        if args.samples else list(C.SAMPLES)
    )

    if args.run:
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            futs = [ex.submit(run_one, s, args.cpu_budget) for s in samples]
            for f in as_completed(futs):
                print("  %s %s" % f.result(), flush=True)

    if not args.analyze:
        return 0

    rows = []
    for smp in samples:
        trk = os.path.join(PROBE, smp["sample"], "prof.LRAA.quant-only.quant.tracking.gz")
        if not os.path.exists(trk):
            continue
        amb, n_reads, n_amb = profile_from_tracking(trk)
        truth = load_truth(smp["truth_quant"])
        mard_pt, p, t = symmetric_mard(amb, truth)
        # A uniform prior over the same support, as the reference point: if the
        # profile is no closer to truth than "spread it evenly", it carries no
        # information the truth cares about.
        uni = {k: 1.0 for k in amb}
        mard_ut, _, _ = symmetric_mard(uni, truth)
        common = sorted(set(amb) & set(truth))
        sp = (
            spearmanr([amb[k] for k in common], [truth[k] for k in common]).statistic
            if len(common) >= 3 else float("nan")
        )
        rows.append(dict(
            sample=smp["sample"], corpus=smp["corpus"], platform=smp["platform"],
            n_reads_assigned=n_reads, n_reads_ambiguous=n_amb,
            frac_reads_ambiguous=n_amb / n_reads if n_reads else float("nan"),
            n_tx_in_profile=len(amb), n_tx_truth=len(truth), n_tx_common=len(common),
            mard_profile_vs_truth=mard_pt,
            mard_uniform_vs_truth=mard_ut,
            profile_beats_uniform=mard_pt < mard_ut,
            spearman_profile_vs_truth=sp,
        ))

    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    pd.set_option("display.width", 220, "display.max_columns", 40)
    print(out.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
