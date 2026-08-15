#!/usr/bin/env python3
"""Collapse every scored arm into one long table.

One row per (sample, arm).  Columns are the scorer's own metrics plus three
dropout counts computed here from the same matched table the metrics were
computed on:

  n_pred_zero      truth-expressed entries the arm predicted at exactly 0.
                   This is the literal false-negative count.
  n_pred_lt_1pct   truth-expressed entries predicted below 1% of their truth
                   value.  On SIRVs nothing is ever exactly 0 -- EM leaves a
                   residue on every reference transcript -- so the literal
                   count is uninformative there and this is the graded form of
                   the same failure.
  n_pred_lt_10pct  same at 10%.

Writes <BASE>/results/grid_metrics.tsv.
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep_config as C  # noqa: E402

METRIC_COLS = [
    "spearman", "pearson_log1p", "ccc_log1p", "kendall",
    "rmse", "nrmse_mean_truth", "mard", "n_scored",
    "n_truth_expressed", "n_truth_unexpressed",
]


def dropouts(path, tool):
    """(n_zero, n_lt_1pct, n_lt_10pct) over truth-expressed rows of a compare table."""
    col = f"{tool}_tpm"
    nz = n1 = n10 = 0
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ref = float(row["ref_tpm"])
            if ref <= 0:
                continue
            pred = float(row[col])
            if pred == 0.0:
                nz += 1
            if pred < 0.01 * ref:
                n1 += 1
            if pred < 0.10 * ref:
                n10 += 1
    return nz, n1, n10


def depth(sample):
    """Primary mapped reads, from the count file LRAA drops beside its outputs."""
    d = os.path.join(C.WORK, sample)
    if not os.path.isdir(d):
        return ""
    for arm in sorted(os.listdir(d)):
        ad = os.path.join(d, arm)
        if not os.path.isdir(ad):
            continue
        for f in os.listdir(ad):
            if f.endswith(".mapped_primary.count"):
                with open(os.path.join(ad, f)) as fh:
                    return fh.read().strip()
    return ""


def cpu_budget(sample, arm):
    """The --cpu_budget the arm actually ran at, read back off its own command line.

    Recorded per arm rather than assumed: budget can in principle change how
    work units are partitioned across workers, and partitioning can change
    float accumulation order.  Two arms compared to each other must have run at
    the same budget, and this column is what makes that checkable after the
    fact instead of remembered.
    """
    log = os.path.join(C.WORK, sample, arm, "run.log")
    if not os.path.exists(log):
        return ""
    with open(log) as fh:
        parts = fh.readline().split()
    return parts[parts.index("--cpu_budget") + 1] if "--cpu_budget" in parts else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(C.BASE, "results", "grid_metrics.tsv"))
    ap.add_argument("--score_dir", default=C.SCORE)
    ap.add_argument("--include_depth", action="store_true",
                    help="also collect the subsampled depth-probe samples")
    args = ap.parse_args()

    rows = []
    for smp in C.SAMPLES + (C.DEPTH_SAMPLES if args.include_depth else []):
        sd = os.path.join(args.score_dir, smp["sample"])
        met = os.path.join(sd, "oarfish_style_metrics.tsv")
        if not os.path.exists(met):
            continue
        nreads = depth(smp["sample"])
        with open(met) as fh:
            for m in csv.DictReader(fh, delimiter="\t"):
                tool = m["name"]
                arm = tool[len("alpha_"):] if tool.startswith("alpha_") else tool
                spec = C.ARMS_BY_NAME.get(arm)
                if spec is None:
                    continue
                cmp_path = os.path.join(sd, f"{tool}.ref_quant_compare.tsv")
                nz, n1, n10 = dropouts(cmp_path, tool) if os.path.exists(cmp_path) else ("", "", "")
                r = dict(
                    sample=smp["sample"], corpus=smp["corpus"], platform=smp["platform"],
                    n_mapped_primary=nreads, cpu_budget=cpu_budget(smp["sample"], arm),
                    arm=arm, alpha=spec["alpha"], threep=int(spec["threep"]),
                    n_pred_zero=nz, n_pred_lt_1pct=n1, n_pred_lt_10pct=n10,
                )
                r.update({k: m.get(k, "") for k in METRIC_COLS})
                rows.append(r)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    cols = (
        ["sample", "corpus", "platform", "n_mapped_primary", "cpu_budget", "arm", "alpha", "threep"]
        + METRIC_COLS + ["n_pred_zero", "n_pred_lt_1pct", "n_pred_lt_10pct"]
    )
    rows.sort(key=lambda r: (r["corpus"], r["sample"], -r["threep"], r["alpha"]))
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    nsamp = len({r["sample"] for r in rows})
    print(f"wrote {args.out}: {len(rows)} arm-rows across {nsamp} samples")


if __name__ == "__main__":
    main()
