#!/usr/bin/env python3
"""Split the alpha response into the transcripts alpha can reach and the rest.

alpha adds `base_alpha * ambiguous_read_count` to a transcript. A transcript
whose ambiguous read count is zero gets a pseudocount of exactly zero at every
alpha, so alpha has no lever on it at all. Those transcripts are the large
majority of the annotation and most of the library's mass, and averaging them
into a library-wide MARD dilutes the measured effect by roughly the reciprocal
of the reachable fraction.

Two questions, answered separately:
  (a) On the reachable partition, how big is the alpha effect really?
  (b) On the unreachable partition, is it indistinguishable from zero -- i.e.
      is moving alpha harmless to the cases that were already fine? This is not
      guaranteed a priori: unreachable transcripts share a normalization with
      reachable ones, so mass moved among the reachable can shift everything
      else by renormalization even with no pseudocount of its own.

Normalization deliberately matches the benchmark: every prediction vector and
the truth are renormalized to 1e6 over the FULL transcript set, then MARD is
averaged within each partition. Library-wide MARD is then exactly the
count-weighted mean of the two partition MARDs, so the decomposition adds up
and neither partition is scored on a private scale.

Note this MARD is computed per TRANSCRIPT, while the benchmark scores per
distinct intron chain (69 -> 60 on SIRVs, since transcripts sharing an intron
chain are indistinguishable and get collapsed). The two are close but not
identical; the full-set column is printed so the difference is visible rather
than assumed away.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep_config as C  # noqa: E402
from ambiguity_profile import PROBE, load_truth, profile_from_tracking  # noqa: E402
from residual_test import expr_path, read_expr  # noqa: E402

EPS = 1e-9


def mard(pred, true, mask=None):
    p = pred / pred.sum() * 1e6
    t = true / true.sum() * 1e6
    v = np.abs(p - t) / (p + t + EPS)
    return float(v.mean() if mask is None else v[mask].mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples")
    ap.add_argument("--threep", type=int, default=1)
    ap.add_argument("--out", default=os.path.join(C.BASE, "results", "partition_analysis.tsv"))
    args = ap.parse_args()

    samples = (
        [C.SAMPLES_BY_NAME[s] for s in args.samples.split(",")]
        if args.samples else list(C.SAMPLES)
    )

    rows = []
    for smp in samples:
        s = smp["sample"]
        trk = os.path.join(PROBE, s, "prof.LRAA.quant-only.quant.tracking.gz")
        if not os.path.exists(trk):
            continue
        amb, _nr, _na = profile_from_tracking(trk)
        truth = load_truth(smp["truth_quant"])

        for a in C.ALPHAS + [C.LIMIT_ALPHA]:
            arm = C.arm_name(a, bool(args.threep))
            path = expr_path(s, arm)
            if not os.path.exists(path):
                continue
            _uniq, allr = read_expr(path)
            keys = sorted(set(allr) | set(truth))
            pred = np.array([allr.get(k, 0.0) for k in keys])
            true = np.array([truth.get(k, 0.0) for k in keys])
            reach = np.array([amb.get(k, 0.0) > 0 for k in keys])
            # Scored only where truth is expressed, matching the benchmark,
            # which never scores a transcript the truth set does not contain.
            scored = true > 0
            r_mask, u_mask = reach & scored, (~reach) & scored
            rows.append(dict(
                sample=s, corpus=smp["corpus"], platform=smp["platform"], alpha=a,
                threep=args.threep,
                n_scored=int(scored.sum()),
                n_reachable=int(r_mask.sum()), n_unreachable=int(u_mask.sum()),
                reachable_pred_mass_frac=float(pred[r_mask].sum() / pred[scored].sum()),
                mard_all=mard(pred, true, scored),
                mard_reachable=mard(pred, true, r_mask) if r_mask.any() else np.nan,
                mard_unreachable=mard(pred, true, u_mask) if u_mask.any() else np.nan,
            ))

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)

    # Effect size on each partition: range across the measured alpha grid, and
    # the paired change from the shipped default to each end of it.
    summ = []
    for s, g in df.groupby("sample"):
        g = g.set_index("alpha").sort_index()
        ref = g.loc[C.DEFAULT_ALPHA]
        row = dict(sample=s, corpus=g.corpus.iloc[0], platform=g.platform.iloc[0],
                   n_reachable=int(ref.n_reachable), n_scored=int(ref.n_scored),
                   frac_tx_reachable=ref.n_reachable / ref.n_scored,
                   reach_mass_frac=ref.reachable_pred_mass_frac)
        for part in ("all", "reachable", "unreachable"):
            c = f"mard_{part}"
            row[f"range_{part}"] = float(g[c].max() - g[c].min())
            row[f"best_alpha_{part}"] = float(g[c].idxmin())
        row["dilution"] = row["range_reachable"] / row["range_all"] if row["range_all"] else np.nan
        summ.append(row)
    su = pd.DataFrame(summ)
    su.to_csv(args.out.replace(".tsv", "_summary.tsv"), sep="\t", index=False)

    pd.set_option("display.width", 240, "display.max_columns", 40)
    print("=== per-sample partition summary (3p %s) ===" % ("on" if args.threep else "off"))
    print(su.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    print("\n=== MARD by partition across alpha ===")
    for s, g in df.groupby("sample"):
        print(f"\n-- {s}")
        print(g.set_index("alpha")[["mard_all", "mard_reachable", "mard_unreachable"]]
              .to_string(float_format=lambda x: f"{x:.6f}"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
