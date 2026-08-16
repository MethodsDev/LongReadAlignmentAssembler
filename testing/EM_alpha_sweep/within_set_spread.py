#!/usr/bin/env python3
"""How similar are the TRUE abundances of the candidates alpha has to choose between?

alpha flattens the abundance distribution WITHIN an ambiguous candidate set: as
it grows, theta over the set approaches the set's ambiguity profile rather than
the likelihood's split. So its sign should be governed by how similar the true
abundances of the candidates in a set actually are -- flattening helps when they
are genuinely comparable and hurts when they are skewed. That predicts SIRV E0
(identical within sets, alpha=1.0 best), E1 (8x range, high alpha), E2 (128x,
alpha=0) and arabidopsis (skewed, alpha=0), all of which are observed. MORF is
the only corpus whose within-set structure has never been measured, and it is
the one that does not fit.

Definition is shared verbatim with the arabidopsis measurement so the two are
comparable; changing it here without changing it there invalidates the whole
comparison.

  source        the 3p-ON arm's quant.tracking.gz at alpha 0.01. Candidate sets
                come from read compatibility, not from theta, so they are
                alpha-independent and the arm only has to be consistent.
  candidate set distinct transcript_ids sharing a read_name, kept when |S| > 1,
                then grouped by frozenset so identical sets collapse and carry a
                read count. Grouped by READ_NAME, not mp_id -- mp ids are unique
                only within a component.
  truth         truth_quant renormalized to TPM over the whole truth set.
                Members absent from truth are decoys at 0.
  variants      EXPRESSED_ONLY drops decoys and keeps sets with >=2 survivors
                (primary: max/min is undefined with a zero, and the hypothesis
                is about flattening among real candidates). ALL_MEMBERS retains
                decoys and reports CV only.
  per set       ratio = max/min of member truth TPM; cv = population std/mean,
                both on the linear TPM scale, so two equal candidates give
                exactly ratio 1.0 and cv 0.0.
  aggregation   READ-weighted (each set contributes its read count -- the mass
                alpha actually acts on) and SET-weighted (each set once, to show
                whether a read-weighted result rides on one huge set).
                Percentiles, not means: ratio is heavy-tailed and one 10,000x
                set would dominate a mean.
"""

import argparse
import os
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

# PROVENANCE: this file is QuantAlpha's testing/EM_alpha_sweep/within_set_spread.py
# at commit b8b8bdf on branch quant-alpha, VERBATIM except for the argument
# plumbing that lets it read an explicit --tracking/--truth_quant pair instead of
# resolving paths through their sweep_config.  The statistic -- candidate_sets(),
# wpct(), and every line that computes a ratio, a CV or a percentile -- is
# untouched, because the arabidopsis and morf2 numbers are only comparable if one
# implementation produced both.  load_truth is imported from their tree rather
# than copied, for the same reason.  Diff sent to QuantAlpha for adoption.
QA = "/home/unix/bhaas/projects/SingleCellOverhaul/LRAA-quant-alpha/testing/EM_alpha_sweep"
sys.path.insert(0, QA)
from ambiguity_profile import load_truth  # noqa: E402
try:                                       # only needed for --samples mode
    import sweep_config as C               # noqa: E402
    from ambiguity_profile import PROBE    # noqa: E402
except Exception:                          # pragma: no cover - path-explicit mode
    C = PROBE = None


def candidate_sets(path):
    """frozenset(transcript_ids) -> number of distinct reads with that set, |S|>1."""
    import csv
    import gzip

    read_to_tx = defaultdict(set)
    with gzip.open(path, "rt") as fh:
        rows = (ln for ln in fh if not ln.startswith("#"))
        for r in csv.DictReader(rows, delimiter="\t"):
            read_to_tx[r["read_name"]].add(r["transcript_id"])
    sets = Counter()
    n_reads = len(read_to_tx)
    for txs in read_to_tx.values():
        if len(txs) > 1:
            sets[frozenset(txs)] += 1
    return sets, n_reads


def wpct(values, weights, qs=(25, 50, 75)):
    """Weighted percentiles. Equivalent to np.percentile when weights are all 1."""
    v = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    if v.size == 0:
        return [float("nan")] * len(qs)
    o = np.argsort(v)
    v, w = v[o], w[o]
    cw = np.cumsum(w)
    cut = cw[-1] * np.asarray(qs, dtype=float) / 100.0
    return [float(v[min(np.searchsorted(cw, c), v.size - 1)]) for c in cut]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", default="morf2_ont,morf2_pacbio")
    ap.add_argument("--tracking", default=None,
                    help="explicit ON-arm quant.tracking.gz; bypasses sweep_config")
    ap.add_argument("--truth_quant", default=None, help="explicit truth quant tsv")
    ap.add_argument("--label", default=None, help="sample label for the output row")
    ap.add_argument("--corpus", default="", help="corpus label for the output row")
    ap.add_argument("--platform", default="", help="platform label for the output row")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    if args.tracking and not args.truth_quant:
        ap.error("--tracking requires --truth_quant")
    if args.out is None:
        if C is None:
            ap.error("--out is required in path-explicit mode")
        args.out = os.path.join(C.BASE, "results", "within_set_spread.tsv")

    rows = []
    if args.tracking:
        targets = [(args.label or "explicit", args.corpus, args.platform,
                    args.tracking, args.truth_quant)]
    else:
        targets = []
        for name in args.samples.split(","):
            smp = C.SAMPLES_BY_NAME[name]
            targets.append((
                name, smp["corpus"], smp["platform"],
                os.path.join(PROBE, smp["sample"], "prof.LRAA.quant-only.quant.tracking.gz"),
                smp["truth_quant"]))
    for name, _corpus, _platform, trk, _truth_path in targets:
        if not os.path.exists(trk):
            print(f"  {name}: no tracking, skipped")
            continue
        sets, n_reads_assigned = candidate_sets(trk)
        truth = load_truth(_truth_path)
        tot = sum(truth.values())
        tpm = {k: v / tot * 1e6 for k, v in truth.items()}

        n_amb_reads = sum(sets.values())
        rec = dict(
            sample=name, corpus=_corpus, platform=_platform,
            n_reads_assigned=n_reads_assigned, n_ambiguous_reads=n_amb_reads,
            frac_reads_ambiguous=n_amb_reads / n_reads_assigned,
            n_distinct_sets=len(sets),
            median_set_size_readwt=wpct([len(s) for s in sets], list(sets.values()), (50,))[0],
            max_set_size=max(len(s) for s in sets) if sets else 0,
        )

        # Decoy exposure: sets containing a member the truth set does not list.
        dec = [(s, n) for s, n in sets.items() if any(t not in tpm for t in s)]
        rec["n_sets_with_decoy"] = len(dec)
        rec["frac_amb_reads_seeing_decoy"] = (
            sum(n for _s, n in dec) / n_amb_reads if n_amb_reads else float("nan")
        )

        for variant in ("EXPRESSED_ONLY", "ALL_MEMBERS"):
            ratios, cvs, wts = [], [], []
            for s, n in sets.items():
                vals = ([tpm[t] for t in s if t in tpm] if variant == "EXPRESSED_ONLY"
                        else [tpm.get(t, 0.0) for t in s])
                if len(vals) < 2:
                    continue
                a = np.asarray(vals, dtype=float)
                cvs.append(float(a.std() / a.mean()) if a.mean() > 0 else float("nan"))
                if variant == "EXPRESSED_ONLY":
                    ratios.append(float(a.max() / a.min()) if a.min() > 0 else float("nan"))
                wts.append(n)
            pre = "exp" if variant == "EXPRESSED_ONLY" else "all"
            rec[f"{pre}_n_sets"] = len(wts)
            rec[f"{pre}_n_reads"] = int(sum(wts))
            for label, w in (("readwt", wts), ("setwt", [1] * len(wts))):
                if ratios:
                    q = wpct(ratios, w)
                    rec[f"{pre}_ratio_p25_{label}"], rec[f"{pre}_ratio_med_{label}"], \
                        rec[f"{pre}_ratio_p75_{label}"] = q
                q = wpct(cvs, w)
                rec[f"{pre}_cv_p25_{label}"], rec[f"{pre}_cv_med_{label}"], \
                    rec[f"{pre}_cv_p75_{label}"] = q
            if ratios:
                r = np.asarray(ratios, dtype=float)
                w = np.asarray(wts, dtype=float)
                ok = np.isfinite(r) & (r > 0)
                rec[f"{pre}_ratio_geomean_readwt"] = float(
                    np.exp(np.average(np.log(r[ok]), weights=w[ok]))
                ) if ok.any() else float("nan")
                rec[f"{pre}_frac_reads_ratio_under_2"] = float(w[ok & (r < 2)].sum() / w[ok].sum())
                rec[f"{pre}_frac_reads_ratio_under_10"] = float(w[ok & (r < 10)].sum() / w[ok].sum())
        rows.append(rec)

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    pd.set_option("display.width", 250, "display.max_columns", 80)
    print(df.T.to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
