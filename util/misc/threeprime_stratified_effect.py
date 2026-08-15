#!/usr/bin/env python3
"""Where does the 3'-agreement weighting actually change accuracy?

Measurement instrumentation. A library-wide MARD delta between the weighting on and off
mixes together transcripts the feature can act on and transcripts it cannot touch,
diluting a real effect with a large inert majority. This splits the delta by whether the
feature had a LEVER on the transcript, and the split is an exact decomposition: the
per-stratum deltas recombine to the library-wide delta, and the script checks that they
do rather than asking to be trusted.

Strata, assigned from the ON arm's quant.tracking (the weights EM actually used):

  LEVER        the transcript received at least one read whose weight vector across its
               candidates was NON-UNIFORM. These are the only transcripts where turning
               the feature off can change an E-step fraction.
  TIED_ONLY    it received ambiguous reads, but every one of them had uniform weights --
               its candidates all sit at the same 3' distance. Uniform weights and
               absent weights give identical E-step fractions, so the feature is inert
               here by construction, not by accident.
  UNAMBIGUOUS  every read assigned to it was compatible with it alone. Weight is forced
               to 1.0 by the sum_weights == 0 branch whether the feature is on or off.
  NO_READS     in the truth set but assigned no reads in the ON arm.

MARD is the scorer's symmetric form, mean(|pred-true| / (pred+true+1e-10)), computed on
the columns of <arm>.ref_quant_compare.tsv, which the benchmark has already renormalized
to 1e6 over the scored universe. Normalization is deliberately NOT redone per stratum:
keeping the library-wide scale is what makes the strata sum back to the whole.

Corpus-agnostic -- it needs one tracking file and two ref_quant_compare.tsv files and
knows nothing about SIRV, MORF or any particular annotation.
"""

import argparse
import gzip
import sys
from collections import defaultdict

EPSILON = 1e-10
SPREAD_EPS = 1e-5
STRATA = ["LEVER", "TIED_ONLY", "UNAMBIGUOUS", "NO_READS"]


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def classify_transcripts(tracking):
    """transcript_id -> 'LEVER' | 'TIED_ONLY' | 'UNAMBIGUOUS', from the ON arm."""
    per_read = defaultdict(list)
    with opener(tracking) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("gene_id\t"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            per_read[f[5]].append((f[1], float(f[7])))

    has_lever, has_tied = set(), set()
    seen = set()
    for members in per_read.values():
        tids = [t for t, _ in members]
        seen.update(tids)
        if len(members) < 2:
            continue
        ws = [w for _, w in members]
        if max(ws) > min(ws) * (1.0 + SPREAD_EPS):
            has_lever.update(tids)
        else:
            has_tied.update(tids)

    out = {}
    for t in seen:
        out[t] = "LEVER" if t in has_lever else ("TIED_ONLY" if t in has_tied
                                                 else "UNAMBIGUOUS")
    return out


def load_compare(path, arm):
    """intronId -> (pred_tpm, ref_tpm, [transcript_ids])."""
    rows = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        try:
            i_pred = hdr.index(f"{arm}_tpm")
        except ValueError:
            sys.exit(f"Error, column {arm}_tpm not in {path}; header={hdr}")
        i_ref, i_tx = hdr.index("ref_tpm"), hdr.index("transcript_ids")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_pred, i_ref, i_tx):
                continue
            rows[f[0]] = (float(f[i_pred]), float(f[i_ref]),
                          [t for t in f[i_tx].split(",") if t])
    return rows


def mard_term(pred, ref):
    return abs(pred - ref) / (pred + ref + EPSILON)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tracking", required=True, help="ON-arm quant.tracking(.gz)")
    ap.add_argument("--on", required=True, help="<on_arm>.ref_quant_compare.tsv")
    ap.add_argument("--off", required=True, help="<off_arm>.ref_quant_compare.tsv")
    ap.add_argument("--on_arm", required=True, help="registry name of the ON arm")
    ap.add_argument("--off_arm", required=True, help="registry name of the OFF arm")
    ap.add_argument("--label", default="")
    args = ap.parse_args()

    tclass = classify_transcripts(args.tracking)
    on = load_compare(args.on, args.on_arm)
    off = load_compare(args.off, args.off_arm)

    shared = set(on) & set(off)
    if not shared:
        sys.exit("Error, the two compare files share no intronId rows")

    # A row can carry several transcript ids sharing an intron chain. Take the strongest
    # claim present: if any of them had a lever, the feature could move this row.
    priority = {"LEVER": 0, "TIED_ONLY": 1, "UNAMBIGUOUS": 2, "NO_READS": 3}
    buckets = defaultdict(list)
    multi_tx_rows = 0
    for iid in shared:
        pred_on, ref, tids = on[iid]
        pred_off = off[iid][0]
        if len(tids) > 1:
            multi_tx_rows += 1
        classes = [tclass.get(t, "NO_READS") for t in tids] or ["NO_READS"]
        cls = min(classes, key=lambda c: priority[c])
        buckets[cls].append((mard_term(pred_on, ref), mard_term(pred_off, ref),
                             pred_on, pred_off, ref))

    n_total = len(shared)
    all_on = sum(b[0] for v in buckets.values() for b in v) / n_total
    all_off = sum(b[1] for v in buckets.values() for b in v) / n_total

    print(f"# {args.label}: {n_total} scored rows, {multi_tx_rows} carrying >1 "
          f"transcript_id; strata from {args.tracking}")
    print(f"# LIBRARY-WIDE  mard_on={all_on:.9f}  mard_off={all_off:.9f}  "
          f"delta={all_off - all_on:+.6e}")
    print("label\tstratum\tn_rows\tpct_rows\tmard_on\tmard_off\tdelta"
          "\tcontrib_to_library_delta\tpct_of_library_delta\tsum_ref_tpm")
    total_delta = all_off - all_on
    contribs = {}
    for cls in STRATA:
        v = buckets.get(cls, [])
        if not v:
            continue
        n = len(v)
        m_on = sum(x[0] for x in v) / n
        m_off = sum(x[1] for x in v) / n
        d = m_off - m_on
        # Each stratum's share of the library-wide delta is its delta weighted by the
        # fraction of rows it holds; these are exactly additive.
        contrib = d * n / n_total
        contribs[cls] = contrib
        pct = 100.0 * contrib / total_delta if total_delta else float("nan")
        print(f"{args.label}\t{cls}\t{n}\t{100.0 * n / n_total:.2f}\t{m_on:.9f}"
              f"\t{m_off:.9f}\t{d:+.6e}\t{contrib:+.6e}\t{pct:.2f}"
              f"\t{sum(x[4] for x in v):.1f}")

    # The decomposition is only worth reading if it reconstructs the whole.
    recon = sum(contribs.values())
    ok = abs(recon - total_delta) < 1e-12
    print(f"# decomposition check: sum(contributions)={recon:+.6e} vs "
          f"library delta={total_delta:+.6e} -> {'EXACT' if ok else 'MISMATCH'}")
    if not ok:
        sys.exit("Error, stratum contributions do not reconstruct the library delta")


if __name__ == "__main__":
    main()
