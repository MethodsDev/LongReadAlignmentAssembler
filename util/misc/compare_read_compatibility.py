#!/usr/bin/env python3
"""Did a config change WHICH transcripts each read is compatible with?

Measurement instrumentation. Two LRAA quant.tracking files, compared per read on the SET
of transcripts the read was assigned to. A config that only reweights an existing
assignment leaves every set identical; a config that also perturbs the splice graph,
polyA sites, or path compatibility changes the sets themselves. Distinguishing those two
is the whole point: they are different mechanisms with different blast radii, and a
metric delta alone cannot tell them apart.

The intended use is isolating max_dist_between_alt_polyA_sites. P is the pseudocount in
the 3'-agreement weight, but it is also a window for polyA-site aggregation, for the
reference-3'-end veto exemption, for terminal-exon trimming, and for snapping a read's
terminus onto a PolyAsite node. Run both arms with the 3' weighting OFF and every weight
is pinned to 1.0, so any difference this script reports cannot be coming through the
weight formula and must be one of the others.

Known-answer validation, on real output: 3'-weighting on vs off, CL_HG002_E2_sirv --
0 of 250,129 compatible sets changed, 277 reads changed fractional assignment. That is
the exact signature of a pure reweighting, and it is pinned in
pylib/test_compare_read_compatibility.py.
"""

import argparse
import gzip
import sys
from collections import defaultdict

FRAC_EPS = 1e-6


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def load(path):
    """(read_name -> frozenset of transcript ids, read_name -> summed frac_assigned)."""
    sets = defaultdict(set)
    frac = defaultdict(float)
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("gene_id\t"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            sets[f[5]].add(f[1])
            frac[f[5]] += float(f[6])
    return {k: frozenset(v) for k, v in sets.items()}, frac


def compare(path_a, path_b):
    """Per-read differences between two tracking files.

    changed_compatible_set counts reads whose candidate transcript set differs, which
    only a change to read/transcript compatibility can cause. changed_total_frac counts
    reads whose assignment merely moved, which reweighting alone explains. Reporting
    both separately is what makes the two mechanisms distinguishable.
    """
    A, fa = load(path_a)
    B, fb = load(path_b)
    if not A or not B:
        raise ValueError(f"empty tracking input: {path_a} / {path_b}")

    shared = set(A) & set(B)
    return {
        "n_reads_A": len(A),
        "n_reads_B": len(B),
        "reads_only_in_A": len(set(A) - set(B)),
        "reads_only_in_B": len(set(B) - set(A)),
        "shared_reads": len(shared),
        "changed_compatible_set": sum(1 for r in shared if A[r] != B[r]),
        "changed_set_size": sum(1 for r in shared if len(A[r]) != len(B[r])),
        "changed_total_frac": sum(1 for r in shared if abs(fa[r] - fb[r]) > FRAC_EPS),
    }


FIELDS = ["n_reads_A", "n_reads_B", "reads_only_in_A", "reads_only_in_B",
          "shared_reads", "changed_compatible_set", "changed_set_size",
          "changed_total_frac"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--label", default="")
    ap.add_argument("--header", action="store_true")
    args = ap.parse_args()

    try:
        r = compare(args.a, args.b)
    except ValueError as e:
        sys.exit(f"Error, {e}")

    shared = r["shared_reads"]
    pct = 100.0 * r["changed_compatible_set"] / shared if shared else float("nan")
    if args.header:
        print("label\t" + "\t".join(FIELDS[:6]) + "\tpct_changed_set\t"
              + "\t".join(FIELDS[6:]))
    print(f"{args.label}\t" + "\t".join(str(r[f]) for f in FIELDS[:6])
          + f"\t{pct:.4f}\t" + "\t".join(str(r[f]) for f in FIELDS[6:]))


if __name__ == "__main__":
    main()
