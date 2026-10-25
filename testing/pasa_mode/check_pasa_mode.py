#!/usr/bin/env python3
"""
Check an LRAA --pasa_mode GTF against fixture/expected.tsv.

Each case's multi-exon models, on its contig and strand, must have exactly the
expected intron chains (as a set): staggered reads sharing introns are chained to
the full chain, and pieces that share no intron are not merged.

usage: check_pasa_mode.py expected.tsv lraa.gtf
"""

import re
import sys
from collections import defaultdict


def main(expected_fn, gtf_fn):
    expected = defaultdict(set)
    for line in open(expected_fn):
        if line.startswith("#"):
            continue
        case, strand, chain = line.rstrip("\n").split("\t")
        expected[(case, strand)].add(chain)

    exons = defaultdict(list)
    where = {}
    for line in open(gtf_fn):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        exons[tid].append((int(f[3]), int(f[4])))
        where[tid] = (f[0], f[6])

    observed = defaultdict(set)
    for tid, ex in exons.items():
        ex.sort()
        if len(ex) > 1:
            chain = ",".join(f"{ex[i][1] + 1}-{ex[i + 1][0] - 1}" for i in range(len(ex) - 1))
            observed[where[tid]].add(chain)

    failures = 0
    for key in sorted(expected):
        ok = observed.get(key, set()) == expected[key]
        failures += not ok
        print(f"{'ok  ' if ok else 'FAIL'} {key[0]}{key[1]}: expected {sorted(expected[key])}, "
              f"got {sorted(observed.get(key, set()))}")
    if failures:
        sys.exit(f"{failures} case(s) failed")
    print("all cases pass")


if __name__ == "__main__":
    main(*sys.argv[1:])
