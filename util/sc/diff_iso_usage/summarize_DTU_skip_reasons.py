#!/usr/bin/env python3
"""
Why the chi-square DTU test skipped each gene x cluster pair, from the per-isoform table
sc_pseudobulk_test_isoform_DiffUsage.py writes with --save_annotated_isoform_details
(<prefix>.diff_iso.annotated_isoforms.tsv, several GB). Every isoform row of a gene x
pair carries the same skip_reason ("." when the pair was tested), so this keeps one row
per gene_symbol, cluster_A, cluster_B.

Output (tsv, gzipped if the name ends in .gz): gene_symbol, cluster_A, cluster_B,
skip_reason -- the --isoform_DTU_skip_reasons input of
util/sc/site_usage/compare_site_vs_isoform_DTU.py.
"""

import argparse
import csv
import gzip
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--annotated_isoforms", required=True, help="<prefix>.diff_iso.annotated_isoforms.tsv[.gz]")
    p.add_argument("--output", required=True, help="output tsv (.gz to compress)")
    args = p.parse_args()

    opener = gzip.open if args.annotated_isoforms.endswith(".gz") else open
    seen = {}
    conflicts = 0
    with opener(args.annotated_isoforms, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            key = (row["gene_symbol"], row["cluster_A"], row["cluster_B"])
            reason = row["skip_reason"]
            prev = seen.get(key)
            if prev is None:
                seen[key] = reason
            elif prev != reason:
                conflicts += 1

    if conflicts:
        print(f"warning: {conflicts} isoform rows disagree with their gene x pair's first skip_reason",
              file=sys.stderr)

    out_open = gzip.open if args.output.endswith(".gz") else open
    with out_open(args.output, "wt") as ofh:
        ofh.write("gene_symbol\tcluster_A\tcluster_B\tskip_reason\n")
        for (g, a, b), r in seen.items():
            ofh.write(f"{g}\t{a}\t{b}\t{r}\n")
    print(f"{len(seen)} gene x cluster pairs", file=sys.stderr)


if __name__ == "__main__":
    main()
