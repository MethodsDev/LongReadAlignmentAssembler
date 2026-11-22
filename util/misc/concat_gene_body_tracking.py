#!/usr/bin/env python3
"""Concatenate per-shard gene-body sidecars (<prefix>.quant.gene_body.tracking.gz).

For a run scattered by chromosome (WDL/LRAA.wdl): each shard is its own LRAA run and
writes its own sidecar. A read belongs to one contig, so the shards' rows are disjoint
and the merge is a concatenation, as for quant.tracking. The first file's comment lines
(LRAA version and command) and header are kept; every other file must carry the same
header, and its comment lines are dropped, except the strand-check lines
(GeneBodyAssigner.STRAND_CHECK_PREFIX), which every shard's are kept for the run total.
"""

import argparse
import gzip
import sys

STRAND_CHECK_PREFIX = "#strand_check\t"  # GeneBodyAssigner.STRAND_CHECK_PREFIX + tab


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def concat(inputs, output):
    header = None
    n_rows = 0
    with gzip.open(output, "wt", compresslevel=1) as ofh:
        for k, path in enumerate(inputs):
            seen_header = False
            with _open(path) as fh:
                for line in fh:
                    if line.startswith("#"):
                        if k == 0 or line.startswith(STRAND_CHECK_PREFIX):
                            ofh.write(line)
                        continue
                    if not seen_header:
                        seen_header = True
                        if header is None:
                            header = line
                            ofh.write(line)
                        elif line != header:
                            sys.exit("gene-body sidecar header differs: {}".format(path))
                        continue
                    ofh.write(line)
                    n_rows += 1
            if not seen_header:
                sys.exit("gene-body sidecar has no header: {}".format(path))
    return n_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--output", required=True, help="output .gz")
    args = ap.parse_args()
    n = concat(args.inputs, args.output)
    print("{} rows from {} sidecars".format(n, len(args.inputs)), file=sys.stderr)


if __name__ == "__main__":
    main()
