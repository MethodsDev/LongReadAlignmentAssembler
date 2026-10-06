#!/usr/bin/env python3

# Fixture for test_parallel_parity: the SIRV bam plus unplaced unmapped reads, so
# the --CPU path's fetch("*") part is exercised, and a reference without SIRV7, so
# one contig has reads but no annotation.

import sys
import pysam

src_bam, out_bam, src_gtf, out_gtf = sys.argv[1:5]

src = pysam.AlignmentFile(src_bam, "rb")
out = pysam.AlignmentFile(out_bam, "wb", template=src)
reads = list(src)
for read in reads:
    out.write(read)
for i, read in enumerate(reads[:20]):
    unplaced = pysam.AlignedSegment(out.header)
    unplaced.query_name = "unplaced_{}".format(i)
    unplaced.query_sequence = read.query_sequence
    unplaced.flag = 4
    unplaced.reference_id = -1
    unplaced.reference_start = -1
    out.write(unplaced)
out.close()
pysam.index(out_bam)

with open(src_gtf) as fh, open(out_gtf, "w") as ofh:
    for line in fh:
        if not line.startswith("SIRV7\t"):
            ofh.write(line)
