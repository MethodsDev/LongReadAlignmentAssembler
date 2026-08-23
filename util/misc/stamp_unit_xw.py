#!/usr/bin/env python3

"""Copy a bam, stamping XW:f:1.0 on every record, preserving the read set exactly.

For the streaming-parity target, which needs a --bam_for_sg that is a MATCHED control:
the same reads as --bam, so that streaming and non-streaming must agree exactly rather
than approximately. LRAA requires a supplied --bam_for_sg to carry XW, because that flag
names already-normalized splice-graph evidence and an untagged thinned bam would be read
as though its survivors were the whole library.

A weight of 1.0 is the honest value here and changes nothing numerically -- it is what
coverage normalization itself writes for a read it kept whole -- so the control declares
"normalized, nothing thinned" without perturbing the comparison.

Running the real normalizer at an unreachable depth target would also stamp 1.0, but it
does not preserve the read set: its strand split discards records whose orientation is
undeterminable (2 of 553 on this corpus), and a control missing reads builds a different
splice graph, which is the one thing this fixture must not do.
"""

import sys

import pysam


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: stamp_unit_xw.py <in.bam> <out.bam>")
    in_bam, out_bam = sys.argv[1], sys.argv[2]

    n = 0
    with pysam.AlignmentFile(in_bam, "rb") as reader:
        with pysam.AlignmentFile(out_bam, "wb", template=reader) as writer:
            for read in reader.fetch(until_eof=True):
                read.set_tag("XW", 1.0, value_type="f")
                writer.write(read)
                n += 1
    pysam.index(out_bam)
    print(f"{out_bam}: stamped XW:f:1.0 on {n} records", file=sys.stderr)


if __name__ == "__main__":
    main()
