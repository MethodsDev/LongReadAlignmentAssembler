#!/usr/bin/env python3

import sys, os, re
import glob
from collections import defaultdict


def main(intron_files_list):

    samplename_to_intron_to_counts = dict()

    samplenames = list()

    for file in intron_files_list:

        sample_name = file
        sample_name = os.path.basename(sample_name).split(".")[0]

        print("-processing: {} -> sample {}".format(file, sample_name), file=sys.stderr)

        samplenames.append(sample_name)

        intron_counts = parse_intron_counts(file)
        samplename_to_intron_to_counts[sample_name] = intron_counts

    print("\t" + "\t".join(samplenames))

    all_introns = set()
    for sample_name, intron_counts in samplename_to_intron_to_counts.items():
        intron_toks = intron_counts.keys()
        all_introns.update(intron_toks)

    for intron in all_introns:
        vals = list()
        vals.append(intron)
        for sample_name in samplenames:
            count = str(samplename_to_intron_to_counts[sample_name][intron])
            vals.append(count)

        print("\t".join(vals))

    sys.exit(0)


def parse_intron_counts(introns_filename):

    intron_counts = defaultdict(int)

    with open(introns_filename) as fh:

        for line in fh:
            """
            0       chr1
            1       14830
            2       14969
            3       -
            4       CT--AC
            5       52
            6       OK
            """

            line = line.rstrip()

            vals = line.split("\t")
            count = vals.pop(5)
            tok = "^".join(vals)
            intron_counts[tok] = count

    return intron_counts


if __name__ == "__main__":

    usage = "usage: build_intron_matrix.py A.introns B.introns ... > introns.matrix"
    if len(sys.argv) < 3:
        exit(usage)

    main(sys.argv[1:])
