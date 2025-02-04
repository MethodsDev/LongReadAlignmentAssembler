#!/usr/bin/env python3

import sys, os, re


def main():

    usage = "\n\tusage: {} sample_id annot.gtf\n\n".format(sys.argv[0])

    if len(sys.argv) < 3:
        exit(usage)

    sample_id = sys.argv[1]
    gtf_filename = sys.argv[2]

    with open(gtf_filename, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            line = line.replace('gene_id "', f'gene_id "{sample_id}:')
            line = line.replace('transcript_id "', f'transcript_id "{sample_id}:')
            print(line)

    sys.exit(0)


if __name__ == "__main__":
    main()
