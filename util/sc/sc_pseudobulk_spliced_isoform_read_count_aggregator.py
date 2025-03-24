#!/usr/bin/env python3

import sys, os, re
import csv
from collections import defaultdict
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    usage = "\n\n\tusage: {} *.expr\n\n".format(sys.argv[0])

    if len(sys.argv) < 2:
        exit(usage)

    files_to_aggregate = sys.argv[1:]

    introns_to_read_count = defaultdict(int)

    for filename in files_to_aggregate:
        logger.info("-processing {}".format(filename))
        with open(filename, "rt") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                introns = row["introns"]
                if introns != "":
                    gene_id = row["gene_id"]
                    chrom = gene_id.split(":")[1]
                    introns = chrom + ":" + introns
                    introns_to_read_count[introns] += float(row["all_reads"])

    print("\t".join(["splice_pattern", "num_reads"]))
    for intron, count in introns_to_read_count.items():
        print("\t".join([intron, str(count)]))

    sys.exit(0)


if __name__ == "__main__":
    main()
