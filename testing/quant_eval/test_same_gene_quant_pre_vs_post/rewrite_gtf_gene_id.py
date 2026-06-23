#!/usr/bin/env python3

import re
import sys


def main():
    if len(sys.argv) != 4:
        raise SystemExit(
            "usage: rewrite_gtf_gene_id.py input.gtf output.gtf replacement_gene_id"
        )

    input_gtf, output_gtf, replacement_gene_id = sys.argv[1:]
    replacement = f'gene_id "{replacement_gene_id}"'

    with open(input_gtf, "rt") as src, open(output_gtf, "wt") as dest:
        for line in src:
            dest.write(re.sub(r'gene_id "[^"]+"', replacement, line))


if __name__ == "__main__":
    main()
