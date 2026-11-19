#!/usr/bin/env python3
"""Genomic span of each gene in a GTF, under the gene id the single-cell matrices use.

LRAA names a gene feature gene_name^gene_id when the GTF gives a gene_name, else
gene_id alone (Transcript.get_output_gene_id). This writes that id with the gene's
span, first exon start to last exon end over all its transcripts:

  gene_id  contig  strand  lend  rend  span

Used by the single-nucleus clustering step (gene_sparseM_to_seurat_clusters_and_umap.R
--gene_spans), whose gene_full counts include intronic reads and so grow with gene
length.
"""

import argparse
import gzip
import os
import re
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PYLIB_DIR = os.path.join(os.path.dirname(os.path.dirname(SCRIPT_DIR)), "pylib")
if PYLIB_DIR not in sys.path:
    sys.path.insert(0, PYLIB_DIR)

from gene_symbol_utils import combine_gene_name_and_id  # noqa: E402

_ATTR = {
    key: re.compile(r'(?:^|;)\s*{} "([^"]*)"'.format(key))
    for key in ("gene_id", "gene_name")
}


def gene_spans(gtf_path):
    spans = {}
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = _ATTR["gene_id"].search(f[8])
            if m is None:
                continue
            m_name = _ATTR["gene_name"].search(f[8])
            gid = combine_gene_name_and_id(m_name.group(1) if m_name else None, m.group(1))
            lend, rend = int(f[3]), int(f[4])
            prev = spans.get(gid)
            if prev is None:
                spans[gid] = [f[0], f[6], lend, rend]
            else:
                prev[2] = min(prev[2], lend)
                prev[3] = max(prev[3], rend)
    return spans


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--gtf", required=True, help="GTF (.gz ok)")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    spans = gene_spans(args.gtf)
    with open(args.output, "wt") as ofh:
        ofh.write("gene_id\tcontig\tstrand\tlend\trend\tspan\n")
        for gid in sorted(spans):
            contig, strand, lend, rend = spans[gid]
            ofh.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(gid, contig, strand, lend, rend, rend - lend + 1))
    print("{} genes".format(len(spans)), file=sys.stderr)


if __name__ == "__main__":
    main()
