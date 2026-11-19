#!/usr/bin/env python3
"""Cell x gene matrices from LRAA's gene-body sidecar (<prefix>.quant.gene_body.tracking.gz).

Written by LRAA (quant-only, or a discovery run's final quant) when
config["gene_body_counts"] is on
(docs/single_nucleus_intronic_gene_counts.md). Each row is one (read, gene) with
frac_assigned and a splice class S / U / A.

Outputs, in the same 10x layout as singlecell_tracking_to_sparse_matrix.py
(features x barcodes, matrix.mtx.gz / features.tsv.gz / barcodes.tsv.gz):

  <prefix>^gene_full-sparseM        every class: the snRNA gene count (intronic included)
  <prefix>^gene_spliced-sparseM     S
  <prefix>^gene_unspliced-sparseM   U
  <prefix>^gene_ambiguous-sparseM   A

All four share one feature list and one barcode list, so S + U + A == gene_full
entry for entry (velocity tools take S and U layers; scVelo/alevin-fry convention
puts A with spliced or splits it, the caller's choice).

--exclude_gene_body keeps only source == "isoform" rows: the existing isoform-assigned
reads, split by class, for comparison with the standard gene matrix.

A prototype: it holds one dict entry per (gene, cell, class), which is fine for a
chr-scale test and bounded by nonzeros, not reads; the genome-scale path would reuse
singlecell_tracking_to_sparse_matrix's chunked LevelAccumulator.
"""

import argparse
import gzip
import logging
import os
import sys
from collections import defaultdict

import numpy as np
from scipy import sparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from singlecell_tracking_to_sparse_matrix import make_sparse_matrix_outputs  # noqa: E402

logger = logging.getLogger(__name__)

CLASSES = {"S": "gene_spliced", "U": "gene_unspliced", "A": "gene_ambiguous"}


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def accumulate(path, exclude_gene_body=False):
    counts = defaultdict(float)  # (gene, barcode, cls) -> sum frac
    with _open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = {name: i for i, name in enumerate(parts)}
                i_g, i_r, i_f = header["gene_id"], header["read_name"], header["frac_assigned"]
                i_c, i_s = header["splice_class"], header["source"]
                continue
            frac = float(parts[i_f])
            if frac == 0.0:
                continue  # gene-ambiguous read under policy "drop"
            if exclude_gene_body and parts[i_s] != "isoform":
                continue
            barcode, tagged, _ = parts[i_r].partition("^")
            if not tagged:
                continue  # no cell tags: no barcode column (as the tracking builder)
            counts[(parts[i_g], barcode, parts[i_c])] += frac
    return counts


def build(counts, output_prefix, gzip_level=1):
    genes = sorted({g for g, _, _ in counts})
    barcodes = sorted({b for _, b, _ in counts})
    gi = {g: k for k, g in enumerate(genes)}
    bi = {b: k for k, b in enumerate(barcodes)}
    layers = {"gene_full": defaultdict(float)}
    for label in CLASSES.values():
        layers[label] = defaultdict(float)
    for (g, b, cls), v in counts.items():
        layers["gene_full"][(gi[g], bi[b])] += v
        layers[CLASSES[cls]][(gi[g], bi[b])] += v
    out = {}
    for label, cells in layers.items():
        if cells:
            rows, cols = zip(*cells.keys())
            vals = np.fromiter(cells.values(), dtype=np.float64, count=len(cells))
        else:
            rows, cols, vals = (), (), np.array([], dtype=np.float64)
        m = sparse.coo_matrix(
            (vals, (np.array(rows, dtype=np.int64), np.array(cols, dtype=np.int64))),
            shape=(len(genes), len(barcodes)),
        ).tocsr()
        make_sparse_matrix_outputs(
            m, genes, barcodes, f"{output_prefix}^{label}-sparseM", gzip_level=gzip_level
        )
        out[label] = m
    return out, genes, barcodes


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--gene_body_tracking", required=True)
    ap.add_argument("--output_prefix", required=True)
    ap.add_argument("--exclude_gene_body", action="store_true")
    ap.add_argument("--gzip_level", type=int, default=1)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    counts = accumulate(args.gene_body_tracking, args.exclude_gene_body)
    mats, genes, barcodes = build(counts, args.output_prefix, args.gzip_level)
    for label, m in mats.items():
        logger.info("%s: %d genes x %d barcodes, total %.1f", label, len(genes), len(barcodes), m.sum())


if __name__ == "__main__":
    main()
