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

Memory follows the nonzeros, not the reads: rows are accumulated as integer
(gene, barcode) triples per class and folded into sparse matrices every
--flush_rows rows, where duplicates sum.
"""

import argparse
import gzip
import logging
from array import array
import os
import sys

import numpy as np
from scipy import sparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from singlecell_tracking_to_sparse_matrix import make_sparse_matrix_outputs  # noqa: E402

logger = logging.getLogger(__name__)

CLASSES = {"S": "gene_spliced", "U": "gene_unspliced", "A": "gene_ambiguous"}


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


class Accumulator:
    """(gene, barcode, class) -> summed frac, as one growing CSR matrix per class."""

    def __init__(self, flush_rows=5_000_000):
        self.flush_rows = flush_rows
        self.gene_ix = {}
        self.bc_ix = {}
        self.mats = {cls: None for cls in CLASSES}
        self._reset()

    def _reset(self):
        self.rows = {cls: array("q") for cls in CLASSES}
        self.cols = {cls: array("q") for cls in CLASSES}
        self.vals = {cls: array("d") for cls in CLASSES}
        self.pending = 0

    def add(self, gene, barcode, cls, frac):
        gi = self.gene_ix.setdefault(gene, len(self.gene_ix))
        bi = self.bc_ix.setdefault(barcode, len(self.bc_ix))
        self.rows[cls].append(gi)
        self.cols[cls].append(bi)
        self.vals[cls].append(frac)
        self.pending += 1
        if self.pending >= self.flush_rows:
            self.flush()

    def flush(self):
        shape = (len(self.gene_ix), len(self.bc_ix))
        for cls in CLASSES:
            part = sparse.coo_matrix(
                (
                    np.frombuffer(self.vals[cls], dtype=np.float64),
                    (
                        np.frombuffer(self.rows[cls], dtype=np.int64),
                        np.frombuffer(self.cols[cls], dtype=np.int64),
                    ),
                ),
                shape=shape,
            ).tocsr()
            prev = self.mats[cls]
            if prev is None:
                self.mats[cls] = part
            else:
                prev.resize(shape)
                self.mats[cls] = prev + part
        self._reset()

    def finish(self):
        """({class: csr}, genes, barcodes), rows and columns in sorted label order."""
        self.flush()
        genes = sorted(self.gene_ix)
        barcodes = sorted(self.bc_ix)
        g_order = np.array([self.gene_ix[g] for g in genes], dtype=np.int64)
        b_order = np.array([self.bc_ix[b] for b in barcodes], dtype=np.int64)
        mats = {}
        for cls, m in self.mats.items():
            m.resize((len(self.gene_ix), len(self.bc_ix)))
            mats[cls] = m[g_order][:, b_order].tocsr()
        return mats, genes, barcodes


def accumulate(path, exclude_gene_body=False, flush_rows=5_000_000):
    acc = Accumulator(flush_rows)
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
            acc.add(parts[i_g], barcode, parts[i_c], frac)
    return acc.finish()


def build(accumulated, output_prefix, gzip_level=1):
    by_class, genes, barcodes = accumulated
    out = {"gene_full": by_class["S"] + by_class["U"] + by_class["A"]}
    for cls, label in CLASSES.items():
        out[label] = by_class[cls]
    for label, m in out.items():
        m = m.tocsr()
        m.eliminate_zeros()
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
    ap.add_argument(
        "--flush_rows",
        type=int,
        default=5_000_000,
        help="rows held as triples before they are folded into the matrices",
    )
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    accumulated = accumulate(args.gene_body_tracking, args.exclude_gene_body, args.flush_rows)
    mats, genes, barcodes = build(accumulated, args.output_prefix, args.gzip_level)
    for label, m in mats.items():
        logger.info("%s: %d genes x %d barcodes, total %.1f", label, len(genes), len(barcodes), m.sum())


if __name__ == "__main__":
    main()
