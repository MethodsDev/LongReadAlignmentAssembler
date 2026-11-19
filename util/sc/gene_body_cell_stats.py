#!/usr/bin/env python3
"""Per-barcode spliced / unspliced / ambiguous totals from the gene-body matrices.

Reads the four matrices gene_body_tracking_to_sparse_matrix.py writes
(<prefix>^gene_full-sparseM and ^gene_{spliced,unspliced,ambiguous}-sparseM) and writes
one row per barcode:

  cell_barcode     the barcode
  gene_full        reads credited to any gene, intronic included (S + U + A)
  spliced          S: the read has a junction and covers no intron of its gene
  unspliced        U: the read covers >= gene_body_min_feature_bp of its gene's introns
  ambiguous        A: unspliced and wholly exonic (mature or pre-mRNA alike)
  n_genes_full     genes with a nonzero gene_full count in the barcode
  spliced_frac     spliced / gene_full
  unspliced_frac   unspliced / gene_full, the intronic fraction

Counts are read counts (fractional where a read was split between genes). A nucleus
typically has a high unspliced_frac; a whole cell or an ambient-RNA droplet a low one,
so the column helps tell nuclei from debris. Rows are in the matrices' barcode order.
"""

import argparse
import gzip
import os
import sys

import numpy as np
from scipy.io import mmread

LEVELS = ["gene_full", "gene_spliced", "gene_unspliced", "gene_ambiguous"]


def _read_lines(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return [line.rstrip("\n") for line in fh]


def _find(dirname, base):
    for name in (base + ".gz", base):
        path = os.path.join(dirname, name)
        if os.path.exists(path):
            return path
    sys.exit("{} not found in {}".format(base, dirname))


def read_matrix(dirname):
    """(csc matrix features x barcodes, barcodes) of one 10x-style directory."""
    m = mmread(_find(dirname, "matrix.mtx")).tocsc()
    barcodes = _read_lines(_find(dirname, "barcodes.tsv"))
    if m.shape[1] != len(barcodes):
        sys.exit("{}: {} columns but {} barcodes".format(dirname, m.shape[1], len(barcodes)))
    return m, barcodes


def cell_stats(dirs):
    """dirs: {level: matrix dir}. Returns (barcodes, {column: np.array})."""
    full, barcodes = read_matrix(dirs["gene_full"])
    col = {b: k for k, b in enumerate(barcodes)}
    out = {
        "gene_full": np.asarray(full.sum(axis=0)).ravel(),
        "n_genes_full": np.diff(full.indptr),
    }
    for level, label in (
        ("gene_spliced", "spliced"),
        ("gene_unspliced", "unspliced"),
        ("gene_ambiguous", "ambiguous"),
    ):
        m, bc = read_matrix(dirs[level])
        sums = np.asarray(m.sum(axis=0)).ravel()
        v = np.zeros(len(barcodes))
        # The builder writes one barcode list for all four, but a matrix from
        # elsewhere need not, so align by name.
        if bc == barcodes:
            v[:] = sums
        else:
            for b, s in zip(bc, sums):
                if b not in col:
                    sys.exit("barcode {} is in {} but not in gene_full".format(b, level))
                v[col[b]] = s
        out[label] = v
    return barcodes, out


def write_table(barcodes, stats, path):
    full = stats["gene_full"]
    with np.errstate(divide="ignore", invalid="ignore"):
        s_frac = np.where(full > 0, stats["spliced"] / full, 0.0)
        u_frac = np.where(full > 0, stats["unspliced"] / full, 0.0)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "wt") as ofh:
        ofh.write(
            "cell_barcode\tgene_full\tspliced\tunspliced\tambiguous\t"
            "n_genes_full\tspliced_frac\tunspliced_frac\n"
        )
        for k, b in enumerate(barcodes):
            ofh.write(
                "{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{}\t{:.4f}\t{:.4f}\n".format(
                    b,
                    full[k],
                    stats["spliced"][k],
                    stats["unspliced"][k],
                    stats["ambiguous"][k],
                    int(stats["n_genes_full"][k]),
                    s_frac[k],
                    u_frac[k],
                )
            )


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument(
        "--sparseM_prefix",
        required=True,
        help="prefix P of P^gene_full-sparseM etc. (the builder's --output_prefix)",
    )
    ap.add_argument("--output", required=True, help="output tsv (.gz to compress)")
    args = ap.parse_args()

    dirs = {level: "{}^{}-sparseM".format(args.sparseM_prefix, level) for level in LEVELS}
    barcodes, stats = cell_stats(dirs)
    write_table(barcodes, stats, args.output)
    full = stats["gene_full"]
    tot = full.sum()
    print(
        "{} barcodes, gene_full {:.1f}, spliced {:.1%}, unspliced {:.1%}, ambiguous {:.1%}".format(
            len(barcodes),
            tot,
            stats["spliced"].sum() / tot if tot else 0,
            stats["unspliced"].sum() / tot if tot else 0,
            stats["ambiguous"].sum() / tot if tot else 0,
        ),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
