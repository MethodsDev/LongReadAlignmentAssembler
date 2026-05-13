#!/usr/bin/env python3

import csv
from collections import defaultdict


def read_tsv(path):
    with open(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def assert_close(actual, expected, tol=1e-6):
    if abs(actual - expected) > tol:
        raise AssertionError(f"expected {expected}, got {actual}")


def main():
    raw_expr = read_tsv("LRAA.pre_cross_gene_em.quant.expr")
    corrected_expr = read_tsv("LRAA.quant.expr")
    raw_tracking = read_tsv("LRAA.pre_cross_gene_em.quant.tracking")
    corrected_tracking = read_tsv("LRAA.quant.tracking")

    raw_expr_by_tx = {row["transcript_id"]: row for row in raw_expr}
    corrected_expr_by_tx = {row["transcript_id"]: row for row in corrected_expr}

    if set(raw_expr_by_tx) != {"txA1", "txA2", "txB1", "txB2"}:
        raise AssertionError(f"unexpected raw transcripts: {sorted(raw_expr_by_tx)}")
    if set(corrected_expr_by_tx) != {"txA1", "txA2", "txB1", "txB2"}:
        raise AssertionError(
            f"unexpected corrected transcripts: {sorted(corrected_expr_by_tx)}"
        )

    assert_close(float(raw_expr_by_tx["txA1"]["all_reads"]), 2.0)
    assert_close(float(raw_expr_by_tx["txA2"]["all_reads"]), 1.0)
    assert_close(float(raw_expr_by_tx["txB1"]["all_reads"]), 2.0)
    assert_close(float(raw_expr_by_tx["txB2"]["all_reads"]), 1.0)

    assert_close(float(corrected_expr_by_tx["txA1"]["all_reads"]), 1.5)
    assert_close(float(corrected_expr_by_tx["txA2"]["all_reads"]), 1.0)
    assert_close(float(corrected_expr_by_tx["txB1"]["all_reads"]), 1.5)
    assert_close(float(corrected_expr_by_tx["txB2"]["all_reads"]), 1.0)

    assert_close(float(corrected_expr_by_tx["txA1"]["isoform_fraction"]), 0.6, tol=1e-3)
    assert_close(float(corrected_expr_by_tx["txA2"]["isoform_fraction"]), 0.4, tol=1e-3)
    assert_close(float(corrected_expr_by_tx["txB1"]["isoform_fraction"]), 0.6, tol=1e-3)
    assert_close(float(corrected_expr_by_tx["txB2"]["isoform_fraction"]), 0.4, tol=1e-3)

    total_tpm = sum(float(row["TPM"]) for row in corrected_expr)
    assert_close(total_tpm, 1_000_000.0, tol=1e-3)

    read_to_genes_raw = defaultdict(set)
    for row in raw_tracking:
        read_to_genes_raw[row["read_name"]].add(row["gene_id"])
    multi_gene_reads_raw = sorted(
        read_name for read_name, genes in read_to_genes_raw.items() if len(genes) > 1
    )
    if multi_gene_reads_raw != ["rShared"]:
        raise AssertionError(f"unexpected multi-gene raw reads: {multi_gene_reads_raw}")

    corrected_shared = [
        row for row in corrected_tracking if row["read_name"] == "rShared"
    ]
    if len(corrected_shared) != 2:
        raise AssertionError(
            f"expected 2 corrected tracking rows for rShared, got {len(corrected_shared)}"
        )
    corrected_by_tx = {row["transcript_id"]: float(row["frac_assigned"]) for row in corrected_shared}
    assert_close(corrected_by_tx["txA1"], 0.5, tol=1e-3)
    assert_close(corrected_by_tx["txB1"], 0.5, tol=1e-3)

    unchanged_reads = [row for row in corrected_tracking if row["read_name"] != "rShared"]
    for row in unchanged_reads:
        assert_close(float(row["frac_assigned"]), 1.0, tol=1e-6)


if __name__ == "__main__":
    main()
