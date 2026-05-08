#!/usr/bin/env python3

import csv
import gzip
import subprocess
import sys
from pathlib import Path


def _write_tsv(path, header, rows):
    with open(path, "wt", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def _read_tsv(path):
    with open(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _read_tsv_gz(path):
    with gzip.open(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_reassign_multigene_tracking_reads_reallocates_only_cross_gene_reads(tmp_path):
    expr_in = tmp_path / "input.quant.expr"
    tracking_in = tmp_path / "input.quant.tracking.gz"
    expr_out = tmp_path / "output.quant.expr"
    tracking_out = tmp_path / "output.quant.tracking.gz"

    expr_header = [
        "gene_id",
        "transcript_id",
        "uniq_reads",
        "all_reads",
        "isoform_fraction",
        "unique_gene_read_fraction",
        "TPM",
        "exons",
        "introns",
        "splice_hash_code",
        "RPM_total_reads",
    ]
    expr_rows = [
        ["geneA", "txA1", "2", "2.0", "0.667", "0.667", "400000.000", "1-10", "", "hA1", "500000.000"],
        ["geneA", "txA2", "1", "1.0", "0.333", "0.333", "200000.000", "11-20", "", "hA2", "250000.000"],
        ["geneB", "txB1", "2", "2.0", "1.000", "1.000", "400000.000", "21-30", "", "hB1", "500000.000"],
    ]
    _write_tsv(expr_in, expr_header, expr_rows)

    tracking_header = [
        "gene_id",
        "transcript_id",
        "transcript_splice_hash_code",
        "mp_id",
        "read_name",
        "frac_assigned",
    ]
    tracking_rows = [
        ["geneA", "txA1", "hA1", "mpA1", "readA_unique", "1.000"],
        ["geneA", "txA2", "hA2", "mpA2", "readA2_unique", "1.000"],
        ["geneB", "txB1", "hB1", "mpB1", "readB_unique", "1.000"],
        ["geneA", "txA1", "hA1", "mpSharedA", "readShared", "1.000"],
        ["geneB", "txB1", "hB1", "mpSharedB", "readShared", "1.000"],
    ]
    with gzip.open(tracking_in, "wt", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(tracking_header)
        writer.writerows(tracking_rows)

    script_path = Path(__file__).resolve().parents[1] / "util" / "reassign_multigene_tracking_reads.py"
    subprocess.run(
        [
            sys.executable,
            str(script_path),
            "--quant_expr",
            str(expr_in),
            "--tracking",
            str(tracking_in),
            "--output_expr",
            str(expr_out),
            "--output_tracking",
            str(tracking_out),
        ],
        check=True,
    )

    expr_rows_out = _read_tsv(expr_out)
    expr_by_tx = {row["transcript_id"]: row for row in expr_rows_out}
    assert expr_by_tx["txA1"]["all_reads"] == "1.5"
    assert expr_by_tx["txA2"]["all_reads"] == "1.0"
    assert expr_by_tx["txB1"]["all_reads"] == "1.5"
    assert expr_by_tx["txA1"]["uniq_reads"] == "1"
    assert expr_by_tx["txA2"]["uniq_reads"] == "1"
    assert expr_by_tx["txB1"]["uniq_reads"] == "1"
    assert expr_by_tx["txA1"]["TPM"] == "375000.000"
    assert expr_by_tx["txA2"]["TPM"] == "250000.000"
    assert expr_by_tx["txB1"]["TPM"] == "375000.000"
    assert expr_by_tx["txA1"]["RPM_total_reads"] == "375000.000"
    assert expr_by_tx["txA2"]["RPM_total_reads"] == "250000.000"
    assert expr_by_tx["txB1"]["RPM_total_reads"] == "375000.000"

    tracking_rows_out = _read_tsv_gz(tracking_out)
    shared_rows = [row for row in tracking_rows_out if row["read_name"] == "readShared"]
    assert len(shared_rows) == 2
    by_tx = {row["transcript_id"]: row["frac_assigned"] for row in shared_rows}
    assert by_tx == {"txA1": "0.500", "txB1": "0.500"}

    unchanged_rows = [row for row in tracking_rows_out if row["read_name"] != "readShared"]
    assert all(row["frac_assigned"] == "1.000" for row in unchanged_rows)
