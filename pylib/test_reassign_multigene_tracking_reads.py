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
        "num_exons",
        "mp_id",
        "read_name",
        "frac_assigned",
    ]
    tracking_rows = [
        ["geneA", "txA1", "hA1", "1", "mpA1", "readA_unique", "1.000"],
        ["geneA", "txA2", "hA2", "1", "mpA2", "readA2_unique", "1.000"],
        ["geneB", "txB1", "hB1", "2", "mpB1", "readB_unique", "1.000"],
        ["geneA", "txA1", "hA1", "1", "mpSharedA", "readShared", "1.000"],
        ["geneA", "txA1", "hA1", "1", "mpSharedA_dup", "readShared", "1.000"],
        ["geneB", "txB1", "hB1", "2", "mpSharedB", "readShared", "1.000"],
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
    assert by_tx == {"txA1": "0.500000", "txB1": "0.500000"}
    assert {row["transcript_id"]: row["num_exons"] for row in shared_rows} == {
        "txA1": "1",
        "txB1": "2",
    }

    unchanged_rows = [row for row in tracking_rows_out if row["read_name"] != "readShared"]
    assert all(row["frac_assigned"] == "1.000000" for row in unchanged_rows)


def _write_two_gene_fixture(tmp_path):
    expr_in = tmp_path / "input.quant.expr"
    tracking_in = tmp_path / "input.quant.tracking.gz"

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
        ["geneA", "txA1", "2", "2.0", "1.000", "1.000", "500000.000", "1-10", "", "hA1", "500000.000"],
        ["geneB", "txB1", "2", "2.0", "1.000", "1.000", "500000.000", "21-30", "", "hB1", "500000.000"],
    ]
    _write_tsv(expr_in, expr_header, expr_rows)

    tracking_header = [
        "gene_id",
        "transcript_id",
        "transcript_splice_hash_code",
        "num_exons",
        "mp_id",
        "read_name",
        "frac_assigned",
    ]
    tracking_rows = [
        ["geneA", "txA1", "hA1", "1", "mpA1", "readA_unique", "1.000"],
        ["geneB", "txB1", "hB1", "1", "mpB1", "readB_unique", "1.000"],
        ["geneA", "txA1", "hA1", "1", "mpSharedA", "readShared", "1.000"],
        ["geneB", "txB1", "hB1", "1", "mpSharedB", "readShared", "1.000"],
    ]
    with gzip.open(tracking_in, "wt", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(tracking_header)
        writer.writerows(tracking_rows)

    return expr_in, tracking_in


def _run_reassign(expr_in, tracking_in, out_dir, extra_args=()):
    expr_out = out_dir / "output.quant.expr"
    tracking_out = out_dir / "output.quant.tracking.gz"
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
            *extra_args,
        ],
        check=True,
    )
    expr_by_tx = {row["transcript_id"]: row for row in _read_tsv(expr_out)}
    return expr_by_tx, _read_tsv_gz(tracking_out)


def test_unique_read_threshold_governs_reported_unique_counts(tmp_path):
    expr_in, tracking_in = _write_two_gene_fixture(tmp_path)

    whole_read_dir = tmp_path / "whole"
    whole_read_dir.mkdir()
    strict, _ = _run_reassign(expr_in, tracking_in, whole_read_dir)

    # the cross-gene read splits evenly, so it is not a whole read for either isoform
    assert strict["txA1"]["all_reads"] == "1.5"
    assert strict["txB1"]["all_reads"] == "1.5"
    assert strict["txA1"]["uniq_reads"] == "1"
    assert strict["txB1"]["uniq_reads"] == "1"

    partial_dir = tmp_path / "partial"
    partial_dir.mkdir()
    lenient, _ = _run_reassign(
        expr_in, tracking_in, partial_dir, extra_args=("--min_uniq_frac", "0.4")
    )

    assert lenient["txA1"]["all_reads"] == "1.5"
    assert lenient["txB1"]["all_reads"] == "1.5"
    assert lenient["txA1"]["uniq_reads"] == "2"
    assert lenient["txB1"]["uniq_reads"] == "2"


def test_candidate_with_no_prior_support_can_regain_reads(tmp_path):
    expr_in = tmp_path / "input.quant.expr"
    tracking_in = tmp_path / "input.quant.tracking.gz"

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
    # txB carries no support going in, so its abundance starts at zero.
    expr_rows = [
        ["geneA", "txA", "2", "2.0", "1.000", "1.000", "666666.667", "1-10", "", "hA", "666666.667"],
        ["geneB", "txB", "0", "0.0", "0.000", "0.000", "0.000", "21-30", "", "hB", "0.000"],
    ]
    _write_tsv(expr_in, expr_header, expr_rows)

    tracking_header = [
        "gene_id",
        "transcript_id",
        "transcript_splice_hash_code",
        "num_exons",
        "mp_id",
        "read_name",
        "frac_assigned",
    ]
    tracking_rows = [
        ["geneA", "txA", "hA", "1", "mpA1", "readA1", "1.000000"],
        ["geneA", "txA", "hA", "1", "mpA2", "readA2", "1.000000"],
        ["geneA", "txA", "hA", "1", "mpSharedA", "readShared", "1.000000"],
        ["geneB", "txB", "hB", "1", "mpSharedB", "readShared", "1.000000"],
    ]
    with gzip.open(tracking_in, "wt", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(tracking_header)
        writer.writerows(tracking_rows)

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    _, tracking_out = _run_reassign(expr_in, tracking_in, out_dir)

    shared = {
        row["transcript_id"]: float(row["frac_assigned"])
        for row in tracking_out
        if row["read_name"] == "readShared"
    }
    # A zero-abundance candidate must not be locked out of the expectation step.
    assert shared["txB"] > 0
    assert shared["txB"] < shared["txA"]
    assert abs(shared["txA"] + shared["txB"] - 1.0) < 1e-9
