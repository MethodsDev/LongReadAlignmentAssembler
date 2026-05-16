import csv
import subprocess
import sys
from pathlib import Path


SCRIPT = Path(__file__).with_name("reassign_multigene_tracking_reads.py")


def write_tsv(path, fieldnames, rows):
    with open(path, "wt", newline="") as ofh:
        writer = csv.DictWriter(ofh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with open(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def run_reassign(tmp_path, expr_rows, tracking_rows):
    expr_fields = [
        "gene_id",
        "transcript_id",
        "uniq_reads",
        "all_reads",
        "isoform_fraction",
        "unique_gene_read_fraction",
        "TPM",
        "RPM_total_reads",
    ]
    tracking_fields = [
        "gene_id",
        "transcript_id",
        "transcript_splice_hash_code",
        "mp_id",
        "read_name",
        "frac_assigned",
        "read_weight",
    ]
    in_expr = tmp_path / "in.quant.expr"
    in_tracking = tmp_path / "in.quant.tracking"
    out_expr = tmp_path / "out.quant.expr"
    out_tracking = tmp_path / "out.quant.tracking"
    write_tsv(in_expr, expr_fields, expr_rows)
    write_tsv(in_tracking, tracking_fields, tracking_rows)

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--quant_expr",
            str(in_expr),
            "--tracking",
            str(in_tracking),
            "--output_expr",
            str(out_expr),
            "--output_tracking",
            str(out_tracking),
            "--tmp_dir",
            str(tmp_path),
        ],
        check=True,
    )
    return read_tsv(out_expr), read_tsv(out_tracking)


def test_no_cross_gene_reads_copies_inputs(tmp_path):
    expr_rows = [
        {
            "gene_id": "geneA",
            "transcript_id": "txA",
            "uniq_reads": "1",
            "all_reads": "1.0",
            "isoform_fraction": "1.000",
            "unique_gene_read_fraction": "1.000",
            "TPM": "500000.000",
            "RPM_total_reads": "500000.000",
        },
        {
            "gene_id": "geneB",
            "transcript_id": "txB",
            "uniq_reads": "1",
            "all_reads": "1.0",
            "isoform_fraction": "1.000",
            "unique_gene_read_fraction": "1.000",
            "TPM": "500000.000",
            "RPM_total_reads": "500000.000",
        },
    ]
    tracking_rows = [
        {
            "gene_id": "geneA",
            "transcript_id": "txA",
            "transcript_splice_hash_code": "txA",
            "mp_id": "MP1",
            "read_name": "rA",
            "frac_assigned": "1.000",
            "read_weight": "1.000",
        },
        {
            "gene_id": "geneB",
            "transcript_id": "txB",
            "transcript_splice_hash_code": "txB",
            "mp_id": "MP2",
            "read_name": "rB",
            "frac_assigned": "1.000",
            "read_weight": "1.000",
        },
    ]

    out_expr, out_tracking = run_reassign(tmp_path, expr_rows, tracking_rows)

    assert out_expr == expr_rows
    assert out_tracking == tracking_rows


def test_cross_gene_reads_are_reassigned_and_duplicate_candidates_collapsed(tmp_path):
    expr_rows = [
        {
            "gene_id": "geneA",
            "transcript_id": "txA",
            "uniq_reads": "1",
            "all_reads": "1.0",
            "isoform_fraction": "1.000",
            "unique_gene_read_fraction": "1.000",
            "TPM": "500000.000",
            "RPM_total_reads": "500000.000",
        },
        {
            "gene_id": "geneB",
            "transcript_id": "txB",
            "uniq_reads": "1",
            "all_reads": "1.0",
            "isoform_fraction": "1.000",
            "unique_gene_read_fraction": "1.000",
            "TPM": "500000.000",
            "RPM_total_reads": "500000.000",
        },
    ]
    tracking_rows = [
        {
            "gene_id": "geneA",
            "transcript_id": "txA",
            "transcript_splice_hash_code": "txA",
            "mp_id": "MP1",
            "read_name": "rA",
            "frac_assigned": "1.000",
            "read_weight": "1.000",
        },
        {
            "gene_id": "geneB",
            "transcript_id": "txB",
            "transcript_splice_hash_code": "txB",
            "mp_id": "MP2",
            "read_name": "rB",
            "frac_assigned": "1.000",
            "read_weight": "1.000",
        },
        {
            "gene_id": "geneA",
            "transcript_id": "txA",
            "transcript_splice_hash_code": "txA",
            "mp_id": "MP3",
            "read_name": "rX",
            "frac_assigned": "1.000",
            "read_weight": "1.000",
        },
        {
            "gene_id": "geneB",
            "transcript_id": "txB",
            "transcript_splice_hash_code": "txB",
            "mp_id": "MP4",
            "read_name": "rX",
            "frac_assigned": "1.000",
            "read_weight": "1.000",
        },
        {
            "gene_id": "geneB",
            "transcript_id": "txB",
            "transcript_splice_hash_code": "txB",
            "mp_id": "MP4_dup",
            "read_name": "rX",
            "frac_assigned": "1.000",
            "read_weight": "1.500",
        },
    ]

    out_expr, out_tracking = run_reassign(tmp_path, expr_rows, tracking_rows)

    assert len(out_tracking) == 4
    assert [row["frac_assigned"] for row in out_tracking if row["read_name"] == "rX"] == [
        "0.355",
        "0.645",
    ]
    assert [row["read_weight"] for row in out_tracking if row["read_name"] == "rX"] == [
        "1.000",
        "1.500",
    ]
    counts_by_tx = {row["transcript_id"]: row["all_reads"] for row in out_expr}
    assert counts_by_tx == {"txA": "1.4", "txB": "1.6"}


def test_component_requant_includes_within_gene_ambiguous_reads(tmp_path):
    expr_rows = [
        {
            "gene_id": "geneA",
            "transcript_id": "txA1",
            "uniq_reads": "1",
            "all_reads": "1.5",
            "isoform_fraction": "0.750",
            "unique_gene_read_fraction": "0.500",
            "TPM": "500000.000",
            "RPM_total_reads": "500000.000",
        },
        {
            "gene_id": "geneA",
            "transcript_id": "txA2",
            "uniq_reads": "0",
            "all_reads": "0.5",
            "isoform_fraction": "0.250",
            "unique_gene_read_fraction": "0.000",
            "TPM": "166666.667",
            "RPM_total_reads": "166666.667",
        },
        {
            "gene_id": "geneB",
            "transcript_id": "txB",
            "uniq_reads": "1",
            "all_reads": "1.0",
            "isoform_fraction": "1.000",
            "unique_gene_read_fraction": "1.000",
            "TPM": "333333.333",
            "RPM_total_reads": "333333.333",
        },
    ]
    tracking_rows = [
        {
            "gene_id": "geneA",
            "transcript_id": "txA1",
            "transcript_splice_hash_code": "txA1",
            "mp_id": "MP_unique_A",
            "read_name": "rA_unique",
            "frac_assigned": "1.000",
            "read_weight": "1.000",
        },
        {
            "gene_id": "geneA",
            "transcript_id": "txA1",
            "transcript_splice_hash_code": "txA1",
            "mp_id": "MP_within_A",
            "read_name": "rA_within",
            "frac_assigned": "0.500",
            "read_weight": "1.000",
        },
        {
            "gene_id": "geneA",
            "transcript_id": "txA2",
            "transcript_splice_hash_code": "txA2",
            "mp_id": "MP_within_A",
            "read_name": "rA_within",
            "frac_assigned": "0.500",
            "read_weight": "1.000",
        },
        {
            "gene_id": "geneA",
            "transcript_id": "txA2",
            "transcript_splice_hash_code": "txA2",
            "mp_id": "MP_cross_A",
            "read_name": "rCross",
            "frac_assigned": "1.000",
            "read_weight": "1.000",
        },
        {
            "gene_id": "geneB",
            "transcript_id": "txB",
            "transcript_splice_hash_code": "txB",
            "mp_id": "MP_cross_B",
            "read_name": "rCross",
            "frac_assigned": "1.000",
            "read_weight": "1.000",
        },
    ]

    out_expr, out_tracking = run_reassign(tmp_path, expr_rows, tracking_rows)

    within_rows = [
        row for row in out_tracking if row["read_name"] == "rA_within"
    ]
    assert len(within_rows) == 2
    assert [row["frac_assigned"] for row in within_rows] != ["0.500", "0.500"]

    assert {row["transcript_id"] for row in out_expr} == {"txA1", "txA2", "txB"}
