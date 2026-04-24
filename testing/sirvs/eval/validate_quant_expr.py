#!/usr/bin/env python3

import argparse
import csv
import math
import sys


def _read_tsv(path):
    with open(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _as_float(row, column):
    try:
        return float(row[column])
    except Exception as exc:
        raise AssertionError(
            "Could not parse {}={!r}".format(column, row.get(column))
        ) from exc


def _assert_close(actual, expected, label, abs_tol=0.001):
    if not math.isclose(actual, expected, rel_tol=0, abs_tol=abs_tol):
        raise AssertionError(
            "{}: observed {}, expected {}".format(label, actual, expected)
        )


def main():
    parser = argparse.ArgumentParser(
        description="Validate SIRV quant output against a semantic golden table."
    )
    parser.add_argument("actual")
    parser.add_argument("expected")
    args = parser.parse_args()

    actual_rows = _read_tsv(args.actual)
    expected_rows = _read_tsv(args.expected)

    if not actual_rows:
        raise AssertionError("actual quant expr is empty")

    required_columns = [
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
    actual_columns = set(actual_rows[0].keys())
    missing_columns = [col for col in required_columns if col not in actual_columns]
    if missing_columns:
        raise AssertionError("actual quant expr missing columns: {}".format(missing_columns))

    actual_by_tx = {row["transcript_id"]: row for row in actual_rows}
    expected_by_tx = {row["transcript_id"]: row for row in expected_rows}

    if set(actual_by_tx) != set(expected_by_tx):
        missing = sorted(set(expected_by_tx) - set(actual_by_tx))
        extra = sorted(set(actual_by_tx) - set(expected_by_tx))
        raise AssertionError(
            "transcript set mismatch; missing={}, extra={}".format(missing, extra)
        )

    exact_columns = ["gene_id", "exons", "introns", "splice_hash_code"]
    numeric_columns = [
        "uniq_reads",
        "all_reads",
        "isoform_fraction",
        "unique_gene_read_fraction",
    ]

    for transcript_id, expected in expected_by_tx.items():
        actual = actual_by_tx[transcript_id]
        for col in exact_columns:
            if actual[col] != expected[col]:
                raise AssertionError(
                    "{} {}: observed {!r}, expected {!r}".format(
                        transcript_id, col, actual[col], expected[col]
                    )
                )
        for col in numeric_columns:
            _assert_close(
                _as_float(actual, col),
                _as_float(expected, col),
                "{} {}".format(transcript_id, col),
            )

    total_reported_reads = sum(_as_float(row, "all_reads") for row in actual_rows)
    tpm_sum = sum(_as_float(row, "TPM") for row in actual_rows)
    _assert_close(tpm_sum, 1e6, "TPM sum", abs_tol=0.1)

    for transcript_id, actual in actual_by_tx.items():
        expected_final_tpm = (
            _as_float(actual, "all_reads") / total_reported_reads * 1e6
            if total_reported_reads > 0
            else 0
        )
        _assert_close(
            _as_float(actual, "TPM"),
            expected_final_tpm,
            "{} final-report TPM".format(transcript_id),
        )
        _assert_close(
            _as_float(actual, "RPM_total_reads"),
            _as_float(expected_by_tx[transcript_id], "TPM"),
            "{} RPM_total_reads".format(transcript_id),
        )

    print(
        "Validated {} SIRV quant rows; TPM sum={:.3f}, total reads={:.1f}".format(
            len(actual_rows), tpm_sum, total_reported_reads
        )
    )


if __name__ == "__main__":
    try:
        main()
    except AssertionError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        sys.exit(1)
