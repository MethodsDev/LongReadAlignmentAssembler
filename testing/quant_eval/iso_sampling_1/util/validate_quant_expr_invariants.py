#!/usr/bin/env python3

import argparse
import csv
import math
import sys


def _read_rows(path):
    with open(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _float(row, column):
    try:
        return float(row[column])
    except Exception as exc:
        raise AssertionError(f"Could not parse {column}={row.get(column)!r}") from exc


def _assert_close(actual, expected, label, abs_tol):
    if not math.isclose(actual, expected, rel_tol=0, abs_tol=abs_tol):
        raise AssertionError(f"{label}: observed {actual}, expected {expected}")


def main():
    parser = argparse.ArgumentParser(description="Validate LRAA quant expr invariants.")
    parser.add_argument("quant_expr")
    parser.add_argument("--tpm-sum-tol", type=float, default=0.1)
    parser.add_argument("--row-tpm-tol", type=float, default=0.001)
    args = parser.parse_args()

    rows = _read_rows(args.quant_expr)
    if not rows:
        raise AssertionError("quant expr has no data rows")

    required_columns = {
        "gene_id",
        "transcript_id",
        "all_reads",
        "TPM",
        "RPM_total_reads",
    }
    missing = sorted(required_columns - set(rows[0].keys()))
    if missing:
        raise AssertionError(f"quant expr missing columns: {missing}")

    total_reported_reads = sum(_float(row, "all_reads") for row in rows)
    tpm_sum = sum(_float(row, "TPM") for row in rows)
    rpm_sum = sum(_float(row, "RPM_total_reads") for row in rows)

    _assert_close(tpm_sum, 1e6, "TPM sum", args.tpm_sum_tol)

    for row in rows:
        expected_tpm = (
            _float(row, "all_reads") / total_reported_reads * 1e6
            if total_reported_reads > 0
            else 0
        )
        _assert_close(
            _float(row, "TPM"),
            expected_tpm,
            f"{row['transcript_id']} TPM",
            args.row_tpm_tol,
        )

    print(
        f"Validated {len(rows)} rows in {args.quant_expr}; "
        f"TPM sum={tpm_sum:.3f}, RPM_total_reads sum={rpm_sum:.3f}"
    )


if __name__ == "__main__":
    try:
        main()
    except AssertionError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
