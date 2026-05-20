#!/usr/bin/env python3

import argparse
import csv
import math
import sys


def _read_tsv(path):
    with open(path, "rt") as fh:
        return list(csv.DictReader(_iter_non_comment_lines(fh), delimiter="\t"))


def _iter_non_comment_lines(fh):
    for line in fh:
        if line.startswith("#"):
            continue
        yield line


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


def _pearson(x_values, y_values):
    if len(x_values) != len(y_values):
        raise AssertionError("Pearson vectors differ in length")
    if len(x_values) < 2:
        raise AssertionError("Need at least two paired values for Pearson correlation")

    x_mean = sum(x_values) / len(x_values)
    y_mean = sum(y_values) / len(y_values)
    numerator = sum((x - x_mean) * (y - y_mean) for x, y in zip(x_values, y_values))
    x_ss = sum((x - x_mean) ** 2 for x in x_values)
    y_ss = sum((y - y_mean) ** 2 for y in y_values)
    if x_ss <= 0 or y_ss <= 0:
        raise AssertionError("Cannot compute Pearson correlation for constant vector")
    return numerator / math.sqrt(x_ss * y_ss)


def main():
    parser = argparse.ArgumentParser(
        description="Validate SIRV quant output by correlation against a golden table."
    )
    parser.add_argument("actual")
    parser.add_argument("expected")
    parser.add_argument(
        "--metric",
        default="TPM",
        choices=["TPM", "RPM_total_reads", "all_reads"],
        help="numeric abundance column to use for Pearson correlation",
    )
    parser.add_argument(
        "--min_pearson",
        type=float,
        default=0.999,
        help="minimum acceptable Pearson correlation",
    )
    parser.add_argument(
        "--max_missing_fraction",
        type=float,
        default=0.001,
        help="maximum fraction of expected transcripts allowed to be absent from actual output",
    )
    parser.add_argument(
        "--allow_extra_transcripts",
        action="store_true",
        help="allow actual output to include transcripts absent from the expected table",
    )
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

    missing = sorted(set(expected_by_tx) - set(actual_by_tx))
    extra = sorted(set(actual_by_tx) - set(expected_by_tx))
    missing_fraction = len(missing) / len(expected_by_tx) if expected_by_tx else 0.0
    if missing_fraction > args.max_missing_fraction:
        raise AssertionError(
            "missing transcript fraction {:.3f} exceeds max {:.3f}; missing={}".format(
                missing_fraction, args.max_missing_fraction, missing
            )
        )
    if extra and not args.allow_extra_transcripts:
        raise AssertionError(
            "actual quant expr includes transcripts absent from expected table: {}".format(
                extra
            )
        )

    expected_values = []
    actual_values = []
    for transcript_id, expected in expected_by_tx.items():
        expected_values.append(_as_float(expected, args.metric))
        actual = actual_by_tx.get(transcript_id)
        actual_values.append(_as_float(actual, args.metric) if actual else 0.0)

    pearson = _pearson(expected_values, actual_values)
    if pearson < args.min_pearson:
        raise AssertionError(
            "{} Pearson {:.6f} is below minimum {:.6f}".format(
                args.metric, pearson, args.min_pearson
            )
        )

    total_reported_reads = sum(_as_float(row, "all_reads") for row in actual_rows)
    tpm_sum = sum(_as_float(row, "TPM") for row in actual_rows)
    _assert_close(tpm_sum, 1e6, "TPM sum", abs_tol=0.1)

    for transcript_id, actual in actual_by_tx.items():
        if transcript_id not in expected_by_tx:
            continue
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

    print(
        "Validated SIRV quant correlation: metric={}, Pearson={:.6f}, expected_rows={}, actual_rows={}, missing={} ({:.3%}), extra={}, TPM sum={:.3f}, total reads={:.1f}".format(
            args.metric,
            pearson,
            len(expected_rows),
            len(actual_rows),
            len(missing),
            missing_fraction,
            len(extra),
            tpm_sum,
            total_reported_reads,
        )
    )


if __name__ == "__main__":
    try:
        main()
    except AssertionError as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        sys.exit(1)
