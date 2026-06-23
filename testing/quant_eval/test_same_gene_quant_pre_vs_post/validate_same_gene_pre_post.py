#!/usr/bin/env python3

import csv
import sys


def iter_non_comment_lines(fh):
    for line in fh:
        if not line.startswith("#"):
            yield line


def read_expr(path):
    with open(path, "rt", newline="") as fh:
        return list(csv.DictReader(iter_non_comment_lines(fh), delimiter="\t"))


def main():
    if len(sys.argv) != 4:
        raise SystemExit(
            "usage: validate_same_gene_pre_post.py run.log pre.quant.expr final.quant.expr"
        )

    run_log, pre_expr_path, final_expr_path = sys.argv[1:]

    with open(run_log, "rt") as fh:
        log_text = fh.read()
    if "found 0 cross-gene reads" not in log_text:
        raise AssertionError("LRAA final cross-gene EM did not report 0 cross-gene reads")

    pre_rows = read_expr(pre_expr_path)
    final_rows = read_expr(final_expr_path)

    if len(pre_rows) <= 50:
        raise AssertionError(f"Expected a complex quant table; found only {len(pre_rows)} rows")

    pre_gene_ids = {row["gene_id"] for row in pre_rows}
    final_gene_ids = {row["gene_id"] for row in final_rows}
    if pre_gene_ids != {"ABLIM1_ALL_ISOFORMS"}:
        raise AssertionError(f"Unexpected pre gene IDs: {sorted(pre_gene_ids)}")
    if final_gene_ids != {"ABLIM1_ALL_ISOFORMS"}:
        raise AssertionError(f"Unexpected final gene IDs: {sorted(final_gene_ids)}")

    if pre_rows != final_rows:
        for idx, (pre_row, final_row) in enumerate(zip(pre_rows, final_rows), start=1):
            if pre_row != final_row:
                raise AssertionError(
                    f"pre/final quant expr mismatch at data row {idx}: "
                    f"pre={pre_row} final={final_row}"
                )
        raise AssertionError(
            f"pre/final quant expr row count mismatch: {len(pre_rows)} vs {len(final_rows)}"
        )

    print(
        f"Validated {len(pre_rows)} same-gene quant rows; pre-cross-gene-EM and final quant expr match exactly."
    )


if __name__ == "__main__":
    try:
        main()
    except AssertionError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
