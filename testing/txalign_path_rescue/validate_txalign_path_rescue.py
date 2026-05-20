#!/usr/bin/env python3

import argparse
import ast
import csv
import sys
from pathlib import Path


def read_tsv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(iter_non_comment_lines(handle), delimiter="\t")
        return list(reader)


def iter_non_comment_lines(handle):
    for line in handle:
        if line.startswith("#"):
            continue
        yield line


def require(condition, message):
    if not condition:
        raise SystemExit(f"ERROR: {message}")


def parse_microexon(spec):
    chrom, coords, strand = spec.split(":")
    start_s, end_s = coords.split("-")
    return chrom, int(start_s), int(end_s), strand


def exon_lists_contain_microexon(exons_text, chrom, start, end):
    require(exons_text.startswith(f"{chrom}:("), f"unexpected exon field format: {exons_text}")
    try:
        exon_lists = ast.literal_eval(exons_text.split(")", 1)[1])
    except Exception as exc:
        raise SystemExit(f"ERROR: could not parse exons field {exons_text!r}: {exc}") from exc
    return [start, end] in exon_lists


def main():
    parser = argparse.ArgumentParser(
        description="Validate a transcriptome-alignment path-rescue execution test."
    )
    parser.add_argument("--case_dir", required=True)
    parser.add_argument("--output_prefix", required=True)
    parser.add_argument("--expected_read_names", required=True)
    parser.add_argument("--expected_total_reads", type=int, required=True)
    parser.add_argument("--expected_rescue_requested", type=int, required=True)
    parser.add_argument("--expected_rescued", type=int, required=True)
    parser.add_argument("--expected_rescue_requested_failed_genome", type=int, default=0)
    parser.add_argument("--expected_rescue_requested_unassigned_quant", type=int)
    parser.add_argument("--expected_microexon")
    args = parser.parse_args()
    expected_unassigned_quant = (
        args.expected_rescue_requested
        if args.expected_rescue_requested_unassigned_quant is None
        else args.expected_rescue_requested_unassigned_quant
    )

    case_dir = Path(args.case_dir)
    summary_path = case_dir / f"{args.output_prefix}.genome_tx_arb.summary.tsv"
    tracking_path = case_dir / f"{args.output_prefix}.quant.tracking"
    expr_path = case_dir / f"{args.output_prefix}.quant.expr"
    expected_names_path = case_dir / args.expected_read_names

    for path in [summary_path, tracking_path, expr_path, expected_names_path]:
        require(path.exists(), f"required file missing: {path}")

    summary_rows = read_tsv(summary_path)
    total = next((row for row in summary_rows if row["row_type"] == "TOTAL"), None)
    require(total is not None, f"TOTAL row missing from {summary_path}")

    checks = {
        "reads_total": args.expected_total_reads,
        "reads_kept_genome": 0,
        "reads_rescue_requested": args.expected_rescue_requested,
        "reads_rescue_rescued": args.expected_rescued,
        "reads_rescue_unrescued": 0,
        "reads_rescue_requested_failed_genome": args.expected_rescue_requested_failed_genome,
        "reads_rescue_requested_unassigned_quant": expected_unassigned_quant,
    }
    for field, expected in checks.items():
        observed = int(total[field])
        require(observed == expected, f"{field}: expected {expected}, observed {observed}")

    expected_names = {
        line.strip()
        for line in expected_names_path.read_text().splitlines()
        if line.strip()
    }
    tracking_rows = read_tsv(tracking_path)
    observed_names = {row["read_name"] for row in tracking_rows}
    require(
        observed_names == expected_names,
        f"tracking read set mismatch: expected {len(expected_names)}, observed {len(observed_names)}",
    )
    require(
        len(tracking_rows) == args.expected_rescued,
        f"tracking row count: expected {args.expected_rescued}, observed {len(tracking_rows)}",
    )

    if args.expected_microexon:
        chrom, start, end, _strand = parse_microexon(args.expected_microexon)
        expr_rows = read_tsv(expr_path)
        containing = [
            row for row in expr_rows
            if exon_lists_contain_microexon(row["exons"], chrom, start, end)
            and int(float(row["uniq_reads"])) == args.expected_rescued
        ]
        require(
            len(containing) == 1,
            f"expected exactly one rescued isoform with microexon {args.expected_microexon}; found {len(containing)}",
        )

    message = (
        "PASS: transcriptome-alignment rescue recovered "
        f"{args.expected_rescued}/{args.expected_rescue_requested} reads"
    )
    if args.expected_microexon:
        message += f" and reconstructed microexon {args.expected_microexon}"
    print(message)


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        sys.exit(1)
