#!/usr/bin/env python3

import argparse
import csv
import logging
import sys


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def _derive_output_fieldnames(fieldnames):
    if fieldnames is None:
        raise RuntimeError("Error, no quant expr header found in input files")

    output_fieldnames = [field for field in fieldnames if field != "RPM_total_reads"]
    output_fieldnames.append("RPM_total_reads")
    return output_fieldnames


def merge_quant_expr_files(quant_files, output_filename):
    rows = list()
    total_reported_read_count = 0.0
    output_fieldnames = None

    for quant_file in quant_files:
        logger.info("Parsing %s", quant_file)
        with open(quant_file, "rt", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if reader.fieldnames is None:
                raise RuntimeError("Error, no header found in {}".format(quant_file))
            if "all_reads" not in reader.fieldnames:
                raise RuntimeError(
                    "Error, required column 'all_reads' not found in {}".format(
                        quant_file
                    )
                )
            if "TPM" not in reader.fieldnames:
                raise RuntimeError(
                    "Error, required column 'TPM' not found in {}".format(quant_file)
                )

            if output_fieldnames is None:
                output_fieldnames = _derive_output_fieldnames(reader.fieldnames)

            for row in reader:
                try:
                    counts = float(row.get("all_reads", 0) or 0)
                except ValueError:
                    raise RuntimeError(
                        "Error, non-numeric all_reads value in {}: {}".format(
                            quant_file, row.get("all_reads")
                        )
                    )
                total_reported_read_count += counts
                rows.append((row, counts))

    output_fieldnames = _derive_output_fieldnames(output_fieldnames)
    logger.info(
        "Writing %s with %d rows and %.3f total reported reads",
        output_filename,
        len(rows),
        total_reported_read_count,
    )

    with open(output_filename, "wt", newline="") as ofh:
        writer = csv.DictWriter(
            ofh, fieldnames=output_fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()

        for row, counts in rows:
            old_total_read_scaled_val = row.get("RPM_total_reads")
            if old_total_read_scaled_val is None:
                old_total_read_scaled_val = row.get("TPM", "0.000")

            out_row = {field: row.get(field, "") for field in output_fieldnames}
            tpm = (
                counts / total_reported_read_count * 1e6
                if total_reported_read_count > 0
                else 0
            )
            out_row["TPM"] = "{:.3f}".format(tpm)
            out_row["RPM_total_reads"] = old_total_read_scaled_val
            writer.writerow(out_row)


def main():
    parser = argparse.ArgumentParser(
        description="Merge LRAA quant expr shards and recompute final-report TPM values",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--output",
        required=True,
        type=str,
        help="output merged .quant.expr filename",
    )
    parser.add_argument(
        "--quant_files",
        required=True,
        nargs="+",
        type=str,
        help="input LRAA .quant.expr files to merge",
    )
    args = parser.parse_args()

    merge_quant_expr_files(args.quant_files, args.output)


if __name__ == "__main__":
    main()
