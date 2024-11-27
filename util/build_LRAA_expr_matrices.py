#!/usr/bin/env python3

import sys, os, re
import csv
from collections import defaultdict
import argparse
import logging


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description="build count and TPM matrices for LRAA.quant.expr files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="output prefix for .count.matrix and .TPM.matrix",
    )

    parser.add_argument(
        "--quant_files",
        type=str,
        nargs="+",
        required=True,
        help="space-delimited list of ${sample_name}.LRAA.quant.expr files to build into matrices (need at least 2)",
    )

    args = parser.parse_args()

    output_prefix = args.output_prefix
    quant_files = args.quant_files

    if len(quant_files) < 2:
        print(
            "Error, need at least 2 quant files specified to build the matrices",
            file=sys.stderr,
        )

    gene_matrix_data = defaultdict(lambda: defaultdict(float))
    transcript_matrix_data = defaultdict(lambda: defaultdict(float))

    sample_to_sum_counts = defaultdict(float)

    gene_ids = set()
    transcript_ids = set()

    for quant_file in quant_files:
        sample_name = os.path.basename(quant_file)
        sample_name = sample_name.replace(".LRAA.quant.expr", "")

        logger.info("Parsing {}".format(quant_file))
        with open(quant_file, "rt") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                gene_id = row["gene_id"]
                transcript_id = row["transcript_id"]
                read_count = float(row["all_reads"])

                transcript_id = gene_id + "^" + transcript_id

                gene_ids.add(gene_id)
                transcript_ids.add(transcript_id)

                sample_to_sum_counts[sample_name] += read_count
                gene_matrix_data[sample_name][gene_id] += read_count
                transcript_matrix_data[sample_name][transcript_id] += read_count

    # write gene count and expr matrices
    gene_count_matrix_filename = output_prefix + ".gene.counts.matrix"
    gene_TPM_matrix_filename = output_prefix + ".gene.TPM.matrix"

    logger.info(
        "-writing gene matrices: {} and {}".format(
            gene_count_matrix_filename, gene_TPM_matrix_filename
        )
    )
    write_matrices(
        gene_matrix_data,
        gene_ids,
        gene_count_matrix_filename,
        gene_TPM_matrix_filename,
        sample_to_sum_counts,
    )

    # write isoform count and expr matrices
    isoform_count_matrix_filename = output_prefix + ".isoform.counts.matrix"
    isoform_TPM_matrix_filename = output_prefix + ".isoform.TPM.matrix"

    logger.info(
        "-writing isoform matrices: {} and {}".format(
            isoform_count_matrix_filename, isoform_TPM_matrix_filename
        )
    )
    write_matrices(
        transcript_matrix_data,
        transcript_ids,
        isoform_count_matrix_filename,
        isoform_TPM_matrix_filename,
        sample_to_sum_counts,
    )

    logger.info("-done")

    sys.exit(0)


def write_matrices(
    matrix_data, ids, counts_matrix_filename, TPM_matrix_filename, sample_to_sum_counts
):

    with open(counts_matrix_filename, "wt") as counts_ofh:
        with open(TPM_matrix_filename, "wt") as TPM_ofh:

            sample_names = matrix_data.keys()

            print("\t" + "\t".join(sample_names), file=counts_ofh)
            print("\t" + "\t".join(sample_names), file=TPM_ofh)

            for id in ids:
                count_vals = [id]
                TPM_vals = [id]
                for sample_name in sample_names:
                    sample_sum_counts = sample_to_sum_counts[sample_name]

                    id_sample_count = matrix_data[sample_name][id]
                    count_vals.append("{:.3f}".format(id_sample_count))

                    id_tpm_val = id_sample_count / sample_sum_counts * 1e6
                    TPM_vals.append("{:.3f}".format(id_tpm_val))

                print("\t".join(count_vals), file=counts_ofh)
                print("\t".join(TPM_vals), file=TPM_ofh)

    return


if __name__ == "__main__":
    main()
