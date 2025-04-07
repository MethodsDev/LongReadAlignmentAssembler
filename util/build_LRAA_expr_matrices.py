#!/usr/bin/env python3

import sys, os, re
import csv
from collections import defaultdict
import argparse
import logging
import hashlib

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

    parser.add_argument(
        "--iso_by_unique_SP",
        action="store_true",
        default=False,
        help="aggregate spliced isoform info by splice pattern",
    )

    parser.add_argument(
        "--iso_unique_read_counts_only",
        action="store_true",
        default=False,
        help="restrict isoform count matrix to unique read assignments",
    )

    args = parser.parse_args()

    output_prefix = args.output_prefix
    quant_files = args.quant_files

    iso_by_unique_SP_flag = args.iso_by_unique_SP
    iso_unique_read_counts_only_flag = args.iso_unique_read_counts_only

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

    splice_pattern_to_hash_code = dict()
    transcript_id_to_splice_pattern = dict()

    def get_splice_pattern_hash_code(splice_pattern):
        if splice_pattern in splice_pattern_to_hash_code:
            return splice_pattern_to_hash_code[splice_pattern]
        else:
            hash_object = hashlib.sha256()
            hash_object.update(splice_pattern.encode("utf-8"))
            hex_digest = hash_object.hexdigest()
            hex_digest = str(hex_digest)
            splice_pattern_to_hash_code[splice_pattern] = hex_digest
            return hex_digest

    for quant_file in quant_files:
        sample_name = os.path.basename(quant_file)
        sample_name = sample_name.replace(".LRAA.quant.expr", "")

        logger.info("Parsing {}".format(quant_file))
        with open(quant_file, "rt") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                gene_id = row["gene_id"]

                if iso_by_unique_SP_flag and row["introns"] != "":
                    transcript_id = get_splice_pattern_hash_code(row["introns"])
                    transcript_id = gene_id + "^" + transcript_id
                    transcript_id_to_splice_pattern[transcript_id] = "\t".join(
                        [gene_id, row["transcript_id"], row["introns"]]
                    )
                else:
                    transcript_id = row["transcript_id"]

                read_count_for_gene = float(row["all_reads"])

                if iso_unique_read_counts_only_flag:
                    read_count_for_transcript = float(row["uniq_reads"])
                else:
                    read_count_for_transcript = float(row["all_reads"])

                # transcript_id = gene_id + "^" + transcript_id

                gene_ids.add(gene_id)
                transcript_ids.add(transcript_id)

                sample_to_sum_counts[sample_name] += read_count_for_gene
                gene_matrix_data[sample_name][gene_id] += read_count_for_gene
                transcript_matrix_data[sample_name][
                    transcript_id
                ] += read_count_for_transcript

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

    isoform_output_prefix = output_prefix
    if iso_by_unique_SP_flag:
        isoform_output_prefix += ".uniqueSP"
    if iso_unique_read_counts_only_flag:
        isoform_output_prefix += ".uniqReadsOnly"

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

    if iso_by_unique_SP_flag:
        # write gene/iso/SP mapping.
        with open(output_prefix + ".isoform_SP_mapping", "wt") as ofh:
            for transcript_id, SP_info in transcript_id_to_splice_pattern.items():
                print("\t".join([transcript_id, SP_info]), file=ofh)

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
