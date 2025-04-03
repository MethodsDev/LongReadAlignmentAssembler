#!/usr/bin/env python3

import sys, os, re
import pysam
import lmdb
import argparse
import csv
import logging
import shutil

FORMAT = (
    "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s:\n\t%(message)s\n"
)

logger = logging.getLogger()
logging.basicConfig(format=FORMAT, level=logging.INFO)


def main():
    parser = argparse.ArgumentParser(
        description="LRAA: Long Read Alignment Assembler",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--bam", type=str, required=True, help="target bam file")
    parser.add_argument(
        "--tracking", type=str, required=True, help="LRAA tracking file"
    )

    args = parser.parse_args()

    bam_file = args.bam
    tracking_file = args.tracking

    quant_read_tracking_lmdb_filename = tracking_file + ".lmdb"

    if os.path.exists(quant_read_tracking_lmdb_filename):
        logger.warning(
            f"{quant_read_tracking_lmdb_filename} exists already and will be removed and rebuilt now."
        )
        shutil.rmtree(quant_read_tracking_lmdb_filename)

    build_read_tracking_lmdb(tracking_file, quant_read_tracking_lmdb_filename)

    output_bam_file = bam_file + ".tagged.bam"
    annotate_bam_with_read_tracking_info(
        bam_file, quant_read_tracking_lmdb_filename, output_bam_file
    )

    sys.exit(0)


def build_read_tracking_lmdb(tracking_file, quant_read_tracking_lmdb_filename):

    logger.info("Building and populating read tracking db")

    env = lmdb.open(
        quant_read_tracking_lmdb_filename,
        map_size=int(1e11),
        # metasync=False,
        # sync=False,
    )

    with env.begin(write=True) as txn:

        record_counter = 0
        with open(tracking_file, "rt") as fh:

            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:

                record_counter += 1
                if record_counter % 1000 == 0:
                    print(f"\r[{record_counter}] ", file=sys.stderr, end="")

                gene_id = row["gene_id"]
                transcript_id = row["transcript_id"]
                read_name = row["read_name"]
                frac_read_assigned = row["frac_assigned"]
                key = read_name.encode("utf-8")
                existing_data = txn.get(key)
                lmdb_data = [
                    gene_id,
                    transcript_id,
                    frac_read_assigned,
                ]
                if existing_data:
                    # Append new data as comma-separated string
                    existing_data = existing_data.decode("utf-8")
                    combined_data = f"{existing_data},{','.join(lmdb_data)}"
                    # print("Storing APPENDED : {} -> {}".format(key, combined_data))
                else:
                    combined_data = ",".join(lmdb_data)
                    # print("Storing init : {} -> {}".format(key, combined_data))

                txn.put(key, combined_data.encode("utf-8"))

                # get_val = txn.get(key)
                # print("retrieved val: {}".format(get_val.decode()))

    logger.info("done building read tracker db.")

    return


def annotate_bam_with_read_tracking_info(
    bam_file, quant_read_tracking_lmdb_filename, output_bam_file
):

    logger.info("-writing bam file with read tracking annotations.")

    # open in read-only mode
    ofh_quant_read_tracking_lmdb = lmdb.open(
        quant_read_tracking_lmdb_filename, readonly=True, lock=False
    )

    with ofh_quant_read_tracking_lmdb.begin(write=False) as txn, pysam.AlignmentFile(
        bam_file, "rb", threads=4
    ) as bam_in, pysam.AlignmentFile(
        output_bam_file, "wb", template=bam_in, threads=4
    ) as bam_out:

        read_counter = 0
        num_reads_annotated = 0

        for read in bam_in.fetch(until_eof=True):

            read_counter += 1
            if read_counter % 1000 == 0:
                print(f"\r[{read_counter}] ", file=sys.stderr, end="")

            read_name = read.query_name
            data = txn.get(read_name.encode("utf-8"))
            if data:
                row = data.decode("utf-8").split(",")
                if read.has_tag("XG") or read.has_tag("XI"):
                    raise ValueError(
                        "Error, read already has XG or XI tag, so cannot set them without overwriting"
                    )

                gene_ids = ",".join(row[i] for i in range(0, len(row), 3))
                transcript_ids = ",".join(row[i + 1] for i in range(0, len(row), 3))
                frac_assigned = ",".join(row[i + 2] for i in range(0, len(row), 3))

                read.set_tag("XG", gene_ids, "Z")
                read.set_tag("XI", transcript_ids, "Z")
                read.set_tag("XF", frac_assigned, "Z")
                num_reads_annotated += 1

            bam_out.write(read)

        logger.info(
            "Number of reads assigned transcript info: {}".format(num_reads_annotated)
        )

        if num_reads_annotated == 0:
            raise RuntimeError(
                "Error, no reads were annotated... check the tracking and bam files match up"
            )

    ofh_quant_read_tracking_lmdb.close()

    return


if __name__ == "__main__":
    main()
