#!/usr/bin/env python3

import sys, os, re
import pysam
import lmdb
import argparse
import csv
import logging

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

    quant_read_tracking_lmdb_filename = os.path.join(tracking_file, ".lmdb")

    build_read_tracking_lmdb(tracking_file)

    output_bam_file = bam_file + ".tagged.bam"
    annotate_bam_with_read_tracking_info(
        bam_file, quant_read_tracking_lmdb_filename, output_bam_file
    )

    sys.exit(0)


def build_read_tracking_lmdb(tracking_file, quant_read_tracking_lmdb_filename):

    ofh_quant_read_tracking_lmdb = lmdb.open(
        quant_read_tracking_lmdb_filename,
        map_size=int(1e11),
        metasync=False,
        sync=False,
    )

    txn = ofh_quant_read_tracking_lmdb.begin(write=True)

    with open(tracking_file, "rt") as fh:

        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene_id = row["gene_id"]
            transcript_id = row["transcript_id"]
            read_name = row["read_name"]
            frac_assigned = row["frac_assigned"]
            key = readname.encode("utf-8")
            existing_data = txn.get(key)
            lmdb_data = [
                gene_id,
                transcript_id,
                "{:.3f}".format(frac_read_assigned),
            ]
            if existing_data:
                # Append new data as comma-separated string
                existing_data = existing_data.decode("utf-8")
                combined_data = f"{existing_data},{','.join(lmdb_data)}"
            else:
                combined_data = ",".join(lmdb_data)

            txn.put(key, combined_data.encode("utf-8"))

    del ofh_quant_read_tracking_lmdb  #  explicitly delete the writeable lmdb environment to close it since there is no .close()

    return


def annotate_bam_with_read_tracking_info(
    bam_file, quant_read_tracking_lmdb_filename, output_bam_file
):

    # open in read-only mode
    ofh_quant_read_tracking_lmdb = lmdb.open(
        quant_read_tracking_lmdb_filename, readonly=True, lock=False
    )

    with ofh_quant_read_tracking_lmdb.begin(write=False) as txn, pysam.AlignmentFile(
        bam_filename, "rb", threads=4
    ) as bam_in, pysam.AlignmentFile(
        output_bam_filename, "wb", template=bam_in, threads=4
    ) as bam_out:

        for read in bam_in.fetch(until_eof=True):
            data = txn.get(read.query_name.encode("utf-8"))
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

            bam_out.write(read)

    ofh_quant_read_tracking_lmdb.close()


if __name__ == "__main__":
    main()
