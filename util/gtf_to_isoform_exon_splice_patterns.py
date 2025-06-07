#!/usr/bin/env python3

import sys, os, re
import logging
import argparse
from collections import defaultdict
import intervaltree as itree
import csv
from hashlib import blake2s

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../pylib"])
)

from Transcript import Transcript, GTF_contig_to_transcripts

FORMAT = (
    "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s:\n\t%(message)s\n"
)

logger = logging.getLogger()
logging.basicConfig(format=FORMAT, level=logging.INFO)


def main():

    parser = argparse.ArgumentParser(
        description="extract exon and intron splice patterns from gtf",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--gtf",
        type=str,
        required=True,
        help="input GTF file ",
    )

    args = parser.parse_args()
    input_gtf_filename = args.gtf

    logger.info("-parsing isoforms from {}".format(input_gtf_filename))
    contig_to_input_transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(
        input_gtf_filename
    )

    tsv_writer = csv.DictWriter(
        sys.stdout,
        fieldnames=["gene_id", "transcript_id", "exons", "introns", "introns_hashcode"],
        delimiter="\t",
        lineterminator="\n",
    )
    tsv_writer.writeheader()

    for contig, transcript_list in contig_to_input_transcripts.items():
        for transcript_obj in transcript_list:
            row = dict()
            row["gene_id"] = transcript_obj.get_gene_id()
            row["transcript_id"] = transcript_obj.get_transcript_id()
            row["exons"] = transcript_obj.get_exons_string()
            row["introns"] = None
            row["introns_hashcode"] = None

            if transcript_obj.has_introns():
                row["introns"] = transcript_obj.get_introns_string()
                row["introns_hashcode"] = get_splice_pattern_hash_code(row["introns"])

            tsv_writer.writerow(row)

    sys.exit(0)


def get_splice_pattern_hash_code(splice_pattern):
    hash_object = blake2s(digest_size=11)
    hash_object.update(splice_pattern.encode("utf-8"))
    hex_digest = hash_object.hexdigest()
    hex_digest = str(hex_digest)
    return hex_digest


if __name__ == "__main__":
    main()
