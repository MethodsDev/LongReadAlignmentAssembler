#!/usr/bin/env python3

import sys, os, re

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../pylib"])
)

from Splice_graph import Splice_graph
from Transcript import Transcript, GTF_contig_to_transcripts
from LRAA import LRAA
import LRAA_Globals
import logging
import argparse
from collections import defaultdict
import Util_funcs

FORMAT = (
    "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s:\n\t%(message)s\n"
)

logger = logging.getLogger()
logging.basicConfig(format=FORMAT, level=logging.INFO)


def main():

    parser = argparse.ArgumentParser(
        description="Merge LRAA gtfs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--gtf",
        type=str,
        required=True,
        nargs="+",
        help="LRAA gtfs to merge",
    )

    parser.add_argument("--genome", type=str, required=True, help="target genome file")

    parser.add_argument(
        "--output_gtf", type=str, default="LRAA", help="prefix for output filenames"
    )

    parser.add_argument(
        "--debug",
        "-d",
        action="store_true",
        default=False,
        help="debug mode, more verbose",
    )

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
        LRAA_Globals.DEBUG = True

    gtf_list = args.gtf
    genome_fasta_file = args.genome
    output_gtf = args.output_gtf

    if len(gtf_list) < 2 and not LRAA_Globals.DEBUG:
        exit("Error, need at least two gtf files to merge")

    ofh = open(output_gtf, "wt")

    contig_strand_to_input_transcripts = defaultdict(list)

    for gtf_file in gtf_list:
        logger.info(f"-capturing input transcripts from gtf {gtf_file}")
        contig_to_input_transcripts = (
            GTF_contig_to_transcripts.parse_GTF_to_Transcripts(gtf_file)
        )
        for contig, transcript_obj_list in contig_to_input_transcripts.items():
            for transcript in transcript_obj_list:
                transcript_strand = transcript.get_strand()
                contig_strand_token = "{}^{}".format(contig, transcript_strand)
                contig_strand_to_input_transcripts[contig_strand_token].append(
                    transcript
                )

    for stranded_contig, transcript_list in contig_strand_to_input_transcripts.items():

        contig_acc, contig_strand = stranded_contig.split("^")

        contig_seq_str = Util_funcs.retrieve_contig_seq_from_fasta_file(
            contig_acc, genome_fasta_file
        )

        ## Build Splice Graph
        logger.info(f"\n// -building splice graph for {contig_acc}")
        sg = Splice_graph()
        sg.build_splice_graph_for_contig(
            contig_acc,
            contig_strand,
            contig_seq_str,
            alignments_bam_file=None,
            region_lend=None,
            region_rend=None,
            input_transcripts=transcript_list,
            quant_mode=False,
        )

        lraa_obj = LRAA(sg)

        lraa_obj.assign_transcripts_paths_in_graph(transcript_list)

        lraa_obj.build_multipath_graph(
            contig_acc,
            contig_strand,
            contig_seq_str,
            bam_file=None,
            input_transcripts=transcript_list,
        )

        # Define merged isoforms
        logger.info(f"\n// -begin merge of isoforms for {contig_acc}")
        transcripts = lraa_obj.reconstruct_isoforms()

        ## report transcripts in GTF format
        logger.info(
            "writing gtf output for {} [{}] containing {} transcripts".format(
                contig_acc, contig_strand, len(transcripts)
            )
        )

        for transcript in transcripts:
            ofh.write(transcript.to_GTF_format(include_TPM=False) + "\n")

    logger.info("Done.")

    sys.exit(0)


if __name__ == "__main__":
    main()
