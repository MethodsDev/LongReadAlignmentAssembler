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
        description="collapse LRAA gtf by splice pattern.  If gtf is annotated with gene_symbol^ prefixes, those will be retained.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--gtf",
        type=str,
        required=True,
        help="LRAA gtf to collapse",
    )

    parser.add_argument(
        "--output_gtf",
        type=str,
        default="LRAA.merged.gtf",
        help="prefix for output filenames",
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

    input_gtf = args.gtf
    output_gtf = args.output_gtf

    ofh = open(output_gtf, "wt")

    logger.info(f"-capturing input transcripts from gtf {input_gtf}")
    contig_to_input_transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(
        input_gtf
    )

    for contig, transcript_obj_list in contig_to_input_transcripts.items():

        transcripts_to_output = list()

        splice_pattern_to_transcripts = defaultdict(list)

        for transcript_obj in transcript_obj_list:

            if transcript_obj.has_introns():
                intron_string = transcript_obj.get_introns_string()
                splice_pattern_code = Util_funcs.get_hash_code(intron_string)
                splice_pattern_to_transcripts[splice_pattern_code].append(
                    transcript_obj
                )
            else:
                transcripts_to_output.append(transcript_obj)

        # attempt splice pattern merge
        for (
            splice_pattern,
            transcripts_same_splice_pattern_list,
        ) in splice_pattern_to_transcripts.items():
            gene_id = transcripts_same_splice_pattern_list[0].get_gene_id()
            gene_symbol = None
            if "^" in gene_id:
                gene_symbol = gene_id.split("^")[0]
                new_transcript_id = "^".join([gene_symbol, splice_pattern])
            else:
                new_transcript_id = splice_pattern

            if len(transcripts_same_splice_pattern_list) == 1:
                transcript_obj = transcripts_same_splice_pattern_list[0]
                transcript_obj.set_transcript_id(new_transcript_id)
                transcripts_to_output.append(transcript_obj)
            else:
                merged_isoform = merge_isoforms(transcripts_same_splice_pattern_list)
                merged_isoform.set_transcript_id(new_transcript_id)
                transcripts_to_output.append(merged_isoform)

        transcripts_to_output = sorted(
            transcripts_to_output, key=lambda x: x._exon_segments[0][0]
        )

        for transcript_obj in transcripts_to_output:
            ofh.write(transcript_obj.to_GTF_format(include_TPM=False) + "\n")

    logger.info("Done.")

    sys.exit(0)


def merge_isoforms(transcript_obj_list):

    first_transcript_obj = transcript_obj_list[0]
    template_exon_coords = first_transcript_obj.get_exon_segments()
    min_lend = template_exon_coords[0][0]
    max_rend = template_exon_coords[-1][1]
    contig_acc = first_transcript_obj.get_contig_acc()
    contig_strand = first_transcript_obj.get_strand()
    gene_id = first_transcript_obj.get_gene_id()
    has_TSS = first_transcript_obj.has_TSS()
    has_PolyA = first_transcript_obj.has_PolyA()

    for transcript_obj in transcript_obj_list[1:]:
        exon_coords = transcript_obj.get_exon_segments()
        lend = exon_coords[0][0]
        rend = exon_coords[-1][1]

        if lend < min_lend:
            min_lend = lend

            if contig_strand == "+":
                has_TSS = transcript_obj.has_TSS()
            else:
                has_PolyA = transcript_obj.has_PolyA()

        if rend > max_rend:
            max_rend = rend
            if contig_strand == "+":
                has_PolyA = transcript_obj.has_PolyA()
            else:
                has_TSS = transcript_obj.has_TSS()

    template_exon_coords[0][0] = min_lend
    template_exon_coords[-1][1] = max_rend

    merged_transcript = Transcript(contig_acc, template_exon_coords, contig_strand)
    merged_transcript.set_gene_id(gene_id)

    merged_transcript._imported_has_TSS = has_TSS
    merged_transcript._imported_has_POLYA = has_PolyA

    return merged_transcript


if __name__ == "__main__":
    main()
