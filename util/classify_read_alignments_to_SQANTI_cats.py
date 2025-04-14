#!/usr/bin/env python3

import sys, os, re
import logging
import argparse
from collections import defaultdict
import intervaltree as itree
import pysam

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../pylib"])
)

from Transcript import Transcript, GTF_contig_to_transcripts
from Pretty_alignment import Pretty_alignment

FORMAT = (
    "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s:\n\t%(message)s\n"
)

logger = logging.getLogger()
logging.basicConfig(format=FORMAT, level=logging.INFO)


def main():

    parser = argparse.ArgumentParser(
        description="Assign reads to sqanti categories",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--gtf",
        type=str,
        required=True,
        help="input GTF",
    )

    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="input bam",
    )

    parser.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="output prefix for bam and tsv files",
    )

    args = parser.parse_args()

    gtf_file = args.gtf
    bam_file = args.bam
    output_prefix = args.output_prefix

    logger.info("-parsing gtf_file: {}".format(gtf_file))
    contig_to_input_transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(
        gtf_file
    )

    stranded_chrom_exon_itrees = defaultdict(lambda: itree.IntervalTree())
    stranded_chrom_intron_itrees = defaultdict(lambda: itree.IntervalTree())
    splice_patterns_to_isoforms = defaultdict(set)
    intron_to_isoforms = defaultdict(set)
    transcript_id_to_obj = dict()
    stranded_splice_sites = set()

    build_isoform_data_structures(
        contig_to_input_transcripts,
        stranded_chrom_exon_itrees,
        stranded_chrom_intron_itrees,
        splice_patterns_to_isoforms,
        intron_to_isoforms,
        stranded_splice_sites,
        transcript_id_to_obj,
    )

    logger.info("Classifying reads from bam: {}".format(bam_file))
    bamfile_reader = pysam.AlignmentFile(bam_file, "rb")

    for read in bamfile_reader:
        chrom = bamfile_reader.get_reference_name(read.reference_id)
        if chrom is None:
            continue

        read_name = read.query_name
        read_strand = "+" if read.is_forward else "-"

        stranded_chrom = "{}:{}".format(chrom, read_strand)

        pretty_alignment = Pretty_alignment.get_pretty_alignment(read)
        align_span_lend, align_span_rend = pretty_alignment.get_alignment_span()

        alignment_segments = pretty_alignment.get_pretty_alignment_segments()
        num_alignment_segments = len(alignment_segments)

        ## Try classify read
        read_class_info = {
            "read_name": read_name,
            "alignment": pretty_alignment.get_pretty_alignment_string(chrom),
            "num_alignment_segments": num_alignment_segments,
            "matching_isoforms": "",
        }

        read_classified = False

        multi_exon_alignment_flag = pretty_alignment.has_introns()

        if multi_exon_alignment_flag:
            pretty_alignment_intron_string = pretty_alignment.get_introns_string(chrom)

            # check FSM
            if pretty_alignment_intron_string in splice_patterns_to_isoforms:
                matching_isoforms = ";".join(
                    sorted(
                        list(
                            splice_patterns_to_isoforms[pretty_alignment_intron_string]
                        )
                    )
                )
                read_class_info["sqanti_cat"] = "FSM"
                read_class_info["matching_isoforms"] = matching_isoforms
                read_classified = True

            # check ISM, NIC, NNIC
            if not read_classified:
                introns_all = None
                introns_any = set()
                introns_none = set()
                found_ref_shared_splice = False
                for intron in pretty_alignment.get_introns():
                    intron_lend, intron_rend = intron
                    if (
                        make_intron_token(chrom, read_strand, intron_lend)
                        in stranded_splice_sites
                        or make_intron_token(chrom, read_strand, intron_rend)
                        in stranded_splice_sites
                    ):
                        found_ref_shared_splice = True

                    intron_tok = make_intron_token(chrom, read_strand, intron)

                    if intron_tok in intron_to_isoforms:
                        isoforms_with_intron = intron_to_isoforms[intron_tok]
                        if introns_all is None:
                            introns_all = set()
                            introns_all.update(isoforms_with_intron)
                        else:
                            introns_all = introns_all & isoforms_with_intron

                        introns_any.update(isoforms_with_intron)
                    else:
                        introns_none.update(intron_tok)
                        if introns_all is not None:
                            introns_all = introns_all.clear()  # ISMs not possible.

                if (
                    introns_all is not None
                    and len(introns_all) > 0
                    and len(introns_none) == 0
                ):
                    read_class_info["sqanti_cat"] = "ISM"
                    read_class_info["matching_isoforms"] = ";".join(
                        sorted(list(introns_all))
                    )
                    read_classified = True
                elif len(introns_any) > 0 and len(introns_none) == 0:
                    read_class_info["sqanti_cat"] = "NIC"
                    read_classified = True
                elif (len(introns_any) > 0 or found_ref_shared_splice) and len(
                    introns_none
                ) > 0:
                    read_class_info["sqanti_cat"] = "NNIC"
                    read_classified = True

        else:
            # single exon mode.
            FSM_candidates = set()
            ISM_candidates = set()
            for stranded_chrom_exon_interval in stranded_chrom_exon_itrees[
                stranded_chrom
            ][align_span_lend : align_span_rend + 1]:
                overlapping_exon_lend = stranded_chrom_exon_interval.begin
                overlapping_exon_rend = stranded_chrom_exon_interval.end
                transcript_id = stranded_chrom_exon_interval.data

                transcript_obj = transcript_id_to_obj[transcript_id]
                transcript_lend, transcript_rend = transcript_obj.get_coords()
                if (
                    overlapping_exon_lend >= transcript_lend
                    and overlapping_exon_rend <= transcript_rend
                ):

                    transcript_id = transcript_obj.get_transcript_id()
                    if transcript_obj.get_num_exon_segments() == 1:
                        FSM_candidates.add(transcript_id)
                    else:
                        ISM_candidates.add(transcript_id)

            if len(FSM_candidates) > 0:
                read_class_info["sqanti_cat"] = "se_FSM"
                read_class_info["matching_isoforms"] = ";".join(
                    sorted(list(FSM_candidates))
                )
                read_classified = True
            elif len(ISM_candidates) > 0:
                read_class_info["sqanti_cat"] = "se_ISM"
                read_class_info["matching_isoforms"] = ";".join(
                    sorted(list(ISM_candidates))
                )
                read_classified = True

        #
        # check genic
        #

        if not read_classified:
            # check for genic - any overlap with exons
            for alignment_segment in alignment_segments:
                align_seg_lend, align_seg_rend = alignment_segment
                overlapping_exon_intervals = stranded_chrom_exon_itrees[stranded_chrom][
                    align_seg_lend : align_seg_rend + 1
                ]
                if len(overlapping_exon_intervals) > 0:
                    read_class_info["sqanti_cat"] = (
                        "genic" if multi_exon_alignment_flag else "se_genic"
                    )
                    read_classified = True
                    break
        #
        # check intronic
        #

        if not read_classified:
            # check for intronic.
            for alignment_segment in alignment_segments:
                align_seg_lend, align_seg_rend = alignment_segment
                overlapping_intron_intervals = stranded_chrom_intron_itrees[
                    stranded_chrom
                ][align_seg_lend : align_seg_rend + 1]
                if len(overlapping_intron_intervals) > 0:
                    read_class_info["sqanti_cat"] = (
                        "intronic" if multi_exon_alignment_flag else "se_intronic"
                    )
                    read_classified = True
                    break

        #
        # check antisense
        #

        if not read_classified:
            # see if overlaps exon from opposite strand
            antisense_strand = "+" if read_strand == "-" else "-"
            antisense_stranded_chrom = "{}:{}".format(chrom, antisense_strand)
            for alignment_segment in alignment_segments:
                align_seg_lend, align_seg_rend = alignment_segment
                overlapping_exon_intervals = stranded_chrom_intron_itrees[
                    antisense_stranded_chrom
                ][align_seg_lend : align_seg_rend + 1]
                if len(overlapping_exon_intervals) > 0:
                    read_class_info["sqanti_cat"] = (
                        "antisense" if multi_exon_alignment_flag else "se_antisense"
                    )
                    read_classified = True
                    break

        #
        # intergenic
        #
        if not read_classified:
            # only thing left is to call it intergenic.
            read_class_info["sqanti_cat"] = (
                "intergenic" if multi_exon_alignment_flag else "se_intergenic"
            )
            read_classified = True

        print(
            "\t".join(
                [
                    read_class_info["read_name"],
                    read_class_info["sqanti_cat"],
                    str(read_class_info["num_alignment_segments"]),
                    read_class_info["alignment"],
                    read_class_info["matching_isoforms"],
                ]
            )
        )

    sys.exit(0)


def make_intron_token(chrom, strand, coord_pair):
    return "{}:{}:{}".format(chrom, strand, coord_pair)


def build_isoform_data_structures(
    contig_to_input_transcripts,
    stranded_chrom_exon_itrees,
    stranded_chrom_intron_itrees,
    splice_patterns_to_isoforms,
    intron_to_isoforms,
    stranded_splice_sites,
    transcript_id_to_obj,
):

    logger.info("-building isoform data structures from parsed gtf")

    for contig, transcript_obj_list in contig_to_input_transcripts.items():

        for transcript_obj in transcript_obj_list:

            transcript_id = transcript_obj.get_transcript_id()
            transcript_id_to_obj[transcript_id] = transcript_obj
            strand = transcript_obj.get_strand()
            stranded_chrom = "{}:{}".format(contig, strand)

            # intron info
            if transcript_obj.has_introns():
                intron_str = transcript_obj.get_introns_string()
                splice_patterns_to_isoforms[intron_str].add(transcript_id)

                for intron in transcript_obj.get_introns():
                    intron_lend, intron_rend = intron
                    stranded_intron_lend_splice = make_intron_token(
                        contig, strand, intron_lend
                    )
                    stranded_intron_rend_splice = make_intron_token(
                        contig, strand, intron_rend
                    )
                    stranded_splice_sites.add(stranded_intron_lend_splice)
                    stranded_splice_sites.add(stranded_intron_rend_splice)

                    intron_token = make_intron_token(contig, strand, intron)
                    intron_to_isoforms[intron_token].add(transcript_id)
                    intron_lend, intron_rend = intron
                    stranded_chrom_intron_itrees[stranded_chrom][
                        intron_lend : intron_rend + 1
                    ] = transcript_id

            # exon info
            for exon in transcript_obj.get_exon_segments():
                exon_lend, exon_rend = exon
                stranded_chrom_exon_itrees[stranded_chrom][
                    exon_lend : exon_rend + 1
                ] = transcript_id

    return


if __name__ == "__main__":
    main()
