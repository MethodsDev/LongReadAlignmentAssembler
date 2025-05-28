#!/usr/bin/env python3

import sys, os, re
import logging
from collections import defaultdict
import intervaltree as itree


from Transcript import Transcript, GTF_contig_to_transcripts
from Pretty_alignment import Pretty_alignment

FORMAT = (
    "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s:\n\t%(message)s\n"
)

logger = logging.getLogger()
logging.basicConfig(format=FORMAT, level=logging.INFO)


######################################################
#
#  SQANTI-like categories assigned
#
#  Spliced features:
#
#      FSM: contains full set of spliced-out introns as a reference isoform
#      ISM:


class SQANTI_like_annotator:

    MIN_FSM_SE_FRAC_OVERLAP = 0.9

    def __init__(self, ref_annot_gtf):

        logger.info("-parsing gtf_file: {}".format(ref_annot_gtf))

        self.contig_to_input_transcripts = (
            GTF_contig_to_transcripts.parse_GTF_to_Transcripts(ref_annot_gtf)
        )

        self.stranded_chrom_exon_itrees = defaultdict(lambda: itree.IntervalTree())
        self.stranded_chrom_intron_itrees = defaultdict(lambda: itree.IntervalTree())
        self.splice_patterns_to_isoforms = defaultdict(set)
        self.intron_to_isoforms = defaultdict(set)
        self.transcript_id_to_obj = dict()
        self.stranded_splice_sites = set()

        self._build_isoform_data_structures()

    def classify_alignment_or_isoform(
        self, feature_chrom, feature_strand, feature_name, transcribed_feature_obj
    ):

        assert type(transcribed_feature_obj) in (
            Transcript,
            Pretty_alignment,
        ), "Error, cannot process transcribed_feature: {}, type {}".format(
            transcribed_feature_obj, type(transcribed_feature_obj)
        )

        stranded_chrom = "{}:{}".format(feature_chrom, feature_strand)
        feature_span_lend, feature_span_rend = (
            transcribed_feature_obj.get_alignment_span()
            if type(transcribed_feature_obj) == Pretty_alignment
            else transcribed_feature_obj.get_coords()
        )

        exon_segments = (
            transcribed_feature_obj.get_pretty_alignment_segments()
            if type(transcribed_feature_obj) == Pretty_alignment
            else transcribed_feature_obj.get_exon_segments()
        )
        num_exon_segments = len(exon_segments)

        ## Try classify read
        feature_class_info = {
            "feature_name": feature_name,
            "structure": (
                transcribed_feature_obj.get_pretty_alignment_string(feature_chrom)
                if type(transcribed_feature_obj) == Pretty_alignment
                else transcribed_feature_obj.get_exons_string()
            ),
            "num_exon_segments": num_exon_segments,
            "matching_isoforms": "",
        }

        feature_classified = False

        multi_exon_alignment_flag = transcribed_feature_obj.has_introns()

        if multi_exon_alignment_flag:
            feature_intron_string = (
                transcribed_feature_obj.get_introns_string(feature_chrom)
                if type(transcribed_feature_obj) == Pretty_alignment
                else transcribed_feature_obj.get_introns_string()
            )

            # check FSM
            if feature_intron_string in self.splice_patterns_to_isoforms:
                matching_isoforms = ",".join(
                    sorted(
                        list(self.splice_patterns_to_isoforms[feature_intron_string])
                    )
                )
                feature_class_info["sqanti_cat"] = "FSM"
                feature_class_info["matching_isoforms"] = matching_isoforms
                feature_classified = True

            # check ISM, NIC, NNIC
            if not feature_classified:
                introns_all = None  # isoforms all contain introns
                introns_any = set()  # any isoform that contains the intron
                introns_none = set()
                found_ref_shared_splice = False
                for intron in transcribed_feature_obj.get_introns():
                    intron_lend, intron_rend = intron
                    if (
                        make_intron_token(feature_chrom, feature_strand, intron_lend)
                        in self.stranded_splice_sites
                        or make_intron_token(feature_chrom, feature_strand, intron_rend)
                        in self.stranded_splice_sites
                    ):
                        found_ref_shared_splice = True

                    intron_tok = make_intron_token(
                        feature_chrom, feature_strand, intron
                    )

                    if intron_tok in self.intron_to_isoforms:
                        isoforms_with_intron = self.intron_to_isoforms[intron_tok]
                        if introns_all is None:
                            # init
                            introns_all = set()
                            introns_all.update(isoforms_with_intron)
                        else:
                            introns_all = introns_all & isoforms_with_intron

                        introns_any.update(isoforms_with_intron)
                    else:
                        # novel intron splice pattern
                        introns_none.update(intron_tok)
                        if introns_all is not None:
                            introns_all = introns_all.clear()  # ISMs not possible.

                if (
                    introns_all is not None
                    and len(introns_all) > 0
                    and len(introns_none) == 0
                ):

                    ordered_splice_matched_isoforms = (
                        self.restrict_splice_matched_isoforms(
                            introns_all, transcribed_feature_obj
                        )
                    )

                    if len(ordered_splice_matched_isoforms) > 0:
                        feature_class_info["sqanti_cat"] = "ISM"
                        feature_class_info["matching_isoforms"] = ",".join(
                            sorted(list(ordered_splice_matched_isoforms))
                        )
                        feature_classified = True

                if (
                    (not feature_classified)
                    and len(introns_any) > 0
                    and len(introns_none) == 0
                ):
                    feature_class_info["sqanti_cat"] = "NIC"
                    feature_classified = True
                elif (len(introns_any) > 0 or found_ref_shared_splice) and len(
                    introns_none
                ) > 0:
                    feature_class_info["sqanti_cat"] = "NNIC"
                    feature_classified = True

        else:
            # single exon mode.
            FSM_candidates = set()
            ISM_candidates = set()
            for stranded_chrom_exon_interval in self.stranded_chrom_exon_itrees[
                stranded_chrom
            ][feature_span_lend : feature_span_rend + 1]:
                overlapping_exon_lend = stranded_chrom_exon_interval.begin
                overlapping_exon_rend = stranded_chrom_exon_interval.end
                transcript_id = stranded_chrom_exon_interval.data

                transcript_obj = self.transcript_id_to_obj[transcript_id]
                transcript_lend, transcript_rend = transcript_obj.get_coords()

                #
                # Only single-exon target transcripts can get FSM or ISM categories with single-exon alignments
                #

                # check amount of overlap

                if transcript_obj.get_num_exon_segments() == 1:

                    frac_transcript_overlap = eval_frac_transcript_overlap(
                        [transcript_lend, transcript_rend],
                        [feature_span_lend, feature_span_rend],
                    )

                    if (
                        frac_transcript_overlap
                        >= SQANTI_like_annotator.MIN_FSM_SE_FRAC_OVERLAP
                    ):
                        FSM_candidates.add(transcript_id)
                    else:
                        ISM_candidates.add(transcript_id)

            if len(FSM_candidates) > 0:
                feature_class_info["sqanti_cat"] = "se_FM"
                feature_class_info["matching_isoforms"] = ",".join(
                    sorted(list(FSM_candidates))
                )
                feature_classified = True
            elif len(ISM_candidates) > 0:
                feature_class_info["sqanti_cat"] = "se_IM"
                feature_class_info["matching_isoforms"] = ",".join(
                    sorted(list(ISM_candidates))
                )
                feature_classified = True

        #
        # check genic
        #

        if not feature_classified:
            # check for genic - any overlap with exons
            for exon_segment in exon_segments:
                exon_seg_lend, exon_seg_rend = exon_segment
                overlapping_exon_intervals = self.stranded_chrom_exon_itrees[
                    stranded_chrom
                ][exon_seg_lend : exon_seg_rend + 1]
                if len(overlapping_exon_intervals) > 0:
                    feature_class_info["sqanti_cat"] = (
                        "genic" if multi_exon_alignment_flag else "se_genic"
                    )
                    feature_classified = True
                elif (
                    len(overlapping_exon_intervals) == 1
                    and not multi_exon_alignment_flag
                ):
                    feature_class_info["sqanti_cat"] = "se_exonic"
                    feature_classified = True

                if feature_classified:
                    break
        #
        # check intronic
        #

        if not feature_classified:
            # check for intronic.
            for exon_segment in exon_segments:
                exon_seg_lend, exon_seg_rend = exon_segment
                overlapping_intron_intervals = self.stranded_chrom_intron_itrees[
                    stranded_chrom
                ][exon_seg_lend : exon_seg_rend + 1]
                if len(overlapping_intron_intervals) > 0:
                    feature_class_info["sqanti_cat"] = (
                        "intronic" if multi_exon_alignment_flag else "se_intronic"
                    )
                    feature_classified = True
                    break

        #
        # check antisense
        #

        if not feature_classified:
            # see if overlaps exon from opposite strand
            antisense_strand = "+" if feature_strand == "-" else "-"
            antisense_stranded_chrom = "{}:{}".format(feature_chrom, antisense_strand)
            for exon_segment in exon_segments:
                exon_seg_lend, exon_seg_rend = exon_segment
                overlapping_exon_intervals = self.stranded_chrom_intron_itrees[
                    antisense_stranded_chrom
                ][exon_seg_lend : exon_seg_rend + 1]
                if len(overlapping_exon_intervals) > 0:
                    feature_class_info["sqanti_cat"] = (
                        "antisense" if multi_exon_alignment_flag else "se_antisense"
                    )
                    feature_classified = True
                    break

        #
        # intergenic
        #
        if not feature_classified:
            # only thing left is to call it intergenic.
            feature_class_info["sqanti_cat"] = (
                "intergenic" if multi_exon_alignment_flag else "se_intergenic"
            )
            feature_classified = True

        return feature_class_info

    def _build_isoform_data_structures(self):

        logger.info("-building isoform data structures from parsed gtf")

        for contig, transcript_obj_list in self.contig_to_input_transcripts.items():

            for transcript_obj in transcript_obj_list:

                transcript_id = transcript_obj.get_transcript_id()
                self.transcript_id_to_obj[transcript_id] = transcript_obj
                strand = transcript_obj.get_strand()
                stranded_chrom = "{}:{}".format(contig, strand)

                # intron info
                if transcript_obj.has_introns():
                    intron_str = transcript_obj.get_introns_string()
                    self.splice_patterns_to_isoforms[intron_str].add(transcript_id)

                    for intron in transcript_obj.get_introns():
                        intron_lend, intron_rend = intron
                        stranded_intron_lend_splice = make_intron_token(
                            contig, strand, intron_lend
                        )
                        stranded_intron_rend_splice = make_intron_token(
                            contig, strand, intron_rend
                        )
                        self.stranded_splice_sites.add(stranded_intron_lend_splice)
                        self.stranded_splice_sites.add(stranded_intron_rend_splice)

                        intron_token = make_intron_token(contig, strand, intron)
                        self.intron_to_isoforms[intron_token].add(transcript_id)
                        intron_lend, intron_rend = intron
                        self.stranded_chrom_intron_itrees[stranded_chrom][
                            intron_lend : intron_rend + 1
                        ] = transcript_id

                # exon info
                for exon in transcript_obj.get_exon_segments():
                    exon_lend, exon_rend = exon
                    self.stranded_chrom_exon_itrees[stranded_chrom][
                        exon_lend : exon_rend + 1
                    ] = transcript_id

        return

    def restrict_splice_matched_isoforms(
        self, isoforms_with_matching_introns, transcribed_feature_obj
    ):

        feature_introns = transcribed_feature_obj.get_introns()

        locally_matching_isoforms = set()

        for isoform_id in isoforms_with_matching_introns:
            isoform_introns = self.transcript_id_to_obj[isoform_id].get_introns()

            idx_start = isoform_introns.index(feature_introns[0])
            isoform_introns = isoform_introns[idx_start:]

            idx_end = isoform_introns.index(feature_introns[-1])
            isoform_introns = isoform_introns[: idx_end + 1]

            if isoform_introns == feature_introns:
                locally_matching_isoforms.add(isoform_id)

        return isoform_id


def make_intron_token(chrom, strand, coord_pair):
    return "{}:{}:{}".format(chrom, strand, coord_pair)


def eval_frac_transcript_overlap(trans_coords, align_coords):
    trans_lend, trans_rend = trans_coords
    align_lend, align_rend = align_coords

    assert (
        trans_lend <= align_rend and trans_rend >= align_lend
    ), "Error, no overlap in coords!"

    overlapping_coords = sorted([trans_lend, trans_rend, align_lend, align_rend])[1:3]
    overlap_len = overlapping_coords[1] - overlapping_coords[0] + 1
    trans_len = trans_rend - trans_lend + 1
    frac_overlap = overlap_len / trans_len

    return frac_overlap
