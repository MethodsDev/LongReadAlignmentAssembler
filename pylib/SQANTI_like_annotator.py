#!/usr/bin/env python3

import sys, os, re
import logging
from collections import defaultdict
import intervaltree as itree
from MockTranscript import MockTranscript

from Transcript import Transcript, GTF_contig_to_transcripts
from Pretty_alignment import Pretty_alignment
import LRAA_Globals

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
#
#  FSM subcategories (two independent columns, filled only for sqanti_cat FSM; empty
#  for every other category, including the mono-exonic se_* ones):
#
#      FSM_full_length_isoforms: the FSM-matched reference isoforms (same intron chain)
#          whose 5' and 3' ends both lie within fsm_full_length_max_end_dist bases of
#          the query's ends. Reads are measured by their alignment span (soft clips
#          excluded). Both ends are compared, so strand does not matter.
#
#      FSM_ISM_of_isoforms: other reference isoforms whose intron chain contains the
#          query's chain as a proper contiguous subset, i.e. the query would be an ISM
#          of them if it did not already match a reference chain exactly. Structural
#          only, like the ISM test itself: any reference isoform on the same contig and
#          strand sharing every query intron in order counts, with no gene_id check.
#          5'/3'/internal fragment direction is not recorded (the ISM category does not
#          record it either).
#
#      FSM_ISM_of_genes: the gene of each FSM_ISM_of_isoforms entry, in the same order,
#          as gene_name^gene_id (gene_id alone when the GTF gives no gene_name).
#
#      FSM_ISM_of_gene_scope: whether those containing isoforms belong to a gene of the
#          FSM-matched isoforms: same_gene, other_gene (none does), same_and_other_gene,
#          or empty when FSM_ISM_of_isoforms is empty.
#
#  A non-empty field is the flag; its value names the reference isoforms.


class SQANTI_like_annotator:

    MIN_FSM_SE_FRAC_OVERLAP = 0.9

    def __init__(self, ref_annot_gtf, fsm_full_length_max_end_dist=None):

        # default from LRAA_Globals; callers that run in worker processes pass it
        # explicitly, since a config override made in the parent is not inherited
        # under the spawn start method.
        if fsm_full_length_max_end_dist is None:
            fsm_full_length_max_end_dist = LRAA_Globals.config[
                "sqlike_fsm_full_length_max_end_dist"
            ]
        self.fsm_full_length_max_end_dist = int(fsm_full_length_max_end_dist)

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

        assert isinstance(
            transcribed_feature_obj, (Transcript, Pretty_alignment, MockTranscript)
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
            "FSM_full_length_isoforms": "",
            "FSM_ISM_of_isoforms": "",
            "FSM_ISM_of_genes": "",
            "FSM_ISM_of_gene_scope": "",
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
                fsm_isoforms = self.splice_patterns_to_isoforms[feature_intron_string]
                matching_isoforms = ",".join(sorted(list(fsm_isoforms)))
                feature_class_info["sqanti_cat"] = "FSM"
                feature_class_info["matching_isoforms"] = matching_isoforms
                feature_classified = True

                # FSM subcategories (do not affect sqanti_cat)
                feature_class_info["FSM_full_length_isoforms"] = ",".join(
                    sorted(
                        self._full_length_isoforms(
                            fsm_isoforms, feature_span_lend, feature_span_rend
                        )
                    )
                )
                containing_isoforms = sorted(
                    self._isoforms_containing_intron_chain(
                        feature_chrom, feature_strand, transcribed_feature_obj
                    )
                    - fsm_isoforms
                )
                if containing_isoforms:
                    feature_class_info.update(
                        self._containment_fields(fsm_isoforms, containing_isoforms)
                    )

            # check ISM, NIC, NNIC
            if not feature_classified:
                introns_all = None  # isoforms all contain introns
                introns_any = set()  # any isoform that contains the intron
                introns_none = set()
                found_ref_shared_splice = False
                found_nonref_splice = False
                for intron in transcribed_feature_obj.get_introns():
                    intron_lend, intron_rend = intron
                    if (
                        make_intron_token(feature_chrom, feature_strand, intron_lend)
                        in self.stranded_splice_sites
                    ):
                        found_ref_shared_splice = True
                    else:
                        found_nonref_splice = True

                    if (
                        make_intron_token(feature_chrom, feature_strand, intron_rend)
                        in self.stranded_splice_sites
                    ):
                        found_ref_shared_splice = True
                    else:
                        found_nonref_splice = True

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
                    and found_ref_shared_splice
                    and (not found_nonref_splice)
                ):
                    feature_class_info["sqanti_cat"] = "NIC"
                    feature_classified = True
                elif found_ref_shared_splice and found_nonref_splice:
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

                overlapping_intron_intervals = self.stranded_chrom_intron_itrees[
                    stranded_chrom
                ][exon_seg_lend : exon_seg_rend + 1]

                if (
                    len(overlapping_exon_intervals) > 0
                    and len(overlapping_intron_intervals) > 0
                ):
                    feature_class_info["sqanti_cat"] = (
                        "genic" if multi_exon_alignment_flag else "se_genic"
                    )
                    feature_classified = True
                elif (
                    len(overlapping_exon_intervals) > 0
                    and len(overlapping_intron_intervals) == 0
                ):
                    feature_class_info["sqanti_cat"] = (
                        "se_exonic" if not multi_exon_alignment_flag else "genic"
                    )
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
                overlapping_exon_intervals = self.stranded_chrom_exon_itrees[
                    antisense_stranded_chrom
                ][exon_seg_lend : exon_seg_rend + 1]
                overlapping_intron_intervals = self.stranded_chrom_intron_itrees[
                    antisense_stranded_chrom
                ][exon_seg_lend : exon_seg_rend + 1]

                if (
                    len(overlapping_exon_intervals) > 0
                    or len(overlapping_intron_intervals) > 0
                ):
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

    def _full_length_isoforms(self, isoform_ids, feature_lend, feature_rend):
        # isoforms whose both ends lie within the tolerance of the feature's ends
        max_dist = self.fsm_full_length_max_end_dist
        full_length = set()
        for isoform_id in isoform_ids:
            iso_lend, iso_rend = self.transcript_id_to_obj[isoform_id].get_coords()
            if (
                abs(feature_lend - iso_lend) <= max_dist
                and abs(feature_rend - iso_rend) <= max_dist
            ):
                full_length.add(isoform_id)
        return full_length

    def _containment_fields(self, fsm_isoforms, containing_isoforms):
        # the containing isoforms, their genes, and whether those genes are the FSM's
        fsm_genes = {
            self.transcript_id_to_obj[tid].get_gene_id() for tid in fsm_isoforms
        }
        containers = [self.transcript_id_to_obj[tid] for tid in containing_isoforms]
        in_fsm_gene = {obj.get_gene_id() in fsm_genes for obj in containers}
        if in_fsm_gene == {True}:
            scope = "same_gene"
        elif in_fsm_gene == {False}:
            scope = "other_gene"
        else:
            scope = "same_and_other_gene"
        return {
            "FSM_ISM_of_isoforms": ",".join(containing_isoforms),
            "FSM_ISM_of_genes": ",".join(
                obj.get_output_gene_id() for obj in containers
            ),
            "FSM_ISM_of_gene_scope": scope,
        }

    def _isoforms_containing_intron_chain(
        self, feature_chrom, feature_strand, transcribed_feature_obj
    ):
        # reference isoforms carrying every feature intron, consecutively and in order
        # (the ISM containment test). Includes isoforms with the identical chain.
        candidates = None
        for intron in transcribed_feature_obj.get_introns():
            intron_tok = make_intron_token(feature_chrom, feature_strand, intron)
            isoforms_with_intron = self.intron_to_isoforms.get(intron_tok)
            if not isoforms_with_intron:
                return set()
            if candidates is None:
                candidates = set(isoforms_with_intron)
            else:
                candidates &= isoforms_with_intron
            if not candidates:
                return set()
        if not candidates:
            return set()
        return self.restrict_splice_matched_isoforms(
            candidates, transcribed_feature_obj
        )

    def restrict_splice_matched_isoforms(
        self, isoforms_with_matching_introns, transcribed_feature_obj
    ):

        feature_introns = transcribed_feature_obj.get_introns()

        # print("Feature introns: {}".format(feature_introns))

        locally_matching_isoforms = set()

        for isoform_id in isoforms_with_matching_introns:

            isoform_introns = self.transcript_id_to_obj[isoform_id].get_introns()
            # print("Isoform {} has introns {}".format(isoform_id, isoform_introns))

            idx_start = isoform_introns.index(feature_introns[0])
            isoform_introns = isoform_introns[idx_start:]

            idx_end = isoform_introns.index(feature_introns[-1])
            isoform_introns = isoform_introns[: idx_end + 1]

            if isoform_introns == feature_introns:
                locally_matching_isoforms.add(isoform_id)

        return locally_matching_isoforms


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
