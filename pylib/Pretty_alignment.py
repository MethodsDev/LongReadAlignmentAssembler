#!/usr/bin/env pythonOA
# encoding: utf-8

import sys, os, re
import subprocess
import logging
import string
import pysam
import LRAA_Globals
import Util_funcs
from collections import defaultdict


logger = logging.getLogger(__name__)


class Pretty_alignment:

    read_aln_gap_merge_int = LRAA_Globals.config["read_aln_gap_merge_int"]
    min_terminal_splice_exon_anchor_length = LRAA_Globals.config[
        "min_terminal_splice_exon_anchor_length"
    ]

    def __init__(self, pysam_alignment, pretty_alignment_segments):

        self._pysam_alignment = pysam_alignment
        self._pretty_alignment_segments = pretty_alignment_segments  # defined in Bam_alignment_extractor //TODO: move logic here.

        # if pysam_alignment.has_tag("RG") and pysam_alignment.get_tag("RG") == "PBLR":
        #    self._read_type = "PBLR"
        # else:
        #    self._read_type = "ILMN"
        self._read_type = "PacBio"

        ## soft clipping before/after polyA trimming.
        # - pre-trimming
        self.orig_left_soft_clipping = None
        self.orig_right_soft_clipping = None
        self.orig_left_soft_clipped_seq = None
        self.orig_right_soft_clipped_seq = None
        # - post-trimming
        self.left_soft_clipping = None
        self.right_soft_clipping = None
        self.left_soft_clipped_seq = None
        self.right_soft_clipped_seq = None

        self._set_read_soft_clipping_info(pysam_alignment)  # sets above

    def __repr__(self):
        return str(self._pretty_alignment_segments)

    def get_read_name(self):
        read = self._pysam_alignment
        return Util_funcs.get_read_name_include_sc_encoding(read)

    def get_pysam_alignment(self):
        return self._pysam_alignment

    def get_pretty_alignment_segments(self):
        return self._pretty_alignment_segments

    def set_pretty_alignment_segments(self, alignment_segments):
        self._pretty_alignment_segments = alignment_segments
        return

    def get_strand(self):
        if self._pysam_alignment.is_forward:
            return "+"
        elif self._pysam_alignment.is_reverse:
            return "-"
        else:
            return "?"

    def get_introns(self):
        intron_coordsets = list()
        exon_segments = self.get_pretty_alignment_segments()
        if len(exon_segments) > 1:
            exon_segments = sorted(exon_segments, key=lambda x: x[0])
            for i in range(1, len(exon_segments)):
                intron_lend = exon_segments[i - 1][1] + 1
                intron_rend = exon_segments[i][0] - 1
                assert intron_lend < intron_rend
                intron_coordsets.append((intron_lend, intron_rend))
        return intron_coordsets

    def get_alignment_span(self):
        lend = self._pretty_alignment_segments[0][0]
        rend = self._pretty_alignment_segments[-1][1]

        return (lend, rend)

    def get_read_type(self):
        return self._read_type

    def _set_read_soft_clipping_info(self, pysam_alignment=None):

        if pysam_alignment is None:
            pysam_alignment = self._pysam_alignment

        cigar_tuples = pysam_alignment.cigartuples

        S = 4  # soft clipping cigar code in pysam

        left_soft_clipping = cigar_tuples[0][1] if cigar_tuples[0][0] == S else 0
        self.orig_left_soft_clipping = left_soft_clipping

        right_soft_clipping = cigar_tuples[-1][1] if cigar_tuples[-1][0] == S else 0
        self.orig_right_soft_clipping = right_soft_clipping

        read_sequence = pysam_alignment.query_sequence
        left_soft_clipped_seq = ""
        if left_soft_clipping > 0:
            left_soft_clipped_seq = read_sequence[0:left_soft_clipping]

        self.orig_left_soft_clipped_seq = left_soft_clipped_seq

        right_soft_clipped_seq = ""
        if right_soft_clipping > 0:
            right_soft_clipped_seq = read_sequence[(-1 * right_soft_clipping) :]

        self.orig_right_soft_clipped_seq = right_soft_clipped_seq

        ## deal with polyA
        min_PolyA_ident_length = LRAA_Globals.config["min_PolyA_ident_length"]

        read_name = Util_funcs.get_read_name_include_sc_encoding(pysam_alignment)

        if pysam_alignment.is_forward and right_soft_clipping >= min_PolyA_ident_length:
            polyA_regex = "A" * min_PolyA_ident_length + "+$"
            right_soft_clippled_seq = re.sub(polyA_regex, "", right_soft_clipped_seq)
            right_soft_clipping = len(right_soft_clipped_seq)
            logger.debug("Stripped polyA from end of read {}".format(read_name))

        elif (
            pysam_alignment.is_reverse and left_soft_clipping >= min_PolyA_ident_length
        ):
            polyT_regex = "^" + "T" * min_PolyA_ident_length + "+"
            left_soft_clipped_seq = re.sub(polyT_regex, "", left_soft_clipped_seq)
            left_soft_clipping = len(left_soft_clipped_seq)
            logger.debug("Stripped polyT from beginning of read {}".format(read_name))

        # set obj vars
        self.left_soft_clipping = left_soft_clipping
        self.right_soft_clipping = right_soft_clipping
        self.left_soft_clipped_seq = left_soft_clipped_seq
        self.right_soft_clipped_seq = right_soft_clipped_seq

    def has_soft_clipping(self):

        assert self.left_soft_clipping is not None
        assert self.right_soft_clipping is not None

        if self.left_soft_clipping > 0 or self.right_soft_clipping > 0:
            return True
        else:
            return False

    @classmethod
    def prune_long_terminal_introns(cls, pretty_alignments, splice_graph):

        logger.info("Removing long terminal introns from pretty alignments.")

        for pretty_alignment in pretty_alignments:

            alignment_segments = pretty_alignment.get_pretty_alignment_segments()
            intron_segments = pretty_alignment.get_introns()

            trimmed = False

            # examine left-most intron:
            if len(intron_segments) > 0:
                leftmost_intron = intron_segments.pop(0)
                leftmost_intron_len = leftmost_intron[1] - leftmost_intron[0] + 1
                if leftmost_intron_len > LRAA_Globals.config["max_intron_length"]:
                    logger.info(
                        "-pruning long intron {} from left of pretty alignment {}".format(
                            leftmost_intron, pretty_alignment
                        )
                    )
                    alignment_segments.pop(0)
                    trimmed = True

            # check right-most intron
            if len(intron_segments) > 0:
                rightmost_intron = intron_segments.pop()
                rightmost_intron_len = rightmost_intron[1] - rightmost_intron[0] + 1
                if rightmost_intron_len > LRAA_Globals.config["max_intron_length"]:
                    logger.info(
                        "-pruning long intron {} from left of pretty alignment {}".format(
                            rightmost_intron, pretty_alignment
                        )
                    )
                    alignment_segments.pop()
                    trimmed = True

            if trimmed is True:
                pretty_alignment.set_pretty_alignment_segments(alignment_segments)

        return

    @classmethod
    def try_correct_alignments(cls, pretty_alignments_list, splice_graph, contig_seq):

        logger.info("Attempting to correct alignments at soft-clips")

        max_softclip_realign_test = LRAA_Globals.config["max_softclip_realign_test"]
        min_softclip_realign_test = LRAA_Globals.config["min_softclip_realign_test"]

        ##################
        # - more extreme verbose setting for debugging this method
        local_debug = False
        ##################

        for pretty_alignment in pretty_alignments_list:

            if not pretty_alignment.has_soft_clipping():
                continue

            alignment_segments = pretty_alignment.get_pretty_alignment_segments()

            left_soft_clipping, right_soft_clipping = (
                pretty_alignment.left_soft_clipping,
                pretty_alignment.right_soft_clipping,
            )
            orig_left_soft_clipping, orig_right_soft_clipping = (
                pretty_alignment.orig_left_soft_clipping,
                pretty_alignment.orig_right_soft_clipping,
            )

            read = pretty_alignment._pysam_alignment

            read_sequence = read.query_sequence
            read_name = read.query_name

            read_adj_lend = orig_left_soft_clipping - left_soft_clipping + 1
            assert read_adj_lend >= 1
            read_adj_rend = (
                len(read_sequence) - orig_right_soft_clipping + right_soft_clipping
            )
            assert read_adj_rend <= len(read_sequence)

            # get mapping of genome pos -> read pos
            aligned_pairs = dict(
                [(y + 1, x + 1) for x, y in read.get_aligned_pairs(matches_only=True)]
            )

            ## examine left-side intron realignment

            # exons   <------------>           <--------------->
            #  read                        XXXXX||||||||||||||||
            # raligned         |||||            ||||||||||||||||

            if (
                left_soft_clipping > 0
                and left_soft_clipping <= max_softclip_realign_test
                and left_soft_clipping >= min_softclip_realign_test
            ):

                if local_debug:
                    print(
                        "evaluating left alignment correction for: {}".format(
                            pretty_alignment
                        )
                    )

                left_alignment_segment = alignment_segments[0]
                exon_seg_lend, exon_seg_rend = left_alignment_segment
                overlapping_introns = list()
                for overlapping_intron in splice_graph.get_overlapping_introns(
                    exon_seg_lend - 1, exon_seg_rend
                ):
                    intron_lend, intron_rend = overlapping_intron.get_coords()
                    if intron_rend + 1 >= exon_seg_lend and intron_rend < exon_seg_rend:
                        overlapping_introns.append(overlapping_intron)

                if local_debug:
                    print(
                        "Got left overlapping introns: {}".format(overlapping_introns)
                    )

                for left_overlapping_intron in overlapping_introns:
                    intron_lend, intron_rend = left_overlapping_intron.get_coords()

                    intron_adjacent_pos = intron_rend + 1
                    if intron_adjacent_pos not in aligned_pairs:
                        continue

                    read_rend = aligned_pairs[intron_adjacent_pos] - 1
                    if read_rend - read_adj_lend + 1 <= max_softclip_realign_test:

                        left_read_seq = read_sequence[read_adj_lend - 1 : read_rend]

                        if local_debug:
                            print(
                                "Checking read realignment for {}".format(left_read_seq)
                            )

                        genomic_rend = intron_lend - 1
                        genomic_lend = genomic_rend - len(left_read_seq) + 1

                        genomic_substr = contig_seq[genomic_lend - 1 : genomic_rend]

                        if local_debug:
                            print(
                                "Comparing to genomic rend seq: {}".format(
                                    genomic_substr
                                )
                            )

                        if left_read_seq.upper() == genomic_substr.upper():

                            if local_debug:
                                print("\tLeft MATCH FOUND")
                            logger.debug(
                                "-left-side alignment correction for read {} repositioning sequence {}".format(
                                    read_name, genomic_substr
                                )
                            )
                            # do reassignment:
                            alignment_segments[0][0] = intron_rend + 1
                            alignment_segments.insert(0, [genomic_lend, genomic_rend])
                            break

            # examine right-side intron realignment

            # exons    <------------>           <--------------->
            #  read    ||||||||||||||XXXXX
            # raligned ||||||||||||||           |||||

            if (
                right_soft_clipping > 0
                and right_soft_clipping <= max_softclip_realign_test
                and right_soft_clipping >= min_softclip_realign_test
            ):

                if local_debug:
                    print(
                        "evaluating right alignment correction for: {}".format(
                            pretty_alignment
                        )
                    )

                right_alignment_segment = alignment_segments[-1]
                exon_seg_lend, exon_seg_rend = right_alignment_segment
                overlapping_introns = list()
                for overlapping_intron in splice_graph.get_overlapping_introns(
                    exon_seg_lend, exon_seg_rend + 1
                ):
                    intron_lend, intron_rend = overlapping_intron.get_coords()
                    if intron_lend > exon_seg_lend and intron_lend <= exon_seg_rend + 1:
                        overlapping_introns.append(overlapping_intron)

                if local_debug:
                    print(
                        "Got right overlapping introns: {}".format(overlapping_introns)
                    )

                for right_overlapping_intron in overlapping_introns:
                    intron_lend, intron_rend = right_overlapping_intron.get_coords()
                    intron_adjacent_pos = intron_lend - 1
                    if intron_adjacent_pos not in aligned_pairs:
                        continue
                    read_lend = aligned_pairs[intron_adjacent_pos]
                    if len(read_sequence) - read_lend + 1 <= max_softclip_realign_test:
                        right_read_seq = read_sequence[read_lend:read_adj_rend]
                        if local_debug:
                            print(
                                "Checking read realignment for {}".format(
                                    right_read_seq
                                )
                            )
                        genomic_lend = intron_rend + 1
                        genomic_rend = genomic_lend + len(right_read_seq) - 1
                        genomic_substr = contig_seq[genomic_lend - 1 : genomic_rend]
                        if local_debug:
                            print(
                                "Comparing to genomic rend seq: {}".format(
                                    genomic_substr
                                )
                            )
                        if right_read_seq.upper() == genomic_substr.upper():
                            if local_debug:
                                print("\tRight MATCH FOUND")
                            logger.debug(
                                "-right-side alignment correction for read {} repositioning sequence {}".format(
                                    read_name, genomic_substr
                                )
                            )
                            # do reassignment:
                            alignment_segments[-1][1] = intron_lend - 1
                            alignment_segments.append([genomic_lend, genomic_rend])
                            break

        return

    @classmethod
    def get_pretty_alignment(cls, pysam_read_alignment):

        alignment_segments = cls.read_to_pretty_alignment_segments(pysam_read_alignment)
        this_pretty_alignment = Pretty_alignment(
            pysam_read_alignment, alignment_segments
        )

        return this_pretty_alignment

    @classmethod
    def get_pretty_alignments(cls, read_alignments_list):
        """
        stores the pysam alignment record along with inferred transcript exon segments,
        where exon segments involve joining nearby alignment segments separated by short indels
        """

        pretty_alignments = [cls.get_pretty_alignment(x) for x in read_alignments_list]

        return pretty_alignments

    @classmethod
    def read_to_pretty_alignment_segments(cls, pysam_read_alignment):

        read_name = pysam_read_alignment.query_name

        aligned_pairs = cls.get_genome_alignment_blocks(pysam_read_alignment)

        # print(aligned_pairs)

        ## merge adjacent blocks within range.
        alignment_segments = list()
        alignment_segments.append(list(aligned_pairs.pop(0)))
        # block coordinates are zero-based, left inclusive, and right-end exclusive
        alignment_segments[0][
            0
        ] += (
            1  # adjust for zero-based.  note, end position doesn't need to be adjusted.
        )

        for aligned_pair in aligned_pairs:
            aligned_pair = list(aligned_pair)
            aligned_pair[0] += 1

            # extend earlier stored segment or append new one
            delta = aligned_pair[0] - alignment_segments[-1][1]
            # logger.debug("comparing {} to {}, delta: {}".format(alignment_segments[-1], aligned_pair, delta))
            if delta < Pretty_alignment.read_aln_gap_merge_int:
                # extend rather than append
                alignment_segments[-1][1] = aligned_pair[1]
            else:
                # append, as too far apart from prev
                alignment_segments.append(list(aligned_pair))

        ##//TODO: make below trimming of short terminal alignment segments an option and configurable
        """

        # trim short terminal segments from each end
        while (len(alignment_segments) > 1 and
            alignment_segments[0][1] - alignment_segments[0][0] + 1 < self._min_terminal_splice_exon_anchor_length):

            alignment_segments.pop(0)

        while (len(alignment_segments) > 1 and
            alignment_segments[len(alignment_segments)-1][1] - alignment_segments[len(alignment_segments)-1][0] + 1 < self._min_terminal_splice_exon_anchor_length):

            alignment_segments.pop()
        """

        logger.debug(
            "read {} pretty alignment segments: {}".format(
                read_name, alignment_segments
            )
        )

        return alignment_segments

    @classmethod
    def get_genome_alignment_blocks(self, read):

        cigartuples = read.cigartuples

        """
        M       BAM_CMATCH      0
        I       BAM_CINS        1
        D       BAM_CDEL        2
        N       BAM_CREF_SKIP   3
        S       BAM_CSOFT_CLIP  4
        H       BAM_CHARD_CLIP  5
        P       BAM_CPAD        6
        =       BAM_CEQUAL      7
        X       BAM_CDIFF       8
        B       BAM_CBACK       9
        """

        read_name = read.query_name

        ref_start = read.reference_start
        read_start = 0

        prev_ref_start = ref_start
        prev_read_start = read_start

        genome_segments = []

        for cigartuple in cigartuples:
            code, val = cigartuple

            token = None

            if code in (0, 7, 8):
                token = "BAM_CMATCH"
                ref_start += val
                read_start += val
                prev_ref_start += 1
                prev_read_start += 1

            elif code == 1:
                token = "BAM_CINS"
                read_start += val
                prev_read_start += 1

            elif code == 2:
                token = "BAM_CDEL"
                ref_start += val
                prev_ref_start += 1

            elif code == 3:
                token = "BAM_CREF_SKIP"
                ref_start += val
                prev_ref_start += 1

            elif code == 4:
                token = "BAM_CSOFT_CLIP"
                read_start += val
                prev_read_start += 1

            elif code == 5:
                token = "BAM_CHARD_CLIP"
                read_start += val
                prev_read_start += 1

            else:
                raise RuntimeError("Not sure how to handle code {}".format(code))

            if token in ["BAM_CMATCH", "BAM_CDEL"]:
                genome_segments.append(
                    [prev_ref_start - 1, ref_start]
                )  # make zero-based left-inclusive right-exclusive for consistency w/ pysam blocks

            prev_read_start = read_start
            prev_ref_start = ref_start

        logger.debug(
            "genome segments from read: {}: {}".format(read_name, genome_segments)
        )

        return genome_segments
