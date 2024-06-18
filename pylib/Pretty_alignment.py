#!/usr/bin/env pythonOA
# encoding: utf-8

import sys, os, re
import subprocess
import logging
import string
import pysam
import PASA_SALRAA_Globals
from collections import defaultdict


logger = logging.getLogger(__name__)


class Pretty_alignment:

    def __init__(self, pysam_alignment, pretty_alignment_segments):

        self._pysam_alignment = pysam_alignment
        self._pretty_alignment_segments = pretty_alignment_segments # defined in Bam_alignment_extractor //TODO: move logic here.

        #if pysam_alignment.has_tag("RG") and pysam_alignment.get_tag("RG") == "PBLR":
        #    self._read_type = "PBLR"
        #else:
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
        
        self._set_read_soft_clipping_info(pysam_alignment) # sets above

        
        
    def __repr__(self):
        return str(self._pretty_alignment_segments)

    def get_read_name(self):
        return self._pysam_alignment.query_name
    
        
    def get_pysam_alignment(self):
        return self._pysam_alignment

    def get_pretty_alignment_segments(self):
        return self._pretty_alignment_segments


    def get_strand(self):
        if self._pysam_alignment.is_forward:
            return('+')
        elif self._pysam_alignment.is_reverse:
            return('-')
        else:
            return("?")

    
    def get_introns(self):
        intron_coordsets = list()
        exon_segments = self.get_pretty_alignment_segments()
        if len(exon_segments) > 1:
            exon_segments = sorted(exon_segments, key=lambda x: x[0])
            for i in range(1, len(exon_segments)):
                intron_lend = exon_segments[i-1][1]+1
                intron_rend = exon_segments[i][0] - 1
                assert intron_lend < intron_rend
                intron_coordsets.append( (intron_lend, intron_rend) )
        return intron_coordsets

    

    def get_alignment_span(self):
        lend = self._pretty_alignment_segments[0][0]
        rend = self._pretty_alignment_segments[-1][1]

        return(lend, rend)
    

    def get_read_type(self):
        return self._read_type

    
    def _set_read_soft_clipping_info(self, pysam_alignment=None):

        if pysam_alignment is None:
            pysam_alignment = self._pysam_alignment
        
        cigar_tuples = pysam_alignment.cigartuples

        S=4 # soft clipping cigar code in pysam
        
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
            right_soft_clipped_seq = read_sequence[ (-1*right_soft_clipping) : ]

        self.orig_right_soft_clipped_seq = right_soft_clipped_seq
        
        ## deal with polyA
        min_PolyA_ident_length = PASA_SALRAA_Globals.config['min_PolyA_ident_length']

        if pysam_alignment.is_forward and right_soft_clipping >= min_PolyA_ident_length :
            polyA_regex = "A" * min_PolyA_ident_length + "+$"
            right_soft_clippled_seq = re.sub(polyA_regex, "", right_soft_clipped_seq)
            right_soft_clipping = len(right_soft_clipped_seq)
            logger.debug("Stripped polyA from end of read {}".format(pysam_alignment.query_name))

        elif pysam_alignment.is_reverse and left_soft_clipping >= min_PolyA_ident_length:
            polyT_regex = "^" + "T"*min_PolyA_ident_length + "+"
            left_soft_clipped_seq = re.sub(polyT_regex, "", left_soft_clipped_seq)
            left_soft_clipping = len(left_soft_clipped_seq)
            logger.debug("Stripped polyT from beginning of read {}".format(pysam_alignment.query_name))
        

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
    def try_correct_alignments(cls, pretty_alignments_list, splice_graph, contig_seq):

        logger.info("Attempting to correct alignments at soft-clips")
        
        max_softclip_realign_test = PASA_SALRAA_Globals.config['max_softclip_realign_test']

        ##################
        # - more extreme verbose setting for debugging this method
        local_debug = False
        ##################
        
        for pretty_alignment in pretty_alignments_list:

            if not pretty_alignment.has_soft_clipping():
                continue
            
            alignment_segments = pretty_alignment.get_pretty_alignment_segments()

            left_soft_clipping, right_soft_clipping = pretty_alignment.left_soft_clipping, pretty_alignment.right_soft_clipping
            orig_left_soft_clipping, orig_right_soft_clipping = pretty_alignment.orig_left_soft_clipping, pretty_alignment.orig_right_soft_clipping

            
            read = pretty_alignment._pysam_alignment
            
            read_sequence = read.query_sequence
            read_name = read.query_name

            read_adj_lend = orig_left_soft_clipping - left_soft_clipping + 1
            assert read_adj_lend >= 1
            read_adj_rend = len(read_sequence) - orig_right_soft_clipping + right_soft_clipping
            assert read_adj_rend <= len(read_sequence)
            
            # get mapping of genome pos -> read pos
            aligned_pairs = dict([ (y+1, x+1) for x,y in read.get_aligned_pairs(matches_only=True) ])



            ## examine left-side intron realignment

            # exons   <------------>           <---------------> 
            #  read                        XXXXX||||||||||||||||
            # raligned         |||||            ||||||||||||||||
            
            if left_soft_clipping > 0 and left_soft_clipping <= max_softclip_realign_test:


                if local_debug:
                    print("evaluating left alignment correction for: {}".format(pretty_alignment))

                left_alignment_segment = alignment_segments[0]
                exon_seg_lend, exon_seg_rend = left_alignment_segment
                overlapping_introns = list()
                for overlapping_intron in splice_graph.get_overlapping_introns(exon_seg_lend-1, exon_seg_rend):
                    intron_lend, intron_rend = overlapping_intron.get_coords()
                    if intron_rend +1 >= exon_seg_lend and intron_rend < exon_seg_rend:
                        overlapping_introns.append(overlapping_intron)

                if local_debug:
                    print("Got left overlapping introns: {}".format(overlapping_introns))

                for left_overlapping_intron in overlapping_introns:
                    intron_lend, intron_rend = left_overlapping_intron.get_coords()

                    intron_adjacent_pos = intron_rend + 1
                    if intron_adjacent_pos not in aligned_pairs:
                        continue
                    
                    read_rend = aligned_pairs[intron_adjacent_pos] - 1
                    if read_rend - read_adj_lend + 1 <=  max_softclip_realign_test:
                        
                        left_read_seq = read_sequence[read_adj_lend -1 : read_rend]

                        if local_debug:
                            print("Checking read realignment for {}".format(left_read_seq))
                        
                        genomic_rend = intron_lend - 1
                        genomic_lend = genomic_rend - len(left_read_seq) + 1
                        
                        genomic_substr = contig_seq[genomic_lend - 1 : genomic_rend]
                        
                        if local_debug:
                            print("Comparing to genomic rend seq: {}".format(genomic_substr))
                        
                        if left_read_seq.upper() == genomic_substr.upper():

                            if local_debug:
                                print("\tLeft MATCH FOUND")
                            logger.debug("-left-side alignment correction for read {} repositioning sequence {}".format(read_name, genomic_substr))
                            # do reassignment:
                            alignment_segments[0][0] = intron_rend + 1
                            alignment_segments.insert(0, [genomic_lend, genomic_rend])
                            break

            

            # examine right-side intron realignment

            # exons    <------------>           <---------------> 
            #  read    ||||||||||||||XXXXX
            # raligned ||||||||||||||           |||||

                        
            if right_soft_clipping > 0 and right_soft_clipping <= max_softclip_realign_test:

                if local_debug:
                    print("evaluating right alignment correction for: {}".format(pretty_alignment))


                right_alignment_segment = alignment_segments[-1]
                exon_seg_lend, exon_seg_rend = right_alignment_segment
                overlapping_introns = list()
                for overlapping_intron in splice_graph.get_overlapping_introns(exon_seg_lend, exon_seg_rend+1):
                    intron_lend, intron_rend = overlapping_intron.get_coords()
                    if intron_lend > exon_seg_lend and intron_lend <= exon_seg_rend+1:
                        overlapping_introns.append(overlapping_intron)

                if local_debug:
                    print("Got right overlapping introns: {}".format(overlapping_introns))

                for right_overlapping_intron in overlapping_introns:
                    intron_lend, intron_rend = right_overlapping_intron.get_coords()
                    intron_adjacent_pos = intron_lend - 1
                    if intron_adjacent_pos not in aligned_pairs:
                        continue
                    read_lend = aligned_pairs[intron_adjacent_pos]
                    if len(read_sequence) - read_lend + 1 <= max_softclip_realign_test:
                        right_read_seq = read_sequence[read_lend:read_adj_rend]
                        if local_debug:
                            print("Checking read realignment for {}".format(right_read_seq))
                        genomic_lend = intron_rend + 1
                        genomic_rend = genomic_lend + len(right_read_seq) - 1
                        genomic_substr = contig_seq[genomic_lend - 1 : genomic_rend]
                        if local_debug:
                            print("Comparing to genomic rend seq: {}".format(genomic_substr))
                        if right_read_seq.upper() == genomic_substr.upper():
                            if local_debug:
                                print("\tRight MATCH FOUND")
                            logger.debug("-right-side alignment correction for read {} repositioning sequence {}".format(read_name, genomic_substr))
                            # do reassignment:
                            alignment_segments[-1][1] = intron_lend -1
                            alignment_segments.append([genomic_lend, genomic_rend])
                            break

        return

    
        
