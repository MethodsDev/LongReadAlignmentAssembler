#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
import subprocess
import logging
import string

from collections import defaultdict

logger = logging.getLogger(__name__)


def retrieve_contig_seq_from_fasta_file(fasta_filename, contig_acc):

    contig_seq_str = subprocess.check_output("samtools faidx {} {}".format(contig_acc, fasta_filename),
                                             shell=True,
                                             encoding="utf-8")
    
    contig_seq_str = contig_seq_str.upper()
    contig_seq_str = contig_seq_str.split("\n")
    contig_seq_str = contig_seq_str[1:]
    contig_seq_str = "".join(contig_seq_str)

    contig_seq_str = re.sub("\s", "", contig_seq_str) # just in case

    return(contig_seq_str)


def coordpairs_overlap(coordset_A, coordset_B):

    A_lend, A_rend = coordset_A

    B_lend, B_rend = coordset_B

    if A_lend <= B_rend and A_rend >= B_lend:
        return True
    else:
        return False

    
def get_num_overlapping_bases(coordset_A, coordset_B):

    if not coordpairs_overlap(coordset_A, coordset_B):
        return 0

    coords = sorted([coordset_A[0], coordset_A[1], coordset_B[0], coordset_B[1]])
    overlap_len = coords[2] - coords[1] + 1

    return overlap_len

                    
