#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
import subprocess
import logging
import string
import time
import pysam
import LRAA_Globals
from hashlib import blake2s

from collections import defaultdict

logger = logging.getLogger(__name__)


def retrieve_contig_seq_from_fasta_file(contig_acc, fasta_filename):

    # samtools faidx <fasta> <region>
    attempts = 3
    contig_seq_str = None
    last_error = None

    for attempt in range(1, attempts + 1):
        try:
            contig_seq_str = subprocess.check_output(
                "samtools faidx {} {}".format(fasta_filename, contig_acc),
                shell=True,
                encoding="utf-8",
            )
            break
        except subprocess.CalledProcessError as exc:
            last_error = exc
            logger.warning(
                "samtools faidx failed for %s (attempt %d/%d): %s",
                contig_acc,
                attempt,
                attempts,
                exc,
            )
            if attempt < attempts:
                time.sleep(1 * attempt)

    if contig_seq_str is None:
        raise last_error

    contig_seq_str = contig_seq_str.upper()
    contig_seq_str = contig_seq_str.split("\n")
    contig_seq_str = contig_seq_str[1:]
    contig_seq_str = "".join(contig_seq_str)

    contig_seq_str = re.sub("\\s", "", contig_seq_str)  # just in case

    return contig_seq_str


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


def get_read_name_include_sc_encoding(pysam_read_alignment):

    cell_barcode_tag = LRAA_Globals.config["cell_barcode_tag"]
    read_umi_tag = LRAA_Globals.config["read_umi_tag"]

    read = pysam_read_alignment
    if read.has_tag(cell_barcode_tag) and read.has_tag(read_umi_tag):
        cell_barcode = read.get_tag(cell_barcode_tag)
        umi = read.get_tag(read_umi_tag)
        read_name = "^".join([cell_barcode, umi, read.query_name])
        return read_name
    else:
        return read.query_name


def frac_base_composition(nuc_seq, nuc_base):

    assert len(nuc_seq) > 0, "Error, nuc_seq is empty"

    counter = 0

    for char in nuc_seq:
        if char.upper() == nuc_base.upper():
            counter += 1

    frac_base = counter / len(nuc_seq)

    return frac_base


def get_hash_code(input_string):
    hash_object = blake2s(digest_size=11)
    hash_object.update(input_string.encode("utf-8"))
    hex_digest = hash_object.hexdigest()
    hex_digest = str(hex_digest)

    return hex_digest
