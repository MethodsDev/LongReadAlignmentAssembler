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


def file_identity_token(path):
    """Short digest of which file this is and whether it has changed since.

    The resolved path separates two inputs that happen to share a basename; size
    and modification time catch one replaced in place. Contents are deliberately
    not hashed: these bams run to gigabytes, and reading one on every startup
    would cost more than the step being cached, to cover a case the stat pair
    already covers.
    """
    resolved = os.path.realpath(path)
    try:
        stat_result = os.stat(resolved)
        identity = "{}:{}:{}".format(
            resolved, stat_result.st_size, stat_result.st_mtime_ns
        )
    except OSError:
        # A missing input is the caller's problem to report; naming it by path
        # alone keeps this from raising first and obscuring that.
        identity = resolved

    return get_hash_code(identity)[:12]


def splice_graph_norm_cache_stem(base_root, normalize_max_cov_level, source_bam):
    """Name for a normalized bam and its work directory.

    Everything that determines the contents belongs in the name, because a hit is
    taken on trust and nothing downstream can audit it.

    The method matters most: reads selected by an algorithm no longer in the tree
    are undetectable once written, since a bam from the read-start-binning era
    carries no weight tag and an untagged read correctly weighs 1.

    The source is identified by content fingerprint rather than by name. base_root
    is only a basename, so without this a second input called sample.bam from a
    different directory, or the same path regenerated, would land on the cache
    built for the first.

    The target depth is here for the same reason, and it also scopes the work
    directory, whose checkpoints would otherwise let a run at one depth reuse the
    strand split and merge performed for another.
    """
    return "{}.norm_{}.{}.{}".format(
        base_root,
        normalize_max_cov_level,
        LRAA_Globals.SPLICE_GRAPH_NORMALIZATION_METHOD,
        file_identity_token(source_bam),
    )
