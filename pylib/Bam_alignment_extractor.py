#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
import subprocess
import logging
import string
import pysam
from collections import defaultdict
from Pretty_alignment import Pretty_alignment
import LRAA_Globals
import Util_funcs

logger = logging.getLogger(__name__)


class Bam_alignment_extractor:

    # ---------------
    # class variables
    # ---------------

    def __init__(self, alignments_bam_filename):

        self._alignments_bam_filename = alignments_bam_filename

        self._pysam_reader = pysam.AlignmentFile(self._alignments_bam_filename, "rb")

        return

    def set_read_aln_gap_merge(self, read_aln_gap_merge_int):

        self._read_aln_gap_merge_int = read_aln_gap_merge_int

        return

    def get_read_alignments(
        self,
        contig_acc,
        contig_strand=None,
        region_lend=None,
        region_rend=None,
        pretty=False,
        per_id_QC_raise_error=False,
        config=LRAA_Globals.config,
    ):

        discarded_read_counter = defaultdict(int)

        read_alignments = list()

        MIN_MAPPING_QUALITY = int(LRAA_Globals.config["min_mapping_quality"])

        # parse read alignments, capture introns and genome coverage info.
        read_fetcher = None
        if region_lend is not None and region_rend is not None:
            if contig_strand is not None:
                logger.debug(
                    "Fetching alignments for {}{}:{}-{}".format(
                        contig_acc, contig_strand, region_lend, region_rend
                    )
                )
            else:
                logger.debug(
                    "Fetching alignments for {}:{}-{}".format(
                        contig_acc, region_lend, region_rend
                    )
                )

            read_fetcher = self._pysam_reader.fetch(
                contig_acc, region_lend, region_rend
            )
        else:
            logger.debug(
                "Fetching all alignments for contig: {} strand {}".format(
                    contig_acc, contig_strand
                )
            )
            read_fetcher = self._pysam_reader.fetch(contig_acc)

        num_alignments_per_id_ok = 0
        num_alignments_per_id_fail = 0

        for read in read_fetcher:

            if contig_strand is not None:
                if read.is_forward and contig_strand != "+":
                    continue
                if read.is_reverse and contig_strand != "-":
                    continue

            if read.mapping_quality < MIN_MAPPING_QUALITY:
                discarded_read_counter["min_mapping_quality"] += 1
                continue

            if read.is_paired and not read.is_proper_pair:
                discarded_read_counter["improper_pair"] += 1
                continue

            if read.is_duplicate:
                discarded_read_counter["duplicate"] += 1
                continue

            if read.is_qcfail:
                discarded_read_counter["qcfail"] += 1
                continue

            if read.is_secondary:
                discarded_read_counter["secondary"] += 1
                continue

            # determine min per_id based on read type:
            min_per_id = LRAA_Globals.config["min_per_id"]

            # check read alignment percent identity
            cigar_stats = read.get_cigar_stats()
            aligned_base_count = cigar_stats[0][0]
            if aligned_base_count == 0:
                aligned_base_count = cigar_stats[0][7] + cigar_stats[0][8]

            mismatch_count = None
            if read.has_tag("NM"):
                mismatch_count = int(read.get_tag("NM"))
            elif read.has_tag("nM"):
                mismatch_count = int(read.get_tag("nM"))
            if mismatch_count is not None:
                per_id = 100 - (mismatch_count / aligned_base_count) * 100
                # logger.info(f"-read per_id: {per_id}")
                if per_id < min_per_id:
                    read_name = Util_funcs.get_read_name_include_sc_encoding(read)
                    logger.debug(
                        "read {} has insufficient per_id {}, < min {} required ".format(
                            read_name, per_id, min_per_id
                        )
                    )
                    discarded_read_counter["low_perID"] += 1
                    # print(read)
                    # print("Cigar_stats: " + str(cigar_stats))
                    num_alignments_per_id_fail += 1
                    continue
                else:
                    num_alignments_per_id_ok += 1

            read_alignments.append(read)

        logger.info(
            "reads kept: {} and discarded: {}".format(
                len(read_alignments), discarded_read_counter
            )
        )

        if (
            num_alignments_per_id_fail + num_alignments_per_id_ok
            >= LRAA_Globals.config["min_total_alignments_engage_frac_per_id_check"]
        ):
            frac_alignments_fail_per_id_check = num_alignments_per_id_fail / (
                num_alignments_per_id_fail + num_alignments_per_id_ok
            )

            if (
                frac_alignments_fail_per_id_check
                < LRAA_Globals.config["min_frac_alignments_pass_per_id_check"]
            ):
                # raise RuntimeError(f"Error, would appear only {frac_alignments_fail_per_id_check} on {contig_acc} have at least {min_per_id} percent identity. Please reevaluate your --min_per_id setting for application of LRAA with these alignments")
                logger.debug(
                    f"Error, would appear only {frac_alignments_fail_per_id_check} on {contig_acc} have at least {min_per_id} percent identity. Please reevaluate your --min_per_id setting for application of LRAA with these alignments"
                )

        if pretty:
            return Pretty_alignment.get_pretty_alignments(read_alignments)
        else:
            return read_alignments
