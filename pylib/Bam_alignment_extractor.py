#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
import time
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

        # Collect either raw reads or Pretty_alignment objects depending on 'pretty'
        read_alignments = [] if not pretty else None
        pretty_alignments = [] if pretty else None

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

        # diagnostics
        def _mem_usage_mb():
            try:
                import psutil  # type: ignore
                return psutil.Process(os.getpid()).memory_info().rss / (1024.0 * 1024.0)
            except Exception:
                try:
                    import resource  # type: ignore
                    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                    if rss > 1e10:
                        return rss / (1024.0 * 1024.0)
                    else:
                        return rss / 1024.0
                except Exception:
                    return None

        LOG_EVERY_N = int(os.environ.get("LRAA_LOG_ALIGN_EVERY_N", "10000"))
        last_log_t = time.time()
        LOG_EVERY_SEC = float(os.environ.get("LRAA_LOG_ALIGN_EVERY_SEC", "10"))
        processed = 0
        kept_so_far = 0

        for read in read_fetcher:
            processed += 1

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

            if pretty:
                # Build Pretty_alignment on the fly to avoid storing all pysam AlignedSegment objects simultaneously
                pretty_alignments.append(Pretty_alignment.get_pretty_alignment(read))
            else:
                read_alignments.append(read)

            kept_so_far += 1

            if (
                (processed % LOG_EVERY_N == 0)
                or (time.time() - last_log_t) >= LOG_EVERY_SEC
            ):
                m = _mem_usage_mb()
                discards = dict(discarded_read_counter)
                try:
                    logger.info(
                        f"progress get_read_alignments: processed={processed}, kept={kept_so_far}, discards={discards}, rss={(m and f'{m:.1f} MB') or '<na>'}"
                    )
                except Exception:
                    logger.info(
                        f"progress get_read_alignments: processed={processed}, kept={kept_so_far}"
                    )
                last_log_t = time.time()

        kept_count = len(pretty_alignments) if pretty else len(read_alignments)
        final_mem = _mem_usage_mb()
        logger.info(
            "reads kept for {} {}: {} and discarded: {} (rss: {})".format(
                contig_acc,
                contig_strand,
                kept_count,
                dict(discarded_read_counter),
                f"{final_mem:.1f} MB" if final_mem is not None else "<na>",
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
            return pretty_alignments
        else:
            return read_alignments
