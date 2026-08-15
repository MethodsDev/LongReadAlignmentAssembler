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
        self._last_discarded_read_names_by_reason = defaultdict(set)

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
        force_lighten_all=False,
    ):

        discarded_read_counter = defaultdict(int)
        self._last_discarded_read_names_by_reason = defaultdict(set)

        # Collect either raw reads or Pretty_alignment objects depending on 'pretty'
        read_alignments = [] if not pretty else None
        pretty_alignments = [] if pretty else None

        MIN_MAPPING_QUALITY = int(LRAA_Globals.config["min_mapping_quality"])
        MAX_INTRON_LENGTH = int(LRAA_Globals.config["max_intron_length"])

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

            # region_lend/region_rend arrive 1-based inclusive, straight from
            # --region contig:lend-rend.  pysam.fetch takes 0-based half-open
            # [start, stop), so passing region_lend unconverted started the
            # window one base late and silently dropped any alignment whose only
            # overlap with the region was its first base.  The right edge needs no
            # adjustment: 1-based region_rend is 0-based region_rend - 1, which is
            # below the exclusive stop and therefore included.
            read_fetcher = self._pysam_reader.fetch(
                contig_acc, region_lend - 1, region_rend
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
        candidates_retained = 0  # number of reads not lightened (eligible for correction)

        for read in read_fetcher:
            processed += 1

            # The retention policy lives in Util_funcs.quant_discard_reason, which
            # this loop is the reference consumer of.  It was inline here, which
            # meant nothing else could ask what this loop keeps: cut selection, the
            # strand split and depth measurement each grew their own approximation,
            # and a filter added here would silently not reach any of them.
            #
            # The reasons are this counter's keys, so the tallies are unchanged.
            reason = Util_funcs.quant_discard_reason(
                read,
                contig_strand,
                max_intron_length=MAX_INTRON_LENGTH,
                min_mapping_quality=MIN_MAPPING_QUALITY,
                min_per_id=LRAA_Globals.config["min_per_id"],
            )

            # A strand mismatch is not a discard: this contig-strand simply is not
            # this read's, and the opposite-strand pass will take it.  Counting it
            # would report half the library as dropped.
            if reason == "wrong_strand" or reason == "unmapped":
                continue

            if reason is not None:
                discarded_read_counter[reason] += 1
                if reason == "low_perID":
                    read_name = Util_funcs.get_read_name_include_sc_encoding(read)
                    logger.debug(
                        "read {} has insufficient per_id, < min {} required ".format(
                            read_name, LRAA_Globals.config["min_per_id"]
                        )
                    )
                    self._last_discarded_read_names_by_reason["low_perID"].add(
                        read_name
                    )
                    num_alignments_per_id_fail += 1
                continue

            # Retained, and its identity was measurable: the QC ratio below is over
            # alignments that carried NM or nM, so a read without either is neither
            # a pass nor a failure.  Two tag probes rather than recomputing the
            # identity the predicate already evaluated.
            if read.has_tag("NM") or read.has_tag("nM"):
                num_alignments_per_id_ok += 1

            if pretty:
                # Build Pretty_alignment on the fly. Immediately lighten non-candidates for
                # soft-clip realignment to avoid retaining heavy pysam objects in memory.
                pa = Pretty_alignment.get_pretty_alignment(read)
                if force_lighten_all:
                    # for oversimplify contigs: never retain the heavy pysam record
                    pa.lighten()
                else:
                    is_candidate = pa.is_softclip_realign_candidate()
                    if not is_candidate:
                        pa.lighten()
                    else:
                        candidates_retained += 1
                pretty_alignments.append(pa)
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
                    frac_cand_kept = (candidates_retained / kept_so_far) if kept_so_far else 0.0
                    frac_cand_proc = (candidates_retained / processed) if processed else 0.0
                    prefix = f"[{contig_acc}{contig_strand}] " if contig_strand is not None else f"[{contig_acc}] "
                    logger.info(
                        prefix
                        + (
                            f"progress get_read_alignments: processed={processed:,}, kept={kept_so_far:,}, candidates={candidates_retained:,} ({frac_cand_kept:.3f} of kept; {frac_cand_proc:.3f} of processed), discards={discards}, rss={(f'{m:.1f} MB' if m is not None else '<na>')}"
                        )
                    )
                except Exception:
                    logger.info(
                        "progress get_read_alignments: processed=%s, kept=%s, candidates=%s",
                        f"{processed:,}",
                        f"{kept_so_far:,}",
                        f"{candidates_retained:,}",
                    )
                last_log_t = time.time()

        kept_count = len(pretty_alignments) if pretty else len(read_alignments)
        final_mem = _mem_usage_mb()
        try:
            frac_cand_final = (candidates_retained / kept_count) if kept_count else 0.0
        except Exception:
            frac_cand_final = 0.0
        logger.info(
            f"[{contig_acc}{contig_strand}] reads kept: {kept_count:,} and discarded: {dict(discarded_read_counter)} | candidates not lightened: {candidates_retained:,} ({frac_cand_final:.3f} of kept) (rss: {f'{final_mem:.1f} MB' if final_mem is not None else '<na>'})"
        )

        checked = num_alignments_per_id_fail + num_alignments_per_id_ok
        if (
            checked
            >= LRAA_Globals.config["min_total_alignments_engage_frac_per_id_check"]
        ):
            frac_pass = num_alignments_per_id_ok / checked

            # Against the pass fraction, which is what the threshold is named
            # for.  This compared the *failure* fraction to it, so it fired
            # whenever fewer than 90% of alignments failed -- that is, on every
            # healthy library -- and the message it then built referenced a
            # local that no longer exists, raising NameError instead of
            # reporting anything.  Not fatal by design: a low pass rate can be a
            # deliberate setting for a noisy library, and refusing to run would
            # leave no way to see what the threshold is doing.
            if (
                frac_pass
                < LRAA_Globals.config["min_frac_alignments_pass_per_id_check"]
            ):
                logger.warning(
                    f"Only {frac_pass:.3f} of {checked:,} identity-bearing alignments on "
                    f"{contig_acc} meet the {LRAA_Globals.config['min_per_id']} percent "
                    f"identity minimum. Please reevaluate --min_per_id for these alignments."
                )

        if pretty:
            return pretty_alignments
        else:
            return read_alignments

    def get_last_discarded_read_names_by_reason(self):
        return {
            reason: set(read_names)
            for reason, read_names in self._last_discarded_read_names_by_reason.items()
        }
