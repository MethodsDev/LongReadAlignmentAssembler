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

    def _fetch_region(
        self, contig_acc, contig_strand=None, region_lend=None, region_rend=None
    ):
        """The pysam iterator every intake in this class reads from.

        region_lend/region_rend arrive 1-based inclusive, straight from
        --region contig:lend-rend.  pysam.fetch takes 0-based half-open
        [start, stop), so passing region_lend unconverted started the window one
        base late and silently dropped any alignment whose only overlap with the
        region was its first base.  The right edge needs no adjustment: 1-based
        region_rend is 0-based region_rend - 1, which is below the exclusive stop
        and therefore included.

        Shared rather than inlined per caller because that conversion is the one
        place a second intake could disagree with get_read_alignments about WHICH
        reads a region holds, and a shard read-count that disagreed would report a
        different denominator than the run it is accounting for.
        """

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
            return self._pysam_reader.fetch(contig_acc, region_lend - 1, region_rend)

        logger.debug(
            "Fetching all alignments for contig: {} strand {}".format(
                contig_acc, contig_strand
            )
        )
        return self._pysam_reader.fetch(contig_acc)

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

        return list(
            self.iter_read_alignments(
                contig_acc,
                contig_strand,
                region_lend,
                region_rend,
                pretty=pretty,
                per_id_QC_raise_error=per_id_QC_raise_error,
                config=config,
                force_lighten_all=force_lighten_all,
            )
        )

    def iter_read_alignments(
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
        """get_read_alignments one read at a time, retaining none of them.

        get_read_alignments is this generator collected into a list, so the two cannot
        disagree about which reads a region holds. A caller that handles each read on
        its own -- the oversimplify paths, which score a read and write its tracking row
        -- iterates this instead and holds one alignment, not the contig's library:
        PBMC chrM is 5.5M records on one strand.

        The summary log line and the per-id QC warning are emitted when the iteration
        is exhausted, and get_last_discarded_read_names_by_reason describes this pass
        only from then on.
        """

        discarded_read_counter = defaultdict(int)
        self._last_discarded_read_names_by_reason = defaultdict(set)

        MIN_MAPPING_QUALITY = int(LRAA_Globals.config["min_mapping_quality"])
        MAX_INTRON_LENGTH = int(LRAA_Globals.config["max_intron_length"])

        read_fetcher = self._fetch_region(
            contig_acc, contig_strand, region_lend, region_rend
        )

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

            yield pa if pretty else read

        kept_count = kept_so_far
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

    def count_retained_reads(
        self,
        contig_acc,
        contig_strand=None,
        region_lend=None,
        region_rend=None,
        include_low_perID=False,
    ):
        """How many distinct reads get_read_alignments would keep over the same window.

        Same fetch window, same Util_funcs.quant_discard_reason policy, same
        Util_funcs.get_read_name_include_sc_encoding naming -- it differs from
        get_read_alignments only in building no Pretty_alignment and retaining no
        pysam record, so a caller that wants the size of a region's read population
        does not pay for the alignments it will not use. With include_low_perID, the
        reads discarded for low_perID count too.

        Exists for LRAA's count-then-report early returns: a work unit that
        quantifies nothing still has to say how many reads it SAW, and the only
        honest number is the one the population that unit would have quantified.
        Reimplementing the retention rule to get it is what produced divergent
        approximations before quant_discard_reason was factored out, which is why
        this asks that predicate rather than repeating it.

        Counted, not collected. The policy discards secondary, supplementary and
        improperly paired alignments, so an unpaired read leaves at most one counted
        record and needs no name to deduplicate; each record takes one outcome, so a
        kept read and a low_perID one are never the same. Only proper-pair mates
        share a name, and only they are deduplicated by it. MEASURED on PBMC Kinnex
        chrM (+): 5,444,358 kept and 1,613 low_perID records, 5,445,971 distinct
        names; the name set this replaces held 0.9 GB.
        """

        min_mapping_quality = int(LRAA_Globals.config["min_mapping_quality"])
        max_intron_length = int(LRAA_Globals.config["max_intron_length"])

        unpaired = {"kept": 0, "low_perID": 0}
        paired_names = {"kept": set(), "low_perID": set()}
        for read in self._fetch_region(
            contig_acc, contig_strand, region_lend, region_rend
        ):
            reason = Util_funcs.quant_discard_reason(
                read,
                contig_strand,
                max_intron_length=max_intron_length,
                min_mapping_quality=min_mapping_quality,
                min_per_id=LRAA_Globals.config["min_per_id"],
            )
            outcome = "kept" if reason is None else reason
            if outcome not in unpaired:
                continue
            if read.is_paired:
                paired_names[outcome].add(
                    Util_funcs.get_read_name_include_sc_encoding(read)
                )
            else:
                unpaired[outcome] += 1

        kept = unpaired["kept"] + len(paired_names["kept"])
        low_perID = unpaired["low_perID"] + len(
            paired_names["low_perID"] - paired_names["kept"]
        )

        logger.info(
            "[%s%s] count_retained_reads: %d read(s) kept, %d discarded for low_perID",
            contig_acc,
            contig_strand,
            kept,
            low_perID,
        )

        return kept + low_perID if include_low_perID else kept

    def get_last_discarded_read_names_by_reason(self):
        return {
            reason: set(read_names)
            for reason, read_names in self._last_discarded_read_names_by_reason.items()
        }
