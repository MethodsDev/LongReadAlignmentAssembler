#!/usr/bin/env python3
# encoding: utf-8

import os
import time
import hashlib
import json
import pysam
import LRAA_Globals
import logging
from Pretty_alignment import Pretty_alignment
import Splice_graph
import intervaltree as itree
from Bam_alignment_extractor import Bam_alignment_extractor
import Util_funcs
import pickle
from typing import Any


# Bump whenever a change alters WHAT A RUN STORES in these pickles for identical
# inputs and identical settings.
#
# The name is deliberately not "extraction": the pickle is not raw extraction
# output. It is extracted, then soft-clip corrected, then terminal-intron pruned,
# then lightened, and this constant is the only element of the cache token that
# identifies ANY of that code. A change to the corrector counts as much as a
# change to the extractor -- without a bump, a pickle written by the previous
# behaviour is a valid hit and the change is silently undone on any resumed or
# --no_cleanup run.
#
#   v1  original
#   v2  region fetch converts 1-based inclusive bounds to pysam's 0-based
#       half-open window, so an alignment ending exactly at region_lend is no
#       longer dropped
CACHED_ALIGNMENT_METHOD_VERSION = 2


logger = logging.getLogger(__name__)


def _contig_seq_digest(contig_seq):
    """Identity of the genome sequence the soft-clip corrector compared reads against.

    Hashed in 1 MB slices rather than with one contig_seq.encode(): a 250 Mb
    chromosome would allocate a second full copy of itself, inside the per-contig
    worker whose peak RSS is already the reason alignments are lightened before they
    are stored. Measured on a 250 Mb string the whole digest costs ~0.8s, against
    minutes of extraction for a contig that size, and it is skipped entirely unless
    correction runs.
    """
    if contig_seq is None:
        return None
    hasher = hashlib.sha256()
    chunk = 1 << 20
    for offset in range(0, len(contig_seq), chunk):
        hasher.update(contig_seq[offset : offset + chunk].encode("ascii", "replace"))
    return f"{len(contig_seq)}:{hasher.hexdigest()[:16]}"


class Pretty_alignment_manager:

    def __init__(self, splice_graph, alignment_cache_dir = "__alignment_cache"):
        self._splice_graph = splice_graph
        self._last_discarded_read_names_by_reason = dict()
        # If caller didn't specify a custom dir and we have a per-worker tmp dir, prefer a structured subdir
        try:
            tmp_root = os.environ.get("LRAA_TMP_DIR")
        except Exception:
            tmp_root = None
        if alignment_cache_dir == "__alignment_cache" and tmp_root:
            self._alignment_cache_dir = os.path.join(tmp_root, "__alignment_cache")
        else:
            self._alignment_cache_dir = alignment_cache_dir

        # Avoid race conditions when multiple processes attempt to create the cache dir
        # Using exist_ok ensures parallel workers don't crash if the dir appears between check and mkdir
        try:
            os.makedirs(self._alignment_cache_dir, exist_ok=True)
        except Exception:
            # Best-effort: if another process created it concurrently, proceed
            if not os.path.isdir(self._alignment_cache_dir):
                raise

    # --- diagnostics helpers ---
    def _mem_usage_mb(self):
        """Return current process RSS in MB if available, else None."""
        try:
            import psutil  # type: ignore
            rss = psutil.Process(os.getpid()).memory_info().rss
            return rss / (1024.0 * 1024.0)
        except Exception:
            try:
                import resource  # type: ignore
                rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                # ru_maxrss units vary: KB on Linux, bytes on macOS. Heuristic convert to MB.
                if rss > 1e10:  # likely bytes
                    return rss / (1024.0 * 1024.0)
                else:  # likely KB
                    return rss / 1024.0
            except Exception:
                return None

    def _log_mem(self, event, extra=None):
        """Log RSS with an event label and optional extras (dict or str)."""
        m = self._mem_usage_mb()
        # Build [contig+strand] prefix when splice graph context is available
        try:
            ca = self._splice_graph.get_contig_acc() if self._splice_graph else None
            cs = self._splice_graph.get_contig_strand() if self._splice_graph else None
            prefix = f"[{ca}{cs}] " if ca and cs else ""
        except Exception:
            prefix = ""
        extra_txt = ""
        if isinstance(extra, dict):
            try:
                extra_txt = ", " + ", ".join(f"{k}={v}" for k, v in extra.items())
            except Exception:
                extra_txt = f", extra={extra}"
        elif isinstance(extra, str) and extra:
            extra_txt = f", {extra}"
        if m is not None:
            logger.info(f"{prefix}[mem] {event}: rss={m:.1f} MB{extra_txt}")
        else:
            logger.info(f"{prefix}[mem] {event}: rss=<unavailable>{extra_txt}")

    def retrieve_pretty_alignments(
        self,
        contig_acc,
        contig_strand,
        contig_seq,
        bam_file,
        region_lend=None,
        region_rend=None,
        use_cache=False,
        restrict_splice_type=None,
        try_correct_alignments=False,
        SE_read_encapsulation_mask=None,
        per_id_QC_raise_error=False,
    ):

        # progress/logging: starting pretty alignment retrieval for this contig/strand (and region if set)
        region_txt = (
            f":{region_lend}-{region_rend}"
            if region_lend is not None and region_rend is not None
            else ""
        )
        try:
            prefix = f"[{contig_acc}{contig_strand}] " if contig_acc and contig_strand else ""
        except Exception:
            prefix = ""
        logger.info(
            f"{prefix}-start: retrieving pretty alignments{region_txt} from {os.path.basename(bam_file)} (restrict_splice_type={restrict_splice_type}, try_correct_alignments={try_correct_alignments})"
        )

        t_start = time.time()
        self._log_mem("start retrieve_pretty_alignments")

        bam_identity = Util_funcs.file_identity_token(bam_file)

        contig_strand_token = f"{contig_acc}^{contig_strand}"
        alignment_cache_dir = self._alignment_cache_dir

        if region_lend is not None and region_rend is not None:
            contig_strand_token = f"{contig_acc}^{contig_strand}:{region_lend}-{region_rend}"

        # Oversimplify is settled here, before the token, because it force-disables
        # correction. While that happened further down -- after the token was already
        # built -- an oversimplify contig wrote UNCORRECTED alignments under a stem that
        # recorded correction as having run, and the next ordinary run over the same bam
        # read them back as corrected. That is a poisoned cache, not a stale one: the
        # entry is wrong rather than merely old, and nothing about it looks wrong.
        #
        # The rule the ordering enforces: every value the token names must be final
        # before the token is built. Anything below this point may read
        # try_correct_alignments but may not change it.
        oversimplify_enabled = bool(LRAA_Globals.config.get("oversimplify_enabled", False))
        oversimplify_contigs = set(LRAA_Globals.config.get("oversimplify_contigs", []) or [])
        oversimplify_this_contig = oversimplify_enabled and (contig_acc in oversimplify_contigs)
        if oversimplify_this_contig and try_correct_alignments:
            logger.info(
                f"[{contig_acc}{contig_strand}] oversimplify enabled for this contig: disabling alignment correction and forcing lightened pretty alignments"
            )
            try_correct_alignments = False

        signature = self._cached_alignment_signature(
            contig_acc,
            contig_strand,
            contig_seq,
            bam_identity,
            region_lend,
            region_rend,
            try_correct_alignments,
            per_id_QC_raise_error,
        )

        # A legible stem plus a digest of the whole signature. The stem is what makes a
        # cache directory readable -- the bam by identity rather than basename, so that
        # sampleA/aligned.bam and sampleB/aligned.bam cannot share an entry -- and the
        # digest is what makes it sound.
        #
        # Every cache defect found in this codebase was a token naming the inputs
        # someone remembered and omitting the rest, and a flat list of f-string
        # fragments is what let that keep happening: nothing about it shows what is
        # missing. The enumeration therefore lives in exactly one place, with a verdict
        # recorded for every input including the ones deliberately absent -- see
        # _cached_alignment_signature -- so a newly found input costs one dict entry.
        #
        # The ME and SE pickles are partitions of one extraction, so they share this
        # stem and its discard-provenance sidecar.
        cache_token = (
            f"{contig_strand_token}.{bam_identity}.pretty_alignments"
            f".v{CACHED_ALIGNMENT_METHOD_VERSION}.{self._signature_digest(signature)}"
        )

        def variant_cache_file(restrict_token):
            return os.path.join(
                alignment_cache_dir, f"{cache_token}.restrict-{restrict_token}.pkl"
            )

        all_alignment_cache_file = variant_cache_file(None)
        ME_alignment_cache_file = variant_cache_file("ME")
        SE_alignment_cache_file = variant_cache_file("SE")

        # Read names dropped during extraction (currently low percent identity) are
        # consumed downstream as transcriptome-rescue candidates, so they must survive a
        # cache hit. Requiring this file for a hit makes caches from earlier builds
        # re-extract once instead of silently reporting no discards.
        discard_cache_file = os.path.join(
            alignment_cache_dir, f"{cache_token}.discards.pkl"
        )

        alignment_cache_file = all_alignment_cache_file

        if restrict_splice_type == "ME":
            alignment_cache_file = ME_alignment_cache_file
        elif restrict_splice_type == "SE":
            # There is deliberately no SE-masked variant: the mask is applied per run
            # to this one. See the mask block at the end of this method.
            alignment_cache_file = SE_alignment_cache_file

        if (
            use_cache
            and self._cache_ready(alignment_cache_file)
            and self._cache_ready(discard_cache_file)
        ):
            logger.info(
                "[%s%s] reusing earlier-generated pretty alignments",
                contig_acc,
                contig_strand,
            )
            pretty_alignments = self._load_pickle_cache(alignment_cache_file)
            self._last_discarded_read_names_by_reason = self._load_discard_cache(
                discard_cache_file
            )
            try:
                cache_sz_mb = os.path.getsize(alignment_cache_file) / (1024.0 * 1024.0)
            except Exception:
                cache_sz_mb = None
            self._log_mem(
                "loaded pretty_alignments from cache",
                extra={
                    "n": len(pretty_alignments) if isinstance(pretty_alignments, list) else "?",
                    "cache_file_mb": f"{cache_sz_mb:.1f}" if cache_sz_mb is not None else "?",
                },
            )

        else:

            self._log_mem("init Bam_alignment_extractor")
            bam_extractor = Bam_alignment_extractor(bam_file)
            self._log_mem("created Bam_alignment_extractor")

            t0 = time.time()
            logger.info("[%s%s] begin get_read_alignments (pretty=True)", contig_acc, contig_strand)
            pretty_alignments = bam_extractor.get_read_alignments(
                contig_acc,
                contig_strand,
                region_lend,
                region_rend,
                pretty=True,
                per_id_QC_raise_error=per_id_QC_raise_error,
                force_lighten_all=oversimplify_this_contig,
            )
            self._last_discarded_read_names_by_reason = (
                bam_extractor.get_last_discarded_read_names_by_reason()
            )
            if use_cache:
                # Written before any alignment pickle so a ready variant cache always
                # has its discard provenance on disk.
                self._write_pickle_cache(
                    discard_cache_file, self._last_discarded_read_names_by_reason
                )
            self._log_mem(
                "completed get_read_alignments",
                extra={
                    "n": (
                        len(pretty_alignments)
                        if isinstance(pretty_alignments, list)
                        else "?"
                    ),
                    "sec": f"{(time.time()-t0):.2f}",
                },
            )

            ## correct alignments containing soft-clips
            if try_correct_alignments:
                # only pass candidates to the corrector to avoid accessing lightened (no pysam) objects
                candidates = [
                    pa for pa in pretty_alignments if pa.is_softclip_realign_candidate()
                ]
                t1 = time.time()
                logger.info(
                    "[%s%s] begin try_correct_alignments on candidates: %d / total: %d",
                    contig_acc,
                    contig_strand,
                    len(candidates),
                    len(pretty_alignments),
                )
                Pretty_alignment.try_correct_alignments(
                    candidates, self._splice_graph, contig_seq
                )
                self._log_mem(
                    "completed try_correct_alignments",
                    extra={
                        "sec": f"{(time.time()-t1):.2f}",
                        "candidates": len(candidates),
                        "total": len(pretty_alignments),
                    },
                )

            Pretty_alignment.prune_long_terminal_introns(
                pretty_alignments, self._splice_graph
            )
            self._log_mem("after prune_long_terminal_introns", extra={"n": len(pretty_alignments)})

            # store for reuse
            # Remove pysam records before caching; if oversimplify, most will already be lightened.
            self._log_mem("before lighten", extra={"n": len(pretty_alignments)})
            pretty_alignments = [x.lighten() for x in pretty_alignments]
            self._log_mem("after lighten", extra={"n": len(pretty_alignments)})
            if use_cache:
                self._write_pickle_cache(all_alignment_cache_file, pretty_alignments)
                logger.info(
                    "[%s%s] Saved %d alignments to cache: %s",
                    contig_acc,
                    contig_strand,
                    len(pretty_alignments),
                    all_alignment_cache_file,
                )
                try:
                    cache_sz_mb = os.path.getsize(all_alignment_cache_file) / (1024.0 * 1024.0)
                    logger.info(
                        "[%s%s] Cache file size (all): %.1f MB",
                        contig_acc,
                        contig_strand,
                        cache_sz_mb,
                    )
                except Exception:
                    pass

            # Define SE and ME alignments and cache them separately when requested
            # To reduce peak memory, only materialize the subset we actually need unless caching all variants
            if restrict_splice_type in ("ME", "SE") or use_cache:
                ME_alignments = []
                SE_alignments = []
                for pa in pretty_alignments:
                    if pa.has_introns():
                        ME_alignments.append(pa)
                    else:
                        SE_alignments.append(pa)

                self._log_mem("partitioned ME/SE", extra={"ME": len(ME_alignments), "SE": len(SE_alignments)})

                if use_cache:
                    # store the ME and SE alignments
                    self._write_pickle_cache(ME_alignment_cache_file, ME_alignments)
                    logger.info(
                        "[%s%s] Saved %d alignments to cache: %s",
                        contig_acc,
                        contig_strand,
                        len(ME_alignments),
                        ME_alignment_cache_file,
                    )
                    try:
                        cache_sz_mb = os.path.getsize(ME_alignment_cache_file) / (1024.0 * 1024.0)
                        logger.info(
                            "[%s%s] Cache file size (ME): %.1f MB",
                            contig_acc,
                            contig_strand,
                            cache_sz_mb,
                        )
                    except Exception:
                        pass

                    self._write_pickle_cache(SE_alignment_cache_file, SE_alignments)
                    logger.info(
                        "[%s%s] Saved %d alignments to cache: %s",
                        contig_acc,
                        contig_strand,
                        len(SE_alignments),
                        SE_alignment_cache_file,
                    )
                    try:
                        cache_sz_mb = os.path.getsize(SE_alignment_cache_file) / (1024.0 * 1024.0)
                        logger.info(
                            "[%s%s] Cache file size (SE): %.1f MB",
                            contig_acc,
                            contig_strand,
                            cache_sz_mb,
                        )
                    except Exception:
                        pass

                if restrict_splice_type == "ME":
                    pretty_alignments = ME_alignments
                elif restrict_splice_type == "SE":
                    pretty_alignments = SE_alignments

        if SE_read_encapsulation_mask is not None:
            assert restrict_splice_type == "SE"

            # Derived per run, never stored.
            #
            # The mask is the multi-exon transcript set this run assembled (LRAA:4515),
            # so two runs over one bam carry different masks whenever their guide
            # annotation or their assembly differs -- and every mask shared the single
            # filename `restrict-SE-masked`, so whichever run went second was served
            # the other one's masked SE reads. Keying a digest of the mask would close
            # that collision; not storing the result removes the question, because the
            # mask is then not an input to anything on disk.
            #
            # Affordable because it is a pure interval-tree filter over already
            # lightened alignments and needs no pysam record. Measured on chr19 2Mb
            # (315,867 alignments, 31,146 of them SE): 0.3-0.5s for an
            # annotation-sized mask of 367 multi-exon transcripts and 1.5s at 5,000,
            # against 17.9s for the extraction the cache is there to avoid. Needing no
            # pysam record is also what makes it derivable at all: soft-clip correction
            # reads one, and a Pretty_alignment still holding a pysam record cannot be
            # pickled at all, which is why lighten() exists.
            self._log_mem(
                "before apply_SE_read_encapsulation_mask",
                extra={"n": len(pretty_alignments)},
            )
            pretty_alignments = self.apply_SE_read_encapsulation_mask(
                pretty_alignments, SE_read_encapsulation_mask
            )
            self._log_mem("after apply_SE_read_encapsulation_mask", extra={"n": len(pretty_alignments)})

        self._log_mem("end retrieve_pretty_alignments", extra={"n": len(pretty_alignments), "sec": f"{(time.time()-t_start):.2f}"})
        return pretty_alignments

    def get_last_discarded_read_names_by_reason(self):
        return {
            reason: set(read_names)
            for reason, read_names in self._last_discarded_read_names_by_reason.items()
        }

    def apply_SE_read_encapsulation_mask(self, pretty_alignments, SE_read_encapsulation_mask):

        # apply the SE read encapsulation mask
        # exclude those SE alignments that have substantial overlap with exons of multi-exon isoforms.
        # Uses percentage-based overlap: if >= min_SE_read_ME_exon_overlap_pct of the SE read length
        # overlaps with any ME exon, the SE read is filtered out.

        import LRAA_Globals

        exon_itree = itree.IntervalTree()

        min_overlap_pct = LRAA_Globals.config["min_SE_read_ME_exon_overlap_pct"]

        non_overlapping_pretty_alignments = list()

        for transcript in SE_read_encapsulation_mask:
            exon_segments = transcript.get_exon_segments()
            for exon in exon_segments:
                exon_lend, exon_rend = exon
                # store each exon as half-open interval [lend, rend+1); data not needed
                exon_itree[exon_lend:exon_rend + 1] = True

        for pretty_alignment in pretty_alignments:
            # only evaluate single-exon (SE) alignments; keep others unchanged
            if hasattr(pretty_alignment, "get_pretty_alignment_segments"):
                segs = pretty_alignment.get_pretty_alignment_segments()
                assert len(segs) == 1, "Error, should only apply mask to SE alignments: {}".format(pretty_alignment)

            align_lend, align_rend = pretty_alignment.get_alignment_span()
            align_len = align_rend - align_lend + 1

            should_filter = False
            # intervaltree expects either a point lookup tree[point] or a slice tree[start:stop]
            # Here we want all exons overlapping the alignment span
            for exon_interval in exon_itree[align_lend:align_rend + 1]:
                exon_lend = exon_interval.begin
                exon_rend = exon_interval.end - 1  # convert from half-open to inclusive

                # Calculate overlap length
                overlap_lend = max(align_lend, exon_lend)
                overlap_rend = min(align_rend, exon_rend)
                overlap_len = max(0, overlap_rend - overlap_lend + 1)

                # Calculate percentage of SE read covered by this ME exon
                overlap_pct = (overlap_len / align_len) * 100.0

                if overlap_pct >= min_overlap_pct:
                    should_filter = True
                    logger.debug("Excluding SE alignment with {:.1f}% overlap with ME exon: {}".format(overlap_pct, pretty_alignment))
                    break

            if not should_filter:
                non_overlapping_pretty_alignments.append(pretty_alignment)

        return(non_overlapping_pretty_alignments)

    # --- cache helpers ---

    def _splice_graph_intron_digest(self):
        """Identity of everything the soft-clip corrector can read from the graph.

        Pretty_alignment.try_correct_alignments touches the splice graph through
        get_overlapping_introns() and through nothing else, and uses only get_coords()
        of what comes back (Pretty_alignment.py:473, 486, 551, 564). The intron
        coordinate set is therefore not a proxy for the graph: it is exactly and
        completely what the graph contributes to a corrected alignment, so two graphs
        sharing this digest correct identically and two differing anywhere that matters
        do not share it.

        Deliberately not LRAA._compute_splice_graph_cache_entry, which digests the
        graph's CONSTRUCTION inputs. That is a broader and weaker proxy for the same
        thing -- two different guide annotations can yield one intron set -- it lives in
        an extension-less driver script this module cannot import, and reusing it would
        couple two cache systems' versioning, so a bump on that side would silently
        invalidate every alignment pickle as well.
        """
        if self._splice_graph is None:
            return None
        intron_coords = self._splice_graph.get_intron_coords()
        hasher = hashlib.sha256()
        for intron_lend, intron_rend in intron_coords:
            hasher.update(f"{intron_lend}-{intron_rend};".encode("ascii"))
        return f"{len(intron_coords)}:{hasher.hexdigest()[:16]}"

    def _cached_alignment_signature(
        self,
        contig_acc,
        contig_strand,
        contig_seq,
        bam_identity,
        region_lend,
        region_rend,
        try_correct_alignments,
        per_id_QC_raise_error,
    ):
        """Every input that decides what a run would store, with a verdict for each.

        The pickle is not raw extraction output, so "what extraction reads" is not the
        boundary of this dict. Anything that changes a field of a stored
        Pretty_alignment, or which alignments are stored at all, must be keyed here or
        else eliminated as an input, as the SE encapsulation mask was. The audit below
        is the point of this method: completeness can be checked against it without
        re-deriving the call graph, and an input that is knowingly absent says so here
        rather than by not appearing.

        KEYED -- which alignments come back at all:
          min_mapping_quality        Bam_alignment_extractor.py:58
          max_intron_length          :59, and again Pretty_alignment.py:330, :341 where
                                     it prunes terminal introns from stored segments
          min_per_id                 :164; --HiFi raises it 80 -> 97

        KEYED -- what each stored alignment holds:
          read_aln_gap_merge_int     merges adjacent aligned blocks into one segment, so
                                     it decides the segments themselves
                                     (Pretty_alignment.py:893). Taken from the class
                                     attribute, not from config: the attribute is bound
                                     once at import (Pretty_alignment.py:35) and is what
                                     :893 consults, so a config value written after
                                     import is not the value that shaped the segments.
                                     Keying the effective value stays honest whichever
                                     way that wiring is later repaired.
          cell_barcode_tag,          compose read_name as barcode^umi^query_name
          read_umi_tag               (Util_funcs.py:83-90 via Pretty_alignment.py:56).
                                     These decide read IDENTITY, so a collision here
                                     does not perturb a fraction of a count, it merges
                                     or splits reads for every downstream stage that
                                     groups or quantifies by name.
          min_PolyA_ident_length,    decide the stored left/right_soft_clipping by
          min_soft_clip_PolyA_       stripping polyA and untemplated G runs
            base_frac_for_           (Pretty_alignment.py:208-266). Those lengths are
            conversion,              what the graph later infers TSS and PolyA sites
          max_untemplated_G_at_TSS   from, and what selects correction candidates.

        KEYED -- correction, which cannot be deferred to a cache hit.
        try_correct_alignments corrects against the splice graph and the contig sequence
        (Pretty_alignment.py:473-505, :551-579) by reading
        pa._pysam_alignment.query_sequence and .get_aligned_pairs(). That record cannot
        be stored at all: pickling a Pretty_alignment still holding one raises
        "TypeError: self._delegate cannot be converted to a Python object for pickling",
        which is why lighten() exists. So caching raw extraction and re-deriving
        correction on every hit is not a cost trade-off, it is impossible, and the graph
        and sequence that determined the correction have to be named instead. Without
        them, two runs sharing the flag and differing in guide annotation reused each
        other's corrected alignments.
          try_correct_alignments     the post-override value; see the ordering rule at
                                     the oversimplify block
          min_softclip_realign_test, decide which alignments are corrected
          max_softclip_realign_test  (Pretty_alignment.py:288-290, :363-364)
          splice_graph_introns       see _splice_graph_intron_digest
          contig_seq                 see _contig_seq_digest
        The last four are None when correction is off, which is truthful -- an
        uncorrected run reads none of them -- and avoids splitting every uncorrected
        run's cache on values that did not reach it.

        ELIMINATED -- no longer an input to anything stored:
          SE_read_encapsulation_mask       applied per run; see the mask block in
          min_SE_read_ME_exon_overlap_pct  retrieve_pretty_alignments, and note that the
                                           pct is reached only by the mask

        KNOWINGLY UNKEYED:
          min_total_alignments_engage_frac_per_id_check,
          min_frac_alignments_pass_per_id_check
                                     read at Bam_alignment_extractor.py:254, :262 and
                                     used only to compose a logger.debug line, so they
                                     change nothing that is stored
          oversimplify_enabled,      their only effect on stored content is to force
          oversimplify_contigs       try_correct_alignments off, which the keyed
                                     post-override value already carries. force_lighten_all
                                     changes nothing else that survives: everything is
                                     lightened before it is written either way.
        """
        config = LRAA_Globals.config
        correcting = bool(try_correct_alignments)
        return {
            "contig": contig_acc,
            "strand": contig_strand,
            "region": [region_lend, region_rend],
            "bam": bam_identity,
            "min_mapping_quality": config["min_mapping_quality"],
            "min_per_id": config["min_per_id"],
            "max_intron_length": config["max_intron_length"],
            "read_aln_gap_merge_int": Pretty_alignment.read_aln_gap_merge_int,
            "cell_barcode_tag": config["cell_barcode_tag"],
            "read_umi_tag": config["read_umi_tag"],
            "min_PolyA_ident_length": config["min_PolyA_ident_length"],
            "min_soft_clip_PolyA_base_frac_for_conversion": config[
                "min_soft_clip_PolyA_base_frac_for_conversion"
            ],
            "max_untemplated_G_at_TSS": config["max_untemplated_G_at_TSS"],
            "try_correct_alignments": correcting,
            "min_softclip_realign_test": (
                config["min_softclip_realign_test"] if correcting else None
            ),
            "max_softclip_realign_test": (
                config["max_softclip_realign_test"] if correcting else None
            ),
            "splice_graph_introns": (
                self._splice_graph_intron_digest() if correcting else None
            ),
            "contig_seq": _contig_seq_digest(contig_seq) if correcting else None,
            # Over-keyed on purpose, and the one entry here that is not an input:
            # per_id_QC_raise_error reaches Bam_alignment_extractor.get_read_alignments
            # and is never read in its body, so it decides nothing that is stored.
            # Over-keying costs a redundant pickle at worst, whereas dropping it would
            # bet that the parameter stays dead.
            "per_id_QC_raise_error": bool(per_id_QC_raise_error),
        }

    @staticmethod
    def _signature_digest(signature) -> str:
        return hashlib.sha256(
            json.dumps(signature, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()[:16]

    def _cache_ready(self, cache_path: str) -> bool:
        """Return True when the pickle cache and its OK marker exist."""
        if not os.path.exists(cache_path):
            return False
        ok_path = self._cache_ok_path(cache_path)
        if not os.path.exists(ok_path):
            logger.warning("Cache file present but OK marker missing: %s", cache_path)
            return False
        return True

    def _load_pickle_cache(self, cache_path: str):
        with open(cache_path, "rb") as handle:
            return pickle.load(handle)

    def _load_discard_cache(self, cache_path: str) -> dict:
        """Return the cached discarded-read-name sets, or empty on an unusable payload."""
        try:
            payload = self._load_pickle_cache(cache_path)
        except Exception:
            logger.warning("Could not read discard provenance cache: %s", cache_path)
            return dict()
        if not isinstance(payload, dict):
            return dict()
        return {reason: set(read_names) for reason, read_names in payload.items()}

    def _write_pickle_cache(self, cache_path: str, payload: Any) -> None:
        """Atomically write a pickle cache and create its OK marker."""
        tmp_path = f"{cache_path}.tmp"
        ok_path = self._cache_ok_path(cache_path)
        try:
            try:
                if os.path.exists(ok_path):
                    os.remove(ok_path)
            except Exception:
                pass
            with open(tmp_path, "wb") as handle:
                pickle.dump(payload, handle)
            os.replace(tmp_path, cache_path)
            with open(ok_path, "w", encoding="utf-8") as marker:
                marker.write("ok\n")
        except Exception:
            try:
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)
            except Exception:
                pass
            try:
                if os.path.exists(ok_path):
                    os.remove(ok_path)
            except Exception:
                pass
            raise

    def _cache_ok_path(self, cache_path: str) -> str:
        return f"{cache_path}.ok"
