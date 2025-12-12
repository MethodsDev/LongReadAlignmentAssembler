#!/usr/bin/env python3
# encoding: utf-8

import os
import time
import pysam
import LRAA_Globals
import logging
from Pretty_alignment import Pretty_alignment
import Splice_graph
import intervaltree as itree
from Bam_alignment_extractor import Bam_alignment_extractor
import pickle
from typing import Any


logger = logging.getLogger(__name__)


class Pretty_alignment_manager:

    def __init__(self, splice_graph, alignment_cache_dir = "__alignment_cache"):
        self._splice_graph = splice_graph
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


    def retrieve_pretty_alignments(self, 
                                    contig_acc, contig_strand, contig_seq, bam_file, 
                                    region_lend=None,
                                    region_rend=None,
                                    use_cache=False,
                                    restrict_splice_type=None,
                                    try_correct_alignments=False,
                                    SE_read_encapsulation_mask=None,
                                    per_id_QC_raise_error=False):
        
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

        bam_file_basename = os.path.basename(bam_file)

        contig_strand_token = f"{contig_acc}^{contig_strand}"
        alignment_cache_dir = self._alignment_cache_dir

        if region_lend is not None and region_rend is not None:
            contig_strand_token = f"{contig_acc}^{contig_strand}:{region_lend}-{region_rend}"

        all_alignment_cache_file = os.path.join(
            alignment_cache_dir,
            f"{contig_strand_token}.{bam_file_basename}.pretty_alignments.restrict-{restrict_splice_type}.corr-{try_correct_alignments}.pkl",
        )

        ME_alignment_cache_file = os.path.join(
            alignment_cache_dir,
            f"{contig_strand_token}.{bam_file_basename}.pretty_alignments.restrict-ME.corr-{try_correct_alignments}.pkl",
        )

        SE_alignment_cache_file = os.path.join(
            alignment_cache_dir,
            f"{contig_strand_token}.{bam_file_basename}.pretty_alignments.restrict-SE.corr-{try_correct_alignments}.pkl",
        )

        SE_masked_alignment_cache_file = os.path.join(
            alignment_cache_dir,
            f"{contig_strand_token}.{bam_file_basename}.pretty_alignments.restrict-SE-masked.corr-{try_correct_alignments}.pkl",
        )

        alignment_cache_file = all_alignment_cache_file

        if restrict_splice_type == "ME":
            alignment_cache_file = ME_alignment_cache_file
        elif restrict_splice_type == "SE":
            if SE_read_encapsulation_mask is None:
                alignment_cache_file = SE_alignment_cache_file
            else:
                alignment_cache_file = SE_masked_alignment_cache_file

        if use_cache and self._cache_ready(alignment_cache_file):
            logger.info(
                "[%s%s] reusing earlier-generated pretty alignments",
                contig_acc,
                contig_strand,
            )
            pretty_alignments = self._load_pickle_cache(alignment_cache_file)
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

            # If oversimplify mode is enabled and this contig is listed, skip error correction
            # and do not retain pysam SAM records within Pretty_alignment objects at all.
            oversimplify_enabled = bool(LRAA_Globals.config.get("oversimplify_enabled", False))
            oversimplify_contigs = set(LRAA_Globals.config.get("oversimplify_contigs", []) or [])
            oversimplify_this_contig = oversimplify_enabled and (contig_acc in oversimplify_contigs)
            # Force-disable correction if oversimplify applies
            if oversimplify_this_contig and try_correct_alignments:
                logger.info(f"[{contig_acc}{contig_strand}] oversimplify enabled for this contig: disabling alignment correction and forcing lightened pretty alignments")
                try_correct_alignments = False

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
            self._log_mem(
                "completed get_read_alignments",
                extra={
                    "n": len(pretty_alignments) if isinstance(pretty_alignments, list) else "?",
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

            if use_cache and self._cache_ready(SE_masked_alignment_cache_file):
                logger.info(
                    "[%s%s] reusing earlier-generated pretty alignments",
                    contig_acc,
                    contig_strand,
                )
                pretty_alignments = self._load_pickle_cache(SE_masked_alignment_cache_file)
                self._log_mem("loaded SE masked alignments from cache", extra={"n": len(pretty_alignments)})
            else:
                # apply ME mask of encapsulated SEs
                self._log_mem("before apply_SE_read_encapsulation_mask", extra={"n": len(pretty_alignments)})
                pretty_alignments = self.apply_SE_read_encapsulation_mask(pretty_alignments, SE_read_encapsulation_mask)
                self._log_mem("after apply_SE_read_encapsulation_mask", extra={"n": len(pretty_alignments)})
                
                if use_cache:
                    self._write_pickle_cache(SE_masked_alignment_cache_file, pretty_alignments)
                    logger.info(
                        "[%s%s] Saved corrected alignments to cache: %s",
                        contig_acc,
                        contig_strand,
                        SE_masked_alignment_cache_file,
                    )

        self._log_mem("end retrieve_pretty_alignments", extra={"n": len(pretty_alignments), "sec": f"{(time.time()-t_start):.2f}"})
        return pretty_alignments                


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
