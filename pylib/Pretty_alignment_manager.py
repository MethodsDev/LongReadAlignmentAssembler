#!/usr/bin/env python3
# encoding: utf-8

import os
import pysam
import LRAA_Globals
import logging
from Pretty_alignment import Pretty_alignment
import Splice_graph
import intervaltree as itree
from Bam_alignment_extractor import Bam_alignment_extractor
import pickle


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
        try:
            region_txt = (
                f":{region_lend}-{region_rend}"
                if region_lend is not None and region_rend is not None
                else ""
            )
            logger.info(
                f"-start: retrieving pretty alignments for {contig_acc}{contig_strand}{region_txt} from {os.path.basename(bam_file)} "
                f"(restrict_splice_type={restrict_splice_type}, try_correct_alignments={try_correct_alignments})"
            )
        except Exception:
            # don't let logging failures interfere with processing
            pass

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

        if use_cache and os.path.exists(alignment_cache_file):
            logger.info(
                "reusing earlier-generated pretty alignments for {}{}".format(
                    contig_acc, contig_strand
                )
            )
            with open(alignment_cache_file, "rb") as f:
                pretty_alignments = pickle.load(f)

        else:

            bam_extractor = Bam_alignment_extractor(bam_file)
            pretty_alignments = bam_extractor.get_read_alignments(
                contig_acc,
                contig_strand,
                region_lend,
                region_rend,
                pretty=True,
                per_id_QC_raise_error=per_id_QC_raise_error,
            )

            ## correct alignments containing soft-clips
            if try_correct_alignments:
                Pretty_alignment.try_correct_alignments(
                    pretty_alignments, self._splice_graph, contig_seq
                )
               

            Pretty_alignment.prune_long_terminal_introns(
                pretty_alignments, self._splice_graph
            )

            # store for reuse
            pretty_alignments = [
                x.lighten() for x in pretty_alignments
            ]  # remove pysam record before storing
            if use_cache:
                with open(all_alignment_cache_file, "wb") as f:
                    pickle.dump(pretty_alignments, f)
                    logger.info(
                        f"Saved {len(pretty_alignments)} alignments to cache: {all_alignment_cache_file}"
                    )


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

                if use_cache:
                    # store the ME and SE alignments
                    with open(ME_alignment_cache_file, "wb") as f:
                        pickle.dump(ME_alignments, f)
                        logger.info(
                            f"Saved {len(ME_alignments)} alignments to cache: {ME_alignment_cache_file}"
                        )

                    with open(SE_alignment_cache_file, "wb") as f:
                        pickle.dump(SE_alignments, f)
                        logger.info(
                            f"Saved {len(SE_alignments)} alignments to cache: {SE_alignment_cache_file}"
                        )

                if restrict_splice_type == "ME":
                    pretty_alignments = ME_alignments
                elif restrict_splice_type == "SE":
                    pretty_alignments = SE_alignments

            
        if SE_read_encapsulation_mask is not None:
            assert restrict_splice_type == "SE"

            if use_cache and os.path.exists(SE_masked_alignment_cache_file):
                logger.info(
                    "reusing earlier-generated pretty alignments for {}{}".format(
                        contig_acc, contig_strand
                    )
                )
                with open(SE_masked_alignment_cache_file, "rb") as f:
                    pretty_alignments = pickle.load(f)
            else:
                # apply ME mask of encapsulated SEs
                pretty_alignments = self.apply_SE_read_encapsulation_mask(pretty_alignments, SE_read_encapsulation_mask)
                
                if use_cache:
                    with open(SE_masked_alignment_cache_file, "wb") as f:
                        pickle.dump(pretty_alignments, f)
                        logger.info(
                            f"Saved corrected alignments to cache: {SE_masked_alignment_cache_file}"
                        )

        return pretty_alignments                


    def apply_SE_read_encapsulation_mask(self, pretty_alignments, SE_read_encapsulation_mask):

        # apply the SE read encapsulation mask
        # exclude those alignments that are fully contained within exons of multi-exon isoforms.

        exon_itree = itree.IntervalTree()

        FUZZDIST = 10 # allow slight extension from exon boundaries

        non_encapsulated_pretty_alignments = list()

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
                assert len(segs) == 1, "Error, should only apply mask to SE alignment introns: {}".format(pretty_alignment)

            align_lend, align_rend = pretty_alignment.get_alignment_span()
            encapsulated = False
            # intervaltree expects either a point lookup tree[point] or a slice tree[start:stop]
            # Here we want all exons overlapping the alignment span
            for exon_interval in exon_itree[align_lend:align_rend + 1]:
                exon_lend = exon_interval.begin
                exon_rend = exon_interval.end - 1  # convert from half-open to inclusive
                if (align_lend >= (exon_lend - FUZZDIST) and
                        align_rend <= (exon_rend + FUZZDIST)):
                    encapsulated = True
                    logger.debug("Excluding SE alignment as encapsulated: {}".format(pretty_alignment))
                    break
            if not encapsulated:
                non_encapsulated_pretty_alignments.append(pretty_alignment)

        return(non_encapsulated_pretty_alignments)
        