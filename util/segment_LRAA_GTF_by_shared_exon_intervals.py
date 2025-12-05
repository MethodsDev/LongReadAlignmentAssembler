#!/usr/bin/env python3

import sys, os, re

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../pylib"])
)

from Splice_graph import Splice_graph
from Transcript import Transcript, GTF_contig_to_transcripts
from LRAA import LRAA
import LRAA_Globals
import logging
import argparse
from collections import defaultdict
import Util_funcs

FORMAT = (
    "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s:\n\t%(message)s\n"
)

logger = logging.getLogger()
logging.basicConfig(format=FORMAT, level=logging.INFO)


def main():

    parser = argparse.ArgumentParser(
        description="collapse LRAA gtf by splice pattern.  If gtf is annotated with gene_symbol^ prefixes, those will be retained.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--gtf",
        type=str,
        required=True,
        help="LRAA gtf to collapse",
    )

    parser.add_argument(
        "--output_gtf",
        type=str,
        required=True,
        help="prefix for output filenames",
    )

    parser.add_argument(
        "--debug",
        "-d",
        action="store_true",
        default=False,
        help="debug mode, more verbose",
    )

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
        LRAA_Globals.DEBUG = True

    input_gtf = args.gtf
    output_gtf = args.output_gtf

    ofh = open(output_gtf, "wt")

    logger.info(f"-capturing input transcripts from gtf {input_gtf}")
    contig_to_input_transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(
        input_gtf
    )

    for contig, transcript_obj_list in contig_to_input_transcripts.items():

        logger.info("-processing {}".format(contig))

        transcripts_to_output = list()

        # Group by strand to segment all overlapping exons on the same strand
        strand_to_transcripts = defaultdict(list)

        for transcript_obj in transcript_obj_list:
            strand_to_transcripts[transcript_obj.get_orient()].append(transcript_obj)

        for transcript_list in strand_to_transcripts.values():

            if len(transcript_list) == 1:
                transcripts_to_output.extend(transcript_list)

            else:
                modified_transcripts = segment_transcript_coordinates(transcript_list)
                transcripts_to_output.extend(modified_transcripts)

        transcripts_to_output = sorted(
            transcripts_to_output, key=lambda x: x._exon_segments[0][0]
        )

        for transcript_obj in transcripts_to_output:
            ofh.write(transcript_obj.to_GTF_format(include_TPM=False) + "\n")

    logger.info("Done.")

    sys.exit(0)


def segment_transcript_coordinates(transcript_list):

    transcript_id_to_obj = dict()
    all_transcript_exon_patterns = dict()
    for transcript_obj in transcript_list:
        transcript_id = transcript_obj.get_transcript_id()
        transcript_id_to_obj[transcript_id] = transcript_obj
        exon_segments = transcript_obj.get_exon_segments()
        # Sort intervals for each transcript for efficient searching
        all_transcript_exon_patterns[transcript_id] = sorted(exon_segments)

    partitioned = partition_intervals(all_transcript_exon_patterns)

    for transcript_id, partitioned_segments in partitioned.items():
        transcript_obj = transcript_id_to_obj[transcript_id]
        transcript_obj._exon_segments = partitioned_segments

    return transcript_list


def flatten_intervals(intervals_by_set):
    """
    Collect all interval boundaries from all sets and flatten into sorted unique positions.
    """
    boundaries = set()
    for intervals in intervals_by_set.values():
        for start, end in intervals:
            boundaries.add(start)
            boundaries.add(end + 1)  # use end+1 to make it half-open for clarity
    return sorted(boundaries)


def build_disjoint_segments(boundaries):
    """
    Create disjoint segments from sorted unique boundaries.
    """
    return [(boundaries[i], boundaries[i + 1] - 1) for i in range(len(boundaries) - 1)]


def segment_membership(segments, intervals_by_set):
    """
    For each disjoint segment, identify which input sets it belongs to.
    Uses binary search for efficiency.
    """
    membership = defaultdict(list)
    
    for set_name, intervals in intervals_by_set.items():
        # intervals are already sorted from segment_transcript_coordinates
        for seg_start, seg_end in segments:
            # Binary search to find the first interval that might contain this segment
            left, right = 0, len(intervals)
            while left < right:
                mid = (left + right) // 2
                if intervals[mid][1] < seg_start:
                    left = mid + 1
                else:
                    right = mid
            
            # Check intervals starting from the found position
            for i in range(left, len(intervals)):
                int_start, int_end = intervals[i]
                if int_start > seg_end:
                    break  # No more intervals can contain this segment
                # segment must be fully contained within the interval
                if int_start <= seg_start and seg_end <= int_end:
                    membership[set_name].append((seg_start, seg_end))
                    break
    
    return membership


def partition_intervals(intervals_by_set):
    boundaries = flatten_intervals(intervals_by_set)
    segments = build_disjoint_segments(boundaries)
    membership = segment_membership(segments, intervals_by_set)

    # Reconstruct new intervals per set - use membership dict directly
    partitioned = defaultdict(list)
    for set_name, seg_list in membership.items():
        # seg_list already contains the segments for this transcript
        partitioned[set_name] = [[start, end] for start, end in seg_list]
    
    return partitioned


if __name__ == "__main__":
    main()
