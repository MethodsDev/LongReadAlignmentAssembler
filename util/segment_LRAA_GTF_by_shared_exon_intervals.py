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
    """
    Fast segmentation using sweep-line algorithm:
    1. Collect all exons from all transcripts
    2. Sort by coordinates and find overlapping regions
    3. Break overlaps into unique segments
    4. Map each transcript's exons to the new segments
    """
    
    # Step 1: Collect all exons with transcript tracking
    all_exons = []
    for transcript_obj in transcript_list:
        transcript_id = transcript_obj.get_transcript_id()
        for exon_start, exon_end in transcript_obj.get_exon_segments():
            all_exons.append((exon_start, exon_end, transcript_id))
    
    # Step 2: Sort exons by start coordinate
    all_exons.sort(key=lambda x: (x[0], x[1]))
    
    # Step 3: Build unique segments from overlapping exons using sweep-line
    segment_boundaries = set()
    for exon_start, exon_end, _ in all_exons:
        segment_boundaries.add(exon_start)
        segment_boundaries.add(exon_end + 1)  # +1 for half-open intervals
    
    segment_boundaries = sorted(segment_boundaries)
    unique_segments = [(segment_boundaries[i], segment_boundaries[i + 1] - 1) 
                       for i in range(len(segment_boundaries) - 1)]
    
    # Step 4: For each transcript, find which segments fall within its exons
    transcript_to_segments = defaultdict(list)
    
    for transcript_obj in transcript_list:
        transcript_id = transcript_obj.get_transcript_id()
        exon_segments = sorted(transcript_obj.get_exon_segments())
        
        # Use two pointers to efficiently match segments to exons
        seg_idx = 0
        for exon_start, exon_end in exon_segments:
            # Skip segments that end before this exon starts
            while seg_idx < len(unique_segments) and unique_segments[seg_idx][1] < exon_start:
                seg_idx += 1
            
            # Collect all segments that fall within this exon
            temp_idx = seg_idx
            while temp_idx < len(unique_segments):
                seg_start, seg_end = unique_segments[temp_idx]
                if seg_start > exon_end:
                    break  # Past this exon
                # Check if segment is fully contained in exon
                if seg_start >= exon_start and seg_end <= exon_end:
                    transcript_to_segments[transcript_id].append([seg_start, seg_end])
                temp_idx += 1
    
    # Step 5: Update transcript objects with new segmented exons
    for transcript_obj in transcript_list:
        transcript_id = transcript_obj.get_transcript_id()
        if transcript_id in transcript_to_segments:
            transcript_obj._exon_segments = transcript_to_segments[transcript_id]
    
    return transcript_list


if __name__ == "__main__":
    main()
