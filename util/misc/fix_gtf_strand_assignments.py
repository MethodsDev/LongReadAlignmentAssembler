#!/usr/bin/env python3

import sys, os, re
import argparse
import logging
import subprocess
import intervaltree as itree
import gzip
from collections import defaultdict

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../../pylib"])
)

import Util_funcs


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description="fix gtf strand info",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--output_gtf", type=str, required=True, help="output gtf file")

    parser.add_argument("--input_gtf", type=str, required=True, help="input gtf file")

    parser.add_argument("--genome", type=str, required=True, help="genome fasta file")

    args = parser.parse_args()

    input_gtf_file = args.input_gtf
    output_gtf_file = args.output_gtf
    genome_fasta = args.genome

    #########
    ### begin

    transcript_structs = parse_transcripts(input_gtf_file)

    prev_chrom = None
    chrom_seq = None

    ofh = open(output_gtf_file, "wt")

    for transcript_struct in transcript_structs:

        transcript_id = transcript_struct["transcript_id"]
        chrom = transcript_struct["chrom"]

        if prev_chrom is None or prev_chrom != chrom:
            prev_chrom = chrom
            chrom_seq = Util_funcs.retrieve_contig_seq_from_fasta_file(
                chrom, genome_fasta
            )

        # first try by dinuc splice sites of spliced introns from read
        strand = infer_spliced_orient(transcript_struct["coords"], chrom_seq)

        if strand == "?":
            logger.info("{} not spliced.".format(transcript_id))
            report_orig(transcript_struct["gtf_lines"], ofh)

        elif strand != transcript_struct["strand"]:
            logger.info("{} MISMATCHED STRAND ASSIGNMENT".format(transcript_id))
            report_fix_strand(transcript_struct["gtf_lines"], strand, ofh)

        else:
            logger.info("{} OK".format(transcript_id))
            report_orig(transcript_struct["gtf_lines"], ofh)

    ofh.close()

    sys.exit(0)


def get_introns(coordsets):

    coordsets = sorted(coordsets, key=lambda x: x[0])

    if len(coordsets) < 2:
        return []

    introns = list()
    for i in range(1, len(coordsets)):
        introns.append([coordsets[i - 1][1] + 1, coordsets[i][0] - 1])

    return introns


def infer_spliced_orient(alignment_segments, contig_seq):

    introns_coordsets = get_introns(alignment_segments)

    if len(introns_coordsets) < 1:
        return "?"

    return majority_vote_intron_orient(introns_coordsets, contig_seq)


def majority_vote_intron_orient(intron_coordsets, contig_seq):

    splice_dinucs_top_strand = {"GTAG", "GCAG", "ATAC"}
    splice_dinucs_bottom_strand = {
        "CTAC",
        "CTGC",
        "GTAT",
    }  # revcomp of top strand dinucs

    # orient_counts = {"+": 0, "-": 0}

    orient_counts = defaultdict(int)

    for intron_coordset in intron_coordsets:
        intron_lend, intron_rend = intron_coordset
        dinuc_left = contig_seq[intron_lend - 1] + contig_seq[intron_lend - 1 + 1]
        dinuc_right = contig_seq[intron_rend - 1 - 1] + contig_seq[intron_rend - 1]
        dinuc_combo = dinuc_left + dinuc_right

        if dinuc_combo in splice_dinucs_top_strand:
            orient_counts["+"] += 1
        elif dinuc_combo in splice_dinucs_bottom_strand:
            orient_counts["-"] += 1

    assert len(orient_counts) == 1, "Error, conflicting orientations: " + str(
        orient_counts
    )

    # check tie or not match
    if orient_counts["+"] == orient_counts["-"]:
        return "?"
    elif orient_counts["+"] > orient_counts["-"]:
        return "+"
    else:
        return "-"


def parse_transcripts(gtf_file):

    logger.info("-parsing transcripts from: " + gtf_file)

    if re.search(gtf_file, "\\.gz"):
        opener = gzip.open
    else:
        opener = open

    prev_transcript_id = None

    transcript_id_to_structs = defaultdict(
        lambda: {
            "transcript_id": None,
            "coords": [],
            "chrom": None,
            "strand": None,
            "gtf_lines": [],
        }
    )

    transcript_structs = list()

    with opener(gtf_file, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue

            line = line.rstrip()
            vals = line.split("\t")
            if len(vals) < 8:
                continue

            if vals[2] != "exon":
                continue

            chrom = vals[0]
            lend = int(vals[3])
            rend = int(vals[4])
            strand = vals[6]

            m = re.search('transcript_id \\"([^\\"]+)\\"', vals[8])
            if m is None:
                raise RuntimeError(
                    "Error, couldn't extract transcript_id from line: {}".format(line)
                )

            transcript_id = m.group(1)

            transcript_struct = transcript_id_to_structs[transcript_id]
            if transcript_id != prev_transcript_id:
                prev_transcript_id = transcript_id
                transcript_struct["transcript_id"] = transcript_id
                transcript_struct["chrom"] = chrom
                transcript_struct["strand"] = strand
                transcript_structs.append(transcript_struct)

            assert (
                transcript_struct["strand"] == strand
            ), "Error, strands don't match for {}".format(transcript_id)

            assert (
                transcript_struct["chrom"] == chrom
            ), "Error, chrom vals don't agree for transcript: {}".format(transcript_id)

            transcript_struct["gtf_lines"].append(line)
            transcript_struct["coords"].append([lend, rend])

    return transcript_structs


def report_orig(gtf_lines, ofh):

    for line in gtf_lines:
        print(line, file=ofh)

    return


def report_fix_strand(gtf_lines, strand, ofh):

    for line in gtf_lines:
        vals = line.split("\t")
        vals[6] = strand
        print("\t".join(vals), file=ofh)

    return


if __name__ == "__main__":
    main()
