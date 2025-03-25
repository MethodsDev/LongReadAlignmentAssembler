#!/usr/bin/env python3

import sys, os, re
import argparse
import pysam
import logging
import subprocess
import intervaltree as itree
import gzip
from collections import defaultdict

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../../pylib"])
)

from Transcript import *
import Util_funcs

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description="annotate splice patterns based on ref annot gtf",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input_matrix",
        type=str,
        required=True,
        help="input splice pattern matrix",
    )

    parser.add_argument(
        "--ref_annot_gtf",
        type=str,
        required=True,
        help="reference annotation used for inferring transcribed orientation of read (required if --infer_read_orient",
    )

    parser.add_argument(
        "--output_matrix",
        type=str,
        required=True,
        help="output matrix including gene annotations",
    )

    args = parser.parse_args()

    input_matrix = args.input_matrix
    ref_annot_gtf = args.ref_annot_gtf
    output_matrix = args.output_matrix

    #########
    ### begin

    splice_pattern_to_gene, intron_to_gene, stranded_chrom_to_itree = (
        build_chrom_itrees(ref_annot_gtf)
    )

    ofh = open(output_matrix, "wt")

    with open(input_matrix, "rt") as fh:
        header = next(fh)
        print(header, end="", file=ofh)

        for line in fh:

            gene_id = None

            vals = line.split("\t")
            splice_pattern = vals[0]

            if splice_pattern in splice_pattern_to_gene:
                gene_id = splice_pattern_to_gene[splice_pattern]
                gene_id = "FSM:" + gene_id
            else:
                # try intron
                contig_acc, strand, intron_coord_strings = parse_splice_pattern(
                    splice_pattern
                )

                for intron_coord_string in intron_coord_strings:
                    if intron_coord_string in intron_to_gene:
                        gene_id = intron_to_gene[intron_coord_string]
                        break

                if gene_id is not None:
                    gene_id = "NIC:" + gene_id
                else:
                    # try coordinate overlap
                    contig_acc, strand, intron_coord_strings = parse_splice_pattern(
                        splice_pattern
                    )
                    gene_id = find_overlapping_gene(
                        contig_acc,
                        strand,
                        intron_coord_strings,
                        stranded_chrom_to_itree,
                    )
                    if gene_id is not None:
                        gene_id = "NNIC:" + gene_id

            if gene_id is not None:
                vals[0] = gene_id + "^" + splice_pattern

            else:
                vals[0] = "novel^" + splice_pattern

            print("\t".join(vals), end="", file=ofh)

    logger.info("Done.")

    sys.exit(0)


def find_overlapping_gene(
    contig_acc, strand, intron_coord_strings, stranded_chrom_to_itree
):

    stranded_chrom = contig_acc + strand

    gene_scorer = defaultdict(int)

    for intron_coord_string in intron_coord_strings:
        m = re.search("\\[\\((\d+), (\d+)\\)\\]$", intron_coord_string)
        if m is None:
            raise RuntimeError(
                "Error, cannot parse intron coords from: {}".format(intron_coord_string)
            )
        lend = int(m.group(1))
        rend = int(m.group(2))

        for interval in stranded_chrom_to_itree[stranded_chrom][lend : rend + 1]:
            gene_id = interval.data
            gene_scorer[gene_id] += 1

    if len(gene_scorer) == 0:
        return None

    scored_genes = sorted(
        gene_scorer.keys(), key=lambda x: gene_scorer[x], reverse=True
    )

    return scored_genes[0]  # top scoring one.


def parse_splice_pattern(splice_pattern):

    m = re.match("^(\\S+):\\(([+-])\\)\\[(.*)\\]$", splice_pattern)

    if m is None:
        raise RuntimeError(
            "Error, cannot parse splice pattern: {}".format(splice_pattern)
        )

    contig = m.group(1)
    strand = m.group(2)
    intron_coords_super_string = m.group(3)

    intron_coord_sub_strings = intron_coords_super_string.split("), ")

    intron_coord_strings = list()
    for intron_coord_sub_string in intron_coord_sub_strings:
        if intron_coord_sub_string[-1] == ")":
            intron_coord_sub_string = intron_coord_sub_string[:-1]

        intron_coord_string = "{}:({})[{})]".format(
            contig, strand, intron_coord_sub_string
        )
        intron_coord_strings.append(intron_coord_string)
        # print(intron_coord_string)

    return contig, strand, intron_coord_strings


def build_chrom_itrees(gtf_file):

    logger.info("-building chrom itrees from: " + gtf_file)

    stranded_chrom_to_itree = defaultdict(lambda: itree.IntervalTree())

    splice_pattern_to_gene = dict()
    intron_to_gene = dict()

    contig_to_transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(gtf_file)

    for contig, transcript_list in contig_to_transcripts.items():

        for transcript in transcript_list:
            gene_id = transcript.get_gene_id()
            gene_name = transcript.get_gene_name()
            if gene_name is not None:
                gene_id = gene_name + "^" + gene_id

            transcript_strand = transcript.get_strand()
            trans_lend, trans_rend = transcript.get_coords()
            if not transcript.is_monoexonic():
                splice_pattern = transcript.get_introns_string()
                splice_pattern_to_gene[splice_pattern] = gene_id

                for intron in transcript.get_introns():
                    intron_lend, intron_rend = intron
                    intron_str = "{}:({})[({}, {})]".format(
                        contig, transcript_strand, intron_lend, intron_rend
                    )
                    intron_to_gene[intron_str] = gene_id

            stranded_chrom = contig + transcript_strand
            stranded_chrom_to_itree[stranded_chrom][
                trans_lend : trans_rend + 1
            ] = gene_id

    return splice_pattern_to_gene, intron_to_gene, stranded_chrom_to_itree


if __name__ == "__main__":
    main()
