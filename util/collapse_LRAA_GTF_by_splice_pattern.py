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
        description="Collapse an LRAA gtf by splice pattern. Isoforms are merged only "
        "when they share BOTH the same gene_id AND the same intron (splice) pattern, "
        "so collapsed identifiers stay consistent with the expression-matrix collapse "
        "(build_LRAA_expr_matrices.py keys on gene_id^splice_hash). gene_symbol^ "
        "prefixes (if present) are retained. Two report files are also written: a merge "
        "report (each collapsed isoform and the isoforms merged into it) and a "
        "gene-conflicts report (intron patterns shared across >1 gene_id, which are "
        "kept separate rather than merged).",
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
        help="output collapsed gtf filename",
    )

    parser.add_argument(
        "--merge_report",
        type=str,
        default=None,
        help="TSV listing each collapsed isoform and the isoforms merged into it "
        "(default: <output_gtf>.isoform_merge_report.tsv)",
    )

    parser.add_argument(
        "--gene_conflicts_report",
        type=str,
        default=None,
        help="TSV listing intron patterns shared across multiple gene_ids that were "
        "kept separate (default: <output_gtf>.gene_conflicts.tsv)",
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
    merge_report_file = args.merge_report or (output_gtf + ".isoform_merge_report.tsv")
    gene_conflicts_file = args.gene_conflicts_report or (
        output_gtf + ".gene_conflicts.tsv"
    )

    ofh = open(output_gtf, "wt")

    merge_ofh = open(merge_report_file, "wt")
    merge_ofh.write(
        "\t".join(
            [
                "collapsed_transcript_id",
                "gene_id",
                "num_isoforms_merged",
                "merged_transcript_ids",
            ]
        )
        + "\n"
    )

    conflict_ofh = open(gene_conflicts_file, "wt")
    conflict_ofh.write(
        "\t".join(
            [
                "contig",
                "splice_pattern_code",
                "num_gene_ids",
                "gene_id",
                "collapsed_transcript_id",
                "num_isoforms",
                "merged_transcript_ids",
            ]
        )
        + "\n"
    )

    logger.info(f"-capturing input transcripts from gtf {input_gtf}")
    contig_to_input_transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(
        input_gtf
    )

    num_conflicting_splice_patterns = 0

    for contig, transcript_obj_list in contig_to_input_transcripts.items():

        transcripts_to_output = list()

        # Group by (gene_id, splice_pattern_code): isoforms are merged only when they
        # share the same gene_id AND the same intron chain. Also track which gene_ids
        # carry each splice pattern, to report cross-gene conflicts.
        gene_splice_to_transcripts = defaultdict(list)
        splice_pattern_to_gene_ids = defaultdict(set)

        for transcript_obj in transcript_obj_list:

            if transcript_obj.has_introns():
                intron_string = transcript_obj.get_introns_string()
                splice_pattern_code = Util_funcs.get_hash_code(intron_string)
                gene_id = transcript_obj.get_gene_id()
                gene_splice_to_transcripts[(gene_id, splice_pattern_code)].append(
                    transcript_obj
                )
                splice_pattern_to_gene_ids[splice_pattern_code].add(gene_id)
            else:
                transcripts_to_output.append(transcript_obj)

        # collapse each (gene_id, splice_pattern) group, recording the mapping.
        group_to_collapsed = dict()  # (gene_id, splice_pattern_code) -> (new_id, member_ids)

        for (
            gene_id,
            splice_pattern,
        ), transcripts_same_group_list in gene_splice_to_transcripts.items():

            gene_symbol = None
            if "^" in gene_id:
                gene_symbol = gene_id.split("^")[0]
                new_transcript_id = "^".join([gene_symbol, splice_pattern])
            else:
                new_transcript_id = splice_pattern

            member_ids = sorted(
                [t.get_transcript_id() for t in transcripts_same_group_list]
            )

            if len(transcripts_same_group_list) == 1:
                transcript_obj = transcripts_same_group_list[0]
                transcript_obj.set_transcript_id(new_transcript_id)
                transcripts_to_output.append(transcript_obj)
            else:
                merged_isoform = merge_isoforms(transcripts_same_group_list)
                merged_isoform.set_transcript_id(new_transcript_id)
                transcripts_to_output.append(merged_isoform)

            group_to_collapsed[(gene_id, splice_pattern)] = (
                new_transcript_id,
                member_ids,
            )

            merge_ofh.write(
                "\t".join(
                    [
                        new_transcript_id,
                        gene_id,
                        str(len(member_ids)),
                        ",".join(member_ids),
                    ]
                )
                + "\n"
            )

        # report intron patterns shared across multiple gene_ids (kept separate).
        for splice_pattern_code, gene_id_set in splice_pattern_to_gene_ids.items():
            if len(gene_id_set) > 1:
                num_conflicting_splice_patterns += 1
                for gene_id in sorted(gene_id_set):
                    new_transcript_id, member_ids = group_to_collapsed[
                        (gene_id, splice_pattern_code)
                    ]
                    conflict_ofh.write(
                        "\t".join(
                            [
                                contig,
                                splice_pattern_code,
                                str(len(gene_id_set)),
                                gene_id,
                                new_transcript_id,
                                str(len(member_ids)),
                                ",".join(member_ids),
                            ]
                        )
                        + "\n"
                    )

        transcripts_to_output = sorted(
            transcripts_to_output, key=lambda x: x._exon_segments[0][0]
        )

        for transcript_obj in transcripts_to_output:
            ofh.write(transcript_obj.to_GTF_format(include_TPM=False) + "\n")

    ofh.close()
    merge_ofh.close()
    conflict_ofh.close()

    logger.info(f"-wrote collapsed gtf: {output_gtf}")
    logger.info(f"-wrote isoform merge report: {merge_report_file}")
    logger.info(
        f"-wrote gene-conflicts report "
        f"({num_conflicting_splice_patterns} splice patterns shared across >1 gene_id): "
        f"{gene_conflicts_file}"
    )
    logger.info("Done.")

    sys.exit(0)


def merge_isoforms(transcript_obj_list):

    first_transcript_obj = transcript_obj_list[0]
    template_exon_coords = first_transcript_obj.get_exon_segments()
    min_lend = template_exon_coords[0][0]
    max_rend = template_exon_coords[-1][1]
    contig_acc = first_transcript_obj.get_contig_acc()
    contig_strand = first_transcript_obj.get_strand()
    gene_id = first_transcript_obj.get_gene_id()
    has_TSS = first_transcript_obj.has_TSS()
    has_PolyA = first_transcript_obj.has_PolyA()

    merged_isoform_ids = [first_transcript_obj.get_transcript_id()]

    candidate_TSS_sites = set()
    candidate_PolyA_sites = set()

    if has_TSS:
        (
            candidate_TSS_sites.add(min_lend)
            if contig_strand == "+"
            else candidate_TSS_sites.add(max_rend)
        )

    if has_PolyA:
        (
            candidate_PolyA_sites.add(max_rend)
            if contig_strand == "+"
            else candidate_PolyA_sites.add(min_lend)
        )

    for transcript_obj in transcript_obj_list[1:]:
        exon_coords = transcript_obj.get_exon_segments()
        lend = exon_coords[0][0]
        rend = exon_coords[-1][1]

        if lend < min_lend:
            min_lend = lend

            if contig_strand == "+":
                has_TSS = transcript_obj.has_TSS()
            else:
                has_PolyA = transcript_obj.has_PolyA()

        if rend > max_rend:
            max_rend = rend
            if contig_strand == "+":
                has_PolyA = transcript_obj.has_PolyA()
            else:
                has_TSS = transcript_obj.has_TSS()

        if transcript_obj.has_TSS():
            (
                candidate_TSS_sites.add(lend)
                if contig_strand == "+"
                else candidate_TSS_sites.add(rend)
            )

        if transcript_obj.has_PolyA():
            (
                candidate_PolyA_sites.add(rend)
                if contig_strand == "+"
                else candidate_PolyA_sites.add(lend)
            )

        merged_isoform_ids.append(transcript_obj.get_transcript_id())

    template_exon_coords[0][0] = min_lend
    template_exon_coords[-1][1] = max_rend

    merged_transcript = Transcript(contig_acc, template_exon_coords, contig_strand)
    merged_transcript.set_gene_id(gene_id)

    merged_transcript._imported_has_TSS = has_TSS
    merged_transcript._imported_has_POLYA = has_PolyA

    merged_transcript.add_meta(
        "merged_isoforms_shared_splice_pattern", ",".join(sorted(merged_isoform_ids))
    )

    if len(candidate_TSS_sites) > 0:
        merged_transcript.add_meta(
            "TSS_sites", ",".join([str(x) for x in sorted(list(candidate_TSS_sites))])
        )

    if len(candidate_PolyA_sites) > 0:
        merged_transcript.add_meta(
            "PolyA_sites",
            ",".join([str(x) for x in sorted(list(candidate_PolyA_sites))]),
        )

    return merged_transcript


if __name__ == "__main__":
    main()
