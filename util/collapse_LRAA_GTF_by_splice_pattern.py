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
        description="Collapse an LRAA gtf by splice pattern. Isoforms are merged when "
        "they share BOTH the same gene_id AND the same intron (splice) pattern, so "
        "collapsed identifiers stay consistent with the expression-matrix collapse "
        "(build_LRAA_expr_matrices.py keys on gene_id^splice_hash). gene_symbol^ "
        "prefixes (if present) are retained. An intron pattern carried by more than "
        "one gene_id is REFUSED before anything is written, because identical splice "
        "patterns are one gene by definition and the collapsed transcript_id is minted "
        "from the pattern, so two such gene_ids would emit duplicate transcript_ids. "
        "Two report files are written: a merge report (each collapsed isoform and the "
        "isoforms merged into it) and a gene-conflicts report, which is what the "
        "refusal names.",
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
        help="TSV listing intron patterns carried by more than one gene_id. Written, "
        "and the run then refused, when any exist (default: "
        "<output_gtf>.gene_conflicts.tsv)",
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

    logger.info(f"-capturing input transcripts from gtf {input_gtf}")
    contig_to_input_transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(
        input_gtf
    )

    # Checked BEFORE the output gtf is opened, so a refused run leaves no partial file
    # to be mistaken for a collapse that worked.
    #
    # Not reconciled here, deliberately. Unifying the gene_ids would be a gene
    # REASSIGNMENT, and this script is not where gene identity is decided: the
    # expression matrices key their features on the gene_id carried by the quant files
    # (build_LRAA_expr_matrices.py:132-135), so a gene_id invented here would appear in
    # the collapsed gtf and nowhere else. Identical splice patterns are made one gene
    # upstream instead, by the clustering that assigns them
    # (GeneCommunityCluster._contract_identical_chains), which keeps the gtf and the
    # matrices describing the same genes.
    conflicts = _find_cross_gene_splice_patterns(contig_to_input_transcripts)
    _write_gene_conflicts_report(gene_conflicts_file, conflicts)
    if conflicts:
        exit(
            "Error, {} intron pattern(s) are carried by more than one gene_id; see {}. "
            "Isoforms with identical splice patterns are one gene by definition, and "
            "the collapsed transcript_id is minted from the splice pattern, so "
            "collapsing these would emit DUPLICATE transcript_ids -- silently, as the "
            "gtf has no uniqueness check. Refusing rather than picking a gene_id, "
            "because that choice belongs to gene clustering and would desynchronize "
            "this gtf from the expression matrices, which take gene_id from the quant "
            "files. This normally indicates stale or externally produced input, or an "
            "upstream clustering invariant that did not hold -- note that "
            "merge_LRAA_GTFs.py logs a warning and proceeds with unreclustered "
            "transcripts if reclustering raises, which bypasses the "
            "guarantee.".format(len(conflicts), gene_conflicts_file)
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

    for contig, transcript_obj_list in contig_to_input_transcripts.items():

        transcripts_to_output = list()

        # Group by (gene_id, splice_pattern_code): isoforms are merged only when they
        # share the same gene_id AND the same intron chain. Cross-gene patterns were
        # refused above, so every pattern here belongs to exactly one gene_id.
        gene_splice_to_transcripts = defaultdict(list)

        for transcript_obj in transcript_obj_list:

            if transcript_obj.has_introns():
                splice_pattern_code = _splice_pattern_code(transcript_obj)
                gene_id = transcript_obj.get_gene_id()
                gene_splice_to_transcripts[(gene_id, splice_pattern_code)].append(
                    transcript_obj
                )
            else:
                transcripts_to_output.append(transcript_obj)

        # collapse each (gene_id, splice_pattern) group.

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

        transcripts_to_output = sorted(
            transcripts_to_output, key=lambda x: x._exon_segments[0][0]
        )

        for transcript_obj in transcripts_to_output:
            ofh.write(transcript_obj.to_GTF_format(include_TPM=False) + "\n")

    ofh.close()
    merge_ofh.close()

    # Cheap, and it pins the property that motivated the refusal above: the collapsed
    # id is derived from the splice pattern, so a duplicate here means two gene_ids
    # shared a pattern and the check missed it. The gtf format has no uniqueness
    # constraint, so nothing else would notice.
    _assert_unique_transcript_ids(output_gtf)

    logger.info(f"-wrote collapsed gtf: {output_gtf}")
    logger.info(f"-wrote isoform merge report: {merge_report_file}")
    logger.info(f"-wrote gene-conflicts report (no conflicts): {gene_conflicts_file}")
    logger.info("Done.")

    sys.exit(0)


def _splice_pattern_code(transcript_obj):
    """The hash the collapsed transcript_id is minted from.

    Its input includes contig and strand (Transcript.get_introns_string), so the code
    identifies a splice pattern genome-wide rather than within a contig.
    """
    return Util_funcs.get_hash_code(transcript_obj.get_introns_string())


def _find_cross_gene_splice_patterns(contig_to_input_transcripts):
    """[(contig, splice_pattern_code, {gene_id: [transcript_id, ...]})] for patterns
    carried by more than one gene_id. Empty list when the input is well formed."""
    found = []
    for contig, transcript_obj_list in contig_to_input_transcripts.items():
        pattern_to_genes = defaultdict(lambda: defaultdict(list))
        for transcript_obj in transcript_obj_list:
            if not transcript_obj.has_introns():
                continue
            pattern_to_genes[_splice_pattern_code(transcript_obj)][
                transcript_obj.get_gene_id()
            ].append(transcript_obj.get_transcript_id())
        for code in sorted(pattern_to_genes):
            genes = pattern_to_genes[code]
            if len(genes) > 1:
                found.append((contig, code, genes))
    return found


def _write_gene_conflicts_report(path, conflicts):
    """Written on every run, conflicts or not: it is a declared output, and a missing
    file is indistinguishable from a run that never checked."""
    with open(path, "wt") as fh:
        fh.write(
            "\t".join(
                [
                    "contig",
                    "splice_pattern_code",
                    "num_gene_ids",
                    "gene_id",
                    "num_isoforms",
                    "transcript_ids",
                ]
            )
            + "\n"
        )
        for contig, code, genes in conflicts:
            for gene_id in sorted(genes):
                transcript_ids = sorted(genes[gene_id])
                fh.write(
                    "\t".join(
                        [
                            contig,
                            code,
                            str(len(genes)),
                            gene_id,
                            str(len(transcript_ids)),
                            ",".join(transcript_ids),
                        ]
                    )
                    + "\n"
                )


def _assert_unique_transcript_ids(gtf_file):
    seen = set()
    with open(gtf_file, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            match = re.search(r'transcript_id\s+"([^"]*)"', cols[8])
            if match is None:
                continue
            tid = match.group(1)
            if tid in seen:
                exit(
                    "Error, collapsed gtf {} contains duplicate transcript_id {}. This "
                    "is a bug in the collapse, not bad input: cross-gene splice "
                    "patterns are refused before writing.".format(gtf_file, tid)
                )
            seen.add(tid)


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
