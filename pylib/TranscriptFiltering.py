#!/usr/bin/env python3

import sys, os, re
import Transcript
import MultiPath
import MultiPathCounter
import Simple_path_utils as SPU
from Quantify import Quantify
from collections import defaultdict
import LRAA_Globals
from LRAA_Globals import SPACER, DEBUG
import logging
from math import log
import intervaltree as itree

logger = logging.getLogger(__name__)


def filter_monoexonic_isoforms_by_TPM_threshold(transcripts, min_TPM):

    transcripts_retained = list()

    for transcript in transcripts:
        tpm = transcript.get_TPM()

        # reftrans logic:
        if (
            transcript.includes_reference_transcript()
            and LRAA_Globals.config["ref_trans_filter_mode"] == "retain_expressed"
            and tpm > 0
        ):
            transcripts_retained.append(transcript)
            continue

        # regular filter logic
        if transcript.is_monoexonic() and tpm < min_TPM:
            # gbye
            pass
        else:
            # keep
            transcripts_retained.append(transcript)

    return transcripts_retained


def filter_isoforms_by_min_isoform_fraction(
    transcripts, min_isoform_fraction, run_EM, max_EM_iterations
):

    min_frac_gene_unique_reads = LRAA_Globals.config["min_frac_gene_unique_reads"]

    logger.info(
        "Filtering transcripts according to min isoform fraction: {}".format(
            min_isoform_fraction
        )
    )

    transcript_id_to_transcript_obj = dict(
        [(x.get_transcript_id(), x) for x in transcripts]
    )

    def get_isoform_unique_assigned_read_count(transcript_id, frac_read_assignments):
        num_unique_reads = 0
        for mp in frac_read_assignments[transcript_id]:
            if (
                frac_read_assignments[transcript_id][mp] >= 0.9999
            ):  # close enough to 1.0
                num_unique_reads += mp.get_read_count()

        return num_unique_reads

    gene_id_to_transcripts = _get_gene_id_to_transcripts(transcripts)

    all_transcripts_retained = list()

    # examine each gene separately:

    for gene_id, isoforms_of_gene in gene_id_to_transcripts.items():

        isoforms_were_filtered = True  # init for loop

        q = Quantify(run_EM, max_EM_iterations, quant_mode="draft")

        filtering_round = 0

        frac_read_assignments = None

        while isoforms_were_filtered:

            filtering_round += 1
            num_total_isoforms = 0
            num_filtered_isoforms = 0
            transcripts_retained = list()
            isoforms_were_filtered = (
                False  # update to True if we do filter an isoform out.
            )

            frac_read_assignments = q._estimate_isoform_read_support(isoforms_of_gene)

            gene_id_to_read_count = Quantify.get_gene_read_counts(
                frac_read_assignments, transcript_id_to_transcript_obj
            )

            gene_read_count = gene_id_to_read_count[gene_id]

            num_isoforms_of_gene_filtered = 0

            isoforms_of_gene = sorted(
                isoforms_of_gene, key=lambda x: x.get_isoform_fraction()
            )
            num_isoforms_of_gene = len(isoforms_of_gene)

            for transcript in isoforms_of_gene:

                num_total_isoforms += 1
                transcript_id = transcript.get_transcript_id()

                transcript_unique_read_count = get_isoform_unique_assigned_read_count(
                    transcript_id, frac_read_assignments
                )

                frac_gene_unique_reads = (
                    transcript_unique_read_count / gene_read_count
                    if gene_read_count > 0
                    else 0
                )

                logger.debug(
                    "Transcript_id: {} has unique read frac of gene total reads: {}".format(
                        transcript_id, frac_gene_unique_reads
                    )
                )

                ## first check to see if we should retain a reftrans
                if (
                    transcript.includes_reference_transcript()
                    and LRAA_Globals.config["ref_trans_filter_mode"]
                    == "retain_expressed"
                    and transcript.get_TPM() > 0
                ):
                    transcripts_retained.append(transcript)
                    continue

                ## if tons of isoforms, allow pruning of multiple in a single round
                if (
                    num_isoforms_of_gene - num_isoforms_of_gene_filtered
                    > LRAA_Globals.config[
                        "min_isoform_count_aggressive_filtering_iso_fraction"
                    ]
                    and frac_gene_unique_reads < min_frac_gene_unique_reads
                    and transcript.get_isoform_fraction() < min_isoform_fraction
                ):
                    isoforms_were_filtered = True
                    num_filtered_isoforms += 1
                    num_isoforms_of_gene_filtered += 1

                # standard isoform fraction based filtering
                elif not isoforms_were_filtered and (
                    frac_gene_unique_reads < min_frac_gene_unique_reads
                    or transcript.get_isoform_fraction() < min_isoform_fraction
                ):

                    isoforms_were_filtered = True
                    num_filtered_isoforms += 1
                    num_isoforms_of_gene_filtered += 1

                    logger.debug(
                        "Filtering out transcript_id {} as low fraction of unique reads: {}".format(
                            transcript_id, frac_gene_unique_reads
                        )
                    )

                elif not isoforms_were_filtered and (
                    transcript.is_novel_isoform() is True
                    and transcript_unique_read_count
                    < LRAA_Globals.config["min_unique_reads_novel_isoform"]
                ):
                    isoforms_were_filtered = True
                    num_filtered_isoforms += 1
                    num_isoforms_of_gene_filtered += 1

                    logger.debug(
                        "Filtering out transcript_id {} as novel isoform with too few unique reads: {}".format(
                            transcript_id, transcript_unique_read_count
                        )
                    )

                else:
                    transcripts_retained.append(transcript)

            logger.debug(
                "gene {} based isoform filtering round {} involved filtering of {} isoforms / {} total isoforms".format(
                    gene_id,
                    filtering_round,
                    num_filtered_isoforms,
                    num_total_isoforms,
                )
            )

            # reset list of transcripts
            isoforms_of_gene = transcripts_retained

        # done with individual gene-based isoform filtering
        all_transcripts_retained.extend(isoforms_of_gene)

    return all_transcripts_retained


def _get_gene_id_to_transcripts(transcripts):

    gene_id_to_transcripts = defaultdict(set)
    for transcript in transcripts:
        gene_id = transcript.get_gene_id()
        gene_id_to_transcripts[gene_id].add(transcript)

    return gene_id_to_transcripts


def prune_likely_degradation_products(transcripts, splice_graph, frac_read_assignments):

    logger.info("Pruning likely degradation products")

    sg = splice_graph

    transcript_id_to_transcript_obj = dict(
        [(x.get_transcript_id(), x) for x in transcripts]
    )
    gene_read_counts = Quantify.get_gene_read_counts(
        frac_read_assignments, transcript_id_to_transcript_obj
    )

    transcripts_ret = list()  # transcripts not pruned and returned.

    # first organize by gene_id
    gene_id_to_transcripts = _get_gene_id_to_transcripts(transcripts)

    for gene_id, transcript_set in gene_id_to_transcripts.items():

        gene_read_count = gene_read_counts[gene_id]

        if gene_read_count == 0:
            # should explore why this is. ref-trans-only?
            # gbye!
            continue

        if len(transcript_set) == 1:
            transcripts_ret.extend(list(transcript_set))
            continue

        # compare isoforms for compatibility ignoring TSS and PolyA sites
        # if an isoform is fully contained by another and has substantially less expression, prune it.

        transcript_list = list(transcript_set)
        contig_strand = transcript_list[0].get_strand()

        """  #TODO: figure out sorting order to speed this up.
        transcript_list = sorted(transcript_list, key=lambda x: (x.get_coords()[0], x.get_coords()[1], x.get_left_boundary_sort_weight(), x.get_right_boundary_sort_weight()))

        if contig_strand == '-':
            transcript_list = list(sorted(transcript_list, key=lambda x: (-1 * x.get_coords()[1],
                                                                          -1 * x.get_coords()[0],
                                                                          -1 * x.get_right_boundary_sort_weight(),
                                                                          -1 * x.get_left_boundary_sort_weight()
                                                                       ) ) ) 
        """

        # sort by desc cdna len
        transcript_list = list(
            reversed(sorted(transcript_list, key=lambda x: x.get_cdna_len()))
        )

        transcript_prune_as_degradation = set()
        for i in range(len(transcript_list)):
            transcript_i = transcript_list[i]
            transcript_i_id = transcript_i.get_transcript_id()

            gene_i_id = transcript_i.get_gene_id()
            assert gene_i_id == gene_id, "Error: gene_i_id {} != gene_id {} ".format(
                gene_i_id, gene_id
            )

            if transcript_i in transcript_prune_as_degradation:
                continue

            transcript_i_simple_path = transcript_i.get_simple_path()
            i_path_trimmed, i_TSS_id, i_polyA_id = SPU.trim_TSS_and_PolyA(
                transcript_i_simple_path, contig_strand
            )
            transcript_i_read_counts_assigned = transcript_i.get_read_counts_assigned()

            frac_gene_expression_i = transcript_i_read_counts_assigned / gene_read_count

            # for j in range(i+1, len(transcript_list)):

            for j in range(len(transcript_list)):

                if i == j:
                    continue

                transcript_j = transcript_list[j]

                if transcript_j in transcript_prune_as_degradation:
                    continue

                transcript_j_id = transcript_j.get_transcript_id()
                gene_j_id = transcript_j.get_gene_id()

                assert gene_j_id == gene_id

                transcript_j_simple_path = transcript_j.get_simple_path()
                j_path_trimmed, j_TSS_id, j_polyA_id = SPU.trim_TSS_and_PolyA(
                    transcript_j_simple_path, contig_strand
                )
                transcript_j_read_counts_assigned = (
                    transcript_j.get_read_counts_assigned()
                )

                if len(j_path_trimmed) > len(i_path_trimmed):
                    # no way can i subsume j
                    continue

                frac_gene_expression_j = (
                    transcript_j_read_counts_assigned / gene_read_count
                )

                logger.debug(
                    "Exploring path_i: {} {} as subsuming path_j {} {}".format(
                        transcript_i_id,
                        transcript_i_simple_path,
                        transcript_j_id,
                        transcript_j_simple_path,
                    )
                )

                frac_gene_expression_j_of_i_and_j = frac_gene_expression_j / (
                    frac_gene_expression_j + frac_gene_expression_i
                )

                ##
                ## Evaluate containments
                ##
                ##
                ## variables of interest: (in LRAA_Globals.config[]
                #    - "collapse_alt_TSS_and_PolyA"
                #  TSS:
                #    - "max_frac_alt_TSS_from_degradation"
                #    - "min_frac_gene_alignments_define_TSS_site"
                #  PolyA:
                #    - "min_frac_alignments_define_polyA_site"

                if SPU.path_A_contains_path_B(i_path_trimmed, j_path_trimmed):

                    logger.debug(
                        "splice compatible & contained transcript_j_id {} has frac gene expression: {}".format(
                            transcript_j_id, frac_gene_expression_j
                        )
                    )

                    subsume_J = False  # init

                    if LRAA_Globals.config["collapse_alt_TSS_and_PolyA"]:
                        logger.debug(
                            "Collapsing compatible path: {} into {}".format(
                                transcript_j, transcript_i
                            )
                        )
                        subsume_J = True

                    else:

                        # no TSS or PolyA on j - prune.
                        if j_TSS_id is None and j_polyA_id is None:
                            logger.debug(
                                "compatible/contained {} being pruned as lacking TSS or polyA annots.".format(
                                    transcript_j_id
                                )
                            )
                            subsume_J = True

                        elif (
                            frac_gene_expression_j_of_i_and_j
                            < LRAA_Globals.config[
                                "max_rel_frac_expr_alt_compat_contained"
                            ]
                        ):
                            # if relative fraction of support for both is below threshold, then prune.
                            logger.debug(
                                "compatible/contained {} being pruned as insufficiently relatively expressed.".format(
                                    transcript_j_id
                                )
                            )
                            subsume_J = True

                    if subsume_J:
                        logger.debug(
                            "Pruning {} as likely degradation product of {}".format(
                                transcript_j, transcript_i
                            )
                        )
                        transcript_prune_as_degradation.add(transcript_j)
                        transcript_i.absorb_other_transcript_multipaths(transcript_j)

        # retain the ones not pruned
        for transcript in transcript_set:
            if transcript not in transcript_prune_as_degradation:
                transcripts_ret.append(transcript)
            else:
                if LRAA_Globals.DEBUG:
                    logger.debug(
                        "FILTERING transcript {} as a likely degradation product".format(
                            transcript
                        )
                    )

    return transcripts_ret


def filter_internally_primed_transcripts(
    transcripts,
    contig_seq_str,
    contig_strand,
    splice_graph,
    restrict_filter_to_monoexonic,
):

    known_polyA_dist_ok_window = LRAA_Globals.config["max_dist_between_alt_polyA_sites"]
    known_polyA_dist_ok_window_half = int(known_polyA_dist_ok_window / 2)

    # build a list of known/acceptable 3' ends that get a free pass
    known_ok_3prime_ends = (
        splice_graph._input_transcript_rend_boundaries
        if contig_strand == "+"
        else splice_graph._input_transcript_lend_boundaries
    )
    known_ok_3prime_ends_itree = itree.IntervalTree()
    for ok_3prime_end in known_ok_3prime_ends:
        known_ok_3prime_ends_itree[
            max(1, ok_3prime_end - known_polyA_dist_ok_window_half) : (
                ok_3prime_end + known_polyA_dist_ok_window_half
            )
        ] = ok_3prime_end

    retained_transcripts = list()

    for transcript in transcripts:

        if restrict_filter_to_monoexonic and not transcript.is_monoexonic():
            # retaining all multi-exonic transcripts when restrict_filter_to_monoexonic set
            retained_transcripts.append(transcript)
            continue

        # evaluate whether transcript looks internally primed.
        transcript_lend, transcript_rend = transcript.get_coords()
        strand = transcript.get_orient()

        looks_internally_primed = _looks_internally_primed(
            transcript_lend, transcript_rend, strand, contig_seq_str
        )

        filter_flag = False

        if looks_internally_primed:

            if transcript.is_monoexonic():
                # gbye
                logger.debug(
                    "FILTERING monoexonic transcript {} as likely internally primed".format(
                        transcript
                    )
                )
                filter_flag = True

            else:
                # different rules for mutlti-exon transcripts.
                end_check = transcript_rend if contig_strand == "+" else transcript_lend

                if len(known_ok_3prime_ends_itree[end_check : end_check + 1]) == 0:
                    # no overlap with known acceptable site.
                    # gbye
                    logger.debug(
                        "FILTERING multiexonic transcript {} as likely internally primed".format(
                            transcript
                        )
                    )
                    filter_flag = True

        if not filter_flag:
            # keep
            retained_transcripts.append(transcript)

    return retained_transcripts


def _looks_internally_primed(
    transcript_lend, transcript_rend, strand, contig_seq_str, check_length=20
):

    if strand not in {"+", "-"}:
        raise ValueError("Strand must be '+' or '-'")

    target_base = "A" if strand == "+" else "T"
    target_polyA_motif = target_base * 6

    contig_length = len(contig_seq_str)

    if strand == "+":
        start = transcript_rend + 1
    else:
        start = transcript_lend - check_length - 1

    end = start + check_length

    # ensure coordinates within contig bounds
    start = max(1, start)
    end = min(end, contig_length)

    extracted_long_sequence = contig_seq_str[start - 1 : end].upper()
    extracted_short_sequence = (
        extracted_long_sequence[-8:] if strand == "-" else extracted_long_sequence[0:8]
    )

    has_flanking_polyA = (
        extracted_long_sequence.count(target_base) >= 12
        or target_polyA_motif in extracted_short_sequence
    )

    return has_flanking_polyA


def evaluate_splice_compatible_alt_isoforms(transcripts):

    transcript_id_to_splice_compatible_containments = defaultdict(set)
    transcript_id_to_splice_compatible_contained_by = defaultdict(set)

    logger.info("-evaluationg splice compatible alt isoforms:")

    if len(transcripts) < 2:
        return (
            transcript_id_to_splice_compatible_containments,
            transcript_id_to_splice_compatible_contained_by,
        )

    if type(transcripts) == set:
        transcripts = list(transcripts)

    transcripts.sort(key=lambda x: x.get_TPM(), reverse=True)

    for i in range(len(transcripts) - 1):

        transcript_i = transcripts[i]
        transcript_i_sp = transcript_i.get_simple_path()
        transcript_i_introns = SPU.get_simple_path_introns(transcript_i_sp)
        transcript_i_id = transcript_i.get_transcript_id()

        if len(transcript_i_introns) == 0:
            continue

        for j in range(i + 1, len(transcripts)):
            transcript_j = transcripts[j]
            transcript_j_sp = transcript_j.get_simple_path()
            transcript_j_introns = SPU.get_simple_path_introns(transcript_j_sp)
            transcript_j_id = transcript_j.get_transcript_id()

            if len(transcript_j_introns) == 0:
                continue

            if len(transcript_j_introns - transcript_i_introns) == 0:

                transcript_id_to_splice_compatible_containments[transcript_i_id].add(
                    transcript_j_id
                )
                transcript_id_to_splice_compatible_contained_by[transcript_j_id].add(
                    transcript_i_id
                )

                logger.debug(
                    "Splice compatible isoforms: {} expr: {} splice compatible with {} expr: {}".format(
                        transcript_i,
                        transcript_i.get_TPM(),
                        transcript_j,
                        transcript_j.get_TPM(),
                    )
                )

            if len(transcript_i_introns - transcript_j_introns) == 0:

                transcript_id_to_splice_compatible_containments[transcript_j_id].add(
                    transcript_i_id
                )
                transcript_id_to_splice_compatible_contained_by[transcript_i_id].add(
                    transcript_j_id
                )

                logger.debug(
                    "Splice compatible isoforms: {} expr: {} splice compatible with {} expr: {}".format(
                        transcript_j,
                        transcript_j.get_TPM(),
                        transcript_i,
                        transcript_i.get_TPM(),
                    )
                )

    return (
        transcript_id_to_splice_compatible_containments,
        transcript_id_to_splice_compatible_contained_by,
    )


def filter_novel_isoforms_by_min_read_support(
    transcripts, min_reads_novel_isoform: int
):

    retained_transcripts = list()

    for transcript in transcripts:
        if transcript.is_novel_isoform() is True:
            if transcript.get_read_counts_assigned() >= min_reads_novel_isoform:
                retained_transcripts.append(transcript)
            else:
                # novel transcript gets pruned as insufficient evidence.
                if LRAA_Globals.DEBUG:
                    logger.debug(
                        "FILTERING {} as insufficient evidence: read_support {} vs. min required {}".format(
                            transcript,
                            transcript.get_read_counts_assigned(),
                            min_reads_novel_isoform,
                        )
                    )
                pass
        else:
            # known transcript, retaining.
            retained_transcripts.append(transcript)

    return retained_transcripts
