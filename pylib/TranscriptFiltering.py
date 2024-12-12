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
import lmdb

logger = logging.getLogger(__name__)


def filter_monoexonic_isoforms_by_TPM_threshold(transcripts, min_TPM):

    transcripts_retained = list()

    for transcript in transcripts:
        tpm = transcript.get_TPM()
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

    def get_idoform_unique_assigned_read_count(transcript_id, frac_read_assignments):
        num_unique_reads = 0
        for read in frac_read_assignments[transcript_id]:
            if (
                frac_read_assignments[transcript_id][read] >= 0.9999
            ):  # close enough to 1.0
                num_unique_reads += 1

        return num_unique_reads

    gene_id_to_transcripts = _get_gene_id_to_transcripts(transcripts)

    all_transcripts_retained = list()

    # examine each gene separately:

    for gene_id, isoforms_of_gene in gene_id_to_transcripts.items():

        isoforms_were_filtered = True  # init for loop

        q = Quantify(run_EM, max_EM_iterations)

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

            for transcript in isoforms_of_gene:
                num_total_isoforms += 1
                transcript_id = transcript.get_transcript_id()

                transcript_unique_read_count = get_idoform_unique_assigned_read_count(
                    transcript_id, frac_read_assignments
                )

                frac_gene_unique_reads = transcript_unique_read_count / gene_read_count

                logger.debug(
                    "Transcript_id: {} has unique read frac of gene total reads: {}".format(
                        transcript_id, frac_gene_unique_reads
                    )
                )

                if not isoforms_were_filtered and (
                    frac_gene_unique_reads < min_frac_gene_unique_reads
                    or transcript.get_isoform_fraction() < min_isoform_fraction
                ):

                    isoforms_were_filtered = True
                    num_filtered_isoforms += 1

                    logger.debug(
                        "Filtering out transcript_id {} as low fraction of unique reads: {}".format(
                            transcript_id, frac_gene_unique_reads
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

            gene_read_count = gene_read_counts[gene_i_id]

            if transcript_i in transcript_prune_as_degradation:
                continue

            transcript_i_simple_path = transcript_i.get_simple_path()
            i_path_trimmed, i_TSS_id, i_polyA_id = SPU.trim_TSS_and_PolyA(
                transcript_i_simple_path, contig_strand
            )
            transcript_i_read_counts_assigned = transcript_i.get_read_counts_assigned()

            # for j in range(i+1, len(transcript_list)):

            for j in range(len(transcript_list)):

                if i == j:
                    continue

                transcript_j = transcript_list[j]

                if transcript_j in transcript_prune_as_degradation:
                    continue

                transcript_j_id = transcript_j.get_transcript_id()
                gene_j_id = transcript_j.get_gene_id()

                assert gene_i_id == gene_j_id

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
                    "Exploring path: {} as subsuming {}".format(
                        transcript_i_simple_path, transcript_j_simple_path
                    )
                )

                if SPU.path_A_contains_path_B(i_path_trimmed, j_path_trimmed):

                    logger.debug(
                        "splice compatible & contained transcript_j_id {} has frac gene expression: {}".format(
                            transcript_j_id, frac_gene_expression_j
                        )
                    )

                    subsume_J = False

                    if LRAA_Globals.config["collapse_alt_TSS_and_PolyA"]:
                        logger.debug(
                            "Collapsing compatible path: {} into {}".format(
                                transcript_j, transcript_i
                            )
                        )
                        subsume_J = True

                    elif j_TSS_id is not None:

                        j_TSS_read_count = sg.get_node_obj_via_id(
                            j_TSS_id
                        ).get_read_support()

                        frac_gene_express_j_TSS = j_TSS_read_count / gene_read_count
                        logger.debug(
                            "frac_gene_express_j_TSS {} {} : {:.3f}".format(
                                transcript_j_id,
                                transcript_j_simple_path,
                                frac_gene_express_j_TSS,
                            )
                        )

                        if (
                            frac_gene_express_j_TSS
                            < LRAA_Globals.config[
                                "min_frac_gene_alignments_define_TSS_site"
                            ]
                        ):
                            logger.debug(
                                "based on j_TSS count frac_gene_expression: {:.3f}, path_i: {} is subsuming path_j: {}".format(
                                    frac_gene_express_j_TSS,
                                    transcript_i_simple_path,
                                    transcript_j_simple_path,
                                )
                            )
                            subsume_J = True

                        elif i_TSS_id is not None:

                            i_TSS_read_count = sg.get_node_obj_via_id(
                                i_TSS_id
                            ).get_read_support()

                            frac_i_TSS = j_TSS_read_count / i_TSS_read_count
                            logger.debug(
                                "frac_i_TSS: {:.3f} of path_j: {} to path_i{}".format(
                                    frac_i_TSS,
                                    transcript_j_simple_path,
                                    transcript_i_simple_path,
                                )
                            )

                            if (
                                frac_i_TSS
                                < LRAA_Globals.config[
                                    "max_frac_alt_TSS_from_degradation"
                                ]
                            ):
                                logger.debug(
                                    "based on frac_i_TSS: {:.3f}, path_i: {} is subsuming path_j: {}".format(
                                        frac_i_TSS,
                                        transcript_i_simple_path,
                                        transcript_j_simple_path,
                                    )
                                )
                                subsume_J = True

                        elif i_TSS_id is None:
                            subsume_J = False

                    elif i_TSS_id is not None and j_TSS_id is None:
                        # no j_TSS but have i_TSS
                        subsume_J = True

                    else:
                        # neither has a TSS assigned.
                        subsume_J = True

                    #####  PolyA check ######
                    ## But dont subsume if they have polyA and they differ
                    if subsume_J:
                        # here we might resurrect it based on polyA status
                        if (
                            frac_gene_expression_j
                            >= LRAA_Globals.config[
                                "min_frac_alignments_define_polyA_site"
                            ]
                        ):

                            if i_polyA_id is not None and j_polyA_id is not None:
                                if i_polyA_id != j_polyA_id:
                                    subsume_J = False
                                    logger.debug(
                                        "resurrecting {} based on alt polyA".format(
                                            transcript_j_id
                                        )
                                    )
                            elif j_polyA_id is not None:
                                logger.debug(
                                    "resurrecting {} based on defined polyA".format(
                                        transcript_j_id
                                    )
                                )
                                subsume_J = False

                    if subsume_J:
                        logger.debug(
                            "Pruning {} as likely degradation product of {}".format(
                                transcript_j, transcript_i
                            )
                        )
                        transcript_prune_as_degradation.add(transcript_j)
                        transcript_i.add_read_names(transcript_j.get_read_names())

        # retain the ones not pruned
        for transcript in transcript_set:
            if transcript not in transcript_prune_as_degradation:
                transcripts_ret.append(transcript)

    return transcripts_ret


def filter_internally_primed_transcripts(
    transcripts, contig_seq_str, restrict_filter_to_monoexonic=True
):

    retained_transcripts = list()

    for transcript in transcripts:

        if restrict_filter_to_monoexonic and not transcript.is_monoexonic():
            retained_transcripts.append(transcript)
            continue

        # evaluate whether transcript looks internally primed.
        transcript_lend, transcript_rend = transcript.get_coords()
        strand = transcript.get_orient()

        if _looks_internally_primed(
            transcript_lend, transcript_rend, strand, contig_seq_str
        ):
            # gbye
            pass
        else:
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

    extracted_sequence = contig_seq_str[start - 1 : end].upper()

    has_flanking_polyA = (
        extracted_sequence.count(target_base) >= 12
        or target_polyA_motif in extracted_sequence
    )

    return has_flanking_polyA
