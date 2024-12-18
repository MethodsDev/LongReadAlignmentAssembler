#!/usr/bin/env python3

import sys, os, re
import Transcript
from collections import defaultdict
import LRAA_Globals
from LRAA_Globals import SPACER, DEBUG
import logging
from math import log


logger = logging.getLogger(__name__)


def run_EM(transcripts, max_EM_iterations=1000):

    read_name_to_transcripts = get_read_name_to_transcripts(transcripts)
    num_mapped_reads = len(read_name_to_transcripts)

    ## Init - split across mapped transcripts according to normalized weights or equally if not using weights..
    (
        transcript_to_read_count,
        transcript_to_fractional_read_assignment,
        transcript_to_expr_val,
    ) = init_read_counts_and_expr_vals(transcripts, read_name_to_transcripts)

    ###################################
    ## go through multiple rounds of EM
    ###################################

    ## compute likelihood.

    prev_log_likelihood = compute_log_likelihood(
        transcript_to_fractional_read_assignment, transcript_to_expr_val
    )
    logger.debug(
        "Log likelihood before starting EM: {:.5E}".format(prev_log_likelihood)
    )

    for EM_round in range(1, max_EM_iterations):

        logger.debug("EM round {}".format(EM_round))

        ## fractionally assign reads based on expr values
        # Estep
        transcript_to_read_count, transcript_to_fractional_read_assignment = E_step(
            transcripts,
            read_name_to_transcripts,
            transcript_to_expr_val,
            transcript_to_fractional_read_assignment,
            EM_round,
        )

        # M-step
        transcript_to_expr_val = M_step(
            transcripts,
            transcript_to_read_count,
            transcript_to_fractional_read_assignment,
            num_mapped_reads,
            EM_round,
        )

        ## Determine log likelihood value and see if reached peak
        log_likelihood = compute_log_likelihood(
            transcript_to_fractional_read_assignment, transcript_to_expr_val
        )

        logger.debug("\tlog_likelihood: {:.5E}".format(log_likelihood))
        if prev_log_likelihood is not None:
            delta = log_likelihood - prev_log_likelihood
            logger.debug(f"\t\tEM_round[{EM_round}] delta = {delta:.5E}")
            if delta < 1e-5:
                logger.debug(
                    f"\t\tEM_round[{EM_round}] delta = {delta:.5E} ...  max likelihood reached. Stopping EM."
                )
                break

        prev_log_likelihood = log_likelihood

    return (
        transcript_to_expr_val,
        transcript_to_fractional_read_assignment,
        transcript_to_read_count,
    )


def E_step(
    transcripts,
    read_name_to_transcripts,
    transcript_to_expr_val,
    transcript_to_fractional_read_assignment,
    EM_round,
):

    local_debug = False  # for highly verbose logging

    ## E-step
    for transcript in transcripts:
        transcript_id = transcript.get_transcript_id()
        transcript_read_count_total = 0
        read_names = transcript.get_read_names()
        transcript_expr = transcript_to_expr_val[transcript_id]

        for read_name in read_names:
            all_transcripts_with_read = read_name_to_transcripts[read_name]

            # get sum of expr values for all isoforms w/ read assigned and weights included.
            sum_denom_expr = 0

            for (
                each_tran_with_read,
                each_tran_read_weight,
            ) in all_transcripts_with_read:

                each_tran_id = each_tran_with_read.get_transcript_id()

                sum_denom_expr_contrib = (
                    each_tran_read_weight * transcript_to_expr_val[each_tran_id]
                )
                sum_denom_expr += sum_denom_expr_contrib

            if local_debug:
                print(
                    f"EM{EM_round} transcript {transcript_id} sum_denom_expr {sum_denom_expr}"
                )

            read_weight = transcript.get_read_weight(read_name)
            if not LRAA_Globals.config["use_weighted_read_assignments"]:
                # not using the weight
                read_weight = 1.0

            frac_read_assignment = (
                read_weight * transcript_expr / sum_denom_expr
                if sum_denom_expr > 0
                else 0.0
            )
            if local_debug:
                print(
                    f"EM{EM_round} transcript {transcript_id} sum_denom_expr {sum_denom_expr} frac_read {read_name} = {frac_read_assignment}"
                )

            # for tracking purposes
            transcript_to_fractional_read_assignment[transcript_id][
                read_name
            ] = frac_read_assignment

    # renormalize fractional read assignments and recount.
    transcript_to_read_count, transcript_to_fractional_read_assignment = (
        renormalize_read_fractions(
            read_name_to_transcripts,
            transcript_to_fractional_read_assignment,
        )
    )

    return transcript_to_read_count, transcript_to_fractional_read_assignment


def renormalize_read_fractions(
    read_name_to_transcripts,
    transcript_to_fractional_read_assignment,
):

    transcript_to_read_count = defaultdict(float)

    # renormalize read fractions so they sum to 1.
    for read_name in read_name_to_transcripts:
        sum_frac_read_assignments = 0
        for transcripts_with_read in read_name_to_transcripts[read_name]:
            each_tran_with_read, read_weight = transcripts_with_read
            t_id = each_tran_with_read.get_transcript_id()
            frac_assignment = transcript_to_fractional_read_assignment[t_id][read_name]
            sum_frac_read_assignments += frac_assignment

        # now normalize read fractions so they sum to 1
        for transcripts_with_read in read_name_to_transcripts[read_name]:
            each_tran_with_read, _ignore_weight_here = transcripts_with_read
            t_id = each_tran_with_read.get_transcript_id()
            frac_assignment = transcript_to_fractional_read_assignment[t_id][read_name]
            transcript_to_fractional_read_assignment[t_id][read_name] = (
                frac_assignment / sum_frac_read_assignments
                if sum_frac_read_assignments > 0
                else 0
            )

            transcript_to_read_count[t_id] += frac_assignment  # for expr val calc later

    return transcript_to_read_count, transcript_to_fractional_read_assignment


def M_step(
    transcripts,
    transcript_to_read_count,
    transcript_to_fractional_read_assignment,
    num_mapped_reads,
    EM_round,
):

    ## M-step
    ## recompute expr_vals
    transcript_to_expr_val = dict()

    sum_transcript_assigned_expr_vals = 0
    for transcript in transcripts:
        transcript_id = transcript.get_transcript_id()
        transcript_read_count = transcript_to_read_count[transcript_id]

        transcript_to_expr_val[transcript_id] = (
            (transcript_read_count / num_mapped_reads) if num_mapped_reads > 0 else 0
        )  # * 1e6
        logger.debug(
            f"-EM round {EM_round} assigning transcript {transcript_id} read count: {transcript_read_count} and expr val {transcript_read_count}/{num_mapped_reads} = {transcript_to_expr_val[transcript_id]}"
        )
        sum_transcript_assigned_expr_vals += transcript_to_expr_val[transcript_id]

    # renormalize fractional expression values
    for transcript in transcripts:
        transcript_to_expr_val[transcript_id] = (
            transcript_to_expr_val[transcript_id] / sum_transcript_assigned_expr_vals
        )

    return transcript_to_expr_val


def compute_log_likelihood(
    transcript_to_fractional_read_assignment, transcript_to_expr_val
):
    # compute log likelihood
    log_likelihood = 0
    for (
        transcript_id,
        reads_to_fracs,
    ) in transcript_to_fractional_read_assignment.items():
        transcript_expr = transcript_to_expr_val[transcript_id]
        if transcript_expr > 0:
            for read, frac in reads_to_fracs.items():
                log_likelihood += frac * log(transcript_expr)

    return log_likelihood


def get_read_name_to_transcripts(transcripts):
    read_name_to_transcripts = dict()

    for transcript in transcripts:
        read_names = transcript.get_read_names()
        for read_name in read_names:
            read_weight = transcript.get_read_weight(read_name)
            if not LRAA_Globals.config["use_weighted_read_assignments"]:
                # not using the read weights.
                read_weight = 1.0

            if read_name not in read_name_to_transcripts:
                read_name_to_transcripts[read_name] = set()

            read_name_to_transcripts[read_name].add((transcript, read_weight))

    return read_name_to_transcripts


def init_read_counts_and_expr_vals(transcripts, read_name_to_transcripts):

    ## Init - split across mapped transcripts according to normalized weights or equally if not using weights..

    transcript_to_fractional_read_assignment = defaultdict(dict)
    transcript_to_read_count = defaultdict(float)
    transcript_to_expr_val = defaultdict(float)

    num_mapped_reads = len(read_name_to_transcripts)

    for transcript in transcripts:

        transcript_id = transcript.get_transcript_id()

        transcript_read_count_total = 0
        read_names = transcript.get_read_names()
        for read_name in read_names:

            frac_read_assignment = 0
            all_transcripts_with_read = read_name_to_transcripts[read_name]

            # split read equally across all copatible reads.
            num_transcripts_with_assigned_read = len(
                read_name_to_transcripts[read_name]
            )
            frac_read_assignment = (
                1 / num_transcripts_with_assigned_read
                if num_transcripts_with_assigned_read > 0
                else 0
            )

            transcript_read_count_total += frac_read_assignment
            transcript_to_fractional_read_assignment[transcript_id][
                read_name
            ] = frac_read_assignment

        # assign initial transcript expr values based on the init read frac assignments.
        transcript_to_read_count[transcript_id] = transcript_read_count_total
        transcript_to_expr_val[transcript_id] = (
            transcript_read_count_total / num_mapped_reads
            if num_mapped_reads > 0
            else 0
        )

        logger.debug(
            f"-INIT: assigning transcript {transcript_id} read count: {transcript_read_count_total} and expr val {transcript_read_count_total}/{num_mapped_reads} = {transcript_to_expr_val[transcript_id]}"
        )

    return (
        transcript_to_read_count,
        transcript_to_fractional_read_assignment,
        transcript_to_expr_val,
    )
