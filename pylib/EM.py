#!/usr/bin/env python3

import numpy as np
from collections import defaultdict
import sys, os, re
from Transcript import Transcript
from collections import defaultdict
import LRAA_Globals
import logging
import time

logger = logging.getLogger(__name__)


def run_EM(
    transcripts: list,
    max_EM_iterations: int = 1000,
):

    if type(transcripts) != list:
        transcripts = list(transcripts)

    local_debug = False

    alpha = LRAA_Globals.config["EM_alpha"]

    num_transcripts = len(transcripts)

    ## assign reads and transcripts to indices
    read_name_to_transcripts_and_weights, transcript_id_to_idx = (
        get_read_to_transcripts_and_weights(transcripts)
    )

    if local_debug:
        print(
            "\n".join(
                [
                    "read name to transcripts and weights:\n",
                    "\n".join(
                        [
                            f"{x}:\t{y}"
                            for (x, y) in read_name_to_transcripts_and_weights.items()
                        ]
                    ),
                    "\ntranscript id to idx:",
                    "\n".join(f"{x}\t{y}" for (x, y) in transcript_id_to_idx.items()),
                ]
            )
            + "\n\n"
        )

    read_names = sorted(list(read_name_to_transcripts_and_weights.keys()))

    num_reads_mapped = len(read_names)

    logger.info(
        f"Running EM for {num_transcripts} transcripts with {num_reads_mapped} mapped reads."
    )

    # inputs to EM
    read_assignments = list()
    read_weights = list()

    # populate inputs to EM
    start_prep_time = time.time()
    for i, read_name in enumerate(read_names):

        indiv_read_to_trans_assignments = list()
        indiv_read_to_trans_weights = list()

        for (
            read_mapped_transcript,
            read_trans_weight_val,
        ) in read_name_to_transcripts_and_weights[read_name]:
            read_mapped_transcript_idx = transcript_id_to_idx[read_mapped_transcript]

            indiv_read_to_trans_assignments.append(read_mapped_transcript_idx)

            if not LRAA_Globals.config["use_weighted_read_assignments"]:
                read_trans_weight_val = 1.0

            indiv_read_to_trans_weights.append(read_trans_weight_val)

        read_assignments.append(indiv_read_to_trans_assignments)
        read_weights.append(indiv_read_to_trans_weights)

    end_prep_time = time.time()

    EM_prep_time = end_prep_time - start_prep_time
    logger.info("EM_prep_time: {:.2f} minutes".format(EM_prep_time / 60))

    if local_debug:

        read_weights_for_printing = list()
        for _list in read_weights:
            read_weights_for_printing.append([f"{x:.3f}" for x in _list])

        trans_assignments_and_weights = zip(read_assignments, read_weights_for_printing)

        print(
            "\n".join(
                [
                    "LRAA interface inputs:",
                    "read_assignments and weights:\n"
                    + "\n".join(
                        [
                            f"r[{x}]-Ts: {y}\nr[{x}]-ws {z}\n"
                            for x, (y, z) in enumerate(trans_assignments_and_weights)
                        ]
                    )
                    + "\n\n",
                ]
            )
        )

    #########
    ## Run EM
    #########

    (
        trans_expr_levels_array,
        transcript_sum_read_counts_array,
        fractional_read_assignments_array,
    ) = em_algorithm_with_weights(
        read_assignments,
        read_weights,
        num_transcripts,
        max_iter=max_EM_iterations,
        base_alpha=alpha,
    )

    if local_debug:

        fractional_read_assignments_array_for_printing = list()
        for _list in fractional_read_assignments_array:
            fractional_read_assignments_array_for_printing.append(
                [f"{x:.4f}" for x in _list]
            )

        print(
            "\n".join(
                [
                    "\n\n# LRAA out from EM method:",
                    "trans_expr_levels_array:\n"
                    + "\n".join([f"{x:.3f}" for x in trans_expr_levels_array]),
                    "\ntranscript_sum_read_counts_array:\n"
                    + "\n".join([str(x) for x in transcript_sum_read_counts_array]),
                    "\nfractional_read_assignments_array: "
                    + "\n".join(
                        [str(x) for x in fractional_read_assignments_array_for_printing]
                    ),
                ]
            )
        )

    # interface results with LRAA data structures.

    # first interface with transcript expr results
    transcript_to_read_count = defaultdict(float)
    transcript_to_expr_val = defaultdict(float)

    for i, trans_expr_frac_val in enumerate(trans_expr_levels_array):
        transcript_id = transcripts[i].get_transcript_id()
        transcript_to_expr_val[transcript_id] = trans_expr_frac_val
        transcript_to_read_count[transcript_id] = transcript_sum_read_counts_array[i]

    # next, interface with fractional read assignment info
    transcript_to_fractional_read_assignment = defaultdict(dict)
    for read_idx, trans_frac_assignments_array in enumerate(
        fractional_read_assignments_array
    ):
        read_name = read_names[read_idx]
        for j, frac_val in enumerate(trans_frac_assignments_array):
            transcript_idx = read_assignments[read_idx][j]
            transcript_id = transcripts[transcript_idx].get_transcript_id()
            transcript_to_fractional_read_assignment[transcript_id][
                read_name
            ] = frac_val

    return (
        transcript_to_expr_val,
        transcript_to_fractional_read_assignment,
        transcript_to_read_count,
    )


def get_read_to_transcripts_and_weights(transcripts):

    read_name_to_transcripts_and_weights = dict()
    transcript_id_to_idx = dict()
    all_read_names = set()

    for i, transcript in enumerate(transcripts):
        transcript_id = transcript.get_transcript_id()
        transcript_id_to_idx[transcript_id] = i

        read_names = transcript.get_read_names()
        for read_name in read_names:

            read_weight = transcript.get_read_weight(read_name)
            if not LRAA_Globals.config["use_weighted_read_assignments"]:
                # not using the read weights.
                read_weight = 1.0

            if read_name not in read_name_to_transcripts_and_weights:
                read_name_to_transcripts_and_weights[read_name] = list()

            read_name_to_transcripts_and_weights[read_name].append(
                (transcript_id, read_weight)
            )

    # order transcripts according to their indices
    for (
        read_name,
        transcript_assigned_list,
    ) in read_name_to_transcripts_and_weights.items():
        transcript_assigned_list.sort(key=lambda x: transcript_id_to_idx[x[0]])

    return read_name_to_transcripts_and_weights, transcript_id_to_idx


def em_algorithm_with_weights(
    read_assignments,
    read_weights,
    num_transcripts,
    max_iter=100,
    tol=1e-6,
    base_alpha=0.1,
):
    """
    Perform the EM algorithm to estimate transcript expression levels with weighted reads.

    Parameters:
        read_assignments (list of lists): Each element is a list of transcript indices (0-based)
                                          to which the read is assigned (unique or ambiguous).
        read_weights (list of lists): A matrix of weights where each sublist corresponds to
                                       the weights of a read for the assigned transcripts.
        num_transcripts (int): Total number of transcripts.
        max_iter (int): Maximum number of iterations.
        tol (float): Convergence tolerance.
        base_alpha (float): Base regularization parameter to scale for each transcript.

    Returns:
        np.ndarray: Estimated expression levels for each transcript.
    """

    time_start = time.time()

    # Initialize expression values uniformly
    transcript_expression_levels = np.ones(num_transcripts) / num_transcripts
    prev_expression_levels = np.zeros(num_transcripts)

    # init fractional read assignments
    fractional_read_assignments = init_fractional_read_assignments(read_assignments)

    # Count ambiguous reads for each transcript
    ambiguous_read_counts = np.zeros(num_transcripts)
    for read in read_assignments:
        if len(read) > 1:  # Ambiguous read
            for trans_id in read:
                ambiguous_read_counts[trans_id] += 1

    # Calculate transcript-specific alpha values based on ambiguous reads
    transcript_alphas = base_alpha * ambiguous_read_counts

    transcript_sum_read_counts = defaultdict(float)

    for iteration in range(max_iter):
        # E-step: Calculate fractional assignments
        transcript_sum_read_counts.clear()
        for read_i, read_mapped_transcripts in enumerate(read_assignments):

            # denominator for fractional assignment is the read-weighted sum of expression for assigned transcripts
            total_weight = sum(
                read_weights[read_i][j] * transcript_expression_levels[trans_id]
                for j, trans_id in enumerate(read_mapped_transcripts)
            )
            for j, trans_id in enumerate(read_mapped_transcripts):
                weight = read_weights[read_i][j]

                # for each transcript this read is assigned,
                # assign a proportion of this read according to its relative expression contribution.
                frac_assignment = (
                    weight * transcript_expression_levels[trans_id] / total_weight
                    if total_weight > 0
                    else 0
                )

                transcript_sum_read_counts[trans_id] += frac_assignment

                fractional_read_assignments[read_i][j] = frac_assignment

        # M-step: Update expression levels
        transcript_expression_levels = np.array(
            [
                transcript_sum_read_counts[trans_id] + transcript_alphas[trans_id]
                for trans_id in range(num_transcripts)
            ]
        )

        # Normalize to ensure expression levels sum to 1
        transcript_expression_levels /= transcript_expression_levels.sum()

        # Check for convergence
        if np.linalg.norm(transcript_expression_levels - prev_expression_levels) < tol:
            runtime = time.time() - time_start
            logger.info(
                f"Converged after {iteration + 1} iterations. (time={runtime:.3f} sec. for {num_transcripts} transcripts)"
            )
            break

        prev_expression_levels = transcript_expression_levels.copy()

    return (
        transcript_expression_levels,
        transcript_sum_read_counts,
        fractional_read_assignments,
    )


def init_fractional_read_assignments(template_list_of_lists):

    init_frac_read_assignments = list()
    for _list in template_list_of_lists:
        zeroed_list = np.zeros(len(_list))
        init_frac_read_assignments.append(zeroed_list)

    return init_frac_read_assignments


# Example usage:
if __name__ == "__main__":

    ######################################################################################################
    # example has 3 transcripts and 5 reads. Some reads map ambiguously and others are unique to isoforms.
    ######################################################################################################

    contig_acc = "fake_contig"
    contig_strand = "+"

    transcript_0_coords = [[100, 200], [300, 400], [500, 900]]
    transcript_0 = Transcript(contig_acc, transcript_0_coords, contig_strand)

    transcript_1_coords = [[100, 200], [250, 400], [500, 750]]
    transcript_1 = Transcript(contig_acc, transcript_1_coords, contig_strand)

    transcript_2_coords = [[100, 200], [500, 900]]
    transcript_2 = Transcript(contig_acc, transcript_2_coords, contig_strand)

    # 5 reads
    read_0 = "read_0"
    read_1 = "read_1"
    read_2 = "read_2"
    read_3 = "read_3"
    read_4 = "read_4"

    transcript_0.add_read_names([read_0, read_2])
    transcript_0.set_read_weights({read_0: 1.0, read_2: 0.5})

    transcript_1.add_read_names([read_1, read_2, read_4])
    transcript_1.set_read_weights({read_1: 1.0, read_2: 0.3, read_4: 0.4})

    transcript_2.add_read_names([read_2, read_3, read_4])
    transcript_2.set_read_weights({read_2: 0.2, read_3: 1.0, read_4: 0.6})

    transcripts = [transcript_0, transcript_1, transcript_2]

    (
        transcript_to_read_count,
        transcript_to_fractional_read_assignment,
        transcript_to_expr_val,
    ) = run_EM(transcripts)

    print("#LRAA-interfaced results:\n\n")
    print(
        "Estimated expression levels, LRAA interfaced:\n",
        "\n".join([f"\t{x}: {y:.4f}" for (x, y) in transcript_to_expr_val.items()]),
        "\n\n",
    )
    print(
        "Transcript read counts:\n",
        "\n".join([f"\t{x}: {y:.1f}" for (x, y) in transcript_to_read_count.items()]),
        "\n\n",
    )
    print("Fractional read assignments:\n")

    for (
        transcript,
        frac_read_assignments,
    ) in transcript_to_fractional_read_assignment.items():
        print(transcript),
        print(
            "\n".join([f"\t{x}: {y:.4f}" for (x, y) in frac_read_assignments.items()])
        )

    # Simulated read assignments
    # Each sublist contains the indices of transcripts a read is assigned to
    read_assignments = [
        [0],  # read_0: Unique assignment to transcript 0
        [1],  # read_1: Unique assignment to transcript 1
        [0, 1, 2],  # read_2:  Ambiguous assignment to transcripts 0, 1, 2
        [2],  # read_3: Unique assignment to transcript 2
        [1, 2],  # read_4: Ambiguous assignment to transcripts 1 and 2
    ]

    # Simulated weights for each read (rows correspond to reads, columns to assigned transcripts)
    read_weights = [
        [1.0],  # read_0: Unique read, full weight to transcript 0
        [1.0],  # read_1: Unique read, full weight to transcript 1
        [0.5, 0.3, 0.2],  # read_2: Ambiguous read, weights for transcripts 0, 1, 2
        [1.0],  # read_3: Unique read, full weight to transcript 2
        [0.4, 0.6],  # read_4: Ambiguous read, weights for transcripts 1, 2
    ]

    # run EM
    num_transcripts = 3
    expression, trans_read_counts, frac_read_assignments = em_algorithm_with_weights(
        read_assignments, read_weights, num_transcripts, base_alpha=0
    )

    print("\n\n" + "##########################\n" + "#:\n")

    print(
        "Estimated expression levels:\n\t"
        + "\n\t".join([f"{x:.4f}" for x in expression]),
    )

    print(
        "\nTrans read counts:\n"
        + "\n".join([f"\t{x}: {y:.1f}" for (x, y) in trans_read_counts.items()]),
    )

    print("\nFrac read assignments:")

    for _list in frac_read_assignments:
        print("\t" + "\t".join([f"{x:.4f}" for x in _list]))
