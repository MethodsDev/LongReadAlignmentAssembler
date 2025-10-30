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
    prefix_str=None,
):

    if type(transcripts) != list:
        transcripts = list(transcripts)

    local_debug = False

    alpha = LRAA_Globals.config["EM_alpha"]

    num_transcripts = len(transcripts)

    ## assign reads and transcripts to indices
    multipath_to_transcripts_and_weights, transcript_id_to_idx = (
        get_multipath_to_transcripts_and_weights(transcripts)
    )

    if local_debug:
        print(
            "\n".join(
                [
                    "mp to transcripts and weights:\n",
                    "\n".join(
                        [
                            f"{x}:\t{y}"
                            for (x, y) in multipath_to_transcripts_and_weights.items()
                        ]
                    ),
                    "\ntranscript id to idx:",
                    "\n".join(f"{x}\t{y}" for (x, y) in transcript_id_to_idx.items()),
                ]
            )
            + "\n\n"
        )

    multipaths = list(multipath_to_transcripts_and_weights.keys())
    num_reads_mapped = 0
    for mp in multipaths:
        num_reads_mapped += mp.get_read_count()

    try:
        if prefix_str:
            logger.info(
                f"{prefix_str}Running EM for {num_transcripts} transcripts with {num_reads_mapped} mapped reads."
            )
        else:
            logger.info(
                f"Running EM for {num_transcripts} transcripts with {num_reads_mapped} mapped reads."
            )
    except Exception:
        logger.info(
            f"Running EM for {num_transcripts} transcripts with {num_reads_mapped} mapped reads."
        )

    # inputs to EM
    mp_assignments = list()  # list of lists of transcript indices
    mp_weights = list()  # list of lists of weights aligned with mp_assignments
    mp_read_counts = list()  # number of reads assigned to each mp

    # populate inputs to EM
    start_prep_time = time.time()
    for i, mp in enumerate(multipaths):

        indiv_mp_to_trans_assignments = list()
        indiv_mp_to_trans_weights = list()

        num_reads_in_mp = mp.get_read_count()
        mp_read_counts.append(num_reads_in_mp)

        for (
            mp_mapped_transcript,
            mp_trans_weight_val,
        ) in multipath_to_transcripts_and_weights[mp]:
            mp_mapped_transcript_idx = transcript_id_to_idx[mp_mapped_transcript]

            indiv_mp_to_trans_assignments.append(mp_mapped_transcript_idx)

            if not LRAA_Globals.config["use_weighted_read_assignments"]:
                mp_trans_weight_val = 1.0

            indiv_mp_to_trans_weights.append(mp_trans_weight_val)

        mp_assignments.append(indiv_mp_to_trans_assignments)
        mp_weights.append(indiv_mp_to_trans_weights)

    end_prep_time = time.time()

    EM_prep_time = end_prep_time - start_prep_time
    try:
        if prefix_str:
            logger.info(f"{prefix_str}EM_prep_time: {EM_prep_time / 60:.2f} minutes")
        else:
            logger.info("EM_prep_time: {:.2f} minutes".format(EM_prep_time / 60))
    except Exception:
        logger.info("EM_prep_time: {:.2f} minutes".format(EM_prep_time / 60))

    if local_debug:
        mp_weights_for_printing = list()
        for _list in mp_weights:
            mp_weights_for_printing.append([f"{x:.3f}" for x in _list])

        trans_assignments_and_weights = zip(mp_assignments, mp_weights_for_printing)

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
        fractional_mp_assignments_array,
    ) = em_algorithm_with_weights(
        mp_assignments,
        mp_weights,
        mp_read_counts,
        num_transcripts,
        max_iter=max_EM_iterations,
        base_alpha=alpha,
        prefix_str=prefix_str,
    )

    if local_debug:

        fractional_mp_assignments_array_for_printing = list()
        for _list in fractional_mp_assignments_array:
            fractional_mp_assignments_array_for_printing.append(
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
                        [str(x) for x in fractional_mp_assignments_array_for_printing]
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
    transcript_to_fractional_mp_assignment = defaultdict(dict)
    for mp_idx, trans_frac_assignments_array in enumerate(
        fractional_mp_assignments_array
    ):
        mp = multipaths[mp_idx]
        for j, frac_val in enumerate(trans_frac_assignments_array):
            transcript_idx = mp_assignments[mp_idx][j]
            transcript_id = transcripts[transcript_idx].get_transcript_id()
            transcript_to_fractional_mp_assignment[transcript_id][mp] = frac_val

    return (
        transcript_to_expr_val,
        transcript_to_fractional_mp_assignment,
        transcript_to_read_count,
    )


def get_multipath_to_transcripts_and_weights(transcripts):

    multipath_to_transcripts_and_weights = dict()
    transcript_id_to_idx = dict()
    all_multipaths = set()

    for i, transcript in enumerate(transcripts):
        transcript_id = transcript.get_transcript_id()
        transcript_id_to_idx[transcript_id] = i

        multipaths = transcript.get_multipaths_evidence_assigned()
        for mp in multipaths:

            mp_weight = transcript.get_multipath_weight(mp)
            if not LRAA_Globals.config["use_weighted_read_assignments"]:
                # not using the read weights.
                mp_weight = 1.0

            if mp not in multipath_to_transcripts_and_weights:
                multipath_to_transcripts_and_weights[mp] = list()

            multipath_to_transcripts_and_weights[mp].append((transcript_id, mp_weight))

    # order transcripts according to their indices
    for (
        mp,
        transcript_assigned_list,
    ) in multipath_to_transcripts_and_weights.items():
        transcript_assigned_list.sort(key=lambda x: transcript_id_to_idx[x[0]])

    return multipath_to_transcripts_and_weights, transcript_id_to_idx


def em_algorithm_with_weights(
    mp_assignments,
    mp_weights,
    mp_read_counts,
    num_transcripts,
    max_iter=100,
    tol=1e-6,
    base_alpha=0.1,
    prefix_str=None,
):
    """
    Perform the EM algorithm to estimate transcript expression levels with weighted multipaths.

    Parameters:
        mp_assignments (list of lists): Each element is a list of transcript indices (0-based)
                                          to which the mp is assigned (unique or ambiguous).
        mp_weights (list of lists): A matrix of weights where each sublist corresponds to
                                       the weights of a mp for the assigned transcripts.
        mp_read_counts (list of ints): each element is the number of reads corresponding to that multipath.
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

    # init fractional mp assignments
    fractional_mp_assignments = init_fractional_mp_assignments(mp_assignments)

    # Count ambiguous mps for each transcript
    ambiguous_read_counts = np.zeros(num_transcripts)
    for i, mp in enumerate(mp_assignments):
        if len(mp) > 1:  # Ambiguous mp
            for trans_id in mp:
                ambiguous_read_counts[trans_id] += mp_read_counts[i]

    # Calculate transcript-specific alpha values based on ambiguous mps
    transcript_alphas = base_alpha * ambiguous_read_counts

    transcript_sum_read_counts = defaultdict(float)

    for iteration in range(max_iter):
        # E-step: Calculate fractional assignments
        transcript_sum_read_counts.clear()
        for mp_i, mp_mapped_transcripts in enumerate(mp_assignments):

            # denominator for fractional assignment is the mp-weighted sum of expression for assigned transcripts
            total_weight = sum(
                mp_weights[mp_i][j] * transcript_expression_levels[trans_id]
                for j, trans_id in enumerate(mp_mapped_transcripts)
            )
            for j, trans_id in enumerate(mp_mapped_transcripts):
                weight = mp_weights[mp_i][j]

                # for each transcript this mp is assigned,
                # assign a proportion of this mp according to its relative expression contribution.
                frac_assignment = (
                    weight * transcript_expression_levels[trans_id] / total_weight
                    if total_weight > 0
                    else 0
                )

                transcript_sum_read_counts[trans_id] += (
                    frac_assignment * mp_read_counts[mp_i]
                )

                fractional_mp_assignments[mp_i][j] = frac_assignment

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
            try:
                if prefix_str:
                    logger.info(
                        f"{prefix_str}Converged after {iteration + 1} iterations. (time={runtime:.3f} sec. for {num_transcripts} transcripts)"
                    )
                else:
                    logger.info(
                        f"Converged after {iteration + 1} iterations. (time={runtime:.3f} sec. for {num_transcripts} transcripts)"
                    )
            except Exception:
                logger.info(
                    f"Converged after {iteration + 1} iterations. (time={runtime:.3f} sec. for {num_transcripts} transcripts)"
                )
            break

        prev_expression_levels = transcript_expression_levels.copy()

    return (
        transcript_expression_levels,
        transcript_sum_read_counts,
        fractional_mp_assignments,
    )


def init_fractional_mp_assignments(template_list_of_lists):

    init_frac_mp_assignments = list()
    for _list in template_list_of_lists:
        zeroed_list = np.zeros(len(_list))
        init_frac_mp_assignments.append(zeroed_list)

    return init_frac_mp_assignments


def test_run_EM():

    ######################################################################################################
    # example has 3 transcripts and 5 mps. Some mps map ambiguously and others are unique to isoforms.
    ######################################################################################################

    contig_acc = "fake_contig"
    contig_strand = "+"

    transcript_0_coords = [[100, 200], [300, 400], [500, 900]]
    transcript_0 = Transcript(contig_acc, transcript_0_coords, contig_strand)

    transcript_1_coords = [[100, 200], [250, 400], [500, 750]]
    transcript_1 = Transcript(contig_acc, transcript_1_coords, contig_strand)

    transcript_2_coords = [[100, 200], [500, 900]]
    transcript_2 = Transcript(contig_acc, transcript_2_coords, contig_strand)

    # 5 mps
    mp_0 = "mp_0"
    mp_1 = "mp_1"
    mp_2 = "mp_2"
    mp_3 = "mp_3"
    mp_4 = "mp_4"

    transcript_0.add_mp_names([mp_0, mp_2])
    transcript_0.set_mp_weights({mp_0: 1.0, mp_2: 0.5})

    transcript_1.add_mp_names([mp_1, mp_2, mp_4])
    transcript_1.set_mp_weights({mp_1: 1.0, mp_2: 0.3, mp_4: 0.4})

    transcript_2.add_mp_names([mp_2, mp_3, mp_4])
    transcript_2.set_mp_weights({mp_2: 0.2, mp_3: 1.0, mp_4: 0.6})

    transcripts = [transcript_0, transcript_1, transcript_2]

    (
        transcript_to_expr_val,
        transcript_to_fractional_mp_assignment,
        transcript_to_mp_count,
    ) = run_EM(transcripts, 1000)

    print("#LRAA-interfaced results:\n\n")
    print(
        "Estimated expression levels, LRAA interfaced:\n",
        "\n".join([f"\t{x}: {y:.4f}" for (x, y) in transcript_to_expr_val.items()]),
        "\n\n",
    )
    print(
        "Transcript mp counts:\n",
        "\n".join([f"\t{x}: {y:.1f}" for (x, y) in transcript_to_mp_count.items()]),
        "\n\n",
    )
    print("Fractional mp assignments:\n")

    for (
        transcript,
        frac_mp_assignments,
    ) in transcript_to_fractional_mp_assignment.items():
        print(transcript),
        print("\n".join([f"\t{x}: {y:.4f}" for (x, y) in frac_mp_assignments.items()]))


def test_em_algorithm_with_weights():

    # Simulated mp assignments
    # Each sublist contains the indices of transcripts a mp is assigned to
    mp_assignments = [
        [0],  # mp_0: Unique assignment to transcript 0
        [1],  # mp_1: Unique assignment to transcript 1
        [0, 1, 2],  # mp_2:  Ambiguous assignment to transcripts 0, 1, 2
        [2],  # mp_3: Unique assignment to transcript 2
        [1, 2],  # mp_4: Ambiguous assignment to transcripts 1 and 2
    ]

    # Simulated weights for each mp (rows correspond to mps, columns to assigned transcripts)
    mp_weights = [
        [1.0],  # mp_0: Unique mp, full weight to transcript 0
        [1.0],  # mp_1: Unique mp, full weight to transcript 1
        [0.5, 0.3, 0.2],  # mp_2: Ambiguous mp, weights for transcripts 0, 1, 2
        [1.0],  # mp_3: Unique mp, full weight to transcript 2
        [0.4, 0.6],  # mp_4: Ambiguous mp, weights for transcripts 1, 2
    ]

    mp_read_counts = [1, 2, 3, 4, 5]

    # run EM
    num_transcripts = 3
    expression, trans_read_counts, frac_mp_assignments = em_algorithm_with_weights(
        mp_assignments, mp_weights, mp_read_counts, num_transcripts, base_alpha=0
    )

    print("\n\n" + "##########################\n" + "#:\n")

    print(
        "Estimated expression levels:\n\t"
        + "\n\t".join([f"{x:.4f}" for x in expression]),
    )

    print(
        "\nTrans mp counts:\n"
        + "\n".join([f"\t{x}: {y:.1f}" for (x, y) in trans_read_counts.items()]),
    )

    print("\nFrac mp assignments:")

    for _list in frac_mp_assignments:
        print("\t" + "\t".join([f"{x:.4f}" for x in _list]))


# Example usage:
if __name__ == "__main__":
    # test_run_EM()
    test_em_algorithm_with_weights()
