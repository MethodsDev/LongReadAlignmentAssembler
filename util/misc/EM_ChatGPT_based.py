#!/usr/bin/env python3

import numpy as np
from collections import defaultdict
import sys, os, re


import logging

logger = logging.getLogger(__name__)


# original from ChatGPT, modified by bhaas


def em_algorithm_with_weights(
    read_assignments, read_weights, num_transcripts, max_iter=100, tol=1e-6
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

    Returns:
        np.ndarray: Estimated expression levels for each transcript.
    """

    # Initialize expression values uniformly
    transcript_expression_levels = np.ones(num_transcripts) / num_transcripts
    prev_expression_levels = np.zeros(num_transcripts)

    # init fractional read assignments
    fractional_read_assignments = init_fractional_read_assignments(read_assignments)

    transcript_sum_read_counts = defaultdict(float)

    for iteration in range(max_iter):
        print("Iteration: {}".format(iteration))
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
                transcript_sum_read_counts[trans_id]
                for trans_id in range(num_transcripts)
            ]
        )

        # Normalize to ensure expression levels sum to 1
        transcript_expression_levels /= transcript_expression_levels.sum()

        # Check for convergence
        if np.linalg.norm(transcript_expression_levels - prev_expression_levels) < tol:
            print(f"Converged after {iteration + 1} iterations.")
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

    # Create 586 reads that map ambiguously to transcript_A and transcript_B
    # transcript_A has index 0
    # transcript_B has index 1

    read_assignments = list()
    read_weights = list()
    for i in range(586):
        read_assignments.append([0, 1])
        read_weights.append([0.5, 0.5])

    # now add one read that maps uniquely to transcript_A
    read_assignments.append([0])
    read_weights.append([1])

    # run the EM
    num_transcripts = 2
    expression, trans_read_counts, frac_read_assignments = em_algorithm_with_weights(
        read_assignments, read_weights, num_transcripts, 250
    )

    print("\n\n" + "##########################\n" + "# ChatGPT-direct interface:\n")

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
