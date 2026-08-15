#!/usr/bin/env python3

"""EM must return fractions and counts computed at the theta it returns.

The loop breaks AFTER its M-step. Without a final E-step, the fractions and counts handed back
came from the previous iteration's abundances while the abundances handed back are the updated
ones -- two different estimates in one return value, with nothing in the output distinguishing
them. Anyone recombining them is mixing estimates: a streaming pass splitting a read it must
resolve itself, or a reader comparing a tracking fraction against the abundance in quant.expr.

The size of the discrepancy is bounded by EM_convergence_tol only when EM converged. On a gene
that exhausts max_iter the two thetas can differ arbitrarily -- and those are the genes with
many ambiguous reads, where the split is hardest and matters most.

These tests pin the invariant rather than any particular number:

    frac[mp][j] == w[mp][j] * theta[t] / sum_over_j(w[mp][j] * theta[t])
    counts[t]   == sum_over_mp(frac[mp][j] * read_count[mp])

both evaluated at the RETURNED theta.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

from EM import em_algorithm_with_weights


def _check_consistency(mp_assignments, mp_weights, mp_read_counts, num_transcripts,
                       max_iter=250, tol=1e-6):
    theta, sums, fracs = em_algorithm_with_weights(
        mp_assignments, mp_weights, mp_read_counts, num_transcripts,
        max_iter=max_iter, tol=tol,
    )

    for mp_i, mapped in enumerate(mp_assignments):
        denom = sum(
            mp_weights[mp_i][j] * theta[t] for j, t in enumerate(mapped)
        )
        for j, t in enumerate(mapped):
            expected = (
                mp_weights[mp_i][j] * theta[t] / denom if denom > 0 else 0.0
            )
            assert fracs[mp_i][j] == pytest.approx(expected, abs=1e-15), (
                f"mp {mp_i} transcript {t}: fraction {fracs[mp_i][j]} was not computed at "
                f"the returned theta (expected {expected})"
            )

    recomputed = {}
    for mp_i, mapped in enumerate(mp_assignments):
        for j, t in enumerate(mapped):
            recomputed[t] = recomputed.get(t, 0.0) + fracs[mp_i][j] * mp_read_counts[mp_i]
    for t, value in recomputed.items():
        assert sums[t] == pytest.approx(value, abs=1e-12), (
            f"transcript {t}: returned count {sums[t]} is not the sum of the returned "
            f"fractions ({value})"
        )
    return theta, sums, fracs


def test_converged_run_returns_fractions_matching_returned_theta():
    """The ordinary case: converges well inside the iteration cap."""
    mp_assignments = [[0], [1], [2], [0, 1], [0, 1, 2]]
    mp_weights = [[1.0], [1.0], [1.0], [1.0, 1.0], [1.0, 1.0, 1.0]]
    counts = [500.0, 120.0, 40.0, 300.0, 900.0]
    theta, _sums, _fracs = _check_consistency(mp_assignments, mp_weights, counts, 3)
    assert theta.sum() == pytest.approx(1.0)


def test_non_converged_run_returns_fractions_matching_returned_theta():
    """The case that motivated this: max_iter exhausted, so the two thetas can differ freely.

    A single iteration guarantees the loop exits without ever meeting the tolerance, which is
    the regime where reusing the previous theta's fractions is unbounded rather than tiny.
    """
    mp_assignments = [[0, 1], [0, 1], [1, 2], [0, 2]]
    mp_weights = [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0]]
    counts = [10.0, 10.0, 10.0, 10.0]
    _check_consistency(mp_assignments, mp_weights, counts, 3, max_iter=1)


def test_weighted_run_returns_fractions_matching_returned_theta():
    """3'-agreement weights are on by default, so the weighted form must hold too."""
    mp_assignments = [[0, 1], [0, 1, 2], [2]]
    mp_weights = [[0.9, 0.1], [0.5, 0.3, 0.2], [1.0]]
    counts = [250.0, 400.0, 75.0]
    _check_consistency(mp_assignments, mp_weights, counts, 3)


def test_zero_abundance_group_yields_zero_fractions_not_a_crash():
    """A multipath whose every candidate has zero abundance keeps the documented behaviour.

    Fractions come out 0.0 rather than raising or redistributing, which is what the streaming
    path deliberately mirrors. Pinned so the final E-step does not quietly change it.
    """
    mp_assignments = [[0], [1]]
    mp_weights = [[1.0], [0.0]]
    counts = [100.0, 0.0]
    theta, sums, fracs = em_algorithm_with_weights(
        mp_assignments, mp_weights, counts, 2, max_iter=50, tol=1e-6
    )
    assert np.isfinite(theta).all(), "theta must stay finite"
    for mp_i, mapped in enumerate(mp_assignments):
        for j, _t in enumerate(mapped):
            assert np.isfinite(fracs[mp_i][j])
            assert fracs[mp_i][j] >= 0.0


def test_final_estep_does_not_move_theta():
    """The final E-step must change only the REPORTED quantities, not the estimator.

    Theta is not re-updated after the loop, so abundances are identical to what the same input
    produced before this change. Pinned with golden values because the consistency tests above
    check fractions against *whatever* theta comes back -- a "fix" that also ran one extra
    M-step would satisfy every one of them while silently shifting everyone's abundances.

    These numbers were captured from the implementation and are asserted tightly on purpose: if
    they move, either the estimator changed or the arithmetic did, and both deserve a failing
    test rather than a quiet diff in every user's quant.expr.
    """
    mp_assignments = [[0], [1], [2], [0, 1], [0, 1, 2]]
    mp_weights = [[1.0], [1.0], [1.0], [1.0, 1.0], [1.0, 1.0, 1.0]]
    counts = [500.0, 120.0, 40.0, 300.0, 900.0]
    theta, _sums, _fracs = em_algorithm_with_weights(
        mp_assignments, mp_weights, counts, 3, max_iter=250, tol=1e-6
    )
    expected = [0.7557989468786793, 0.19485563493769187, 0.04934541818362883]
    for i, want in enumerate(expected):
        assert theta[i] == pytest.approx(want, rel=1e-12), (
            f"transcript {i}: theta moved from {want} to {theta[i]}; the final E-step is "
            "supposed to leave abundances alone"
        )


def test_read_mass_is_conserved_regardless_of_theta():
    """Counts must sum to the input read total, whatever theta the fractions were computed at.

    Per multipath the fractions sum to 1 whenever the denominator is positive, so this holds
    for any abundance vector -- which is exactly why the final E-step cannot change total read
    mass, and why the only movement observed in quant.expr on real data was a single value
    crossing a one-decimal rounding boundary.
    """
    mp_assignments = [[0, 1], [1, 2], [0, 1, 2], [2]]
    mp_weights = [[0.7, 0.3], [1.0, 1.0], [0.5, 0.25, 0.25], [1.0]]
    counts = [321.0, 88.0, 654.0, 17.0]
    _theta, sums, _fracs = em_algorithm_with_weights(
        mp_assignments, mp_weights, counts, 3, max_iter=250, tol=1e-6
    )
    assert sum(sums.values()) == pytest.approx(sum(counts), rel=1e-12), (
        f"read mass not conserved: {sum(sums.values())} vs input {sum(counts)}"
    )
