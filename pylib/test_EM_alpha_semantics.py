#!/usr/bin/env python3

"""What EM_alpha actually is, pinned so the sweep's recommendation stays legible.

`EM_alpha` reads like a sparsity penalty and is documented as a "regularization
factor". It is the opposite. `em_algorithm_with_weights` forms

    transcript_alphas = base_alpha * ambiguous_read_counts        (EM.py)
    theta_t          proportional to  counts_t + transcript_alphas[t]

so alpha ADDS mass, and it adds it in proportion to how ambiguous a transcript
is. Raising alpha spreads abundance across indistinguishable isoforms; lowering
it lets EM collapse toward winner-take-all. Anyone who reads "regularization" as
L1 will tune it backwards.

Three properties are pinned here, all of them load bearing for the tuning work:

1. The pseudocount is proportional to AMBIGUOUS support only. A transcript
   supported solely by unique reads gets no pseudocount at all, so alpha cannot
   inflate it. A flat-pseudocount bug would break this.

2. Raising alpha moves indistinguishable isoforms TOWARD each other, and alpha=0
   lets the winner take the most. This is the direction of the effect, and a
   sign error in the M-step would break it while leaving EM converging happily.

3. alpha is DIMENSIONLESS: scaling every read count by k scales both the counts
   term and the alpha term by k, so the normalized abundances are unchanged.
   This is what licenses transferring one alpha default across libraries of
   different depth. It is asserted here on exact powers of two, where the
   scaling is exact in floating point, and to tolerance on a non-power of two.

These are properties of the estimator, not of any corpus, so they hold whatever
the SIRV/MORF sweep concludes about the best numeric value.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

from EM import em_algorithm_with_weights  # noqa: E402


def _two_indistinguishable_plus_one_unique(n_ambig=100.0, n_uniq_a=10.0, n_uniq_c=50.0):
    """Transcripts 0 and 1 share an ambiguous multipath; 2 is unique-only.

    mp0: ambiguous over {0, 1}      -- the mass alpha redistributes
    mp1: unique to 0                -- gives 0 a real edge over 1
    mp2: unique to 2                -- untouched by alpha, the control
    """
    mp_assignments = [[0, 1], [0], [2]]
    mp_weights = [[1.0, 1.0], [1.0], [1.0]]
    mp_read_counts = [n_ambig, n_uniq_a, n_uniq_c]
    return mp_assignments, mp_weights, mp_read_counts, 3


def _run(alpha, counts=None, **kw):
    mp_assignments, mp_weights, mp_read_counts, n = _two_indistinguishable_plus_one_unique(**kw)
    if counts is not None:
        mp_read_counts = counts
    theta, sums, _fracs = em_algorithm_with_weights(
        mp_assignments, mp_weights, mp_read_counts, n, base_alpha=alpha
    )
    return theta, sums


def test_pseudocount_is_proportional_to_ambiguous_support_only():
    """A transcript with no ambiguous reads is untouched by alpha.

    Transcript 2 carries only a unique multipath, so its ambiguous_read_count is
    0 and its pseudocount is 0 * alpha == 0 at every alpha. Its EM read count
    must therefore be exactly its unique support, whatever alpha does elsewhere.
    """
    for alpha in (0.0, 0.01, 1.0, 100.0):
        _theta, sums = _run(alpha)
        assert sums[2] == pytest.approx(50.0, rel=1e-12), f"alpha={alpha}"


def test_raising_alpha_spreads_mass_across_indistinguishable_isoforms():
    """alpha is anti-sparsity: more of it, less separation between 0 and 1.

    Transcript 0 wins on unique support, so it should always hold more of the
    ambiguous mass than 1. The gap is what alpha closes.
    """
    gaps = []
    for alpha in (0.0, 0.001, 0.01, 0.1, 1.0):
        theta, _sums = _run(alpha)
        assert theta[0] > theta[1], f"unique support must still win at alpha={alpha}"
        gaps.append(theta[0] - theta[1])

    assert gaps == sorted(gaps, reverse=True), f"gap must shrink monotonically: {gaps}"
    # And the effect is large enough to matter, not a rounding artifact.
    assert gaps[0] > 2.0 * gaps[-1]


def test_alpha_zero_is_the_winner_take_most_boundary():
    """alpha=0 is the natural lower boundary, not an arbitrary grid edge.

    Below zero the pseudocount becomes a SUBTRACTION applied hardest to the most
    ambiguous transcripts, which drives the unnormalized M-step vector negative
    and normalizes garbage. Pinning the resulting sign flip documents why the
    alpha sweep stops at 0 rather than extending below it.
    """
    theta_zero, _ = _run(0.0)
    gap_zero = theta_zero[0] - theta_zero[1]
    theta_small, _ = _run(1e-6)
    assert gap_zero >= theta_small[0] - theta_small[1]

    theta_neg, _ = _run(-0.01)
    assert np.any(theta_neg < 0.0), (
        "a negative alpha must be visibly broken, not quietly plausible; "
        f"got {theta_neg}"
    )


@pytest.mark.parametrize("k", [2.0, 1024.0])
def test_alpha_is_dimensionless_exact_under_power_of_two_scaling(k):
    """Scaling depth by k leaves the normalized abundances bit-identical.

    counts_t and alpha * ambiguous_t both scale by k, so the M-step vector
    scales by k and normalization divides the factor straight back out. On a
    power of two the scaling is exact in binary floating point, so this is an
    equality, not an approximation.
    """
    for alpha in (0.0, 0.01, 1.0):
        base_counts = [100.0, 10.0, 50.0]
        theta_1, _ = _run(alpha, counts=base_counts)
        theta_k, _ = _run(alpha, counts=[c * k for c in base_counts])
        assert np.array_equal(theta_1, theta_k), f"alpha={alpha}, k={k}"


def test_alpha_is_dimensionless_to_tolerance_under_arbitrary_scaling():
    """Same claim at a non-power-of-two factor, where only rounding differs."""
    for alpha in (0.0, 0.003, 0.01, 0.3):
        base_counts = [137.0, 11.0, 53.0]
        theta_1, _ = _run(alpha, counts=base_counts)
        theta_k, _ = _run(alpha, counts=[c * 10.0 for c in base_counts])
        np.testing.assert_allclose(theta_1, theta_k, rtol=1e-9, atol=0.0)


def test_read_counts_scale_linearly_while_abundances_do_not_move():
    """The counts EM hands back are on the read scale; the abundances are not.

    Both must be true at once, and only one of them scales. Getting this
    backwards would make alpha look depth-dependent in exactly the way the sweep
    is checking for.
    """
    base_counts = [100.0, 10.0, 50.0]
    theta_1, sums_1 = _run(0.01, counts=base_counts)
    theta_4, sums_4 = _run(0.01, counts=[c * 4.0 for c in base_counts])

    assert np.array_equal(theta_1, theta_4)
    for t in range(3):
        assert sums_4[t] == pytest.approx(4.0 * sums_1[t], rel=1e-12)
