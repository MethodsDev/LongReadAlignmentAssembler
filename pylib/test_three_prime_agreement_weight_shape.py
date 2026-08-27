#!/usr/bin/env python3

"""Pin the SHAPE of the 3'-end-agreement read weighting.

config['weight_reads_by_3prime_agreement'] apportions an ambiguous read across the N
transcripts it is compatible with, using each transcript's 3'-end distance d_i and the
pseudocount P = config['max_dist_between_alt_polyA_sites']:

    w_i = 1 - (d_i + P) / (sum_j d_j + N*P)     then renormalized to sum to 1

Everything downstream -- EM's E-step denominator, the read_weight column of
quant.tracking -- consumes the renormalized w_i, so these are the numbers that decide
where an ambiguous read's mass lands.

The formula has three properties that are easy to assume and worth pinning, because a
reshaping would break them silently and the effect on any single MARD number is far
below what a benchmark can resolve:

  1. N=1 falls through the sum_weights == 0 branch to weight 1.0. It never reaches the
     ratio at all, so no unambiguous read is ever rescaled.
  2. Equal distances give exactly equal weights at every N. This is not a corner case:
     on clean data most multi-compatible reads are equidistant from all their candidate
     3' ends, so the weighting is identically uniform for them.
  3. Relative spread is bounded by 1 + (d_max - d_min)/(N*(dbar + P)) and therefore
     COLLAPSES toward 1 as N grows with the distance spread held fixed. The weighting
     is structurally least able to discriminate exactly where ambiguity is worst.

Property 3 is the one that decides whether the weighting is correctly shaped, so it is
pinned as a monotone claim over N rather than as a single magic number.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import LRAA_Globals
from Quantify import Quantify
from GenomeFeature import Exon, Intron


class _FakeSplicegraph:
    def __init__(self, strand="+"):
        self._strand = strand

    def get_contig_acc(self):
        return "chrFake"

    def get_contig_strand(self):
        return self._strand


class _FakeMultipath:
    def get_simple_path(self):
        return ["fake_mp_simple_path"]


class _FakeTranscript:
    def __init__(self, transcript_id):
        self._id = transcript_id

    def get_transcript_id(self):
        return self._id

    def get_simple_path(self):
        return [f"fake_sp_{self._id}"]


def _weights_for(distances, P=None, monkeypatch=None, strand="+"):
    """Renormalized 3'-agreement weights for a read at these 3'-end distances."""
    q = Quantify(run_EM=False, max_EM_iterations=1)
    transcripts = [_FakeTranscript(f"t{i}") for i in range(len(distances))]
    dist_by_id = {t.get_transcript_id(): d for t, d in zip(transcripts, distances)}

    def fake_dist(self, splice_graph, source_sp, target_sp, from_which_end):
        # target_sp is the transcript's simple path, built as fake_sp_<id> above.
        return dist_by_id[target_sp[0].removeprefix("fake_sp_")]

    monkeypatch.setattr(Quantify, "_get_simple_path_dist_to_termini", fake_dist)
    if P is not None:
        monkeypatch.setitem(LRAA_Globals.config, "max_dist_between_alt_polyA_sites", P)

    weights = q._assign_read_weights_based_on_read_end_agreement(
        _FakeSplicegraph(strand), _FakeMultipath(), transcripts
    )
    return [weights[t.get_transcript_id()] for t in transcripts]


def _three_prime_keys_used(strand, monkeypatch):
    """Which simple-path coordinate the weighting measured, for one contig strand."""
    seen = []

    def record(self, splice_graph, source_sp, target_sp, from_which_end):
        seen.append(from_which_end)
        return 0

    monkeypatch.setattr(Quantify, "_get_simple_path_dist_to_termini", record)
    Quantify(run_EM=False, max_EM_iterations=1)._assign_read_weights_based_on_read_end_agreement(
        _FakeSplicegraph(strand), _FakeMultipath(), [_FakeTranscript("t0")]
    )
    return seen


def test_single_compatible_transcript_gets_weight_one(monkeypatch):
    """N=1 must never be rescaled, whatever its 3'-end distance is.

    adj_sum_dists == d + P exactly, so the raw weight is 0 and the sum_weights == 0
    branch supplies 1.0. A reshaping that removed that branch would send every
    unambiguous read to weight 0 and delete the unambiguous evidence EM depends on.
    """
    for d in (0, 1, 50, 5000):
        assert _weights_for([d], monkeypatch=monkeypatch) == [1.0]


@pytest.mark.parametrize("n", [2, 3, 5, 10])
def test_equal_distances_give_uniform_weights(n, monkeypatch):
    """Equidistant candidates are indistinguishable, and the weighting says so."""
    w = _weights_for([37] * n, monkeypatch=monkeypatch)
    assert w == pytest.approx([1.0 / n] * n)


@pytest.mark.parametrize("n", [1, 2, 3, 5, 10])
def test_weights_always_sum_to_one(n, monkeypatch):
    """EM's E-step divides by a weighted sum; unnormalized weights would rescale reads."""
    w = _weights_for([10 * i for i in range(n)], monkeypatch=monkeypatch)
    assert sum(w) == pytest.approx(1.0)


def test_closer_three_prime_end_gets_more_weight(monkeypatch):
    """Direction check: the point of the weighting is to favour 3'-end agreement."""
    w = _weights_for([0, 200], monkeypatch=monkeypatch)
    assert w[0] > w[1]


def test_two_transcript_closed_form(monkeypatch):
    """At N=2 the renormalized weight is exactly (d_other + P) / (d_1 + d_2 + 2P).

    Pinning the closed form makes the role of P explicit: it is a pseudocount added to
    BOTH distances, so it damps the ratio toward 1 and cannot be tuned away.
    """
    d1, d2, P = 10, 400, 50
    w = _weights_for([d1, d2], P=P, monkeypatch=monkeypatch)
    denom = d1 + d2 + 2 * P
    assert w == pytest.approx([(d2 + P) / denom, (d1 + P) / denom])


def test_pseudocount_damps_the_weight_ratio(monkeypatch):
    """Larger P flattens the weights; smaller P sharpens them. Monotone, at fixed N.

    At N=2 the ratio is (d_2 + P)/(d_1 + P), so P interpolates the whole scale between
    the raw distance ratio and no preference at all. With d = (50, 300) the undamped
    ratio is 6.0; the default P=50 already cuts it to 3.5, and P=500 to 1.45.
    """
    ratios = []
    for P in (5, 15, 50, 150, 500):
        w = _weights_for([50, 300], P=P, monkeypatch=monkeypatch)
        ratios.append(max(w) / min(w))
    assert ratios == sorted(ratios, reverse=True), ratios
    assert ratios == pytest.approx([305 / 55, 315 / 65, 350 / 100, 450 / 200, 800 / 550])


def test_relative_spread_collapses_as_N_grows(monkeypatch):
    """The structural defect: adding equidistant competitors mutes a real 3' difference.

    Two transcripts 0 and 300 nt from the read's 3' end are held fixed while ever more
    candidates at the mean distance are added. The evidence distinguishing the original
    pair has not changed, but the weight ratio between them decays toward 1, so EM sees
    progressively less of it.
    """
    ratios = []
    for n_extra in (0, 1, 3, 8, 18):
        w = _weights_for([0, 300] + [150] * n_extra, monkeypatch=monkeypatch)
        ratios.append(w[0] / w[1])
    assert ratios == sorted(ratios, reverse=True), ratios
    assert ratios[0] > 1.9   # N=2:  a real 2x preference
    assert ratios[-1] < 1.2  # N=20: the same evidence, now nearly uniform


def test_spread_collapse_is_not_rescued_by_shrinking_P(monkeypatch):
    """P cannot undo the N-flattening: the N*P term grows with N too.

    If sharpening P were a fix, the recommendation would be "retune P". It is not --
    even at P=5 the same fixed distance gap is nearly erased by N=20, because sum_j d_j
    grows with N whether or not P does.
    """
    r2 = _weights_for([0, 300], P=5, monkeypatch=monkeypatch)
    r20 = _weights_for([0, 300] + [150] * 18, P=5, monkeypatch=monkeypatch)
    assert max(r2) / min(r2) > 15.0
    assert max(r20) / min(r20) < 1.25


def test_large_P_drives_weights_to_uniform(monkeypatch):
    """P -> infinity sends every weight to 1/N regardless of the distances.

    (d_i + P)/(sum_j d_j + N*P) -> 1/N as P dominates, so the raw weight tends to
    1 - 1/N for every i and renormalizes to 1/N. Distance information is erased.

    Convergence is O(d/(N*P)), NOT exact at any finite P: at P=1e6 with a 1200nt spread
    the residual is still ~3e-5. That matters for interpreting a large-P run against the
    off arm -- they converge, they do not coincide.
    """
    d = [0, 300, 1200, 75]
    residuals = []
    for P in (10**4, 10**5, 10**6):
        w = _weights_for(d, P=P, monkeypatch=monkeypatch)
        residuals.append(max(abs(x - 0.25) for x in w))
    assert residuals == sorted(residuals, reverse=True), residuals
    assert residuals[-1] < 1e-4
    assert residuals[-1] > 0.0


def test_uniform_weights_match_no_weighting_in_EM_to_rounding(monkeypatch):
    """Uniform weights and all-ones weights give EM the same answer, to float rounding.

    This is why the large-P limit is not a new regime but the OFF arm: the E-step forms
    frac = w_i*theta_i / sum_j(w_j*theta_j), which is invariant under scaling every
    weight of a multipath by the same constant, and disabling the feature simply sets
    every weight to 1.0. So sweeping P upward interpolates between the shipped
    behaviour and switching the feature off -- it cannot reach anything beyond that,
    which puts a computable ceiling on how much large P can cost.

    The equivalence is algebraic but NOT bit-exact in floating point: dividing by N and
    renormalizing does not reassociate exactly, and this asserts 1-ULP agreement rather
    than identity. Measured at 1 ULP on theta here. Pinned on EM directly rather than
    argued, because the whole P analysis rests on it.
    """
    from EM import em_algorithm_with_weights

    mp_assignments = [[0, 1], [1, 2], [0], [0, 1, 2], [2]]
    mp_read_counts = [120.0, 45.0, 300.0, 17.0, 88.0]
    ones = [[1.0] * len(mp) for mp in mp_assignments]
    uniform = [[1.0 / len(mp)] * len(mp) for mp in mp_assignments]

    theta_a, sums_a, fracs_a = em_algorithm_with_weights(
        mp_assignments, ones, mp_read_counts, 3, max_iter=200, tol=1e-10
    )
    theta_b, sums_b, fracs_b = em_algorithm_with_weights(
        mp_assignments, uniform, mp_read_counts, 3, max_iter=200, tol=1e-10
    )

    assert list(theta_b) == pytest.approx(list(theta_a), rel=1e-12)
    # fracs is a ragged list-of-lists, one row per multipath; flatten to compare.
    flat_a = [x for row in fracs_a for x in row]
    flat_b = [x for row in fracs_b for x in row]
    assert flat_b == pytest.approx(flat_a, rel=1e-12)
    for t in sums_a:
        assert sums_b[t] == pytest.approx(sums_a[t], rel=1e-12)


def test_nonuniform_weights_do_change_EM(monkeypatch):
    """Guard for the test above: it must not pass because EM ignores weights entirely."""
    from EM import em_algorithm_with_weights

    mp_assignments = [[0, 1], [1, 2], [0], [0, 1, 2], [2]]
    mp_read_counts = [120.0, 45.0, 300.0, 17.0, 88.0]
    ones = [[1.0] * len(mp) for mp in mp_assignments]
    skewed = [[0.9, 0.1], [0.5, 0.5], [1.0], [0.7, 0.2, 0.1], [1.0]]

    theta_a, _, _ = em_algorithm_with_weights(
        mp_assignments, ones, mp_read_counts, 3, max_iter=200, tol=1e-10
    )
    theta_b, _, _ = em_algorithm_with_weights(
        mp_assignments, skewed, mp_read_counts, 3, max_iter=200, tol=1e-10
    )
    assert list(theta_a) != list(theta_b)


# ------------------------------------------------ which end the weighting measures
#
# The weighting compares the read's 3' end to each candidate's 3' end. Simple paths
# are ordered by genomic coordinate, so which COORDINATE that is depends on the
# strand, and getting it wrong measures the 5' end instead -- producing weights that
# still sum to 1, still look uniform on equidistant candidates, and are reversed.
#
# The hazard these pin is not a wrong branch but a missing strand. Splice_graph's
# _contig_strand defaults to None (Splice_graph.py:62) and is assigned only by
# build_splice_graph_for_contig (Splice_graph.py:289), and a bare
# `"rend" if strand == "+" else "lend"` sends None down the MINUS branch. Every
# production caller is downstream of that builder and it only ever receives "+" or
# "-" (LRAA:2572 iterates exactly those), so this is a refusal for a state that
# should be unreachable -- which is the point: unreachable is not the same as safe,
# and the failure mode if it is ever reached is silent.


def test_plus_strand_measures_the_rend(monkeypatch):
    assert _three_prime_keys_used("+", monkeypatch) == ["rend"]


def test_minus_strand_measures_the_lend(monkeypatch):
    """Not a symmetry to assume: on '-' the 3' end IS the lower coordinate."""
    assert _three_prime_keys_used("-", monkeypatch) == ["lend"]


def test_an_unset_contig_strand_is_refused(monkeypatch):
    """A graph constructed but never built reports strand None.

    The refusal must NAME the contig and the value: on a whole-genome run "some
    splice graph had no strand" is not actionable, and the value is what says whether
    the graph was unbuilt (None) or carries something unexpected.
    """

    with pytest.raises(RuntimeError) as err:
        _weights_for([10, 200], monkeypatch=monkeypatch, strand=None)

    message = str(err.value)
    assert "chrFake" in message
    assert "None" in message
    assert "3'" in message


def test_a_stringified_none_strand_is_refused(monkeypatch):
    """"None" the STRING, not the object, because CpuBudget.WorkUnit stringifies.

    WorkUnit.__new__ stores str(contig_strand) (pylib/CpuBudget.py:109), so a strand
    that went missing upstream arrives here as the four-character string "None" and
    an `is None` check would wave it through into the minus-strand branch.
    """

    with pytest.raises(RuntimeError) as err:
        _weights_for([10, 200], monkeypatch=monkeypatch, strand="None")
    assert "'None'" in str(err.value)


def test_a_present_strand_still_weights(monkeypatch):
    """Guard for the three above: the refusal must not have swallowed the happy path."""
    assert sum(_weights_for([10, 200], monkeypatch=monkeypatch)) == pytest.approx(1.0)


# --------------------------------------- what counts as exonic in the 3' distance


class _AnnotatedExon(Exon):
    """A subclass of Exon.

    None exists in the tree today (GenomeFeature.py:147 has no subclasses), so the
    exact-type comparison this replaced was latent rather than live. It is declared
    here because that is the only way to test the distinction, and because the day
    someone adds a real one is the day the distance silently shortens.
    """


class _NodeLookupGraph:
    def __init__(self, nodes):
        self._nodes = {node.get_id(): node for node in nodes}

    def get_contig_acc(self):
        return "chrFake"

    def get_contig_strand(self):
        return "+"

    def get_node_obj_via_id(self, node_id):
        return self._nodes[node_id]


def _dist_over(nodes, shared):
    graph = _NodeLookupGraph(list(nodes) + [shared])
    return Quantify(
        run_EM=False, max_EM_iterations=1
    )._get_simple_path_dist_to_termini(
        graph,
        [shared.get_id()],
        [node.get_id() for node in nodes] + [shared.get_id()],
        "lend",
    )


def test_an_exon_subclass_counts_toward_the_three_prime_distance():
    """The distance is exonic bases, and a subclass of Exon is still exonic.

    An exact type comparison skips it, shortening the distance and so raising that
    candidate's weight -- with no error, no log line, and a weight vector that still
    sums to 1. Here the subclassed exon carries 50 of the 150 bases, so skipping it
    understates the distance by a third.
    """

    plain = Exon("chrFake", 101, 200, "+", 10.0)          # 100 bases
    subclassed = _AnnotatedExon("chrFake", 301, 350, "+", 10.0)  # 50 bases
    shared = Exon("chrFake", 401, 500, "+", 10.0)

    assert _dist_over([plain, subclassed], shared) == 150


def test_introns_still_do_not_count():
    """Guard on the widening: Intron is a SIBLING of Exon, not a subclass.

    The distance is exonic bases only, so relaxing the test all the way to
    GenomeFeature would add intron length to every spliced candidate. isinstance(Exon)
    admits Exon and its subclasses and nothing else, which is what this pins.
    """

    plain = Exon("chrFake", 101, 200, "+", 10.0)
    intron = Intron("chrFake", 201, 300, "+", 5)
    shared = Exon("chrFake", 401, 500, "+", 10.0)

    assert _dist_over([plain, intron], shared) == 100


# ------------------------------------- the E-step's per-read renormalization


def test_flat_weights_conserve_the_read_total():
    """What disabling 3' weighting relies on, asserted rather than assumed.

    EM replaces every candidate weight with a flat 1.0 when
    weight_reads_by_3prime_agreement is off (EM.py:117, :250). That is equivalent to an
    even split ONLY because the E-step renormalizes each multipath across its own
    candidates: total_weight sums w*theta over exactly that mp's candidates
    (EM.py:333-336) and frac_assignment divides by it (EM.py:342-344), so each mp's
    fractions sum to 1 and the per-transcript counts sum back to the input support.
    Were it not renormalized, turning the feature off would rescale totals rather than
    only flattening ratios. Asserted for both weight vectors so the invariant is shown
    to be the E-step's and not a property of uniform input.
    """

    from EM import em_algorithm_with_weights

    mp_assignments = [[0, 1], [1, 2], [0], [0, 1, 2], [2]]
    support = [120.0, 45.0, 300.0, 17.0, 88.0]
    flat = [[1.0] * len(mp) for mp in mp_assignments]
    skewed = [[0.9, 0.1], [0.5, 0.5], [1.0], [0.7, 0.2, 0.1], [1.0]]

    for weights in (flat, skewed):
        _theta, sums, fracs = em_algorithm_with_weights(
            mp_assignments, weights, support, 3, max_iter=200, tol=1e-12
        )
        for row in fracs:
            assert sum(row) == pytest.approx(1.0)
        assert sum(sums.values()) == pytest.approx(sum(support))
