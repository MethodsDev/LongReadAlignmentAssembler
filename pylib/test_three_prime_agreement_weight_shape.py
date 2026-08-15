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


class _FakeSplicegraph:
    def __init__(self, strand="+"):
        self._strand = strand

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


def _weights_for(distances, P=None, monkeypatch=None):
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
        _FakeSplicegraph(), _FakeMultipath(), transcripts
    )
    return [weights[t.get_transcript_id()] for t in transcripts]


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
