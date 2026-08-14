#!/usr/bin/env python3

"""AssignmentTable.build must assign each read exactly once per gene.

The tracking file's contract is that a read's fractions over one gene's isoforms sum to
1.0: the read is fully accounted for, none invented and none lost. Two ways to break that
without producing anything that looks broken:

  over-assignment -- emit a row per candidate-list entry rather than per transcript, so a
  duplicated transcript contributes its fraction twice and the read totals 1.75

  skewed split    -- divide by the candidate-list length rather than the transcript count,
  so two isoforms split 2/3 and 1/3 instead of evenly, which still totals 1.0

Both are invisible in the output. These tests pin the sums and the row cardinality
directly, so neither can appear without a failure.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import LRAA_Globals
import StreamingQuant
from StreamingQuant import AssignmentTable


class _FakeSpliceGraph:
    def canonical_simple_path(self, simple_path):
        return "canon:" + ",".join(simple_path)


class _FakeMP:
    def __init__(self, path, mp_id, weights=None):
        self._path = path
        self._id = mp_id
        self._weights = weights or {}

    def get_simple_path(self):
        return self._path

    def get_id(self):
        return self._id


class _FakeTranscript:
    def __init__(self, tid, gene_id, num_exons=2, weight=1.0):
        self._tid = tid
        self._gene = gene_id
        self._num_exons = num_exons
        self._weight = weight

    def get_transcript_id(self):
        return self._tid

    def get_gene_id(self):
        return self._gene

    def get_output_gene_id(self):
        return self._gene

    def get_num_exon_segments(self):
        return self._num_exons

    def get_introns_string(self):
        return "i:" + self._tid

    def get_multipath_weight(self, mp):
        return self._weight


FRAC_IDX = 5
TX_IDX = 1
GENE_IDX = 0


def _sum_fracs_by_gene(rows):
    sums = {}
    for r in rows:
        sums[r[GENE_IDX]] = sums.get(r[GENE_IDX], 0.0) + r[FRAC_IDX]
    return sums


def test_theta_split_assigns_each_read_exactly_once():
    """Fractions over one gene's isoforms must sum to 1.0 and follow theta."""
    sg = _FakeSpliceGraph()
    mp = _FakeMP(["E:1", "E:2"], "MP1")
    t1 = _FakeTranscript("t1", "g1")
    t2 = _FakeTranscript("t2", "g1")
    theta = {"t1": 0.75, "t2": 0.25}

    old = LRAA_Globals.config["weight_reads_by_3prime_agreement"]
    LRAA_Globals.config["weight_reads_by_3prime_agreement"] = False
    try:
        table = AssignmentTable.build(sg, {mp: [t1, t2]}, theta, True, ())
    finally:
        LRAA_Globals.config["weight_reads_by_3prime_agreement"] = old

    rows = table.lookup("canon:E:1,E:2")
    assert rows is not None and len(rows) == 2
    got = {r[TX_IDX]: r[FRAC_IDX] for r in rows}
    assert got["t1"] == pytest.approx(0.75)
    assert got["t2"] == pytest.approx(0.25)
    assert _sum_fracs_by_gene(rows)["g1"] == pytest.approx(1.0)


def test_equal_split_without_theta_divides_by_transcript_count():
    """No EM: an even split over the gene's isoforms, not over list entries."""
    sg = _FakeSpliceGraph()
    mp = _FakeMP(["E:1", "E:2"], "MP1")
    t1 = _FakeTranscript("t1", "g1")
    t2 = _FakeTranscript("t2", "g1")

    old = LRAA_Globals.config["weight_reads_by_3prime_agreement"]
    LRAA_Globals.config["weight_reads_by_3prime_agreement"] = False
    try:
        table = AssignmentTable.build(sg, {mp: [t1, t2]}, {}, False, ())
    finally:
        LRAA_Globals.config["weight_reads_by_3prime_agreement"] = old

    rows = table.lookup("canon:E:1,E:2")
    got = {r[TX_IDX]: r[FRAC_IDX] for r in rows}
    assert got["t1"] == pytest.approx(0.5), "even split expected over 2 isoforms"
    assert got["t2"] == pytest.approx(0.5)
    assert _sum_fracs_by_gene(rows)["g1"] == pytest.approx(1.0)


def test_two_overlapping_genes_each_sum_to_one():
    """A path compatible with two genes gets one full split per gene, not one pooled.

    Theta is normalized within a gene, so pooling across genes would divide numbers that
    were never normalized together.
    """
    sg = _FakeSpliceGraph()
    mp = _FakeMP(["E:1", "E:2"], "MP1")
    a1 = _FakeTranscript("a1", "gA")
    a2 = _FakeTranscript("a2", "gA")
    b1 = _FakeTranscript("b1", "gB")
    theta = {"a1": 0.6, "a2": 0.4, "b1": 1.0}

    old = LRAA_Globals.config["weight_reads_by_3prime_agreement"]
    LRAA_Globals.config["weight_reads_by_3prime_agreement"] = False
    try:
        table = AssignmentTable.build(sg, {mp: [a1, a2, b1]}, theta, True, ())
    finally:
        LRAA_Globals.config["weight_reads_by_3prime_agreement"] = old

    rows = table.lookup("canon:E:1,E:2")
    sums = _sum_fracs_by_gene(rows)
    assert sums["gA"] == pytest.approx(1.0)
    assert sums["gB"] == pytest.approx(1.0)


def test_duplicate_transcript_in_group_is_refused():
    """The tripwire: duplicates cannot arise today, and must not pass quietly if they do."""
    sg = _FakeSpliceGraph()
    mp = _FakeMP(["E:1", "E:2"], "MP1")
    t1 = _FakeTranscript("t1", "g1")
    t1_again = _FakeTranscript("t1", "g1")
    t2 = _FakeTranscript("t2", "g1")

    with pytest.raises(RuntimeError, match="duplicate transcript ids"):
        AssignmentTable.build(sg, {mp: [t1, t1_again, t2]}, {"t1": 0.75, "t2": 0.25}, True, ())


def test_same_transcript_under_two_genes_is_refused():
    """The violation the per-gene check cannot see.

    Grouping is by gene id, so one transcript id under two gene ids puts one copy in each
    group and every group still looks internally unique. The intra-group check is blind to
    it. It matters because the candidate list is merged across genes with dict update(),
    which keeps the last object -- so the real consequence is an isoform silently dropped
    from whichever gene lost, not a duplicated one.
    """
    sg = _FakeSpliceGraph()
    mp = _FakeMP(["E:1", "E:2"], "MP1")
    t_a = _FakeTranscript("shared", "gA")
    t_b = _FakeTranscript("shared", "gB")

    with pytest.raises(RuntimeError, match="two gene ids"):
        AssignmentTable.build(sg, {mp: [t_a, t_b]}, {"shared": 1.0}, True, ())


def test_zero_denominator_matches_default_path_rather_than_summing_to_one():
    """Parity, not arithmetic elegance.

    All-zero theta gives a zero denominator, and the fractions come out 0.0, so this read's
    per-gene fractions sum to 0 instead of 1. That is what the default path's E-step does
    (EM.py:309-313), and reproducing the default path is the whole point of this mode, so
    an equal-split "fix" here would be a divergence bug rather than an improvement.

    Pinned so that nobody later reads the zero as an oversight and helpfully repairs it.
    """
    sg = _FakeSpliceGraph()
    mp = _FakeMP(["E:1", "E:2"], "MP1")
    t1 = _FakeTranscript("t1", "g1")
    t2 = _FakeTranscript("t2", "g1")

    old = LRAA_Globals.config["weight_reads_by_3prime_agreement"]
    LRAA_Globals.config["weight_reads_by_3prime_agreement"] = False
    try:
        table = AssignmentTable.build(sg, {mp: [t1, t2]}, {"t1": 0.0, "t2": 0.0}, True, ())
    finally:
        LRAA_Globals.config["weight_reads_by_3prime_agreement"] = old

    rows = table.lookup("canon:E:1,E:2")
    assert len(rows) == 2, "rows are still emitted; only the fractions are zero"
    assert all(r[FRAC_IDX] == 0.0 for r in rows)
    assert _sum_fracs_by_gene(rows)["g1"] == pytest.approx(0.0)


def test_missing_theta_is_refused_rather_than_assigned_zero():
    """A transcript with no theta means the state and theta are from different runs."""
    sg = _FakeSpliceGraph()
    mp = _FakeMP(["E:1", "E:2"], "MP1")
    t1 = _FakeTranscript("t1", "g1")
    t2 = _FakeTranscript("t2", "g1")

    with pytest.raises(RuntimeError, match="no EM theta"):
        AssignmentTable.build(sg, {mp: [t1, t2]}, {"t1": 0.75}, True, ())


def test_empty_theta_with_em_run_is_refused_not_equal_split():
    """The ambiguity of an empty theta must be resolved by the caller, not guessed.

    An empty theta means either run_EM=False or EM ran and the state was reset in between.
    Inferring from emptiness silently takes the equal-split branch in the second case, which
    produces a complete, internally consistent table built from no abundance information --
    and nothing in the output distinguishes it from a legitimate no-EM run.
    """
    sg = _FakeSpliceGraph()
    mp = _FakeMP(["E:1", "E:2"], "MP1")
    t1 = _FakeTranscript("t1", "g1")
    t2 = _FakeTranscript("t2", "g1")

    with pytest.raises(RuntimeError, match="EM ran but no theta"):
        AssignmentTable.build(sg, {mp: [t1, t2]}, {}, True, ())


def test_two_multipaths_on_one_canonical_path_are_refused():
    """Both row sets would land under one key and the read would be written twice."""
    sg = _FakeSpliceGraph()
    mp_a = _FakeMP(["E:1", "E:2"], "MP1")
    mp_b = _FakeMP(["E:1", "E:2"], "MP2")
    t1 = _FakeTranscript("t1", "g1")

    with pytest.raises(RuntimeError, match="two multipaths map to canonical path"):
        AssignmentTable.build(sg, {mp_a: [t1], mp_b: [t1]}, {"t1": 1.0}, True, ())
