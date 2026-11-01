#!/usr/bin/env python3
"""Weak-junction guards of the SE pass.

A single-exon gene carrying a rare splice (TUSCO RNF26, SEPHS2, C6orf47, NCBP2AS2 in
BT474: ~1,000 unspliced full-length reads, a 2-read GT-AG junction) was reported only
as the spliced model, because the spliced model's exons masked the unspliced reads from
the SE graph and its intron's "background" floor -- the read-through pile itself -- was
subtracted from the SE coverage. Both steps now weigh the junction's spliced-read
support J against the unspliced reads through the intron.
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA_Globals

# Splice_graph must be imported first: it and Pretty_alignment_manager import each other.
from Splice_graph import Splice_graph
from Pretty_alignment_manager import (
    Pretty_alignment_manager,
    count_readthrough_reads,
    weak_junction_exons_for_SE_mask,
)


class _SE_read:
    def __init__(self, lend, rend, weight=1.0):
        self._span = (lend, rend)
        self._weight = weight

    def get_alignment_span(self):
        return self._span

    def get_normalization_weight(self):
        return self._weight


class _Model:
    def __init__(self, exons, support=None):
        self._exons = exons
        self._support = support

    def get_exon_segments(self):
        return self._exons

    def get_intron_read_support(self):
        return self._support


# a 2-exon model, intron 201-300
EXONS = [(101, 200), (301, 400)]
INTRON = (201, 300)


def _readthrough(n, weight=1.0):
    # identical spans on purpose: an IntervalTree would collapse them
    return [_SE_read(110, 390, weight) for _ in range(n)]


def test_readthrough_counts_identical_spans_and_weights():
    reads = _readthrough(5, weight=2.0) + [_SE_read(150, 250), _SE_read(250, 350)]
    assert count_readthrough_reads(reads, {INTRON}) == {INTRON: 10.0}


def test_readthrough_requires_both_flanking_bases():
    reads = [_SE_read(200, 301), _SE_read(201, 301), _SE_read(200, 300)]
    assert count_readthrough_reads(reads, {INTRON}) == {INTRON: 1}


@pytest.fixture
def config(monkeypatch):
    def set_(**kv):
        for k, v in kv.items():
            monkeypatch.setitem(LRAA_Globals.config, k, v)

    set_(SE_weak_junction_stats_path=None)
    return set_


def test_weak_junction_exempts_flanking_exons(config):
    config(SE_mask_max_readthrough_to_junction_ratio=10)
    model = _Model(EXONS, {INTRON: 2})
    weak = weak_junction_exons_for_SE_mask(_readthrough(100), [model])
    assert weak == {(id(model), 101, 200), (id(model), 301, 400)}


def test_well_supported_junction_still_masks(config):
    config(SE_mask_max_readthrough_to_junction_ratio=10)
    model = _Model(EXONS, {INTRON: 50})
    assert weak_junction_exons_for_SE_mask(_readthrough(100), [model]) == set()


def test_unknown_support_and_disabled_guard_never_exempt(config):
    config(SE_mask_max_readthrough_to_junction_ratio=10)
    assert weak_junction_exons_for_SE_mask(_readthrough(100), [_Model(EXONS, None)]) == set()
    assert weak_junction_exons_for_SE_mask(_readthrough(100), [_Model(EXONS, {INTRON: None})]) == set()
    config(SE_mask_max_readthrough_to_junction_ratio=None)
    assert weak_junction_exons_for_SE_mask(_readthrough(100), [_Model(EXONS, {INTRON: 2})]) == set()


def test_only_exons_flanking_the_weak_intron_are_exempt(config):
    config(SE_mask_max_readthrough_to_junction_ratio=10)
    exons = [(101, 200), (301, 400), (501, 600)]
    model = _Model(exons, {(201, 300): 2, (401, 500): 500})
    weak = weak_junction_exons_for_SE_mask(_readthrough(100), [model])
    assert weak == {(id(model), 101, 200), (id(model), 301, 400)}


def test_mask_keeps_readthrough_reads_of_a_weak_model_only(config):
    # SEPHS2 shape: a long exon holds most of each unspliced read, a short intron
    exons, intron = [(101, 400), (421, 500)], (401, 420)
    mgr = Pretty_alignment_manager.__new__(Pretty_alignment_manager)
    reads = [_SE_read(110, 450) for _ in range(100)]  # 85% inside the first exon
    config(SE_mask_max_readthrough_to_junction_ratio=10, min_SE_read_ME_exon_overlap_pct=50)
    kept = mgr.apply_SE_read_encapsulation_mask(reads, [_Model(exons, {intron: 2})])
    assert len(kept) == 100
    kept = mgr.apply_SE_read_encapsulation_mask(reads, [_Model(exons, {intron: 50})])
    assert len(kept) == 0
    config(SE_mask_max_readthrough_to_junction_ratio=None)
    kept = mgr.apply_SE_read_encapsulation_mask(reads, [_Model(exons, {intron: 2})])
    assert len(kept) == 0


def _coverage_graph(depth_exon, depth_intron):
    sg = Splice_graph.__new__(Splice_graph)
    cov = [0] * 501
    for i in range(101, 401):
        cov[i] = depth_exon
    for i in range(INTRON[0], INTRON[1] + 1):
        cov[i] = depth_intron
    sg._contig_base_cov = cov
    sg._contig_seq_len = 500
    return sg


def test_floor_far_above_junction_support_is_not_subtracted(config):
    config(SE_intronic_background_max_floor_to_junction_ratio=10,
           SE_intronic_background_intron_length_pctile=80)
    sg = _coverage_graph(0, 900)
    sg._subtract_intronic_background([_Model(EXONS, {INTRON: 2})], "chr1", "+")
    assert sg._contig_base_cov[250] == 900


def test_floor_below_junction_support_is_subtracted(config):
    config(SE_intronic_background_max_floor_to_junction_ratio=10,
           SE_intronic_background_intron_length_pctile=80)
    sg = _coverage_graph(0, 3)
    sg._subtract_intronic_background([_Model(EXONS, {INTRON: 200})], "chr1", "+")
    assert sg._contig_base_cov[250] == 0


def test_subtraction_guard_disabled_or_support_unknown_subtracts(config):
    config(SE_intronic_background_intron_length_pctile=80)
    for ratio, support in ((None, {INTRON: 2}), (10, None), (10, {INTRON: None})):
        config(SE_intronic_background_max_floor_to_junction_ratio=ratio)
        sg = _coverage_graph(0, 900)
        sg._subtract_intronic_background([_Model(EXONS, support)], "chr1", "+")
        assert sg._contig_base_cov[250] == 0, (ratio, support)
