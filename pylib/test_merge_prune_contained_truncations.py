#!/usr/bin/env python3

"""Structure-only containment prune applied after a GTF merge.

Each per-cluster run absorbs a contained model with no terminal feature of its own into
its container (prune_likely_degradation_products). A cluster where the container is weak or
absent keeps the truncation, so the merge sees it next to the container rebuilt in another
cluster. prune_contained_truncations_structural re-applies the structural cases on the
merged set. A merge carries no reads, so the function must not depend on read counts.
"""

import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if os.path.join(REPO_ROOT, "pylib") not in sys.path:
    sys.path.insert(0, os.path.join(REPO_ROOT, "pylib"))

from Transcript import Transcript
import TranscriptFiltering


def _model(tid, exons, path, strand="+", gene_id="g1"):
    t = Transcript("chr1", exons, strand)
    t.set_transcript_id(tid)
    t.set_gene_id(gene_id)
    t.set_simple_path(path)
    return t


# + strand container L: TSS, E1, I1, E2a, E2b, I2, E3, PolyA
# E2 is split at the short model's start (E2a | E2b), as the merge graph splits it.
L_EXONS = [[100, 200], [300, 450], [600, 700]]
L_PATH = ["TSS:1", "E:1", "I:1", "E:2a", "E:2b", "I:2", "E:3", "POLYA:1"]


def test_5prime_truncation_sharing_polyA_is_absorbed():
    L = _model("L", L_EXONS, L_PATH)
    S = _model("S", [[380, 450], [600, 700]], ["E:2b", "I:2", "E:3", "POLYA:1"])

    kept, pruned = TranscriptFiltering.prune_contained_truncations_structural(
        [S, L], "+"
    )

    assert [t.get_transcript_id() for t in kept] == ["L"]
    assert pruned == [("S", "L")]


def test_own_TSS_is_kept():
    # an internal-promoter model with its own TSS node is not a truncation
    L = _model("L", L_EXONS, L_PATH)
    S = _model(
        "S", [[380, 450], [600, 700]], ["TSS:2", "E:2b", "I:2", "E:3", "POLYA:1"]
    )

    kept, pruned = TranscriptFiltering.prune_contained_truncations_structural(
        [L, S], "+"
    )

    assert {t.get_transcript_id() for t in kept} == {"L", "S"}
    assert pruned == []


def test_distinct_polyA_is_kept():
    # no TSS, but its own PolyA: possible alternative polyadenylation, left alone
    L = _model("L", L_EXONS, L_PATH)
    S = _model("S", [[380, 450], [600, 700]], ["E:2b", "I:2", "E:3", "POLYA:2"])

    kept, _ = TranscriptFiltering.prune_contained_truncations_structural([L, S], "+")

    assert {t.get_transcript_id() for t in kept} == {"L", "S"}


def test_not_contained_is_kept():
    # S's chain is not inside L (different intron)
    L = _model("L", L_EXONS, L_PATH)
    S = _model("S", [[380, 450], [650, 700]], ["E:2b", "I:9", "E:3b", "POLYA:1"])

    kept, _ = TranscriptFiltering.prune_contained_truncations_structural([L, S], "+")

    assert {t.get_transcript_id() for t in kept} == {"L", "S"}


def test_other_gene_is_not_a_container():
    L = _model("L", L_EXONS, L_PATH, gene_id="g1")
    S = _model(
        "S", [[380, 450], [600, 700]], ["E:2b", "I:2", "E:3", "POLYA:1"], gene_id="g2"
    )

    kept, _ = TranscriptFiltering.prune_contained_truncations_structural([L, S], "+")

    assert {t.get_transcript_id() for t in kept} == {"L", "S"}


def test_minus_strand_and_spare_predicate():
    # - strand paths run low to high coordinate: PolyA first, TSS last
    L = _model(
        "L",
        L_EXONS,
        ["POLYA:1", "E:1", "I:1", "E:2a", "E:2b", "I:2", "E:3", "TSS:1"],
        strand="-",
    )
    S = _model(
        "S", [[100, 200], [300, 400]], ["POLYA:1", "E:1", "I:1", "E:2a"], strand="-"
    )

    kept, pruned = TranscriptFiltering.prune_contained_truncations_structural(
        [L, S], "-"
    )
    assert [t.get_transcript_id() for t in kept] == ["L"]
    assert pruned == [("S", "L")]

    L2 = _model(
        "L",
        L_EXONS,
        ["POLYA:1", "E:1", "I:1", "E:2a", "E:2b", "I:2", "E:3", "TSS:1"],
        strand="-",
    )
    S2 = _model(
        "S", [[100, 200], [300, 400]], ["POLYA:1", "E:1", "I:1", "E:2a"], strand="-"
    )
    kept, pruned = TranscriptFiltering.prune_contained_truncations_structural(
        [L2, S2], "-", spare=lambda t: t.get_transcript_id() == "S"
    )
    assert {t.get_transcript_id() for t in kept} == {"L", "S"}
    assert pruned == []


def test_tolerant_end_overhang_is_absorbed():
    # S overhangs L's 3' end by 6 bases (different exon segment), no TSS, no PolyA
    L = _model("L", L_EXONS, L_PATH)
    S = _model("S", [[380, 450], [600, 706]], ["E:2b", "I:2", "E:3", "E:3x"])

    kept, _ = TranscriptFiltering.prune_contained_truncations_structural([L, S], "+")
    assert {t.get_transcript_id() for t in kept} == {"L", "S"}

    kept, pruned = TranscriptFiltering.prune_contained_truncations_structural(
        [L, S], "+", end_tolerance=50
    )
    assert [t.get_transcript_id() for t in kept] == ["L"]
    assert pruned == [("S", "L")]


def test_tolerant_rule_keeps_far_end_and_own_TSS():
    L = _model("L", L_EXONS, L_PATH)
    far = _model("far", [[380, 450], [600, 900]], ["E:2b", "I:2", "E:3", "E:3y"])
    own_tss = _model(
        "own", [[380, 450], [600, 706]], ["TSS:2", "E:2b", "I:2", "E:3", "E:3x"]
    )
    kept, _ = TranscriptFiltering.prune_contained_truncations_structural(
        [L, far, own_tss], "+", end_tolerance=50
    )
    assert {t.get_transcript_id() for t in kept} == {"L", "far", "own"}


def test_tolerant_rule_minus_strand():
    # - strand: 3' end is the low coordinate; S lacks L's 5'-most (high) intron
    L = _model(
        "L",
        L_EXONS,
        ["POLYA:1", "E:1", "I:1", "E:2a", "E:2b", "I:2", "E:3", "TSS:1"],
        strand="-",
    )
    S = _model(
        "S", [[95, 200], [300, 400]], ["POLYA:2", "E:0", "E:1", "I:1", "E:2a"], strand="-"
    )
    kept, pruned = TranscriptFiltering.prune_contained_truncations_structural(
        [L, S], "-", end_tolerance=50
    )
    assert [t.get_transcript_id() for t in kept] == ["L"]

    # S's 5' exon extends past L's exon: not inside the container, kept
    S2 = _model(
        "S2", [[100, 200], [300, 500]], ["POLYA:1", "E:1", "I:1", "E:2a", "E:2b", "E:2c"], strand="-"
    )
    L2 = _model(
        "L",
        L_EXONS,
        ["POLYA:1", "E:1", "I:1", "E:2a", "E:2b", "I:2", "E:3", "TSS:1"],
        strand="-",
    )
    kept, _ = TranscriptFiltering.prune_contained_truncations_structural(
        [L2, S2], "-", end_tolerance=50
    )
    assert {t.get_transcript_id() for t in kept} == {"L", "S2"}


def test_tolerant_rule_absorbs_own_polyA_within_tolerance():
    # Documents a deliberate departure from the per-cluster rule, which needs the SAME
    # PolyA node: with end_tolerance set, a no-TSS 5' truncation ending within the
    # tolerance of its container's 3' end is absorbed even with a PolyA node of its own.
    L = _model("L", L_EXONS, L_PATH)
    S = _model("S", [[380, 450], [600, 700]], ["E:2b", "I:2", "E:3", "POLYA:2"])

    kept, _ = TranscriptFiltering.prune_contained_truncations_structural(
        [L, S], "+", end_tolerance=50
    )
    assert [t.get_transcript_id() for t in kept] == ["L"]


def test_prune_3prime_false_keeps_TSS_carrying_3prime_truncation():
    # S shares L's TSS and lacks a PolyA node (a 3' truncation): absorbed by default,
    # kept with prune_3prime=False; a TSS-less 5' truncation is absorbed either way.
    def models():
        L = _model("L", L_EXONS, L_PATH)
        S3 = _model("S3", [[100, 200], [300, 400]], ["TSS:1", "E:1", "I:1", "E:2a"])
        S5 = _model("S5", [[380, 450], [600, 700]], ["E:2b", "I:2", "E:3", "POLYA:1"])
        return [L, S3, S5]

    kept, _ = TranscriptFiltering.prune_contained_truncations_structural(models(), "+")
    assert [t.get_transcript_id() for t in kept] == ["L"]

    kept, pruned = TranscriptFiltering.prune_contained_truncations_structural(
        models(), "+", prune_3prime=False
    )
    assert {t.get_transcript_id() for t in kept} == {"L", "S3"}
    assert pruned == [("S5", "L")]
