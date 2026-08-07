#!/usr/bin/env python3

import TranscriptFiltering
from Transcript import Transcript


class _DummySpliceGraph:
    def get_contig_acc(self):
        return "chr1"

    def get_contig_strand(self):
        return "-"


class _DummyMP:
    def __init__(self, read_count):
        self._read_count = read_count

    def get_read_count(self):
        return self._read_count


def _build_transcript(transcript_id, exons, simple_path, read_count):
    transcript = Transcript("chr1", exons, "-")
    transcript.set_gene_id("g1")
    transcript.set_transcript_id(transcript_id)
    transcript.set_simple_path(simple_path)
    transcript.set_read_counts_assigned(read_count)
    return transcript


def test_prune_degradation_products_handles_terminal_split_exon_paths():
    transcript_with_boundaries = _build_transcript(
        "t.keep",
        [[100, 150], [200, 250]],
        ["POLYA:1", "E:1", "E:2", "I:1", "E:3", "TSS:1"],
        10.0,
    )
    transcript_terminal_trunc = _build_transcript(
        "t.prune",
        [[101, 150], [200, 260]],
        ["E:2", "I:1", "E:3", "E:4"],
        3.0,
    )

    frac_read_assignments = {
        "t.keep": {_DummyMP(10): 1.0},
        "t.prune": {_DummyMP(3): 1.0},
    }

    retained = TranscriptFiltering.prune_likely_degradation_products(
        [transcript_with_boundaries, transcript_terminal_trunc],
        _DummySpliceGraph(),
        frac_read_assignments,
    )

    retained_ids = {transcript.get_transcript_id() for transcript in retained}
    assert retained_ids == {"t.keep"}


def test_prune_degradation_products_tolerates_zero_support_pair():
    # A gene with expression elsewhere can still contain two models that carry no
    # assigned reads; their relative support is undefined rather than zero.
    expressed = _build_transcript(
        "t.expressed",
        [[100, 150], [200, 250]],
        ["POLYA:1", "E:1", "E:2", "I:1", "E:3", "TSS:1"],
        12.0,
    )
    zero_long = _build_transcript(
        "t.zero_long",
        [[300, 350], [400, 450]],
        ["POLYA:2", "E:5", "I:2", "E:6", "TSS:2"],
        0.0,
    )
    zero_contained = _build_transcript(
        "t.zero_contained",
        [[300, 350], [400, 440]],
        ["POLYA:2", "E:5", "I:2", "E:6", "TSS:2"],
        0.0,
    )

    frac_read_assignments = {
        "t.expressed": {_DummyMP(12): 1.0},
        "t.zero_long": {},
        "t.zero_contained": {},
    }

    retained = TranscriptFiltering.prune_likely_degradation_products(
        [expressed, zero_long, zero_contained],
        _DummySpliceGraph(),
        frac_read_assignments,
    )

    retained_ids = {transcript.get_transcript_id() for transcript in retained}
    assert "t.expressed" in retained_ids
    # neither zero-support model may be pruned on relative expression grounds
    assert {"t.zero_long", "t.zero_contained"} <= retained_ids


def test_prune_degradation_products_prefers_better_supported_terminal_twin():
    # Terminal variants of one intron chain are the same length and each qualifies as
    # the other's degradation product, so the structural tests cannot separate them.
    # The better-supported model must win, and the winner must not depend on the order
    # the models happen to arrive in.
    def build_pair():
        more = _build_transcript(
            "t.more", [[102, 152], [202, 252]], ["E:2", "I:1", "E:3", "E:4"], 5.0
        )
        less = _build_transcript(
            "t.less", [[100, 150], [200, 250]], ["E:1", "E:2", "I:1", "E:3"], 2.0
        )
        assert more.get_cdna_len() == less.get_cdna_len()
        return more, less

    frac_read_assignments = {
        "t.more": {_DummyMP(5): 1.0},
        "t.less": {_DummyMP(2): 1.0},
    }

    for ordering in ((0, 1), (1, 0)):
        pair = build_pair()
        retained = TranscriptFiltering.prune_likely_degradation_products(
            [pair[i] for i in ordering],
            _DummySpliceGraph(),
            frac_read_assignments,
        )
        retained_ids = {transcript.get_transcript_id() for transcript in retained}
        assert retained_ids == {"t.more"}
