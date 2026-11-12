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

    def get_read_weight(self):
        # no read here carried a normalization weight, so weight is the count
        return float(self._read_count)


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


class _PlusSpliceGraph(_DummySpliceGraph):
    def get_contig_strand(self):
        return "+"


def _plus_transcript(transcript_id, exons, simple_path, read_count):
    transcript = Transcript("chr1", exons, "+")
    transcript.set_gene_id("g1")
    transcript.set_transcript_id(transcript_id)
    transcript.set_simple_path(simple_path)
    transcript.set_read_counts_assigned(read_count)
    return transcript


def _prune_pair(short_end, longer_path, enabled=True):
    """A longer model ending at 400 (path `longer_path`) and an abundant 5'-truncated
    model of the same chain, no TSS call, PolyA call, ending at `short_end`."""
    import LRAA_Globals

    saved = LRAA_Globals.config.get("contained_PolyA_match_by_end_coords")
    LRAA_Globals.config["contained_PolyA_match_by_end_coords"] = enabled
    try:
        longer = _plus_transcript(
            "t.long", [[100, 150], [200, 250], [300, 400]], longer_path, 2.0
        )
        shorter = _plus_transcript(
            "t.short",
            [[210, 250], [300, short_end]],
            ["E:2", "I:2", "E:3", "POLYA:1"],
            8.0,
        )
        frac = {"t.long": {_DummyMP(2): 1.0}, "t.short": {_DummyMP(8): 1.0}}
        retained = TranscriptFiltering.prune_likely_degradation_products(
            [longer, shorter], _PlusSpliceGraph(), frac
        )
        return {t.get_transcript_id(): t for t in retained}
    finally:
        LRAA_Globals.config["contained_PolyA_match_by_end_coords"] = saved


NO_POLYA = ["E:1", "I:1", "E:2", "I:2", "E:3"]


def test_polyA_less_container_absorbs_and_adopts_coincident_3prime_end():
    # The container ends 1 bp from the contained model's PolyA site but never reached
    # the site vertex: identity sees a distinguishing PolyA, end coordinates do not.
    retained = _prune_pair(401, NO_POLYA)
    assert set(retained) == {"t.long"}
    survivor = retained["t.long"]
    assert survivor.get_simple_path()[-1] == "POLYA:1"
    assert survivor.has_PolyA()
    assert survivor.get_coords() == (100, 401)

    assert set(_prune_pair(401, NO_POLYA, enabled=False)) == {"t.long", "t.short"}


def test_container_with_its_own_polyA_is_judged_by_identity_only():
    # The container carries a DIFFERENT PolyA call within 50 bp: the coordinate rule
    # does not apply and the abundant contained model survives, as before.
    retained = _prune_pair(401, NO_POLYA + ["POLYA:2"])
    assert set(retained) == {"t.long", "t.short"}
    assert retained["t.long"].get_simple_path()[-1] == "POLYA:2"


def test_nested_alternative_polyA_is_not_absorbed_by_end_coords():
    # Ends farther apart than max_dist_between_alt_polyA_sites: a nested APA isoform,
    # left to the expression test, which it passes.
    import LRAA_Globals

    assert 400 - 330 > LRAA_Globals.config["max_dist_between_alt_polyA_sites"]
    retained = _prune_pair(330, NO_POLYA)
    assert set(retained) == {"t.long", "t.short"}
    assert not retained["t.long"].has_PolyA()


def _prune_tss_pair(short_start, longer_path, enabled=True):
    """The 5' mirror: a longer model starting at 100 (path `longer_path`) and an
    abundant 3'-truncated model of the same chain, TSS call, no PolyA call, starting
    at `short_start`."""
    import LRAA_Globals

    saved = LRAA_Globals.config.get("contained_TSS_match_by_end_coords")
    LRAA_Globals.config["contained_TSS_match_by_end_coords"] = enabled
    try:
        longer = _plus_transcript(
            "t.long", [[100, 250], [300, 350], [400, 500]], longer_path, 2.0
        )
        shorter = _plus_transcript(
            "t.short",
            [[short_start, 250], [300, 320]],
            ["TSS:1", "E:1", "I:1", "E:2"],
            8.0,
        )
        frac = {"t.long": {_DummyMP(2): 1.0}, "t.short": {_DummyMP(8): 1.0}}
        retained = TranscriptFiltering.prune_likely_degradation_products(
            [longer, shorter], _PlusSpliceGraph(), frac
        )
        return {t.get_transcript_id(): t for t in retained}
    finally:
        LRAA_Globals.config["contained_TSS_match_by_end_coords"] = saved


NO_TSS = ["E:1", "I:1", "E:2", "I:2", "E:3"]


def test_TSS_less_container_absorbs_and_adopts_coincident_5prime_end():
    retained = _prune_tss_pair(99, NO_TSS)
    assert set(retained) == {"t.long"}
    survivor = retained["t.long"]
    assert survivor.get_simple_path()[0] == "TSS:1"
    assert survivor.has_TSS()
    assert survivor.get_coords() == (99, 500)

    assert set(_prune_tss_pair(99, NO_TSS, enabled=False)) == {"t.long", "t.short"}


def test_container_with_its_own_TSS_is_judged_by_identity_only():
    retained = _prune_tss_pair(99, ["TSS:2"] + NO_TSS)
    assert set(retained) == {"t.long", "t.short"}
    assert retained["t.long"].get_simple_path()[0] == "TSS:2"


def test_nested_alternative_TSS_is_not_absorbed_by_end_coords():
    import LRAA_Globals

    assert 190 - 100 > LRAA_Globals.config["max_dist_between_alt_TSS_sites"]
    retained = _prune_tss_pair(190, NO_TSS)
    assert set(retained) == {"t.long", "t.short"}
    assert not retained["t.long"].has_TSS()
