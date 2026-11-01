"""terminal_boundary_coord: where an isoform end without a TSS / PolyA site goes.

Default (terminal_boundary_method_low_count = None): with fewer than
min_reads_for_terminal_adjustment reads the end stays at the terminal exon's edge.
With "median", 1 <= n < min_reads places the end at the median read end; n >= min_reads
follows terminal_boundary_method unchanged.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))

import LRAA_Globals  # noqa: E402
from LRAA import (  # noqa: E402
    own_vote_min_peak,
    extension_lacks_peak,
    read_start_vote,
    terminal_boundary_coord,
    terminal_low_count_methods,
    tss_end_by_rule,
)


def test_default_key_is_none():
    assert LRAA_Globals.config["terminal_boundary_method_low_count"] is None
    assert LRAA_Globals.config["terminal_boundary_method"] == "percentile"
    assert LRAA_Globals.config["min_reads_for_terminal_adjustment"] == 7


@pytest.mark.parametrize("low", [None, "median"])
def test_enough_reads_uses_percentile_either_way(low):
    starts = [100, 110, 120, 130, 140, 150, 160, 170, 180, 190]
    # statistics.quantiles(n=100) exclusive: 10th pct of 10 values -> 101
    assert terminal_boundary_coord(starts, "left", 50, "percentile", 90, 7, low) == 101
    ends = [1000, 1010, 1020, 1030, 1040, 1050, 1060, 1070, 1080, 1090]
    assert terminal_boundary_coord(ends, "right", 2000, "percentile", 90, 7, low) == 1089


def test_few_reads_default_keeps_exon_edge():
    assert terminal_boundary_coord([120, 130, 200], "left", 50, "percentile", 90, 7) == 50
    assert terminal_boundary_coord([900, 950], "right", 2000, "percentile", 90, 7) == 2000


def test_few_reads_median_both_sides():
    assert terminal_boundary_coord([120, 130, 200], "left", 50, "percentile", 90, 7, "median") == 130
    assert terminal_boundary_coord([900, 950, 990], "right", 2000, "percentile", 90, 7, "median") == 950
    # one read
    assert terminal_boundary_coord([140], "left", 50, "percentile", 90, 7, "median") == 140
    assert terminal_boundary_coord([1500], "right", 2000, "percentile", 90, 7, "median") == 1500


def test_few_reads_median_even_count_floors():
    assert terminal_boundary_coord([100, 105], "left", 50, "percentile", 90, 7, "median") == 102
    assert terminal_boundary_coord([1000, 1005], "right", 2000, "percentile", 90, 7, "median") == 1002


def test_median_is_clamped_to_exon_edge():
    assert terminal_boundary_coord([10, 20, 30], "left", 25, "percentile", 90, 7, "median") == 25
    assert terminal_boundary_coord([10, 20, 30], "right", 15, "percentile", 90, 7, "median") == 15


def test_no_reads_keeps_edge():
    assert terminal_boundary_coord([], "left", 50, "percentile", 90, 7, "median") == 50
    assert terminal_boundary_coord([], "right", 2000, "percentile", 90, 7, "median") == 2000


def test_extreme_method_unaffected():
    assert terminal_boundary_coord([120, 130], "left", 50, "extreme", 90, 7, "median") == 120
    assert terminal_boundary_coord([120, 130], "right", 2000, "extreme", 90, 7, "median") == 130


def test_median_everywhere_config_A():
    # terminal_boundary_method=median, min_reads=1: median at any n >= 1
    assert terminal_boundary_coord([120, 130, 200], "left", 50, "median", 90, 1) == 130
    assert terminal_boundary_coord([120], "right", 2000, "median", 90, 1) == 120


def test_unknown_low_count_value_rejected():
    with pytest.raises(ValueError):
        terminal_boundary_coord([1, 2], "left", 0, "percentile", 90, 7, "mean")


def test_matches_previous_inline_rules():
    """Byte-for-byte the rules the inline code applied before the refactor."""
    import statistics

    starts = [300, 310, 305, 400, 320, 330, 350, 301]
    lim = 290
    assert terminal_boundary_coord(starts, "left", lim, "mean", 90, 7) == max(lim, int(statistics.mean(starts)))
    assert terminal_boundary_coord(starts, "left", lim, "quartile", 90, 7) == max(
        lim, int(statistics.quantiles(starts, n=4)[0])
    )
    assert terminal_boundary_coord(starts, "right", 1000, "quartile", 90, 7) == min(
        1000, int(statistics.quantiles(starts, n=4)[2])
    )
    assert terminal_boundary_coord(starts, "left", lim, "percentile", 80, 7) == max(
        lim, int(statistics.quantiles(sorted(starts), n=100)[19])
    )


@pytest.mark.parametrize(
    "setting,strand,expected",
    [
        (None, "+", (None, None)),
        (None, "-", (None, None)),
        ("median", "+", ("median", "median")),
        ("median", "-", ("median", "median")),
        ("median_5prime", "+", ("median", None)),
        ("median_5prime", "-", (None, "median")),
    ],
)
def test_low_count_setting_per_strand(setting, strand, expected):
    assert terminal_low_count_methods(setting, strand) == expected


def test_low_count_setting_rejects_unknown_value_and_strand():
    with pytest.raises(ValueError):
        terminal_low_count_methods("mean", "+")
    with pytest.raises(ValueError):
        terminal_low_count_methods("median_5prime", ".")


def test_median_5prime_leaves_the_3prime_end_at_the_edge():
    starts, ends = [120, 130, 200], [900, 950, 990]
    left, right = terminal_low_count_methods("median_5prime", "+")
    assert terminal_boundary_coord(starts, "left", 50, "percentile", 90, 7, left) == 130
    assert terminal_boundary_coord(ends, "right", 2000, "percentile", 90, 7, right) == 2000
    left, right = terminal_low_count_methods("median_5prime", "-")
    assert terminal_boundary_coord(starts, "left", 50, "percentile", 90, 7, left) == 50
    assert terminal_boundary_coord(ends, "right", 2000, "percentile", 90, 7, right) == 950


def test_vote_needs_two_agreeing_reads():
    assert read_start_vote([]) is None
    assert read_start_vote([100]) is None
    assert read_start_vote([100, 103, 200]) is None
    assert read_start_vote([100, 102, 200]) is not None


def test_vote_window_then_exact_then_downstream():
    # 100-102 all have 3 reads in the window (2 at 200): most downstream of the tie
    assert read_start_vote([100, 101, 102, 200, 200]) == 102
    # equal windows: 200 has two exact reads, 100/101 one each
    assert read_start_vote([100, 101, 200, 200, 300]) == 200
    # equal window and exact: the more downstream position
    assert read_start_vote([100, 101, 300]) == 101


def _rule(rule, pos, own, side, edge, fallback, min_reads=7):
    return tss_end_by_rule(rule, pos, own, side, edge, fallback, min_reads)


def test_rules_defer_to_fallback_at_high_count_when_low_only():
    pos = list(range(100, 170, 10))  # n = 7
    for r in ("low_outermost_own", "low_median_own", "low_vote"):
        assert _rule(r, pos, pos, "left", 50, 111) == 111


def test_low_own_rules_use_own_reads_only():
    pos, own = [100, 150, 160], [150, 160]
    assert _rule("low_outermost_own", pos, own, "left", 50, 50) == 150
    assert _rule("low_median_own", pos, own, "left", 50, 50) == 155
    # minus strand: the 5' end is the right end, outermost = largest
    pos, own = [900, 950, 990], [900, 950]
    assert _rule("low_outermost_own", pos, own, "right", 2000, 2000) == 950


def test_no_own_reads_gives_fallback():
    assert _rule("low_median_own", [100], [], "left", 50, 50) == 50
    assert _rule("own_vote", [100, 100], [], "right", 2000, 2000) == 2000


def test_low_vote_falls_back_when_no_agreement():
    assert _rule("low_vote", [100, 150], [], "left", 50, 50) == 50
    assert _rule("low_vote", [100, 101], [], "left", 50, 50) == 101


def test_vote_is_downstream_on_minus_strand():
    # right end: reads at 1000 and 999 tie, the more downstream is the smaller coordinate
    assert _rule("vote", [1000, 999, 800], [], "right", 2000, 2000) == 999


def test_outer30_vote_then_20th_percentile():
    # outermost 3 of 10 are 100, 101, 150: vote 101
    pos = [100, 101, 150, 160, 170, 180, 190, 200, 210, 220]
    assert _rule("outer30_vote", pos, [], "left", 50, 50) == 101
    # no vote in the outermost 3: 20th percentile (lower) of the 10 reads, index 1
    pos = [100, 110, 120, 160, 170, 180, 190, 200, 200, 220]
    assert _rule("outer30_vote", pos, [], "left", 50, 50) == 110


def test_own_vote_then_20th_percentile_of_own():
    pos = [90, 91, 100, 120, 140]
    own = [100, 120, 140]
    assert _rule("own_vote", pos, own, "left", 50, 50) == 100
    assert _rule("own_vote", pos, [100, 100, 140], "left", 50, 50) == 100
    # right end, n = 3, index 0 from the 5' (largest) side
    assert _rule("own_vote", pos, [100, 120, 140], "right", 2000, 2000) == 140


def test_own_vote_min_peak_size():
    # min(3 reads, 20% of n), never below 2
    assert own_vote_min_peak(5, 3, 0.20) == 2
    assert own_vote_min_peak(17, 3, 0.20) == 3
    assert own_vote_min_peak(17, None, None) == 2
    # 17 own reads, a chance 2-read pile-up far inside: no longer a peak, 20th percentile
    own = [100, 130, 160, 190, 220, 250, 280, 310, 340, 370, 400, 430, 460, 490, 2500, 2501, 520]
    assert tss_end_by_rule("own_vote", own, own, "left", 50, 50) == 2501
    assert (
        tss_end_by_rule("own_vote", own, own, "left", 50, 50, min_peak_reads=3, min_peak_frac=0.2)
        == 190
    )
    # a 3-read peak still wins
    own3 = own + [2502]
    assert (
        tss_end_by_rule("own_vote", own3, own3, "left", 50, 50, min_peak_reads=3, min_peak_frac=0.2)
        == 2502
    )


def test_own_vote_shared_first_intron_fallback():
    # minus strand (5' = right end): own reads 5'-truncated, no peak; the sibling's capped
    # reads at 1661 vote before the own 20th percentile (1584)
    own = [1702, 1584, 1500, 1450, 1400, 1300, 1200, 1100]
    shared = own + [1661, 1661, 1661, 1660, 1500, 1300]
    assert tss_end_by_rule("own_vote", own, own, "right", 1716, 1716) == 1584
    kw = dict(min_peak_reads=3, min_peak_frac=0.2)
    assert tss_end_by_rule("own_vote", own, own, "right", 1716, 1716, **kw) == 1584
    assert (
        tss_end_by_rule("own_vote", own, own, "right", 1716, 1716, shared_positions=shared, **kw)
        == 1661
    )
    # an own peak keeps priority over the shared reads
    own_pk = own + [1450, 1451]
    assert (
        tss_end_by_rule("own_vote", own_pk, own_pk, "right", 1716, 1716, shared_positions=shared, **kw)
        == 1450
    )
    # shared reads with no peak: own 20th percentile as before
    assert (
        tss_end_by_rule("own_vote", own, own, "right", 1716, 1716, shared_positions=[1690, 1600], **kw)
        == 1584
    )


def test_v048_own_vote_peak_defaults():
    assert LRAA_Globals.config["terminal_end_rule_own_vote_min_peak_reads"] == 3
    assert LRAA_Globals.config["terminal_end_rule_own_vote_min_peak_frac"] == 0.20
    assert LRAA_Globals.config["terminal_end_rule_own_vote_shared_first_intron"] is True


def test_rule_clamped_to_exon_and_unknown_rejected():
    assert _rule("vote", [10, 10], [], "left", 25, 25) == 25
    with pytest.raises(ValueError):
        _rule("nope", [1], [], "left", 0, 0)


def test_v048_end_rule_defaults():
    # v0.48.0: own_vote at the 5' end with extension and the 200 bp no-peak cap; 3' unchanged
    assert LRAA_Globals.config["terminal_boundary_percentile_TSS"] is None
    assert LRAA_Globals.config["terminal_TSS_end_rule"] == "own_vote"
    assert LRAA_Globals.config["terminal_end_rule_extend"] is True
    assert LRAA_Globals.config["terminal_end_rule_extend_unvoted_max"] == 200


def test_percentile_tss_key_mirrors_per_side():
    starts = list(range(100, 200, 10))
    assert terminal_boundary_coord(starts, "left", 0, "percentile", 80, 7) == 112
    assert terminal_boundary_coord(starts, "right", 1000, "percentile", 80, 7) == 178


def test_polyA_key_default_off():
    assert LRAA_Globals.config["terminal_PolyA_end_rule"] is None


def test_own_vote_mirrored_at_3prime_end():
    # + strand 3' end is the right end: outer = largest; no vote -> 20% from the outer end
    own = [100, 200, 300, 400, 500, 600]
    assert _rule("own_vote", own, own, "right", 2000, 2000) == 500
    # one read: that read
    assert _rule("own_vote", [700], [700], "right", 2000, 2000) == 700
    # peak wins over the percentile
    assert _rule("own_vote", [100, 300, 301, 600], [100, 300, 301, 600], "right", 2000, 2000) == 300


def test_extend_goes_past_the_exon_edge_only_when_asked():
    # minus-strand 5' end: right end, exon edge 1000, three own reads reach 1219
    own = [1219, 1219, 1210]
    assert tss_end_by_rule("own_vote", own, own, "right", 1000, 1000) == 1000
    assert tss_end_by_rule("own_vote", own, own, "right", 1000, 1000, extend=True) == 1219
    # left end mirrored
    own = [80, 80, 300]
    assert tss_end_by_rule("own_vote", own, own, "left", 100, 100, extend=True) == 80


def test_extend_stray_read_does_not_set_the_end_when_a_peak_exists():
    own = [1500, 950, 951, 930, 920]
    assert tss_end_by_rule("own_vote", own, own, "right", 1000, 1000, extend=True) == 950


class _FakeIntron:
    def __init__(self, c):
        self._c = c

    def get_coords(self):
        return self._c


class _FakeMP:
    def __init__(self, introns, names):
        self._i, self._n = [_FakeIntron(c) for c in introns], names

    def get_introns(self):
        return self._i

    def get_read_names(self):
        return self._n


class _FakeTranscript:
    def __init__(self, mps):
        self._m = mps

    def get_multipaths_evidence_assigned(self):
        return self._m

    def get_transcript_id(self):
        return "t"


def _fake_names(src, tid, span, pooled=0, k=2):
    return [
        f"fake_for_merge|src={src}|tid={tid}|flags=-|n={k}|uid={i}|span={span}|pooled={pooled}"
        for i in range(1, k + 1)
    ]


def test_source_ends_count_each_model_once_and_pooled_as_residual():
    from LRAA import own_chain_source_model_ends

    chain = [(200, 300)]
    t = _FakeTranscript(
        [
            _FakeMP(chain, _fake_names("a.gtf", "t1", "100-500")),
            _FakeMP(chain, _fake_names("b.gtf", "t2", "150-500")),
            _FakeMP(chain, _fake_names("p.gtf", "t3", "101-700", pooled=1)),
            _FakeMP([(200, 300), (350, 400)], _fake_names("c.gtf", "t4", "10-900")),
        ]
    )
    ends = own_chain_source_model_ends(t, chain)
    assert sorted(ends, key=str) == sorted([(100, 500), (150, 500), (None, 700)], key=str)


def test_merge_end_keys_default_on():
    for key in (
        "merge_end_rule_from_sources",
        "merge_trim_from_model_end",
        "merge_repath_before_end_adjustment",
        "merge_end_stop_at_called_site",
        "merge_dedupe_identical_models",
        "merge_site_end_requires_source_support",
    ):
        assert LRAA_Globals.config[key] is True


def test_source_ends_any_chain_when_no_chain_given():
    from LRAA import own_chain_source_model_ends

    t = _FakeTranscript(
        [
            _FakeMP([(200, 300)], _fake_names("a.gtf", "t1", "100-500")),
            _FakeMP([(200, 300), (350, 400)], _fake_names("c.gtf", "t4", "10-900")),
        ]
    )
    assert sorted(own_chain_source_model_ends(t)) == [(10, 900), (100, 500)]


@pytest.mark.parametrize(
    "end,new_end,want",
    [
        (100, 80, 90),  # outward, stops on the site
        (100, 120, 110),  # inward, stops on the site
        (100, 60, 90),  # the first site met, not the farthest
        (100, 90, 90),  # landing on a site is not crossing it
        (100, 95, 95),  # no site between
        (90, 80, 80),  # a site at the end itself does not hold it
        (100, 100, 100),
        (200, 300, 300),  # past the last site
    ],
)
def test_stop_end_at_site(end, new_end, want):
    from LRAA import stop_end_at_site

    assert stop_end_at_site(end, new_end, [70, 90, 110, 130]) == want


class _Node:
    def __init__(self, coords):
        self._coords = coords

    def get_coords(self):
        return self._coords


class _Graph:
    def __init__(self, strand, nodes):
        self._strand, self._nodes = strand, nodes

    def get_contig_strand(self):
        return self._strand

    def get_node_obj_via_id(self, node_id):
        return _Node(self._nodes[node_id])


@pytest.mark.parametrize("site,window", [("TSS", "max_dist_between_alt_TSS_sites"), ("POLYA", "max_dist_between_alt_polyA_sites")])
def test_full_window_trim_measures_from_the_model_end(site, window):
    """MEASURED chr22: a model 42 bp from a site, in a node starting 63 bp from it."""
    import Simple_path_utils as SPU

    w = LRAA_Globals.config[window]
    edge_dist = w + 13
    own_dist = w - 8
    pos = 1000
    if site == "TSS":  # + strand: the model's left end
        graph = _Graph("+", {"E:1": (pos - edge_dist, pos - 1), "TSS:1": (pos, pos), "E:2": (pos, pos + 100)})
        path = ["E:1", "TSS:1", "E:2"]
        refine, trimmed, own = SPU.refine_TSS_simple_path, ["TSS:1", "E:2"], (pos - own_dist, pos + 100)
    else:  # + strand: the model's right end
        graph = _Graph("+", {"E:1": (pos - 100, pos), "POLYA:1": (pos, pos), "E:2": (pos + 1, pos + edge_dist)})
        path = ["E:1", "POLYA:1", "E:2"]
        refine, trimmed, own = SPU.refine_PolyA_simple_path, ["E:1", "POLYA:1"], (pos - 100, pos + own_dist)
    assert refine(graph, list(path), full_window=True) == [n for n in path if n in ("E:1", "E:2")]
    assert refine(graph, list(path), full_window=True, own_end=own) == trimmed


def test_collapse_identical_models_keeps_the_evidence():
    from LRAA import collapse_identical_models

    class T:
        def __init__(self, exons, mps):
            self.exons, self.mps = exons, list(mps)

        def get_exon_segments(self):
            return self.exons

        def get_multipaths_evidence_assigned(self):
            return self.mps

        def add_multipaths_evidence_assigned(self, mps):
            self.mps.extend(mps)

    a, b, c = T([[1, 10], [20, 30]], ["m1"]), T([[1, 10], [20, 30]], ["m2"]), T([[2, 10], [20, 30]], ["m3"])
    assert collapse_identical_models([a, b, c]) == [a, c]
    assert a.mps == ["m1", "m2"]


def test_source_end_coords_include_pooled_and_every_chain():
    from LRAA import source_model_end_coords

    t = _FakeTranscript(
        [
            _FakeMP([(200, 300)], _fake_names("a.gtf", "t1", "100-500")),
            _FakeMP([(200, 300)], _fake_names("p.gtf", "t2", "101-700", pooled=1)),
            _FakeMP([(350, 400)], _fake_names("b.gtf", "t3", "10-900")),
            _FakeMP([(350, 400)], ["reftranscript:x", "read1"]),
        ]
    )
    assert source_model_end_coords(t) == ({100, 101, 10}, {500, 700, 900})


@pytest.mark.parametrize(
    "site,ends,want",
    [
        (1000, {1000}, True),
        (1000, {1050}, True),
        (1000, {950}, True),
        (1000, {1051}, False),
        (1000, {949, 2000}, False),
        (1000, set(), False),
    ],
)
def test_site_end_supported_window(site, ends, want):
    from LRAA import site_end_supported

    assert site_end_supported(site, ends, 50) is want


def _gated_path(path, nodes, spans, remap=None, strand="-"):
    """Runs the check on `path`; remap: what mapping without the trim returns."""
    from LRAA import LRAA

    lraa = LRAA.__new__(LRAA)
    lraa._splice_graph = _Graph(strand, nodes)
    calls = []

    def map_read(segments, refine_TSS_simple_path, refine_PolyA_simple_path, full_window_boundary_trim):
        calls.append((refine_TSS_simple_path, refine_PolyA_simple_path))
        return list(path) if remap is None else remap

    lraa._map_read_to_graph = map_read
    mps = [_FakeMP([], _fake_names("a.gtf", f"t{i}", s)) for i, s in enumerate(spans)]
    return lraa._path_without_unsupported_site_ends(_FakeTranscript(mps), [[1, 2]], list(path), True), calls


TIGD7_LIKE = {
    "E:1": (100, 200),
    "E:2": (300, 728),
    "TSS:1": (728, 728),
    "POLYA:1": (100, 100),
}


def test_site_end_past_all_sources_comes_off_the_path():
    path, calls = _gated_path(["E:1", "E:2", "TSS:1"], TIGD7_LIKE, ["100-431", "100-410"])
    assert path == ["E:1", "E:2"]
    assert calls == [(False, True)]


def test_site_end_with_a_source_in_the_window_stays():
    path, calls = _gated_path(["E:1", "E:2", "TSS:1"], TIGD7_LIKE, ["100-431", "100-678"])
    assert path == ["E:1", "E:2", "TSS:1"]
    assert calls == []


def test_polyA_end_gated_with_its_own_side_and_tss_trim_kept():
    path, calls = _gated_path(["POLYA:1", "E:1", "E:2"], TIGD7_LIKE, ["400-431", "160-800"])
    assert path == ["E:1", "E:2"]
    assert calls == [(True, False)]


def test_trim_onto_unsupported_site_is_undone_without_moving_the_model():
    # the trimmed path ended on the site; untrimmed the site is interior, so no flag
    nodes = dict(TIGD7_LIKE, **{"E:3": (729, 750)})
    path, _ = _gated_path(
        ["E:1", "E:2", "TSS:1"], nodes, ["100-431"], remap=["E:1", "E:2", "TSS:1", "E:3"]
    )
    assert path == ["E:1", "E:2", "TSS:1", "E:3"]


def test_model_without_source_names_is_left_alone():
    path, calls = _gated_path(["E:1", "E:2", "TSS:1"], TIGD7_LIKE, [])
    assert path == ["E:1", "E:2", "TSS:1"]
    assert calls == []


def test_strip_needs_an_exon_left():
    path, _ = _gated_path(["TSS:1", "POLYA:1"], TIGD7_LIKE, ["400-431"])
    assert path == ["TSS:1", "POLYA:1"]


def test_stop_does_not_stop_on_the_site_the_end_leaves():
    from LRAA import stop_end_at_site

    assert stop_end_at_site(728, 431, [728]) == 431


def test_long_extension_needs_a_read_peak():
    assert extension_lacks_peak(700, 1000, [700, 950, 990], 200)
    assert not extension_lacks_peak(700, 1000, [700, 702, 990], 200)
    assert not extension_lacks_peak(900, 1000, [900], 200)
    assert extension_lacks_peak(1300, 1000, [1300], 200)


@pytest.mark.parametrize(
    "key",
    [
        "terminal_PolyA_outward_pileup",
        "terminal_end_rule_own_chain_strict",
        "terminal_trim_from_model_end",
    ],
)
def test_rejected_experiment_keys_are_gone(key):
    # tested on the v0.48.0 candidate and not adopted; the code is removed with the key
    assert key not in LRAA_Globals.config


# SPACER revert (terminal_end_rule_extend guard). Model exon edges 1000 (first exon lend)
# and 5000 (last exon rend); the rule put the ends at 400 (extended) and 5600 (extended),
# its in-exon coordinates being 1050 and 4950.


def _revert(new_first, new_last, mappable, adjust_left=True, adjust_right=True, **kw):
    from LRAA import revert_unmappable_extension

    calls = []

    def _m(lend, rend):
        calls.append((lend, rend))
        return mappable(lend, rend)

    got = revert_unmappable_extension(
        new_first, new_last, (1000, 5000), (1050, 4950), adjust_left, adjust_right, _m, **kw
    )
    return got, calls


def test_spacer_revert_keeps_a_mappable_extension():
    got, calls = _revert(400, 5000, lambda l, r: True)
    assert got == (400, 5000)
    assert calls == [(400, 5000)]


def test_spacer_revert_undoes_an_unmappable_left_extension():
    # the extended first exon reaches a segment the exon has no edge to
    got, _ = _revert(400, 5000, lambda l, r: l >= 1000)
    assert got == (1050, 5000)


def test_spacer_revert_never_leaves_the_end_past_the_edge():
    # the rule's unextended coordinate is the edge itself or beyond: clamp to the edge
    from LRAA import revert_unmappable_extension

    got = revert_unmappable_extension(
        400, 5000, (1000, 5000), (900, 5100), True, False, lambda l, r: l >= 1000
    )
    assert got == (1000, 5000)


def test_spacer_revert_right_end_only_when_still_unmappable():
    # left revert alone fixes it: the right extension stays
    got, _ = _revert(400, 5600, lambda l, r: l >= 1000)
    assert got == (1050, 5600)
    # both ends unmappable: both revert
    got, _ = _revert(400, 5600, lambda l, r: l >= 1000 and r <= 5000)
    assert got == (1050, 4950)


def test_spacer_revert_ignores_ends_that_did_not_extend():
    # inward moves are never tested or changed, even when the model is unmappable
    got, calls = _revert(1050, 4950, lambda l, r: False)
    assert got == (1050, 4950)
    assert calls == []
    # an extension on a side that is not adjusted (site-anchored end) is not ours to revert
    got, calls = _revert(400, 5000, lambda l, r: False, adjust_left=False)
    assert got == (400, 5000) and calls == []


def test_spacer_revert_stops_on_a_called_site():
    # merge_end_stop_at_called_site: the reverted end still may not cross a site
    got, _ = _revert(400, 5000, lambda l, r: l >= 1000, left_sites=[1020])
    assert got == (1020, 5000)
