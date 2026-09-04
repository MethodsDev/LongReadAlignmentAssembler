#!/usr/bin/env python3

"""The near-adjacent aggregation branch must weigh support, and its bound is INCLUSIVE.

`Splice_graph._prune_spurious_introns_shared_boundary` removes alternative introns that
share a donor (or acceptor) with a better-supported one. Two branches purge:

  1. support: alt/total or alt/top below `min_alt_splice_freq` (0.03)
  2. distance: the other boundary within `aggregate_splice_boundary_dist` (5 bp)

Branch 2 was support-BLIND, so a 464-read junction three bases from a 465-read junction
was deleted despite branch 1 having just passed it. `aggregate_splice_boundary_max_rel_support`
adds the missing condition: collapse only when alt/top support is at or below the bound.

The bound must be INCLUSIVE. Ties are explicitly permitted upstream -- the function
asserts `alt_support <= most_supported_support` -- so an equal-support alternative has
alt/top == 1.0 exactly. Under a strict `<`, the documented no-op default of 1.0 would
SPARE those ties while legacy purged them, silently changing behaviour at the setting
that is supposed to change nothing. Since alt support can never exceed the top's,
`<= 1.0` is always true and the default reproduces legacy exactly.

One property surfaces from writing this: on an exact tie the support sort has no
tie-break, so which member of the pair is treated as "top" -- and therefore which one
survives -- follows the iteration order of `_intron_objs` rather than anything about the
data. The tests below assert that a tie collapses to ONE junction, not which one.
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA_Globals
from GenomeFeature import Intron
from Splice_graph import Splice_graph

CONTIG = "chr1"
SHARED_DONOR = 1000
TOP_ACCEPTOR = 2000
ALT_ACCEPTOR = 2003          # 3 bp away: inside aggregate_splice_boundary_dist (5)


@pytest.fixture(autouse=True)
def _restore_config():
    keys = (
        "aggregate_adjacent_splice_boundaries",
        "aggregate_splice_boundary_dist",
        "aggregate_splice_boundary_max_rel_support",
        "min_alt_splice_freq",
        "max_intron_length",
    )
    saved = {k: LRAA_Globals.config[k] for k in keys}
    LRAA_Globals.config["aggregate_adjacent_splice_boundaries"] = True
    LRAA_Globals.config["aggregate_splice_boundary_dist"] = 5
    LRAA_Globals.config["min_alt_splice_freq"] = 0.03
    try:
        yield
    finally:
        LRAA_Globals.config.update(saved)


def _graph_with_pair(top_support, alt_support):
    """Two introns sharing a donor, acceptors 3 bp apart, with the given support."""
    sg = Splice_graph()
    sg._contig_acc = CONTIG
    sg._contig_strand = "+"
    Splice_graph._min_alt_splice_freq = LRAA_Globals.config["min_alt_splice_freq"]

    def _intron(acceptor, support):
        i = Intron(CONTIG, SHARED_DONOR, acceptor, "+", support)
        sg._intron_objs["{}:{}".format(SHARED_DONOR, acceptor)] = i
        return i

    top = _intron(TOP_ACCEPTOR, top_support)
    alt = _intron(ALT_ACCEPTOR, alt_support)
    return sg, top, alt


def _surviving_acceptors(sg):
    return {i.get_coords()[1] for i in sg._intron_objs.values()}


def test_equal_support_tie_is_purged_at_the_legacy_default():
    """alt/top == 1.0 exactly. Legacy purged it on distance alone; 1.0 must still purge.

    This is the assertion that fails under a strict `<`: 1.0 < 1.0 is false, the alt
    survives, and the "1.0 reproduces the old behaviour" contract is broken.
    """
    LRAA_Globals.config["aggregate_splice_boundary_max_rel_support"] = 1.0
    sg, _top, _alt = _graph_with_pair(top_support=465, alt_support=465)
    sg._prune_spurious_introns_shared_boundary("left")
    survivors = _surviving_acceptors(sg)
    # WHICH member survives is not asserted, and deliberately so: the support sort has
    # no tie-break, so `intron_list[-1]` on equal support is whichever iteration order
    # of _intron_objs put last. The contract legacy provided is that a tied
    # near-adjacent pair COLLAPSES TO ONE, and that is what 1.0 must preserve.
    assert len(survivors) == 1, (
        "an equal-support alternative survived at the legacy default of 1.0 "
        f"(survivors: {sorted(survivors)}): the ratio bound is exclusive where it "
        "must be inclusive, so the documented no-op setting silently spares ties the "
        "old distance-only branch purged"
    )
    assert survivors <= {TOP_ACCEPTOR, ALT_ACCEPTOR}


def test_near_tie_is_spared_at_the_new_default():
    """464 vs 465 is the case the rule exists for: 0.9978 is above 0.2, so it stays."""
    LRAA_Globals.config["aggregate_splice_boundary_max_rel_support"] = 0.2
    sg, _top, _alt = _graph_with_pair(top_support=465, alt_support=464)
    sg._prune_spurious_introns_shared_boundary("left")
    survivors = _surviving_acceptors(sg)
    assert {TOP_ACCEPTOR, ALT_ACCEPTOR} <= survivors, (
        "a near-tie was collapsed at 0.2: the distance branch is still support-blind"
    )


def test_a_genuinely_weak_alternative_is_still_collapsed():
    """The rule must not become a blanket exemption: 40/465 = 0.086 is below 0.2.

    Chosen to sit ABOVE min_alt_splice_freq (0.086 of top, 0.079 of total, both > 0.03)
    so branch 1 cannot claim it and the removal is attributable to branch 2.
    """
    LRAA_Globals.config["aggregate_splice_boundary_max_rel_support"] = 0.2
    sg, _top, _alt = _graph_with_pair(top_support=465, alt_support=40)
    sg._prune_spurious_introns_shared_boundary("left")
    survivors = _surviving_acceptors(sg)
    assert TOP_ACCEPTOR in survivors
    assert ALT_ACCEPTOR not in survivors, (
        "a 40-vs-465 alternative survived at 0.2, so near-adjacent aggregation has "
        "stopped collapsing genuinely lopsided pairs"
    )


def test_the_bound_itself_is_inclusive():
    """alt/top exactly equal to the configured bound must collapse, not survive."""
    LRAA_Globals.config["aggregate_splice_boundary_max_rel_support"] = 0.2
    sg, _top, _alt = _graph_with_pair(top_support=500, alt_support=100)  # exactly 0.2
    sg._prune_spurious_introns_shared_boundary("left")
    assert ALT_ACCEPTOR not in _surviving_acceptors(sg), (
        "alt/top exactly at the bound survived: the comparison is exclusive"
    )
