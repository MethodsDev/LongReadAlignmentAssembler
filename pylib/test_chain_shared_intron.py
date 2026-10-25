#!/usr/bin/env python3

"""chain_requires_shared_intron: chaining links two read paths only if they share an intron.

End-to-end behavior (staggered reads chained, exon-only overlaps kept apart) is covered
by testing/pasa_mode; these check the predicate and the config defaults.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA  # noqa: F401  (resolves the pylib import cycle before MultiPathGraphNode)
import LRAA_Globals
from MultiPathGraphNode import MultiPathGraphNode


def _node(simple_path):
    node = MultiPathGraphNode.__new__(MultiPathGraphNode)
    node.get_simple_path = lambda: simple_path
    return node


def test_reads_sharing_an_intron_share_it():
    a = _node(["E:1", "I:1", "E:2", "I:2", "E:3"])
    b = _node(["E:2", "I:2", "E:3", "I:3", "E:4"])
    assert a.shares_intron(b) and b.shares_intron(a)


def test_overlap_within_an_exon_is_not_a_shared_intron():
    a = _node(["E:1", "I:1", "E:2"])
    b = _node(["E:2", "I:2", "E:3"])
    assert not a.shares_intron(b)


def test_unspliced_reads_never_share_an_intron():
    se = _node(["E:2"])
    me = _node(["E:1", "I:1", "E:2"])
    assert not se.shares_intron(me) and not me.shares_intron(se)
    assert not se.shares_intron(_node(["E:2"]))


def test_defaults_keep_stock_behavior():
    assert LRAA_Globals.config["chain_requires_shared_intron"] is False
    assert LRAA_Globals.config["restrict_asm_to_collapse"] is True
    assert LRAA_Globals.config["snap_terminal_overhang_fuzz"] == 0


def test_pasa_mode_chains_on_shared_introns():
    preset = LRAA_Globals.PASA_MODE_CONFIG
    assert preset["restrict_asm_to_collapse"] is False
    assert preset["chain_requires_shared_intron"] is True
    assert set(preset) <= set(LRAA_Globals.config)
