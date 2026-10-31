#!/usr/bin/env python3

"""Trellis extension at score 0 (pasa_mode).

Nodes supported only by input-transcript templates score 0. Stock LRAA keeps an extension
only if it scores above 0; pasa_mode keeps the best compatible one even at 0, so two such
nodes chain. Without it, at AT4G14030 (Gen 1 models + Trinity) the template node with an
extra 5' exon could not join the other, got a path of its own, and the final validation
stopped the run ("overlap and are compatible - should have gotten assembled together").
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA  # noqa: F401  (resolves the pylib import cycle before MultiPathGraphNode)
import LRAA_Globals
from MultiPathGraphNode import MultiPathGraphNode
from Vertex import Vertex


class _Path:
    def __init__(self, score, incompatible=False):
        self._score, self._incompatible = score, incompatible

    def incompatibility_detected(self, mpgn):
        return self._incompatible

    def create_scored_path_extension(self, mpgn):
        return _Path(self._score)

    def get_score(self):
        return self._score

    def get_simple_path(self):
        return ["E:1"]


class _Prev:
    def __init__(self, paths):
        self._paths = paths

    def get_fromPaths(self):
        return list(self._paths)


def _vertex():
    v = Vertex.__new__(Vertex)
    node = MultiPathGraphNode.__new__(MultiPathGraphNode)
    node.get_simple_path = lambda: ["E:1"]
    v._multipath_graph_node = node
    v._fromPaths = [_Path(0)]
    return v


def _extend(prev_paths, pasa_mode):
    saved = LRAA_Globals.config.get("pasa_mode", False)
    LRAA_Globals.config["pasa_mode"] = pasa_mode
    try:
        v = _vertex()
        v.add_highest_scoring_path_extension(_Prev(prev_paths))
        return v.get_fromPaths()[1:]
    finally:
        LRAA_Globals.config["pasa_mode"] = saved


def test_zero_score_extension_kept_in_pasa_mode():
    assert len(_extend([_Path(0)], pasa_mode=True)) == 1


def test_zero_score_extension_dropped_by_default():
    assert _extend([_Path(0)], pasa_mode=False) == []


def test_incompatible_extension_never_kept():
    assert _extend([_Path(0, incompatible=True)], pasa_mode=True) == []


def test_highest_score_still_wins():
    ext = _extend([_Path(0), _Path(3), _Path(1)], pasa_mode=True)
    assert [p.get_score() for p in ext] == [3]
