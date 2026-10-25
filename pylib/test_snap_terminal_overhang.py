#!/usr/bin/env python3

"""snap_terminal_overhang_fuzz: PASA's terminal fuzz applied to a read before mapping.

A read end that runs at most `fuzz` bases past a splice site into the intron, while the
same terminal exon covers the exon beside that splice site, is trimmed to the splice
site. Anything else is left alone.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

from LRAA import LRAA


class _Intron:
    def __init__(self, lend, rend):
        self._coords = (lend, rend)

    def get_coords(self):
        return self._coords


class _Graph:
    def __init__(self, introns):
        self._introns = [_Intron(*i) for i in introns]

    def get_overlapping_introns(self, lend, rend):
        return [i for i in self._introns if i.get_coords()[0] <= rend and i.get_coords()[1] >= lend]


def _snap(segments, introns, fuzz=20):
    obj = LRAA.__new__(LRAA)
    obj._splice_graph = _Graph(introns)
    return obj._snap_terminal_overhangs_into_introns(segments, fuzz)


INTRON = (201, 299)  # exons ...-200 and 300-...


def test_right_end_overhanging_into_intron_is_trimmed():
    assert _snap([[150, 215]], [INTRON]) == [[150, 200]]


def test_left_end_overhanging_into_intron_is_trimmed():
    assert _snap([[285, 350]], [INTRON]) == [[300, 350]]


def test_overhang_beyond_fuzz_is_kept():
    segs = [[150, 221]]
    assert _snap(segs, [INTRON]) is segs


def test_end_not_past_a_splice_site_is_kept():
    # ends inside the intron but does not cover the exon beside it
    segs = [[210, 230]]
    assert _snap(segs, [INTRON]) is segs


def test_internal_boundaries_are_untouched():
    segs = [[100, 200], [300, 400]]
    assert _snap(segs, [INTRON, (401, 500)]) == segs


def test_input_is_not_modified():
    segs = [[150, 215]]
    _snap(segs, [INTRON])
    assert segs == [[150, 215]]
