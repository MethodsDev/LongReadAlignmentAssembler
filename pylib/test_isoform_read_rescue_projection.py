#!/usr/bin/env python3

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

from IsoformReadRescue import _merge_contiguous_genomic_segments


def test_merge_contiguous_genomic_segments_collapses_split_exons():
    segments = [
        (100, 110),
        (111, 125),
        (200, 210),
        (211, 211),
        (300, 350),
    ]

    assert _merge_contiguous_genomic_segments(segments) == [
        (100, 125),
        (200, 211),
        (300, 350),
    ]


def test_merge_contiguous_genomic_segments_preserves_true_introns():
    segments = [
        (100, 150),
        (250, 300),
    ]

    assert _merge_contiguous_genomic_segments(segments) == segments
