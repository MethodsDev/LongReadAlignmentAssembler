#!/usr/bin/env python3

"""merge_LRAA_GTFs.py --cpu: the groups run in parallel, the output order does not move.

The (contig, strand) groups of a merge are independent -- each reads only its own input
transcripts and its contig's sequence -- so ``--cpu`` runs them in forked workers. What
that must not change is what comes out: the parent writes results in the original group
order, whichever worker finishes first, and a failure in a worker is a failure of the
merge rather than a silent gap in the output.

The group body itself (splice graph, isoform reconstruction) is replaced here by a stub
that takes longer for LARGER groups, because the pool submits the largest first: the
completion order is then the reverse of the order the results have to come back in,
which is the case an ordering bug would show up in.
"""

import sys
import time
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "util"))

import merge_LRAA_GTFs as merge  # noqa: E402


@pytest.fixture()
def groups():
    # sizes 1..6 in key order, so submission order (largest first) is the exact reverse
    return {"chr{}^+".format(i): list(range(i)) for i in range(1, 7)}


def _stub(key, transcripts, genome):
    time.sleep(0.04 * len(transcripts))
    return ["gtf:{}:{}".format(key, len(transcripts))], [{"key": key}]


@pytest.mark.parametrize("workers", [1, 2, 3, 8])
def test_results_come_back_in_group_order_whatever_finishes_first(
    monkeypatch, groups, workers
):
    monkeypatch.setattr(merge, "_merge_one_group", _stub)

    out = list(merge._iter_group_results(groups, "genome.fa", workers))

    assert [key for key, _ in out] == list(groups)
    for key, (gtf_lines, rows) in out:
        assert gtf_lines == ["gtf:{}:{}".format(key, len(groups[key]))]
        assert rows == [{"key": key}]


def test_the_pool_returns_what_the_serial_path_returns_for_a_deterministic_group(
    monkeypatch, groups
):
    """Plumbing only: the stub group is a pure function of its input.

    This does NOT claim a real merge is byte-identical between --cpu 1 and --cpu N. It is
    not: the serial pass carries state from earlier contigs into later ones, and on a real
    14-input merge that moved two tied isoforms of one gene (8 of 1.08 M gtf lines), where
    the pooled result matches a fresh single-group run exactly. The pool is also identical
    across worker counts (6 vs 8, byte for byte), which is what a caller can rely on.
    """

    monkeypatch.setattr(merge, "_merge_one_group", _stub)

    serial = list(merge._iter_group_results(groups, "genome.fa", 1))
    parallel = list(merge._iter_group_results(groups, "genome.fa", 4))

    assert serial == parallel


def test_a_failing_group_fails_the_merge(monkeypatch, groups):
    def boom(key, transcripts, genome):
        if key == "chr3^+":
            raise RuntimeError("splice graph failed for chr3")
        return _stub(key, transcripts, genome)

    monkeypatch.setattr(merge, "_merge_one_group", boom)

    with pytest.raises(RuntimeError, match="chr3"):
        list(merge._iter_group_results(groups, "genome.fa", 3))


def test_cpu_defaults_to_the_serial_pass():
    """Unset must mean the historical single pass, so existing callers are unchanged."""

    import argparse
    import inspect

    src = inspect.getsource(merge.main)
    assert '"--cpu"' in src
    assert "default=1" in src[src.index('"--cpu"') : src.index('"--cpu"') + 200]
    assert "max(1, args.cpu)" in src
    # the task reads the cores of the machine it picked, not the number it asked for
    wdl = (REPO / "WDL" / "LRAA-cell_cluster_guided.wdl").read_text()
    assert "--cpu ~{c3d_effective_cpu}" in wdl
    del argparse
