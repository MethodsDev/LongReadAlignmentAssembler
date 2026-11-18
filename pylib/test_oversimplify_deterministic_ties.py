#!/usr/bin/env python3
"""--oversimplify breaks best-overlap ties deterministically.

A read whose exonic overlap ties between several models used to go to random.choice
(unseeded), so reruns of the same input differed in quant.expr / quant.tracking, and
an on/off comparison of any config key differed on chrM. The tie is now broken by
_oversimplify_tie_break, which depends only on the read name and the tied ids:

  1. the same read and the same tied set always give the same model, whatever order
     the ids arrive in (and in a fresh process: no salted hash());
  2. a sole candidate is returned as is;
  3. tied reads are spread over the tied models, not all sent to one;
  4. end to end, two runs (models supplied in opposite orders) write identical
     tracking and quant rows.
"""

import importlib.util
import subprocess
import sys
from importlib.machinery import SourceFileLoader
from io import StringIO
from pathlib import Path

import LRAA_Globals
from Transcript import Transcript


LRAA_PATH = Path(__file__).resolve().parents[1] / "LRAA"


def _load_lraa_module():
    loader = SourceFileLoader("lraa_oversimplify_ties_test", str(LRAA_PATH))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def _tx(contig, exons, gene, tid, strand="+"):
    t = Transcript(contig, [list(e) for e in exons], strand)
    t.set_gene_id(gene)
    t.set_transcript_id(tid)
    return t


def test_tie_break_is_order_independent_and_stable():
    lraa = _load_lraa_module()
    tids = ["t3", "t1", "t2"]
    for name in ("bc^umi^r1", "bc^umi^r2", "x"):
        pick = lraa._oversimplify_tie_break(name, tids)
        assert pick in tids
        assert lraa._oversimplify_tie_break(name, list(reversed(tids))) == pick
        assert lraa._oversimplify_tie_break(name, sorted(tids)) == pick
    assert lraa._oversimplify_tie_break("anything", ["only"]) == "only"


def test_tie_break_is_stable_across_processes():
    # Python's str hash() is salted per process; the pick must not depend on it.
    code = (
        "import sys; from importlib.machinery import SourceFileLoader; "
        "import importlib.util; "
        f"l = SourceFileLoader('m', {str(LRAA_PATH)!r}); "
        "s = importlib.util.spec_from_loader(l.name, l); "
        "m = importlib.util.module_from_spec(s); l.exec_module(m); "
        "print(','.join(m._oversimplify_tie_break(f'r{i}', ['a','b','c']) for i in range(50)))"
    )
    pylib = str(Path(__file__).resolve().parent)
    outs = set()
    for seed in ("1", "2"):
        env = {"PYTHONHASHSEED": seed, "PYTHONPATH": pylib, "PATH": "/usr/bin:/bin"}
        out = subprocess.run(
            [sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True
        ).stdout.strip().splitlines()[-1]
        outs.add(out)
    assert len(outs) == 1


def test_tied_reads_are_spread_over_models():
    lraa = _load_lraa_module()
    picks = [lraa._oversimplify_tie_break(f"bc^u^read{i}", ["a", "b"]) for i in range(400)]
    assert 120 < picks.count("a") < 280


def _run(lraa, monkeypatch, transcripts):
    class _Alignment:
        def __init__(self, name, segs):
            self._name, self._segs = name, segs

        def get_pretty_alignment_segments(self):
            return self._segs

        def get_read_name(self):
            return self._name

    def _alignments(*args, **kwargs):
        # every read lies on the shared exon, so it ties between t_a and t_b
        for i in range(40):
            yield _Alignment(f"bc^u^r{i}", [[20, 90]])

    monkeypatch.setattr(lraa, "_oversimplify_alignments", _alignments)
    quant, track = StringIO(), StringIO()
    lraa._run_oversimplify_best_overlap(
        "chrM", "+", "A" * 3000, "unused.bam", None, None, transcripts, quant, track
    )
    return quant.getvalue(), track.getvalue()


def test_oversimplify_runs_are_identical(monkeypatch, tmp_path):
    lraa = _load_lraa_module()
    monkeypatch.setenv("LRAA_TMP_DIR", str(tmp_path))
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 1)

    def models():
        return [
            _tx("chrM", [(1, 100), (300, 400)], "MT-X", "t_a"),
            _tx("chrM", [(1, 100), (500, 600)], "MT-X", "t_b"),
        ]

    q1, t1 = _run(lraa, monkeypatch, models())
    q2, t2 = _run(lraa, monkeypatch, list(reversed(models())))
    assert t1 == t2
    assert q1 == q2
    tids = [row.split("\t")[1] for row in t1.splitlines()]
    assert set(tids) == {"t_a", "t_b"}  # spread, not all to one
    # tied reads are never marked unique
    assert {row.split("\t")[8] for row in t1.splitlines()} == {"0"}
