#!/usr/bin/env python3

"""The two subtle correctness cases for Quantify.quantify(prior=True).

run_transcript_assembly (LRAA) calls one draft-mode Quantify object's quantify()
three times against the same mp_counter, over a shrinking transcript set. Calls 2
and 3 pass prior=True so that a read-sharing component untouched by the intervening
filter step is carried forward from the prior call instead of being reassigned and
re-EM'd from scratch. Two things have to hold for that to be correct, and both are
easy to get wrong silently (see Quantify.quantify's prior-capture comments):

(a) an isoform removed with ZERO assigned multipaths must not trigger its
    component's rerun at all -- affected_mps (and therefore dirty_gene_ids) must
    stay empty, which this pins with a spy on the cascade and the estimator, not
    just an output comparison an accidentally-identical rerun could also produce.

(b) a carried-forward transcript's multipaths_evidence_assigned /
    _multipaths_evidence_weights / read_counts_assigned must survive an unrelated
    component's rerun untouched. Transcript.init_quant_info() resets those three
    fields on the transcript OBJECT ITSELF, unconditionally, if quantify() calls it
    without checking dirty_gene_ids first -- and a carried-forward component's
    multipaths are never reprocessed to repopulate them.
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA_Globals
from GenomeFeature import Exon
from MultiPath import MultiPath
from MultiPathCounter import MultiPathCounter
from Quantify import Quantify
from Splice_graph import Splice_graph


EXON_COORDS = [(100, 150), (200, 250), (300, 350), (400, 450), (500, 550), (600, 650)]


@pytest.fixture(autouse=True)
def _quiet_and_deterministic_config():
    saved = {
        key: LRAA_Globals.config[key]
        for key in (
            "aggressively_assign_reads",
            "show_progress_quant_assign",
            "fraction_read_align_overlap",
        )
    }
    LRAA_Globals.config["aggressively_assign_reads"] = False
    LRAA_Globals.config["show_progress_quant_assign"] = False
    LRAA_Globals.config["fraction_read_align_overlap"] = 0.4
    try:
        yield
    finally:
        LRAA_Globals.config.update(saved)


def _graph():
    """A splice graph of 6 disjoint exon "islands", far enough apart that no two
    transcripts built from different pairs ever share a splice-graph node."""

    Exon.reset_counter()
    splice_graph = Splice_graph()
    exons = [Exon("chr1", lend, rend, "+", 1) for lend, rend in EXON_COORDS]
    for exon in exons:
        splice_graph._node_id_to_node[exon.get_id()] = exon
    splice_graph._contig_acc = "chr1"
    splice_graph._contig_strand = "+"
    node_ids = [exon.get_id() for exon in exons]
    return splice_graph, node_ids


def _transcript(splice_graph, node_ids, exon_idxs, gene_id, transcript_id):
    transcript = MultiPath(
        splice_graph, [[node_ids[i] for i in exon_idxs]]
    ).toTranscript()
    transcript.set_gene_id(gene_id)
    transcript.set_transcript_id(transcript_id)
    return transcript


def _counter(splice_graph, node_ids, reads):
    """reads: iterable of (read_name, exon_idxs) pairs, one multipath each."""

    mp_counter = MultiPathCounter()
    for read_name, exon_idxs in reads:
        mp_counter.add(
            MultiPath(
                splice_graph,
                [[node_ids[i] for i in exon_idxs]],
                read_types={"PacBio"},
                read_names={read_name},
                read_count=1,
            )
        )
    return mp_counter


def test_zero_support_isoform_removal_skips_its_components_rerun(monkeypatch):
    """(a) removing a zero-support isoform must not trigger ANY component's rerun."""

    sg, node_ids = _graph()
    t1a = _transcript(sg, node_ids, [0, 1], "g1", "t1a")
    t1b = _transcript(sg, node_ids, [2, 3], "g1", "t1b")
    t2 = _transcript(sg, node_ids, [4, 5], "g2", "t2")

    # readA is only compatible with t1a; nothing is compatible with t1b at all, so
    # t1b carries zero assigned multipaths after call 1 -- exactly the routine
    # filter_novel_isoforms_by_min_read_support case the brief calls out.
    counter = _counter(sg, node_ids, [("readA", [0, 1]), ("readB", [4, 5])])

    q = Quantify(run_EM=True, max_EM_iterations=10)
    q.quantify(sg, [t1a, t1b, t2], counter)

    assert t1b.get_read_counts_assigned() == 0, "fixture must give t1b zero support"
    assert t1b.get_multipaths_evidence_assigned() == [], (
        "fixture must give t1b zero assigned multipaths"
    )

    cascade_calls = []
    real_resolve = Quantify.resolve_mp_to_transcripts

    def spy_resolve(self, *args, **kwargs):
        cascade_calls.append(args)
        return real_resolve(self, *args, **kwargs)

    monkeypatch.setattr(Quantify, "resolve_mp_to_transcripts", spy_resolve)

    estimator_calls = []
    real_estimator = Quantify._estimate_isoform_read_support

    def spy_estimator(self, transcripts, **kwargs):
        estimator_calls.append(sorted(t.get_transcript_id() for t in transcripts))
        return real_estimator(self, transcripts, **kwargs)

    monkeypatch.setattr(Quantify, "_estimate_isoform_read_support", spy_estimator)

    # t1b dropped; t1a and t2 (and every multipath) are untouched.
    frac_read_assignments = q.quantify(sg, [t1a, t2], counter, prior=True)

    assert cascade_calls == [], (
        "removing a zero-support isoform must not run the compatibility cascade for "
        "ANY multipath -- affected_mps has to be empty, not merely small"
    )
    assert estimator_calls == [], (
        "removing a zero-support isoform must not run EM/read-support estimation for "
        "ANY component -- an unaffected component must be carried forward rather than "
        "recomputed to a coincidentally-identical answer"
    )

    # The carried-forward answer is still the right one.
    assert frac_read_assignments["t1a"], "t1a's carried-forward assignment is empty"
    assert frac_read_assignments["t2"], "t2's carried-forward assignment is empty"
    assert "t1b" not in frac_read_assignments


def test_carried_forward_transcript_state_survives_a_sibling_components_rerun():
    """(b) an unrelated component's rerun must not corrupt a carried-forward transcript."""

    sg, node_ids = _graph()
    t1a = _transcript(sg, node_ids, [0, 1], "g1", "t1a")
    t1b = _transcript(sg, node_ids, [2, 3], "g1", "t1b")
    t2 = _transcript(sg, node_ids, [4, 5], "g2", "t2")

    # This time t1b DOES carry real support, so removing it genuinely dirties g1's
    # own (singleton) component -- proving the guard is scoped per component, not
    # merely "nothing was removed so nothing resets" as in the zero-support case.
    counter = _counter(
        sg, node_ids, [("readA", [0, 1]), ("readB1", [2, 3]), ("readB2", [4, 5])]
    )

    q = Quantify(run_EM=True, max_EM_iterations=10)
    q.quantify(sg, [t1a, t1b, t2], counter)

    t2_read_counts_before = t2.get_read_counts_assigned()
    t2_mps_before = list(t2.get_multipaths_evidence_assigned())
    t2_weights_before = {mp: t2.get_multipath_weight(mp) for mp in t2_mps_before}
    assert t2_read_counts_before, "fixture must give t2 nonzero read support"
    assert t2_mps_before, "fixture must give t2 at least one assigned multipath"

    # g2 (t2) never shares a multipath with g1, so it is a separate component; g1's
    # component is the only one that should rerun here.
    q.quantify(sg, [t1a, t2], counter, prior=True)

    assert t2.get_multipaths_evidence_assigned() == t2_mps_before, (
        "a carried-forward transcript's multipaths_evidence_assigned must survive an "
        "unrelated component's rerun untouched, not be wiped by an unguarded "
        "init_quant_info() call"
    )
    for mp in t2_mps_before:
        assert t2.get_multipath_weight(mp) == t2_weights_before[mp], (
            "a carried-forward transcript's _multipaths_evidence_weights must survive"
        )
    assert t2.get_read_counts_assigned() == t2_read_counts_before, (
        "a carried-forward transcript's read_counts_assigned must survive"
    )

    # g1's own (genuinely affected) component was in fact recomputed, not silently
    # skipped too: t1a keeps its own real support.
    assert t1a.get_read_counts_assigned() > 0
    assert [mp.get_id() for mp in t1a.get_multipaths_evidence_assigned()]
