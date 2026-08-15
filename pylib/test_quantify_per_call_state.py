#!/usr/bin/env python3

"""Quantify's per-call state must not survive into the next call.

`_mp_to_transcripts`, `_unassigned_mp_count_pairs`, `_read_name_to_multipath` and
`_transcript_id_to_theta` are all populated by the per-gene loops inside `quantify()`.
Accumulating across genes within one call is intended; accumulating across calls is not.

Discovery quantifies several times over a shrinking transcript set, and a streaming
assignment pass reads `get_mp_to_transcripts()` together with `get_transcript_id_to_theta()`
to build its lookup table. If the multipaths outlive a call while theta is refreshed, that
table pairs multipaths from an earlier pass with this pass's abundances -- a plausible table
built from two different runs, which is the shape of error this codebase produces most
easily.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import LRAA_Globals
from GenomeFeature import Exon
from MultiPath import MultiPath
from MultiPathCounter import MultiPathCounter
from Quantify import Quantify
from Splice_graph import Splice_graph
from Transcript import Transcript

CONTIG = "chr_test"


def _graph():
    sg = Splice_graph()
    Exon.reset_counter()
    sg._contig_acc = CONTIG
    sg._contig_strand = "+"
    exons = [
        Exon(CONTIG, 100, 200, "+", 10.0),
        Exon(CONTIG, 300, 400, "+", 10.0),
        Exon(CONTIG, 500, 600, "+", 10.0),
    ]
    for e in exons:
        sg._node_id_to_node[e.get_id()] = e
    return sg, exons


def _transcript(sg, exons, transcript_id, gene_id):
    segments = [list(e.get_coords()) for e in exons]
    t = Transcript(CONTIG, segments, "+")
    t.set_transcript_id(transcript_id)
    t.set_gene_id(gene_id)
    t._simplepath = [e.get_id() for e in exons]
    t.set_splice_graph(sg) if hasattr(t, "set_splice_graph") else None
    return t


def _counter(sg, exons, read_names):
    counter = MultiPathCounter()
    for rn in read_names:
        mp = MultiPath(
            sg,
            [[e.get_id() for e in exons]],
            read_types={"PacBio"},
            read_names={rn},
            read_count=1,
        )
        counter.add(mp)
    return counter


def test_repeated_quantify_does_not_carry_state_forward():
    """Two quantify() calls on one object: the second must see only its own multipaths."""
    sg, exons = _graph()
    t1 = _transcript(sg, exons, "t1", "g1")
    q = Quantify(run_EM=False, max_EM_iterations=10)

    q.quantify(sg, [t1], _counter(sg, exons, ["readA", "readB"]))
    first_mps = set(q.get_mp_to_transcripts().keys())
    assert first_mps, "first pass produced no multipath assignments; fixture is inert"

    q.quantify(sg, [t1], _counter(sg, exons, ["readC"]))
    second_mps = set(q.get_mp_to_transcripts().keys())

    assert not (first_mps & second_mps), (
        "multipaths from the first quantify() survived into the second; a streaming "
        "assignment table built from these would mix passes"
    )
    assert len(second_mps) == 1, (
        f"second pass should hold only its own multipath, holds {len(second_mps)}"
    )


def test_repeated_quantify_clears_unassigned_and_read_map():
    sg, exons = _graph()
    t1 = _transcript(sg, exons, "t1", "g1")
    q = Quantify(run_EM=False, max_EM_iterations=10)

    q.quantify(sg, [t1], _counter(sg, exons, ["readA", "readB"]))
    n_first = len(q.get_unassigned_mp_count_pairs())

    q.quantify(sg, [t1], _counter(sg, exons, ["readC"]))
    n_second = len(q.get_unassigned_mp_count_pairs())

    # whatever the fixture's assignment outcome, the second call must not report the
    # first call's unassigned pairs on top of its own
    assert n_second <= max(n_first, 1), (
        f"unassigned pairs accumulated across calls: {n_first} then {n_second}"
    )
    assert len(q._read_name_to_multipath) <= 1, (
        "read-name -> multipath map retained names from the previous quantify()"
    )


def test_second_quantify_cannot_assign_to_a_dropped_transcript():
    """The sharper case: pass 2's reads must not reach pass 1's isoforms.

    `_assign_path_nodes_to_gene` only ever adds, so a transcript dropped between calls
    stays reachable through `_path_node_id_to_gene_ids` and
    `_gene_id_to_transcript_objs`, and `_get_all_genes_with_node_matches_to_simplepath`
    will still offer its gene as a candidate. That is not a stale leftover in a results
    container -- it is a live assignment to an isoform the current run does not contain,
    which then reaches the tracking file.

    Uses disjoint node sets per pass so any assignment across the two can only come from
    a retained index, never from genuine shared structure.
    """
    sg = Splice_graph()
    Exon.reset_counter()
    sg._contig_acc = CONTIG
    sg._contig_strand = "+"

    left = [Exon(CONTIG, 100, 200, "+", 10.0), Exon(CONTIG, 300, 400, "+", 10.0)]
    right = [Exon(CONTIG, 5000, 5100, "+", 10.0), Exon(CONTIG, 5300, 5400, "+", 10.0)]
    for e in left + right:
        sg._node_id_to_node[e.get_id()] = e

    t_old = _transcript(sg, left, "t_old", "g_old")
    t_new = _transcript(sg, right, "t_new", "g_new")

    q = Quantify(run_EM=False, max_EM_iterations=10)

    # pass 1: only the left locus exists
    q.quantify(sg, [t_old], _counter(sg, left, ["readL"]))
    assert "g_old" in q._gene_id_to_transcript_objs, "pass 1 did not index its gene"

    # pass 2: a disjoint locus, and t_old is gone from the transcript set
    q.quantify(sg, [t_new], _counter(sg, right, ["readR"]))

    assert "g_old" not in q._gene_id_to_transcript_objs, (
        "gene from the previous quantify() is still indexed, so its dropped isoforms "
        "remain candidates for this pass's reads"
    )
    for node_id, gene_ids in q._path_node_id_to_gene_ids.items():
        assert "g_old" not in gene_ids, (
            f"node {node_id} still maps to the previous pass's gene"
        )

    assigned_tids = {
        t.get_transcript_id()
        for transcripts in q.get_mp_to_transcripts().values()
        for t in transcripts
    }
    assert "t_old" not in assigned_tids, (
        "a read from pass 2 was assigned to a transcript pass 2 does not contain"
    )


def test_indexing_refuses_a_transcript_under_two_genes():
    """Caught at indexing time, which is the only place the check is complete.

    _assign_path_nodes_to_gene visits every transcript exactly once, so it detects a
    malformed annotation globally and in O(transcripts). The alternative -- checking inside
    candidate_isoforms_for_genes -- runs once per multipath and can only notice a collision
    on those paths where both offending genes happen to be candidate genes, so a malformed
    input whose two genes never compete for the same read would pass unnoticed.

    It must also fire when the two registrations arrive in separate calls, which is why the
    backing map outlives a single call and is reset per quantify() instead.
    """
    sg, exons = _graph()
    q = Quantify(run_EM=False, max_EM_iterations=10)

    shared_a = _transcript(sg, exons, "shared_tx", "gene_A")
    shared_b = _transcript(sg, exons, "shared_tx", "gene_B")

    with pytest.raises(RuntimeError, match="two gene ids"):
        q._assign_path_nodes_to_gene([shared_a, shared_b])

    # and across separate registrations
    q2 = Quantify(run_EM=False, max_EM_iterations=10)
    q2._assign_path_nodes_to_gene([shared_a])
    with pytest.raises(RuntimeError, match="two gene ids"):
        q2._assign_path_nodes_to_gene([shared_b])


def test_candidate_isoforms_returns_one_object_per_transcript():
    """The ordinary case still collapses overlapping genes to unique transcript ids."""
    sg, exons = _graph()
    q = Quantify(run_EM=False, max_EM_iterations=10)

    t1 = _transcript(sg, exons, "t1", "gene_A")
    t2 = _transcript(sg, exons, "t2", "gene_A")
    t3 = _transcript(sg, exons, "t3", "gene_B")
    q._assign_path_nodes_to_gene([t1, t2, t3])

    got = q.candidate_isoforms_for_genes(["gene_A", "gene_B"])
    tids = [t.get_transcript_id() for t in got]
    assert tids == ["t1", "t2", "t3"], "insertion order must follow the gene sequence"
    assert len(set(tids)) == len(tids)


def test_theta_is_empty_without_EM_and_scoped_with_it():
    """Theta must never describe a transcript the current call did not quantify."""
    sg, exons = _graph()
    t1 = _transcript(sg, exons, "t1", "g1")

    q = Quantify(run_EM=False, max_EM_iterations=10)
    q.quantify(sg, [t1], _counter(sg, exons, ["readA"]))
    assert q.get_transcript_id_to_theta() == {}, (
        "run_EM=False must leave theta empty; a non-empty theta here would be read as "
        "EM abundances by a consumer that cannot tell the difference"
    )


def test_direct_estimator_calls_do_not_publish_theta():
    """The aggregate is reset at quantify() entry, so only quantify() may publish into it.

    _estimate_isoform_read_support has two kinds of caller. quantify() drives it once per EM
    unit and wants them to accumulate; it resets the aggregate first. TranscriptFiltering
    (TranscriptFiltering.py:609) and LRAA's post-filter report (LRAA:5649) drive it DIRECTLY on
    an object that never enters quantify(), in a loop -- so no reset ever runs for them, and an
    unconditional publish would leave theta growing across filtering rounds and describing a
    transcript set no single call ever saw.

    Nothing reads theta off those objects today, so this is a latent defect rather than a live
    one; pinning it is what stops the next consumer from inheriting it. Clearing at the
    estimator's entry instead would be actively wrong -- quantify() calls it per unit, so the
    aggregate would keep only the last unit.
    """
    sg, exons = _graph()
    t1 = _transcript(sg, exons, "t1", "g1")
    counter = _counter(sg, exons, ["readA", "readB"])

    q = Quantify(run_EM=True, max_EM_iterations=10)
    q._assign_path_nodes_to_gene([t1])
    q._assign_reads_to_transcripts(sg, counter)

    # drive the estimator the way the filtering loop does: repeatedly, no quantify()
    for _round in range(3):
        q._estimate_isoform_read_support([t1])

    with pytest.raises(RuntimeError, match="no valid EM theta"):
        q.get_transcript_id_to_theta()

    # and the opt-in path does publish
    q2 = Quantify(run_EM=True, max_EM_iterations=10)
    q2._assign_path_nodes_to_gene([t1])
    q2._assign_reads_to_transcripts(sg, counter)
    q2.quantify(sg, [t1], _counter(sg, exons, ["readA", "readB"]))
    assert q2.get_transcript_id_to_theta(), (
        "a completed quantify() must publish, or streaming gets an empty aggregate"
    )


def test_theta_from_a_completed_run_is_invalidated_by_a_later_direct_call():
    """The sequence the fresh-object test cannot reach: publish, then drive the estimator.

    The accumulate_theta gate stops NEW leakage but does not invalidate what is already there.
    An object that completes quantify() holds a valid aggregate; if it is then driven through
    _estimate_isoform_read_support directly -- the filtering path's entry point -- that
    aggregate describes a transcript set the object has since moved on from, and the accessor
    would keep serving it.

    Invalidating on the non-accumulating branch is safe precisely because quantify() is the
    only caller that passes True, so nothing needing the per-unit accumulation reaches it.
    """
    sg, exons = _graph()
    t1 = _transcript(sg, exons, "t1", "g1")

    q = Quantify(run_EM=True, max_EM_iterations=10)
    q.quantify(sg, [t1], _counter(sg, exons, ["readA", "readB"]))
    published = q.get_transcript_id_to_theta()
    assert published, "fixture produced no theta; the test would prove nothing"

    # now drive the estimator the way TranscriptFiltering does, on the SAME object
    q._estimate_isoform_read_support([t1])

    with pytest.raises(RuntimeError, match="no valid EM theta"):
        q.get_transcript_id_to_theta()


def test_a_failed_second_quantify_invalidates_the_previous_theta():
    """Validity must not be inherited across calls, or a raise mid-run serves a partial theta.

    Sequence the other tests miss: one call succeeds and publishes; a second call begins, empties
    the aggregate, then raises before finishing. With validity only ever SET on success, the flag
    still says True from the first call while the aggregate is now empty -- so the accessor hands
    out an empty theta stamped valid, which a consumer reads as "EM ran and found nothing" rather
    than "this object has no answer".

    The second call is made to raise through a real in-band failure -- two transcripts sharing an
    id under different genes, which quantify() refuses after the reset -- rather than by patching,
    so the test exercises the ordering as shipped.
    """
    sg, exons = _graph()
    t1 = _transcript(sg, exons, "t1", "g1")

    q = Quantify(run_EM=True, max_EM_iterations=10)
    q.quantify(sg, [t1], _counter(sg, exons, ["readA", "readB"]))
    assert q.get_transcript_id_to_theta(), "first call must publish; fixture is inert otherwise"

    collide_a = _transcript(sg, exons, "shared_tx", "gene_A")
    collide_b = _transcript(sg, exons, "shared_tx", "gene_B")
    with pytest.raises(RuntimeError, match="two gene ids"):
        q.quantify(sg, [collide_a, collide_b], _counter(sg, exons, ["readC"]))

    with pytest.raises(RuntimeError, match="no valid EM theta"):
        q.get_transcript_id_to_theta()


def test_a_direct_estimator_call_whose_EM_raises_invalidates_the_previous_theta():
    """The same failed-call hole, one level down from quantify().

    Invalidating after EM.run_EM returns is not enough: a direct estimator call following a
    valid quantify() will, if EM raises, leave the earlier run's theta both present and marked
    valid -- and the accessor serves it. That is the identical defect as inheriting validity
    across quantify() calls, so it needs the identical remedy: invalidate at entry, before any
    work that can fail.

    EM is patched to raise because the point under test is the ORDERING, not EM. Everything
    else about the sequence is real.
    """
    import EM as EM_module

    sg, exons = _graph()
    t1 = _transcript(sg, exons, "t1", "g1")

    q = Quantify(run_EM=True, max_EM_iterations=10)
    q.quantify(sg, [t1], _counter(sg, exons, ["readA", "readB"]))
    assert q.get_transcript_id_to_theta(), "first call must publish; fixture is inert otherwise"

    real_run_EM = EM_module.run_EM

    def exploding_run_EM(*_a, **_k):
        raise RuntimeError("synthetic EM failure")

    EM_module.run_EM = exploding_run_EM
    try:
        with pytest.raises(RuntimeError, match="synthetic EM failure"):
            q._estimate_isoform_read_support([t1])
    finally:
        EM_module.run_EM = real_run_EM

    with pytest.raises(RuntimeError, match="no valid EM theta"):
        q.get_transcript_id_to_theta()


def test_theta_is_invalidated_even_when_argument_validation_rejects_the_call():
    """The guard must precede the asserts, not follow them.

    quantify() validates its arguments by indexing transcripts[0]. With the invalidation placed
    after that, quantify(sg, [], mp_counter) raises before reaching it and the previous call's
    theta stays both present and marked valid -- so the accessor answers with abundances for a
    transcript set the caller has since abandoned.

    Found by reading the other branch's equivalent guard, which was placed above its own
    asserts for exactly this reason.
    """
    sg, exons = _graph()
    t1 = _transcript(sg, exons, "t1", "g1")

    q = Quantify(run_EM=True, max_EM_iterations=10)
    q.quantify(sg, [t1], _counter(sg, exons, ["readA", "readB"]))
    assert q.get_transcript_id_to_theta(), "first call must publish; fixture is inert otherwise"

    # rejected by argument validation, before any quantification work
    with pytest.raises((IndexError, AssertionError)):
        q.quantify(sg, [], _counter(sg, exons, ["readC"]))

    with pytest.raises(RuntimeError, match="no valid EM theta"):
        q.get_transcript_id_to_theta()
