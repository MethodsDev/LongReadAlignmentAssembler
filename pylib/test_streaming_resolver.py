#!/usr/bin/env python3

"""The resolver must reproduce the table, not merely resemble it.

A streaming pass looks each read's splice-graph path up in a table precomputed from a
finished quantification, and resolves the paths that quantification never saw. Those two
routes must produce identical rows for identical input, because the output file gives no
indication of which route produced a row: a read whose path survived coverage normalization
and one whose path did not would be assigned by different code, and any disagreement would
read as biology.

The check here is the one that can actually fail: take a path the table already holds, run
the resolver on it, and require the rows to match exactly. Both routes share
`rows_for_multipath`, so this is really testing that the resolver feeds it the same
anchoring, the same candidate set, the same cascade result and the same weights that
`build()` did.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import LRAA_Globals
import StreamingQuant
from GenomeFeature import Exon
from MultiPath import MultiPath
from MultiPathCounter import MultiPathCounter
from Quantify import Quantify
from Splice_graph import Splice_graph
from StreamingQuant import AssignmentTable, make_path_resolver
from Transcript import Transcript

CONTIG = "chr_res"


def _graph_with_two_isoforms():
    sg = Splice_graph()
    Exon.reset_counter()
    sg._contig_acc = CONTIG
    sg._contig_strand = "+"
    e1 = Exon(CONTIG, 100, 200, "+", 10.0)
    e2 = Exon(CONTIG, 300, 400, "+", 10.0)
    e3 = Exon(CONTIG, 500, 600, "+", 10.0)
    for e in (e1, e2, e3):
        sg._node_id_to_node[e.get_id()] = e
    return sg, (e1, e2, e3)


def _transcript(sg, exons, transcript_id, gene_id):
    t = Transcript(CONTIG, [list(e.get_coords()) for e in exons], "+")
    t.set_transcript_id(transcript_id)
    t.set_gene_id(gene_id)
    t.set_simple_path([e.get_id() for e in exons])
    return t


def _counter(sg, path, read_names):
    counter = MultiPathCounter()
    for rn in read_names:
        counter.add(
            MultiPath(sg, [path], read_types={"PacBio"}, read_names={rn}, read_count=1)
        )
    return counter


def test_resolver_reproduces_the_table_for_a_path_the_table_holds():
    sg, (e1, e2, e3) = _graph_with_two_isoforms()
    long_tx = _transcript(sg, (e1, e2, e3), "long", "g1")
    short_tx = _transcript(sg, (e1, e2), "short", "g1")
    path = [e1.get_id(), e2.get_id()]

    q = Quantify(run_EM=False, max_EM_iterations=10)
    q.quantify(sg, [long_tx, short_tx], _counter(sg, path, ["r1", "r2"]))

    mp_to_transcripts = q.get_mp_to_transcripts()
    if not mp_to_transcripts:
        pytest.skip("fixture produced no assignments; nothing to compare")

    # The real map from the quantification, not a stand-in: the table and the resolver must
    # group by the same components or their splits differ for reasons this test exists to
    # rule out.
    components = q.get_transcript_id_to_component_id()

    table = AssignmentTable.build(
        sg,
        mp_to_transcripts,
        q.get_transcript_id_to_theta(),
        q.em_was_run(),
        q.get_unassigned_mp_count_pairs(),
        components,
    )

    resolver = make_path_resolver(
        sg, q, q.get_transcript_id_to_theta(), q.em_was_run(), components
    )

    for mp in mp_to_transcripts:
        simple_path = mp.get_simple_path()
        canon = sg.canonical_simple_path(simple_path)
        from_table = table.lookup(canon)
        from_resolver = resolver(simple_path)

        assert from_table is not None, "path missing from the table it was built from"
        # mp id is the one field that legitimately differs: the table carries the first
        # pass's multipath id, the resolver mints its own for a synthetic single-read
        # multipath. Everything that affects a count must match.
        strip = lambda rows: [(r[0], r[1], r[2], r[3], r[5], r[6]) for r in rows]
        assert strip(from_resolver) == strip(from_table), (
            "resolver and table disagree for path {}:\n  table    {}\n  resolver {}".format(
                canon, strip(from_table), strip(from_resolver)
            )
        )


def test_resolver_returns_empty_for_a_path_anchoring_to_no_gene():
    """An unanchorable path yields [], which is a cached answer rather than a dropped read."""
    sg, (e1, e2, e3) = _graph_with_two_isoforms()
    tx = _transcript(sg, (e1, e2), "short", "g1")

    q = Quantify(run_EM=False, max_EM_iterations=10)
    q._assign_path_nodes_to_gene([tx])

    resolver = make_path_resolver(sg, q, {}, False, {})
    orphan = Exon(CONTIG, 9000, 9100, "+", 5.0)
    sg._node_id_to_node[orphan.get_id()] = orphan

    assert resolver([orphan.get_id()]) == []


def test_resolver_refuses_when_em_ran_without_theta():
    """Same refusal as the table: an empty theta must not silently become an equal split."""
    sg, (e1, e2, e3) = _graph_with_two_isoforms()
    tx = _transcript(sg, (e1, e2), "short", "g1")

    q = Quantify(run_EM=True, max_EM_iterations=10)
    q._assign_path_nodes_to_gene([tx])

    # A real component map, so the refusal under test is the theta one rather than the
    # component-coverage check that would otherwise fire first.
    resolver = make_path_resolver(sg, q, {}, True, {"short": "comp0"})
    with pytest.raises(RuntimeError, match="EM ran but no theta"):
        resolver([e1.get_id(), e2.get_id()])
