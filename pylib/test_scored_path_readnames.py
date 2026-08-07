#!/usr/bin/env python3

import os
import pytest

import LRAA_Globals as LG
from ReadNameStore import ReadNameStore

from Splice_graph import Splice_graph
from GenomeFeature import Exon
from MultiPath import MultiPath
from MultiPathGraphNode import MultiPathGraphNode
from MultiPathCounter import MultiPathCounter
from Scored_path import Scored_path


def build_minimal_sg_with_exons():
    sg = Splice_graph()
    # Reset exon counter to deterministic IDs if needed
    Exon.reset_counter()
    # Create three simple exons in order
    e1 = Exon("chr1", 100, 150, "+", 1)
    e2 = Exon("chr1", 200, 250, "+", 1)
    e3 = Exon("chr1", 300, 350, "+", 1)
    # Register them in the splice graph so MultiPath can resolve IDs
    sg._node_id_to_node[e1.get_id()] = e1
    sg._node_id_to_node[e2.get_id()] = e2
    sg._node_id_to_node[e3.get_id()] = e3
    return sg, e1, e2, e3


def test_scored_path_reads_without_containments(tmp_path):
    sg, e1, e2, _ = build_minimal_sg_with_exons()
    # initialize name store so string names are mapped to IDs
    base = tmp_path / "spnames1"
    os.makedirs(base, exist_ok=True)
    LG.READ_NAME_STORE = ReadNameStore(str(base / "names"))

    class DummyMPG:
        def __init__(self, sg):
            self._sg = sg

        def get_splice_graph(self):
            return self._sg

    # path over two exons; embed read names directly to avoid external store dependency
    mp = MultiPath(
        sg, [[e1.get_id(), e2.get_id()]], read_names={"r1", "r2"}, read_count=2
    )
    lend, rend = mp.get_coords()
    mpgn = MultiPathGraphNode(mp, count=2, lend=lend, rend=rend, mpg=DummyMPG(sg))

    sp = Scored_path([mpgn])
    names = sp.get_all_represented_reads()

    assert names == {"r1", "r2"}

    # cleanup
    if LG.READ_NAME_STORE is not None:
        LG.READ_NAME_STORE.close()
        LG.READ_NAME_STORE = None


def test_scored_path_reads_with_containments(tmp_path):
    sg, e1, e2, e3 = build_minimal_sg_with_exons()
    base = tmp_path / "spnames2"
    os.makedirs(base, exist_ok=True)
    LG.READ_NAME_STORE = ReadNameStore(str(base / "names"))

    class DummyMPG:
        def __init__(self, sg):
            self._sg = sg

        def get_splice_graph(self):
            return self._sg

    # parent path e1-e3, contained path e1-e2
    mp_parent = MultiPath(
        sg, [[e1.get_id(), e3.get_id()]], read_names={"p1"}, read_count=1
    )
    mp_child = MultiPath(
        sg, [[e1.get_id(), e2.get_id()]], read_names={"c1", "c2"}, read_count=2
    )

    lend_p, rend_p = mp_parent.get_coords()
    lend_c, rend_c = mp_child.get_coords()

    node_parent = MultiPathGraphNode(
        mp_parent, count=1, lend=lend_p, rend=rend_p, mpg=DummyMPG(sg)
    )
    node_child = MultiPathGraphNode(
        mp_child, count=2, lend=lend_c, rend=rend_c, mpg=DummyMPG(sg)
    )

    # manually register containment
    node_parent.add_containment(node_child)

    sp = Scored_path([node_parent])
    names = sp.get_all_represented_reads()

    # union of parent and contained names
    assert names == {"p1", "c1", "c2"}

    if LG.READ_NAME_STORE is not None:
        LG.READ_NAME_STORE.close()
        LG.READ_NAME_STORE = None


class _DummyMPG:
    def __init__(self, sg):
        self._sg = sg

    def get_splice_graph(self):
        return self._sg


def _scored_path(sg, exon_ids, read_names, count):
    mp = MultiPath(sg, [exon_ids], read_names=read_names, read_count=count)
    lend, rend = mp.get_coords()
    mpgn = MultiPathGraphNode(mp, count=count, lend=lend, rend=rend, mpg=_DummyMPG(sg))
    return mp, Scored_path([mpgn])


@pytest.fixture
def name_store(tmp_path):
    base = tmp_path / "synscore"
    os.makedirs(base, exist_ok=True)
    LG.READ_NAME_STORE = ReadNameStore(str(base / "names"))
    preserved = set(LG.SYNTHETIC_READ_IDS)
    LG.SYNTHETIC_READ_IDS.clear()
    yield
    LG.SYNTHETIC_READ_IDS.clear()
    LG.SYNTHETIC_READ_IDS.update(preserved)
    if LG.READ_NAME_STORE is not None:
        LG.READ_NAME_STORE.close()
        LG.READ_NAME_STORE = None


def test_path_supported_only_by_a_synthetic_read_scores_zero(name_store):
    """A structure no read supports must not be selectable as a candidate."""
    sg, e1, e2, _ = build_minimal_sg_with_exons()
    mp, sp = _scored_path(sg, [e1.get_id(), e2.get_id()], {"reftranscript:t1"}, 1)
    assert sp.compute_path_score() == 1  # before registration, it counts

    LG.SYNTHETIC_READ_IDS.update(mp.get_read_ids())
    assert sp.compute_path_score() == 0
    # must not fall through to the node-count fallback, which still holds the read
    assert sp._all_represented_read_ids


def test_synthetic_read_does_not_inflate_a_real_supported_path(name_store):
    sg, e1, e2, _ = build_minimal_sg_with_exons()
    mp, sp = _scored_path(
        sg, [e1.get_id(), e2.get_id()], {"reftranscript:t1", "r1", "r2"}, 3
    )
    synthetic = {
        rid
        for rid in mp.get_read_ids()
        if str(LG.READ_NAME_STORE.get_name(int(rid))).startswith("reftranscript:")
    }
    assert len(synthetic) == 1
    LG.SYNTHETIC_READ_IDS.update(synthetic)
    # scores the two real reads only
    assert sp.compute_path_score() == 2


def test_fallback_still_applies_when_no_read_ids_exist(name_store):
    """The node-count fallback covers missing id data, not zero support."""
    sg, e1, e2, _ = build_minimal_sg_with_exons()
    _, sp = _scored_path(sg, [e1.get_id(), e2.get_id()], set(), 5)
    assert not sp._all_represented_read_ids
    assert sp.compute_path_score() == 5


def _mpgn(sg, exon_ids, read_names, count):
    mp = MultiPath(sg, [exon_ids], read_names=read_names, read_count=count)
    lend, rend = mp.get_coords()
    return mp, MultiPathGraphNode(
        mp, count=count, lend=lend, rend=rend, mpg=_DummyMPG(sg)
    )


def test_unrepresented_nodes_separate_templates_from_lost_reads(name_store):
    from LRAA import LRAA

    sg, e1, e2, e3 = build_minimal_sg_with_exons()
    # a node whose only support is the synthetic read standing in for an annotation
    tmpl_mp, tmpl = _mpgn(sg, [e1.get_id(), e2.get_id()], {"reftranscript:t1"}, 1)
    LG.SYNTHETIC_READ_IDS.update(tmpl_mp.get_read_ids())
    # a node carrying real reads that no reported transcript covered
    _, real = _mpgn(sg, [e2.get_id(), e3.get_id()], {"r1", "r2", "r3"}, 3)

    template_only, nodes_with_reads, real_reads = LRAA._summarize_unrepresented_mpgns(
        [tmpl, real]
    )
    assert (template_only, nodes_with_reads, real_reads) == (1, 1, 3)

    # a template node that also accumulated real reads counts as a loss, not a template
    mixed_mp, mixed = _mpgn(
        sg, [e1.get_id(), e3.get_id()], {"reftranscript:t2", "r9"}, 2
    )
    synthetic = {
        rid
        for rid in mixed_mp.get_read_ids()
        if str(LG.READ_NAME_STORE.get_name(int(rid))).startswith("reftranscript:")
    }
    LG.SYNTHETIC_READ_IDS.update(synthetic)
    assert LRAA._summarize_unrepresented_mpgns([mixed]) == (0, 1, 1)


def test_merge_mode_fake_reads_are_not_discounted_from_scoring(tmp_path):
    """Discovery templates are discounted from path scoring; merge-mode reads are not.

    In discovery a `reftranscript:` read stands in for reads that may or may not exist,
    so a path supported by nothing else must score zero. Merge mode has no real reads
    at all -- its fake reads carry each input model's weight and are the only evidence
    present -- so discounting them there scores every path zero and merges nothing.
    """
    from LRAA import LRAA
    from Transcript import Transcript

    sg, e1, e2, _e3 = build_minimal_sg_with_exons()
    sg._contig_acc = "chr1"
    sg._contig_strand = "+"

    base = tmp_path / "mergenames"
    os.makedirs(base, exist_ok=True)
    LG.READ_NAME_STORE = ReadNameStore(str(base / "names"))
    preserved = set(LG.SYNTHETIC_READ_IDS)
    LG.SYNTHETIC_READ_IDS.clear()
    try:
        transcript = Transcript("chr1", [[100, 150], [200, 250]], "+")
        transcript.set_gene_id("g1")
        transcript.set_transcript_id("t1")
        transcript.set_simple_path([e1.get_id(), e2.get_id()])

        lraa = LRAA(sg)

        # merge mode: no BAM, so the fake reads are the only evidence there is
        merge_counter = MultiPathCounter()
        lraa._incorporate_transcripts_into_mp_counter(
            merge_counter, [transcript], bam_file=None
        )
        assert LG.SYNTHETIC_READ_IDS == set()

        # discovery: the template is registered and therefore discounted
        discovery_counter = MultiPathCounter()
        lraa._incorporate_transcripts_into_mp_counter(
            discovery_counter, [transcript], bam_file="placeholder.bam"
        )
        assert LG.SYNTHETIC_READ_IDS != set()
    finally:
        LG.SYNTHETIC_READ_IDS.clear()
        LG.SYNTHETIC_READ_IDS.update(preserved)
        LG.READ_NAME_STORE = None


def test_template_breaks_a_score_tie_but_never_creates_support():
    from LRAA import _path_selection_sort_key

    class _Stub:
        def __init__(self, score, ids):
            self._score = score
            self._ids = set(ids)

        def get_score(self):
            return self._score

        def get_all_represented_read_ids(self):
            return self._ids

    preserved = set(LG.SYNTHETIC_READ_IDS)
    LG.SYNTHETIC_READ_IDS.clear()
    LG.SYNTHETIC_READ_IDS.update({901})
    try:
        plain = _Stub(5, {1, 2, 3})
        templated = _Stub(5, {1, 2, 3, 901})
        # paths are popped from the end, so the templated one is taken first on a tie
        assert sorted([plain, templated], key=_path_selection_sort_key)[-1] is templated

        # a template cannot outrank real support
        stronger = _Stub(6, {1, 2, 3, 4})
        assert (
            sorted([templated, stronger], key=_path_selection_sort_key)[-1] is stronger
        )
    finally:
        LG.SYNTHETIC_READ_IDS.clear()
        LG.SYNTHETIC_READ_IDS.update(preserved)
