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
    """A structure no read supports must never be reconstructed.

    The template alone is not evidence, so such a path scores zero and is not
    selected. Asserted on the SIRV suite too, via test_unsupported_ref_model.
    """
    sg, e1, e2, _ = build_minimal_sg_with_exons()
    mp, sp = _scored_path(sg, [e1.get_id(), e2.get_id()], {"reftranscript:t1"}, 1)
    assert sp.compute_path_score() == 1  # before registration, it counts

    LG.SYNTHETIC_READ_IDS.update(mp.get_read_ids())
    assert sp.compute_path_score() == 0
    # must not fall through to the node-count fallback, which still holds the read
    assert sp._all_represented_read_ids


def test_synthetic_read_floors_a_path_whose_real_reads_were_claimed(name_store):
    """A supplied model with real evidence stays selectable once that evidence is taken.

    This is the case the greedy loop otherwise loses: the reads are real but a dominant
    overlapping isoform got them first, so without the template floor the path rescores
    to zero and never reaches quantification.
    """
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

    # unexcluded: the template counts alongside the two real reads
    assert sp.compute_path_score() == 3

    # every real read claimed by a competitor: the floor keeps it selectable
    real = set(mp.get_read_ids()) - synthetic
    assert sp.compute_path_score(exclude_read_ids=real) == 1


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


def test_detach_leaves_a_rescued_read_on_only_its_rescued_path(tmp_path):
    """A rescued read must not keep supporting the path the rescue improved on.

    Otherwise one read backs two competing structures, and a path whose reads have all
    moved elsewhere keeps their support.
    """
    from MultiPathCounter import MultiPathCounter

    sg, e1, e2, e3 = build_minimal_sg_with_exons()
    base = tmp_path / "detachnames"
    os.makedirs(base, exist_ok=True)
    LG.READ_NAME_STORE = ReadNameStore(str(base / "names"))
    try:
        original = MultiPath(
            sg, [[e1.get_id(), e2.get_id()]], read_names={"shared", "own"}
        )
        rescued = MultiPath(
            sg, [[e1.get_id(), e2.get_id(), e3.get_id()]], read_names={"shared"}
        )
        counter = MultiPathCounter()
        counter.add(original)
        counter.add(rescued)

        detached = counter.detach_reads_from_other_paths(
            rescued.get_read_ids(), {rescued.get_simple_path_tuple()}
        )
        assert detached == 1

        by_path = {
            mp.get_simple_path_tuple(): count
            for mp, count in (
                pair.get_multipath_and_count()
                for pair in counter.get_all_MultiPathCountPairs()
            )
        }
        # the rescued path keeps the read; the original keeps only its own
        assert by_path[rescued.get_simple_path_tuple()] == 1
        assert by_path[original.get_simple_path_tuple()] == 1
        assert "shared" not in {
            LG.READ_NAME_STORE.get_name(rid) for rid in original.get_read_ids()
        }
    finally:
        LG.READ_NAME_STORE = None


def _paths_to_counts(counter):
    return {
        mp.get_simple_path_tuple(): count
        for mp, count in (
            pair.get_multipath_and_count()
            for pair in counter.get_all_MultiPathCountPairs()
        )
    }


def _names_on(multipath):
    return {LG.READ_NAME_STORE.get_name(rid) for rid in multipath.get_read_ids()}


def test_detach_keeps_the_unrescued_read_on_the_original_path(tmp_path):
    """A read that was not rescued keeps its genome path attachment.

    Detachment is justified only by an accepted rescue having explained at least as
    much of that read. It says nothing about the other reads sharing the path, so
    they must survive by name, not merely leave the count non-zero.
    """
    sg, e1, e2, e3 = build_minimal_sg_with_exons()
    base = tmp_path / "detach_unrescued"
    os.makedirs(base, exist_ok=True)
    LG.READ_NAME_STORE = ReadNameStore(str(base / "names"))
    try:
        original = MultiPath(
            sg, [[e1.get_id(), e2.get_id()]], read_names={"moved", "stayed"}
        )
        rescued = MultiPath(
            sg, [[e1.get_id(), e2.get_id(), e3.get_id()]], read_names={"moved"}
        )
        counter = MultiPathCounter()
        counter.add(original)
        counter.add(rescued)

        counter.detach_reads_from_other_paths(
            rescued.get_read_ids(), {rescued.get_simple_path_tuple()}
        )

        assert _names_on(original) == {"stayed"}
        assert original.get_simple_path_tuple() in _paths_to_counts(counter)
    finally:
        LG.READ_NAME_STORE = None


def test_detach_leaves_a_wholly_unrescued_path_untouched(tmp_path):
    """A path sharing no reads with the rescue is not touched at all."""
    sg, e1, e2, e3 = build_minimal_sg_with_exons()
    base = tmp_path / "detach_bystander"
    os.makedirs(base, exist_ok=True)
    LG.READ_NAME_STORE = ReadNameStore(str(base / "names"))
    try:
        bystander = MultiPath(sg, [[e1.get_id(), e2.get_id()]], read_names={"solo"})
        rescued = MultiPath(
            sg, [[e1.get_id(), e2.get_id(), e3.get_id()]], read_names={"moved"}
        )
        counter = MultiPathCounter()
        counter.add(bystander)
        counter.add(rescued)

        detached = counter.detach_reads_from_other_paths(
            rescued.get_read_ids(), {rescued.get_simple_path_tuple()}
        )

        assert detached == 0
        assert _paths_to_counts(counter)[bystander.get_simple_path_tuple()] == 1
        assert _names_on(bystander) == {"solo"}
    finally:
        LG.READ_NAME_STORE = None


def test_no_accepted_rescue_detaches_nothing(tmp_path):
    """Offering a read to rescue must not cost it its genome path.

    Rescue considers reads that already have a graph path, and an alignment can be
    declined by the aligned-length floor or the percent-identity floor. A declined
    read is not a rescued read, so its single-read path must survive: it is exactly
    the case where dropping the attachment would delete the model outright.
    """
    from LRAA import _detach_rescued_reads_from_original_paths

    sg, e1, e2, _ = build_minimal_sg_with_exons()
    base = tmp_path / "detach_declined"
    os.makedirs(base, exist_ok=True)
    LG.READ_NAME_STORE = ReadNameStore(str(base / "names"))
    try:
        original = MultiPath(sg, [[e1.get_id(), e2.get_id()]], read_names={"candidate"})
        counter = MultiPathCounter()
        counter.add(original)

        # rescue ran and accepted nothing
        detached = _detach_rescued_reads_from_original_paths(counter, [])

        assert detached == 0
        assert _paths_to_counts(counter)[original.get_simple_path_tuple()] == 1
        assert _names_on(original) == {"candidate"}
    finally:
        LG.READ_NAME_STORE = None
