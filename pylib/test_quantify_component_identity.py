#!/usr/bin/env python3

"""The component identity Quantify exposes alongside theta.

theta is normalized over a read-sharing component of genes, not over a gene, so a
consumer that recomputes a per-read isoform split has to renormalize over the
component.  Renormalizing over a gene gives a read compatible with transcripts of
two genes in one component a split summing to 2, which is arithmetically plausible
and produces inflated counts silently.  These tests pin the mapping that lets a
consumer avoid that, and pin that it never carries an entry between calls.
"""

import logging
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


EXON_COORDS = [(100, 150), (200, 250), (300, 350)]


def _build(read_exon_idxs, transcript_specs):
    """One splice graph, transcripts named by (gene_id, transcript_id), one read."""

    Exon.reset_counter()
    splice_graph = Splice_graph()
    exons = [Exon("chr1", lend, rend, "+", 1) for lend, rend in EXON_COORDS]
    for exon in exons:
        splice_graph._node_id_to_node[exon.get_id()] = exon
    splice_graph._contig_acc = "chr1"
    splice_graph._contig_strand = "+"

    node_ids = [exon.get_id() for exon in exons]

    transcripts = list()
    for gene_id, transcript_id, exon_idxs in transcript_specs:
        transcript = MultiPath(
            splice_graph, [[node_ids[i] for i in exon_idxs]]
        ).toTranscript()
        transcript.set_gene_id(gene_id)
        transcript.set_transcript_id(transcript_id)
        transcripts.append(transcript)

    mp_counter = MultiPathCounter()
    mp_counter.add(
        MultiPath(
            splice_graph,
            [[node_ids[i] for i in read_exon_idxs]],
            read_types={"PacBio"},
            read_names={"r1"},
            read_count=1,
        )
    )

    return splice_graph, transcripts, mp_counter


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
    logging.disable(logging.CRITICAL)
    try:
        yield
    finally:
        logging.disable(logging.NOTSET)
        LRAA_Globals.config.update(saved)


# Read spans all three exons, so it is compatible with a transcript of each gene.
_SHARED_READ = ([0, 1, 2], [("g1", "t1", [0, 1]), ("g2", "t2", [1, 2])])

# Read covers the first exon only, so only g1's transcript takes it.
_UNSHARED_READ = ([0], [("g1", "t1", [0, 1]), ("g2", "t2", [1, 2])])


def test_two_genes_sharing_a_read_are_one_component():
    """The case a per-gene denominator gets wrong: one read, two genes, one EM."""

    quantify = Quantify(False, 1)
    quantify.quantify(*_build(*_SHARED_READ))

    assigned = list(quantify.get_mp_to_transcripts().values())
    assert len(assigned) == 1
    assert sorted(t.get_transcript_id() for t in assigned[0]) == ["t1", "t2"]

    tx_to_component = quantify.get_transcript_id_to_component_id()
    assert tx_to_component["t1"] == tx_to_component["t2"], (
        "transcripts competing for one read are quantified by one joint EM, so a "
        "consumer must renormalize over both and not over each gene"
    )

    component_id = tx_to_component["t1"]
    assert quantify.get_component_id_to_gene_ids() == {component_id: ["g1", "g2"]}


def test_genes_sharing_no_read_stay_separate_components():
    quantify = Quantify(False, 1)
    quantify.quantify(*_build(*_UNSHARED_READ))

    tx_to_component = quantify.get_transcript_id_to_component_id()
    assert tx_to_component["t1"] != tx_to_component["t2"]
    assert quantify.get_component_id_to_gene_ids() == {"g1": ["g1"], "g2": ["g2"]}


def test_a_transcript_taking_no_read_still_has_a_component():
    """Completeness: a consumer never needs a fallback for the singleton case."""

    quantify = Quantify(False, 1)
    splice_graph, transcripts, mp_counter = _build(*_UNSHARED_READ)
    quantify.quantify(splice_graph, transcripts, mp_counter)

    assigned_transcript_ids = {
        t.get_transcript_id()
        for transcripts_assigned in quantify.get_mp_to_transcripts().values()
        for t in transcripts_assigned
    }
    assert "t2" not in assigned_transcript_ids, "t2 takes no read in this fixture"

    tx_to_component = quantify.get_transcript_id_to_component_id()
    assert set(tx_to_component) == {t.get_transcript_id() for t in transcripts}


def test_component_identity_does_not_survive_into_a_later_call():
    """A stale entry pairs this run's theta with a previous run's component.

    quantify() is called repeatedly on one object over a shrinking transcript set,
    so a map that accumulated would keep answering for transcripts the current run
    does not contain.
    """

    quantify = Quantify(False, 1)
    quantify.quantify(*_build(*_SHARED_READ))
    assert set(quantify.get_transcript_id_to_component_id()) == {"t1", "t2"}

    quantify.quantify(*_build([2], [("g9", "t9", [2])]))

    assert set(quantify.get_transcript_id_to_component_id()) == {"t9"}, (
        "transcripts absent from this call must not retain a component id"
    )
    assert quantify.get_component_id_to_gene_ids() == {"g9": ["g9"]}


def test_a_call_that_fails_before_component_construction_refuses_rather_than_answers():
    """A failed call must not serve the previous call's components.

    Invalidated at entry, before any work that can fail, so a raise anywhere in
    quantify() leaves the accessors refusing.  Invalidating where the maps are
    filled would expose the earlier call's components in the quietest possible
    form: a complete, plausible mapping belonging to a transcript set the object
    no longer holds, with no error to notice.
    """

    quantify = Quantify(False, 1)
    quantify.quantify(*_build(*_SHARED_READ))
    assert set(quantify.get_transcript_id_to_component_id()) == {"t1", "t2"}

    def _fail(*args, **kwargs):
        raise RuntimeError("read assignment failed")

    quantify._assign_reads_to_transcripts = _fail

    with pytest.raises(RuntimeError, match="read assignment failed"):
        quantify.quantify(*_build([2], [("g9", "t9", [2])]))

    with pytest.raises(RuntimeError, match="quantify\\(\\) has not completed"):
        quantify.get_transcript_id_to_component_id()
    with pytest.raises(RuntimeError, match="quantify\\(\\) has not completed"):
        quantify.get_component_id_to_gene_ids()


def test_the_accessors_refuse_before_the_first_call():
    """An empty mapping reads as "no components"; a refusal cannot be misread.

    A consumer treating empty as a real answer would renormalize per gene and
    inflate every read compatible with two genes in one component.
    """

    quantify = Quantify(False, 1)
    with pytest.raises(RuntimeError, match="quantify\\(\\) has not completed"):
        quantify.get_transcript_id_to_component_id()
    with pytest.raises(RuntimeError, match="quantify\\(\\) has not completed"):
        quantify.get_component_id_to_gene_ids()


def test_a_call_rejected_by_argument_validation_still_invalidates():
    """The narrowest window: a raise ABOVE the body, in the asserts themselves.

    quantify() indexes transcripts[0], so an empty list raises before any work
    runs at all.  An invalidation placed after argument validation -- the obvious
    place, since nothing has happened yet -- leaves the previous call's components
    readable and marked valid, which is the same stale answer as any other failure
    path and the easiest one to overlook.
    """

    quantify = Quantify(False, 1)
    splice_graph, transcripts, mp_counter = _build(*_SHARED_READ)
    quantify.quantify(splice_graph, transcripts, mp_counter)
    assert set(quantify.get_transcript_id_to_component_id()) == {"t1", "t2"}

    with pytest.raises((AssertionError, IndexError)):
        quantify.quantify(splice_graph, [], mp_counter)

    with pytest.raises(RuntimeError, match="quantify\\(\\) has not completed"):
        quantify.get_transcript_id_to_component_id()
    with pytest.raises(RuntimeError, match="quantify\\(\\) has not completed"):
        quantify.get_component_id_to_gene_ids()


def test_a_direct_estimator_call_invalidates_a_published_map():
    """The non-publishing entry point, which is also a real one.

    run_quant_only completes quantify() and returns its Quantify (LRAA:4433,
    :4469); a later caller drives _estimate_isoform_read_support directly
    (LRAA:5079) over a transcript set that may have been filtered since.  That
    updates theta and the assignments while the component map still describes the
    set quantify() saw, so a consumer reading both gets this run's theta against
    an earlier grouping.

    Publishing the map only from quantify() does not prevent this: an
    already-valid map stays serveable until something invalidates it.
    """

    quantify = Quantify(False, 1)
    splice_graph, transcripts, mp_counter = _build(*_SHARED_READ)
    quantify.quantify(splice_graph, transcripts, mp_counter)
    assert set(quantify.get_transcript_id_to_component_id()) == {"t1", "t2"}

    # exactly what LRAA:5079 does with the object run_quant_only handed back
    quantify._estimate_isoform_read_support(transcripts)

    with pytest.raises(RuntimeError, match="quantify\\(\\) has not completed"):
        quantify.get_transcript_id_to_component_id()
    with pytest.raises(RuntimeError, match="quantify\\(\\) has not completed"):
        quantify.get_component_id_to_gene_ids()


def test_a_direct_read_assignment_call_invalidates_a_published_map():
    """The third mutating entry point: it rewrites what components derive FROM.

    _assign_reads_to_transcripts writes _mp_to_transcripts, which is the input
    _build_read_sharing_gene_components unions over.  A map built before that call
    no longer describes the object's assignment state, so it must stop answering.

    Every caller today reaches this on a fresh object, so nothing depends on the
    guard yet.  It is held anyway: "no caller reaches it with a valid map" is a
    claim about callers, and the equivalent claim has already been wrong twice.
    """

    quantify = Quantify(False, 1)
    splice_graph, transcripts, mp_counter = _build(*_SHARED_READ)
    quantify.quantify(splice_graph, transcripts, mp_counter)
    assert set(quantify.get_transcript_id_to_component_id()) == {"t1", "t2"}

    _, _, other_counter = _build(*_UNSHARED_READ)
    quantify._assign_reads_to_transcripts(splice_graph, other_counter)

    with pytest.raises(RuntimeError, match="quantify\\(\\) has not completed"):
        quantify.get_transcript_id_to_component_id()
    with pytest.raises(RuntimeError, match="quantify\\(\\) has not completed"):
        quantify.get_component_id_to_gene_ids()


def test_the_accessors_hand_back_copies():
    quantify = Quantify(False, 1)
    quantify.quantify(*_build(*_SHARED_READ))

    quantify.get_transcript_id_to_component_id()["t1"] = "mutated"
    quantify.get_component_id_to_gene_ids()["g1"].append("mutated")

    assert quantify.get_transcript_id_to_component_id()["t1"] == "g1"
    assert quantify.get_component_id_to_gene_ids()["g1"] == ["g1", "g2"]


def test_a_bridge_from_an_earlier_call_does_not_fuse_two_surviving_genes():
    """The stale entry that changes an answer rather than merely lingering.

    Component construction unions the genes of every retained _mp_to_transcripts
    entry, and its guard only skips genes absent from the current call.  Two genes
    that both keep isoforms are therefore still eligible to be fused by a multipath
    from a previous call, which is reachable: de novo runs quantify() three times on
    one object, after degradation pruning and after isoform-fraction filtering, and
    a bridging transcript removed by either one leaves its multipath behind.

    The fusion only ever adds genes to a component, so theta is normalized across
    genes that no longer share a read and every fraction in both is depressed.
    Nothing in the output marks it: the component is well formed, just not this
    call's.
    """

    quantify = Quantify(False, 1)

    # call 1 genuinely bridges g1 and g2
    quantify.quantify(*_build(*_SHARED_READ))
    fused = quantify.get_transcript_id_to_component_id()
    assert fused["t1"] == fused["t2"]

    # call 2 keeps both genes and both isoforms; only the bridging read is gone
    quantify.quantify(*_build(*_UNSHARED_READ))

    tx_to_component = quantify.get_transcript_id_to_component_id()
    assert tx_to_component["t1"] != tx_to_component["t2"], (
        "a multipath from an earlier call must not fuse genes this call's reads "
        "leave separate"
    )
    assert quantify.get_component_id_to_gene_ids() == {"g1": ["g1"], "g2": ["g2"]}
