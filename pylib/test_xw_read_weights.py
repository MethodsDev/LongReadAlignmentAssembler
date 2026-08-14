#!/usr/bin/env python3

"""Weighted support must track the reads a multipath stands for, not the reads kept.

Coverage normalization accepts a read with probability p and records 1/p in XW, so a
multipath's support is the sum of its reads' weights. The invariants that matter:

  - with nothing registered, weight is indistinguishable from the read count, so an
    unnormalized bam quantifies exactly as before;
  - aggregating the same read twice must not double it, in weight or in count;
  - EM is scale-free in these quantities, so a locus thinned uniformly recovers the
    same isoform fractions it would have had unthinned.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import LRAA_Globals
from MultiPath import MultiPath
from MultiPathCounter import MultiPathCounter
from Splice_graph import Splice_graph
from GenomeFeature import Exon
from EM import em_algorithm_with_weights


CONTIG = "chr_test"


@pytest.fixture(autouse=True)
def clean_registry():
    LRAA_Globals.reset_read_weight_registry()
    yield
    LRAA_Globals.reset_read_weight_registry()


def _splice_graph():
    sg = Splice_graph()
    Exon.reset_counter()
    sg._contig_acc = CONTIG
    sg._contig_strand = "+"
    exons = [
        Exon(CONTIG, 100, 200, "+", 1.0),
        Exon(CONTIG, 300, 400, "+", 1.0),
        Exon(CONTIG, 500, 600, "+", 1.0),
    ]
    for e in exons:
        e._id = e.get_id()
        sg._node_id_to_node[e.get_id()] = e
    return sg, exons


def _mp(sg, exons, read_name):
    return MultiPath(
        sg,
        [[e.get_id() for e in exons]],
        read_types={"PacBio"},
        read_names={read_name},
        read_count=1,
    )


def _rid(read_name):
    return MultiPath._coerce_read_identifier(None, read_name)


def test_unregistered_reads_weigh_the_same_as_they_count():
    """An unnormalized bam must be quantified exactly as it was before weighting."""
    sg, exons = _splice_graph()
    counter = MultiPathCounter()
    for i in range(5):
        counter.add(_mp(sg, exons, f"read{i}"))
    (mp, count), = [p.get_multipath_and_count() for p in counter.get_all_MultiPathCountPairs()]
    assert count == 5
    assert mp.get_read_weight() == pytest.approx(5.0)


def test_weight_sums_registered_values_while_count_stays_a_read_tally():
    """The two quantities answer different questions and must not be conflated."""
    sg, exons = _splice_graph()
    for name, w in (("a", 1.0), ("b", 4.0), ("c", 10.0)):
        LRAA_Globals.register_read_weight(_rid(name), w)

    counter = MultiPathCounter()
    for name in ("a", "b", "c"):
        counter.add(_mp(sg, exons, name))
    (mp, count), = [p.get_multipath_and_count() for p in counter.get_all_MultiPathCountPairs()]

    assert count == 3, "three reads survived thinning"
    assert mp.get_read_weight() == pytest.approx(15.0), "they stand for fifteen"


def test_reaggregating_a_read_changes_neither_count_nor_weight():
    """Uniqueness is enforced on read IDs; weight must respect the same dedupe."""
    sg, exons = _splice_graph()
    LRAA_Globals.register_read_weight(_rid("dup"), 7.0)

    counter = MultiPathCounter()
    counter.add(_mp(sg, exons, "dup"))
    counter.add(_mp(sg, exons, "dup"))

    (mp, count), = [p.get_multipath_and_count() for p in counter.get_all_MultiPathCountPairs()]
    assert count == 1
    assert mp.get_read_weight() == pytest.approx(7.0)


def test_removing_a_read_removes_exactly_its_weight():
    sg, exons = _splice_graph()
    LRAA_Globals.register_read_weight(_rid("keep"), 2.0)
    LRAA_Globals.register_read_weight(_rid("drop"), 9.0)

    mp = _mp(sg, exons, "keep")
    mp.include_read_name("drop")
    assert mp.get_read_weight() == pytest.approx(11.0)

    mp.remove_read_name("drop")
    assert mp.get_read_weight() == pytest.approx(2.0)
    assert mp.get_read_count() == 1


def test_a_rebuilt_multipath_recovers_the_weight_without_being_told():
    """Transcriptome rescue reconstructs multipaths after minimap2 drops the tag.

    The registry is keyed on the read name's compact ID, so reconstruction anywhere
    finds the weight; without this, the rescue cohort would silently weigh 1.
    """
    sg, exons = _splice_graph()
    LRAA_Globals.register_read_weight(_rid("rescued"), 6.0)
    rebuilt = _mp(sg, exons, "rescued")
    assert rebuilt.get_read_weight() == pytest.approx(6.0)


def test_an_explicit_zero_count_is_not_an_observation_however_it_is_named():
    """A reference template must not become evidence by way of its weight.

    Input transcripts are placed in the graph as multipaths so a structure has a
    foothold, carrying a synthetic 'reftranscript:' name but a zero count. Deriving
    weight from that name would hand EM support no read stands behind.
    """
    sg, exons = _splice_graph()
    template = MultiPath(
        sg,
        [[e.get_id() for e in exons]],
        read_types={"reftranscript"},
        read_names={"reftranscript:t1"},
        read_count=0,
    )
    assert template.get_read_count() == 0
    assert template.get_read_weight() == pytest.approx(0.0)


def test_pruning_a_reference_template_removes_its_weight_too():
    """Weight must follow the IDs, not lag them.

    prune_reftranscript_as_evidence drops the template's ID and decrements the count;
    a weight left behind would keep contributing support for the rest of the run, and
    every abundance consumer now reads the weight rather than the count.
    """
    from ReadNameStore import ReadNameStore
    import tempfile

    sg, exons = _splice_graph()
    tmpdir = tempfile.mkdtemp()
    prior_store = LRAA_Globals.READ_NAME_STORE
    LRAA_Globals.READ_NAME_STORE = ReadNameStore(os.path.join(tmpdir, "names"))
    try:
        real_id = LRAA_Globals.READ_NAME_STORE.get_or_add("read_real")
        tmpl_id = LRAA_Globals.READ_NAME_STORE.get_or_add("reftranscript:t1")
        LRAA_Globals.register_read_weight(real_id, 5.0)
        LRAA_Globals.register_read_weight(tmpl_id, 1.0)

        mp = MultiPath(
            sg,
            [[e.get_id() for e in exons]],
            read_types={"PacBio"},
            read_names=set(),
            read_count=0,
        )
        mp.include_read_id(real_id)
        mp.include_read_id(tmpl_id)
        assert mp.get_read_count() == 2
        assert mp.get_read_weight() == pytest.approx(6.0)

        mp.prune_reftranscript_as_evidence()

        assert mp.get_read_count() == 1, "the template no longer counts"
        assert mp.get_read_weight() == pytest.approx(5.0), "nor does it weigh"
    finally:
        try:
            LRAA_Globals.READ_NAME_STORE.close()
        except Exception:
            pass
        LRAA_Globals.READ_NAME_STORE = prior_store


def test_registration_is_authoritative_not_cumulative():
    """A read has one weight; a later authoritative write replaces the provisional one."""
    LRAA_Globals.register_read_weight(_rid("r"), 3.0)
    LRAA_Globals.register_read_weight(_rid("r"), 5.0)
    assert LRAA_Globals.read_weight_for_id(_rid("r")) == pytest.approx(5.0)


def test_reset_scopes_weights_to_one_bam():
    """A disabled or control pass must not inherit a weighted pass's registry."""
    LRAA_Globals.register_read_weight(_rid("stale"), 8.0)
    LRAA_Globals.reset_read_weight_registry()
    assert LRAA_Globals.read_weight_for_id(_rid("stale")) == pytest.approx(1.0)


def test_uniform_thinning_recovers_the_unthinned_isoform_fractions():
    """Why weighting is worth the plumbing.

    A locus thinned uniformly at rate p keeps 1/p-weighted reads. EM's abundances are
    scale-free in the support vector, so the weighted thinned locus must reproduce the
    unthinned isoform fractions exactly, while the unweighted survivor counts do not
    reproduce the counts.
    """
    assign = [[0], [1], [0, 1]]
    w = [[1.0], [1.0], [1.0, 1.0]]
    full = [600.0, 200.0, 400.0]
    p = 0.1
    weighted = [c for c in full]              # survivors x 1/p, in expectation
    survivors = [c * p for c in full]         # what counting kept reads would give

    e_full, c_full, _ = em_algorithm_with_weights(assign, w, full, 2, max_iter=2000, tol=1e-12)
    e_w, c_w, _ = em_algorithm_with_weights(assign, w, weighted, 2, max_iter=2000, tol=1e-12)
    e_s, c_s, _ = em_algorithm_with_weights(assign, w, survivors, 2, max_iter=2000, tol=1e-12)

    assert e_w == pytest.approx(e_full, abs=1e-12)
    assert e_s == pytest.approx(e_full, abs=1e-12), "abundances are scale-free either way"
    for i in range(2):
        assert c_w[i] == pytest.approx(c_full[i], rel=1e-9), "weighted counts recover scale"
        assert c_s[i] < c_full[i] * 0.5, "survivor counts do not"
