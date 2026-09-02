"""Shared splice junctions as gene-clustering evidence, and the must-link that is not
a threshold.

Two separate things live here, and conflating them was the original defect:

  RULE 1  an EDGE, when the pair shares at least `min_recluster_shared_introns`
          introns AND that count exceeds `min_recluster_shared_intron_frac` of the
          smaller isoform's own intron count. An edge only makes a cut expensive.

  RULE 2  a MUST-LINK, when two isoforms have the SAME intron chain. Identical splice
          patterns are one gene by definition, so this cannot be expressed as an edge
          -- Leiden optimizes an objective and may still cut an edge -- and it cannot
          be derived from rule 1 either, because any threshold strict enough to
          suppress single-junction bridges also rejects an identical 1-intron pair.
          It is enforced structurally: contraction to one vertex under Leiden, an
          unconditional union under the DSU fallback.

The empty-chain case is the trap worth a permanent test. A single-exon transcript has
no introns, so "identical intron chain" is trivially true of every single-exon model
against every other, and treating that as a must-link collapses all of them into one
gene.
"""
import pytest

import LRAA_Globals
from GeneCommunityCluster import (
    _contract_identical_chains,
    _intron_sets,
    _shares_introns,
    partition_with_dsu,
)


class FakeTranscript:
    """Exon segments are all the clustering criteria read."""

    def __init__(self, exons):
        self._exons = sorted(exons)

    def get_exon_segments(self):
        return list(self._exons)

    def get_introns(self):
        return [
            (a[1] + 1, b[0] - 1) for a, b in zip(self._exons, self._exons[1:])
        ]

    def get_cdna_len(self):
        return sum(b - a + 1 for a, b in self._exons)

    def get_lend(self):
        return self._exons[0][0]

    def get_rend(self):
        return self._exons[-1][1]


def _chain(n_introns, start=1000, step=1000):
    """A transcript with n_introns introns at predictable coordinates."""
    exons = []
    pos = start
    for _ in range(n_introns + 1):
        exons.append((pos, pos + 99))
        pos += step
    return FakeTranscript(exons)


def _introns_of(t):
    return frozenset(t.get_introns())


@pytest.fixture(autouse=True)
def restore_thresholds(monkeypatch):
    """Several tests below drive these to hostile values.

    LRAA_Globals.config is process-global, so without this every later test in the
    session would inherit whichever threshold ran last -- and a floor of 99 silently
    disables the shared-intron rule everywhere.
    """
    for key in (
        "min_recluster_shared_intron_frac",
        "min_recluster_shared_introns",
    ):
        monkeypatch.setitem(LRAA_Globals.config, key, LRAA_Globals.config[key])


# --------------------------------------------------------------- rule 1, the fraction


def test_exactly_at_the_fraction_does_not_edge():
    """The fraction is a strict >, so a pair sitting exactly on it is not linked.

    Pinned because flipping this to >= silently widens every locus by one step, and
    nothing in the output would say so.
    """
    ten = frozenset((i * 10, i * 10 + 5) for i in range(10))
    two_shared = frozenset([(0, 5), (10, 15)]) | frozenset(
        (i * 1000 + 90000, i * 1000 + 90005) for i in range(8)
    )
    assert len(ten & two_shared) / min(len(ten), len(two_shared)) == pytest.approx(0.2)
    assert _shares_introns(ten, two_shared, 0.20, 2) is False


def test_above_the_fraction_edges():
    ten = frozenset((i * 10, i * 10 + 5) for i in range(10))
    three_shared = frozenset([(0, 5), (10, 15), (20, 25)]) | frozenset(
        (i * 1000 + 90000, i * 1000 + 90005) for i in range(7)
    )
    assert _shares_introns(ten, three_shared, 0.20, 2) is True


def test_denominator_is_the_smaller_intron_count():
    """A short model is judged on ITS OWN splice pattern, not the long one's.

    Two of a 2-intron model's introns inside a 20-intron model is 1.00, not 0.10.
    Normalizing by the longer isoform would reject exactly the fragmentary models the
    rule exists to keep attached.
    """
    small = frozenset([(1000, 1999), (3000, 3999)])
    large = small | frozenset((i * 10000 + 50000, i * 10000 + 50999) for i in range(18))
    assert len(small & large) / min(len(small), len(large)) == 1.0
    assert len(small & large) / max(len(small), len(large)) == pytest.approx(0.1)
    assert _shares_introns(small, large, 0.20, 2) is True


# ------------------------------------------------------------------ rule 1, the floor


def test_one_shared_junction_does_not_edge_a_small_model():
    """The floor exists because the fraction alone cannot reject this.

    One shared intron satisfies >20% whenever the smaller model has <= 4 introns, which
    is enough for a single fragment to bind two neighbouring genes into one gene.
    """
    a = frozenset([(1000, 1999), (3000, 3999)])
    b = frozenset([(1000, 1999), (9000, 9999)])
    assert len(a & b) == 1
    assert 1 / min(len(a), len(b)) == 0.5  # passes the fraction
    assert _shares_introns(a, b, 0.20, 2) is False  # rejected by the floor
    assert _shares_introns(a, b, 0.20, 1) is True  # and only by the floor


# ------------------------------------------------------------------- the empty chain


def test_single_exon_never_shares_introns():
    """No introns means no junction evidence, not universal agreement."""
    empty = frozenset()
    spliced = frozenset([(1000, 1999), (3000, 3999)])
    assert _shares_introns(empty, spliced, 0.20, 2) is False
    assert _shares_introns(empty, empty, 0.20, 2) is False


def test_single_exon_transcripts_are_never_contracted_together():
    """The trap: every single-exon model has the same (empty) chain as every other.

    Contracting on that would make one gene of all of them, and because contraction is
    a hard must-link no later step could undo it.
    """
    singles = [
        FakeTranscript([(10, 20)]),
        FakeTranscript([(30, 40)]),
        FakeTranscript([(50, 60)]),
    ]
    vertices = _contract_identical_chains(_intron_sets(singles))
    assert len(set(vertices)) == 3


# ----------------------------------------------------- rule 2, structural and absolute


def test_identical_chains_contract_to_one_vertex():
    """Same chain, different boundaries: one vertex, so no partition can separate."""
    a = FakeTranscript([(1000, 1100), (2000, 2100)])
    b = FakeTranscript([(1050, 1100), (2000, 2500)])  # same intron, other boundaries
    far = FakeTranscript([(900000, 900100), (901000, 901100)])
    assert _introns_of(a) == _introns_of(b)
    vertices = _contract_identical_chains(_intron_sets([a, b, far]))
    assert vertices[0] == vertices[1]
    assert vertices[2] != vertices[0]


@pytest.mark.parametrize("frac,floor", [(0.20, 2), (1.0, 2), (0.20, 99), (2.0, 99)])
def test_contraction_ignores_both_thresholds(monkeypatch, frac, floor):
    """The must-link is not a tunable, so no threshold may switch it off."""
    monkeypatch.setitem(LRAA_Globals.config, "min_recluster_shared_intron_frac", frac)
    monkeypatch.setitem(LRAA_Globals.config, "min_recluster_shared_introns", floor)
    a = FakeTranscript([(1000, 1100), (2000, 2100)])
    b = FakeTranscript([(1050, 1100), (2000, 2500)])
    assert _contract_identical_chains(_intron_sets([a, b])) == [0, 0]


@pytest.mark.parametrize(
    "frac,floor",
    [(0.20, 2), (1.0, 2), (0.20, 99), (2.0, 99)],
)
def test_dsu_must_link_survives_thresholds_that_disable_rule_one(
    monkeypatch, frac, floor
):
    """The DSU fallback must enforce rule 2 on its own.

    Regression: this branch once relied on an identical pair scoring 1.0 through the
    fraction rule. That fails for a fraction >= 1.0, and it fails for ANY floor above
    the chain length -- including the shipped floor of 2 against a 1-intron chain,
    which is this fixture. The union is therefore unconditional and taken before the
    sweep, so neither threshold nor the sweep's span window can drop it.
    """
    monkeypatch.setitem(LRAA_Globals.config, "min_recluster_shared_intron_frac", frac)
    monkeypatch.setitem(LRAA_Globals.config, "min_recluster_shared_introns", floor)

    a = FakeTranscript([(1000, 1100), (2000, 2100)])
    b = FakeTranscript([(1050, 1100), (2000, 2500)])
    far = FakeTranscript([(900000, 900100), (901000, 901100)])

    membership = partition_with_dsu([a, b, far], "chrTEST", "+")
    assert membership[0] == membership[1], "identical chains must be one gene"
    assert membership[2] != membership[0], "an unrelated locus must stay separate"


def test_identical_one_intron_chain_is_held_only_by_rule_two(monkeypatch):
    """Shows rule 2 is not redundant with rule 1 under the shipped settings.

    A 1-intron identical pair shares one intron, so the floor of 2 rejects it as an
    edge. If rule 2 were derived from rule 1 rather than enforced structurally, this
    pair would be free to land in two genes.
    """
    monkeypatch.setitem(LRAA_Globals.config, "min_recluster_shared_intron_frac", 0.20)
    monkeypatch.setitem(LRAA_Globals.config, "min_recluster_shared_introns", 2)

    a = FakeTranscript([(1000, 1100), (2000, 2100)])
    b = FakeTranscript([(1050, 1100), (2000, 2500)])
    assert _shares_introns(_introns_of(a), _introns_of(b), 0.20, 2) is False
    assert _contract_identical_chains(_intron_sets([a, b])) == [0, 0]


def test_dsu_intron_rule_is_evaluated_before_the_exon_overlap_prefilter(monkeypatch):
    """The intron OR must not sit behind the cDNA-overlap gate.

    That gate needs span overlap >= 20% of the LONGER cDNA. A short model whose whole
    span sits inside a long gene fails it however well their junctions agree, so an
    intron rule evaluated after it would never fire for exactly the fragmentary models
    it is for.

    The fixture is built so the pair CANNOT be linked any other way: 2 shared introns
    out of the smaller model's 2, exon overlap gated out by the prefilter, and
    different chains so the identical-chain must-link does not apply. Raising the floor
    out of reach must therefore separate them -- which is what proves the union came
    from the intron rule rather than from the exon branch.
    """
    monkeypatch.setitem(LRAA_Globals.config, "min_recluster_shared_intron_frac", 0.20)
    monkeypatch.setitem(LRAA_Globals.config, "min_recluster_shared_introns", 2)

    short = FakeTranscript([(990, 1010), (2000, 2010), (3000, 3010)])
    long_ = FakeTranscript(
        [(1000, 1010), (2000, 2010), (3000, 3010), (100000, 200000)]
    )
    far = FakeTranscript([(900000, 900100), (901000, 901100)])

    shared = _introns_of(short) & _introns_of(long_)
    assert len(shared) == 2, "fixture must share two introns"
    assert _introns_of(short) != _introns_of(long_), "must not be an identical chain"

    # The prefilter the intron rule has to precede.
    span_overlap = 3010 - 1000 + 1
    shorter_cdna, longer_cdna = short.get_cdna_len(), long_.get_cdna_len()
    min_required = max(int(0.50 * shorter_cdna), int(0.20 * longer_cdna))
    assert span_overlap < min_required, (
        "fixture must FAIL the exon-overlap prefilter, else it proves nothing"
    )

    membership = partition_with_dsu([short, long_, far], "chrTEST", "+")
    assert membership[0] == membership[1], (
        "two shared introns must link the pair despite failing the cDNA-overlap gate"
    )
    assert membership[2] != membership[0]

    # Floor out of reach: rule 1 cannot fire, and nothing else can link them.
    monkeypatch.setitem(LRAA_Globals.config, "min_recluster_shared_introns", 99)
    blocked = partition_with_dsu([short, long_, far], "chrTEST", "+")
    assert blocked[0] != blocked[1], (
        "with rule 1 disabled the pair must separate; if not, the exon branch linked "
        "them and this test no longer isolates the intron rule"
    )


def test_production_graph_carries_multiplicity_and_self_loops(monkeypatch):
    """Assert on the graph partition_with_leiden actually builds.

    The objective-equivalence test below reasons about a graph it constructs itself, so
    it would keep passing if the production path stopped weighting. This captures the
    real igraph object handed to find_partition.
    """
    pytest.importorskip("igraph")
    leidenalg = pytest.importorskip("leidenalg")
    import GeneCommunityCluster as G

    monkeypatch.setitem(LRAA_Globals.config, "min_recluster_shared_intron_frac", 0.20)
    monkeypatch.setitem(LRAA_Globals.config, "min_recluster_shared_introns", 2)

    captured = {}
    real_find_partition = leidenalg.find_partition

    def recording_find_partition(graph, partition_type, **kwargs):
        captured["edges"] = graph.get_edgelist()
        captured["weights"] = list(graph.es["weight"]) if graph.ecount() else []
        captured["weights_kwarg"] = kwargs.get("weights")
        captured["vcount"] = graph.vcount()
        return real_find_partition(graph, partition_type, **kwargs)

    monkeypatch.setattr(G.leidenalg, "find_partition", recording_find_partition)

    # three transcripts sharing ONE chain -> one vertex, plus a neighbour linked to all
    # three by two shared introns.
    chain = [(1000, 1099), (2000, 2099), (3000, 3099), (4000, 4099)]
    t1 = FakeTranscript(chain)
    t2 = FakeTranscript([(900, 1099)] + chain[1:])
    t3 = FakeTranscript(chain[:-1] + [(4000, 4200)])
    neighbour = FakeTranscript([(1000, 1099), (2000, 2099), (3000, 3099)])

    assert _introns_of(t1) == _introns_of(t2) == _introns_of(t3)
    assert len(_introns_of(neighbour) & _introns_of(t1)) == 2

    G.partition_with_leiden([t1, t2, t3, neighbour], "chrTEST", "+")

    assert captured["weights_kwarg"] == "weight", "weights must be passed to Leiden"
    assert captured["vcount"] == 2, "the identical-chain class must be one vertex"

    weight_of = dict(zip(captured["edges"], captured["weights"]))
    normalized = {
        (min(a, b), max(a, b)): w for (a, b), w in weight_of.items()
    }
    assert normalized.get((0, 0)) == 3, (
        "the class's 3 internal pairs must survive as a self-loop of weight 3; "
        "dropping it strips the vertex of the degree its members had"
    )
    assert normalized.get((0, 1)) == 3, (
        "all 3 transcript-level edges to the neighbour must survive as weight 3, "
        "not be deduplicated to 1"
    )


# ------------------------------------------------- contraction preserves the objective


def test_weighted_contraction_preserves_the_rbconfiguration_objective():
    """Contraction must impose the must-link WITHOUT changing what Leiden optimizes.

    RBConfiguration's null model is degree-based, so collapsing an identical-chain
    class to one vertex and deduplicating its edges strips the class of most of its
    degree: Leiden then maximizes a different function. Keeping edge multiplicity as
    weights and internal edges as self-loops makes the contracted quality EQUAL to the
    full-graph quality for every must-link-respecting partition, so the argmax is the
    original one subject to the constraint.

    Asserted as a constant difference, since equal argmax is all that is required.
    """
    igraph = pytest.importorskip("igraph")
    leidenalg = pytest.importorskip("leidenalg")

    import itertools

    klass = [0, 1, 2, 3]  # one identical-chain class
    edges = [(m, nb) for m in klass for nb in (4, 5, 6)]
    edges += list(itertools.combinations(klass, 2))  # internal to the class
    edges += [(4, 5), (5, 6), (6, 7), (7, 8), (4, 8)]
    n = 9
    vertex_of = [0, 0, 0, 0, 1, 2, 3, 4, 5]
    gamma = 0.2

    def contracted(dedup):
        weights = {}
        for i, j in edges:
            a, b = vertex_of[i], vertex_of[j]
            key = (a, b) if a <= b else (b, a)
            if dedup:
                if a == b:
                    continue
                weights[key] = 1
            else:
                weights[key] = weights.get(key, 0) + 1
        keys = sorted(weights)
        g = igraph.Graph(n=max(vertex_of) + 1, edges=keys, directed=False)
        g.es["weight"] = [weights[k] for k in keys]
        return g

    full = igraph.Graph(n=n, edges=edges, directed=False)
    weighted, deduped = contracted(False), contracted(True)

    vertex_partitions = [
        [0, 0, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 1],
        [0, 0, 1, 1, 2, 2],
        [0, 1, 2, 0, 1, 2],
        [0, 0, 0, 1, 1, 1],
    ]

    def quality(graph, membership, weights=None):
        return leidenalg.RBConfigurationVertexPartition(
            graph,
            initial_membership=membership,
            weights=weights,
            resolution_parameter=gamma,
        ).quality()

    weighted_gaps, deduped_gaps = [], []
    for vm in vertex_partitions:
        expanded = [vm[vertex_of[i]] for i in range(n)]
        q_full = quality(full, expanded)
        weighted_gaps.append(q_full - quality(weighted, vm, "weight"))
        deduped_gaps.append(q_full - quality(deduped, vm, "weight"))

    assert max(weighted_gaps) - min(weighted_gaps) < 1e-9, (
        "weighted contraction must preserve the objective up to a constant"
    )
    assert max(deduped_gaps) - min(deduped_gaps) > 1.0, (
        "deduplicated contraction is expected to distort the objective; if this no "
        "longer holds, the fixture stopped exercising edge multiplicity"
    )
