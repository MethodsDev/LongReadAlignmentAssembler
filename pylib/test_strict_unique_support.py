#!/usr/bin/env python3

"""Unique support must mean EXCLUSIVELY ASSIGNABLE, not "EM gave it almost all the mass".

`get_isoform_unique_support` counts a multipath only when exactly one isoform of the
gene can be assigned it (`n_compat[mp] == 1`). It previously counted any multipath whose
fractional assignment cleared `unique_read_filter_min_frac` (0.9995), which is a
statement about how EM happened to distribute mass, not about exclusivity: a read
compatible with two isoforms still clears it whenever EM strongly favours one.

The two definitions coincide on any corpus where no read is both multi-compatible and
strongly allocated, and such a corpus cannot exercise the change at all: the outputs are
identical either way. Where they differ it is on a small minority of decisions. This
file builds the divergent population directly instead of hoping a corpus contains it.

The fixture is one multipath reachable by two isoforms. Rather than relying on EM to
drive its fraction above 0.9995 -- which depends on convergence dynamics and, measured
here, does not happen with 50:1 evidence -- the test lowers
`unique_read_filter_min_frac` so that an equal split clears it. That isolates the
semantics under test: whether a threshold on ALLOCATION can admit a read, versus
EXCLUSIVITY refusing it. The threshold's value is immaterial to that distinction, and
the key is now a diagnostic anyway.

Both magnitudes are asserted, because one membership test now feeds both: the LITERAL
count consumed by the novel-isoform floor and the rescue, and the WEIGHTED support
consumed by frac_gene_unique_reads. The shared read is given a large weight so that
including it in the weighted figure is unmistakable.
"""

import os
import sys
import tempfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA_Globals
from GenomeFeature import Exon
from MultiPath import MultiPath
from Splice_graph import Splice_graph

CONTIG = "chr1"
#           0            1            2            3
EXONS = [(100, 150), (200, 250), (400, 450), (600, 650)]
N_EXCLUSIVE = 50          # exclusive reads on the dominant isoform, weight 1 each
SHARED_READ_WEIGHT = 100.0  # so a weighted figure that wrongly includes it is obvious


@pytest.fixture(autouse=True)
def _isolate():
    keys = (
        "min_frac_gene_unique_reads",
        "min_FSM_reads_gate",
        "min_FSM_reads_retain_isoform",
        "min_unique_reads_novel_isoform",
        "unique_read_filter_min_frac",
        "aggregate_splice_boundary_max_rel_support",
        "aggressively_assign_reads",
        "show_progress_quant_assign",
        "fraction_read_align_overlap",
    )
    saved = {k: LRAA_Globals.config[k] for k in keys}
    LRAA_Globals.config["min_frac_gene_unique_reads"] = 0.0
    LRAA_Globals.config["min_FSM_reads_gate"] = 0
    LRAA_Globals.config["min_FSM_reads_retain_isoform"] = 0
    LRAA_Globals.config["min_unique_reads_novel_isoform"] = 0
    # An equal two-way split is 0.5, so this makes the shared read clear the OLD
    # allocation gate deterministically. Under exclusivity it must still be refused.
    LRAA_Globals.config["unique_read_filter_min_frac"] = 0.4
    LRAA_Globals.config["aggressively_assign_reads"] = False
    LRAA_Globals.config["show_progress_quant_assign"] = False
    LRAA_Globals.config["fraction_read_align_overlap"] = 0.4
    LRAA_Globals.reset_read_weight_registry()
    try:
        yield
    finally:
        LRAA_Globals.config.update(saved)
        LRAA_Globals.reset_read_weight_registry()


def _graph():
    Exon.reset_counter()
    sg = Splice_graph()
    exons = [Exon(CONTIG, lend, rend, "+", 1) for lend, rend in EXONS]
    for e in exons:
        sg._node_id_to_node[e.get_id()] = e
    sg._contig_acc = CONTIG
    sg._contig_strand = "+"
    return sg, [e.get_id() for e in exons]


def _transcript(sg, node_ids, idxs, tid):
    t = MultiPath(sg, [[node_ids[i] for i in idxs]]).toTranscript()
    t.set_gene_id("g1")
    t.set_transcript_id(tid)
    return t


def _mp(sg, node_ids, idxs, read_name, weight=None):
    if weight is not None:
        LRAA_Globals.register_read_weight(
            MultiPath._coerce_read_identifier(None, read_name), weight
        )
    return MultiPath(
        sg,
        [[node_ids[i] for i in idxs]],
        read_types={"PacBio"},
        read_names={read_name},
        read_count=1,
    )


def _fixture():
    """dominant = exons 0,1,2 with 50 exclusive reads; scarce = exons 0,1,3 with nothing
    of its own. ONE multipath over exons 0,1 is handed to BOTH, so it is compatible with
    two isoforms while EM concentrates its mass on the dominant one."""
    sg, node_ids = _graph()
    dominant = _transcript(sg, node_ids, [0, 1, 2], "t_dominant")
    scarce = _transcript(sg, node_ids, [0, 1, 3], "t_scarce")

    exclusive = [_mp(sg, node_ids, [0, 1, 2], f"ex{i}") for i in range(N_EXCLUSIVE)]
    dominant.add_multipaths_evidence_assigned(exclusive)

    shared = _mp(sg, node_ids, [0, 1], "shared", weight=SHARED_READ_WEIGHT)
    dominant.add_multipaths_evidence_assigned([shared])
    scarce.add_multipaths_evidence_assigned([shared])

    return sg, dominant, scarce, shared


def _decisions(transcripts):
    """Drive the real filter with decision logging on and return rows by transcript_id."""
    import importlib

    with tempfile.TemporaryDirectory() as td:
        os.environ["LRAA_FILTER_DECISION_LOG"] = os.path.join(td, "f")
        import TranscriptFiltering

        importlib.reload(TranscriptFiltering)
        TranscriptFiltering.filter_isoforms_by_min_isoform_fraction(
            transcripts, min_isoform_fraction=0.0, run_EM=False, max_EM_iterations=1
        )
        rows = {}
        for name in os.listdir(td):
            if not name.endswith(".tsv") or "prefilter" in name:
                continue
            with open(os.path.join(td, name)) as fh:
                hdr = fh.readline().rstrip("\n").split("\t")
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    rows[f[hdr.index("transcript_id")]] = dict(zip(hdr, f))
        os.environ.pop("LRAA_FILTER_DECISION_LOG", None)
        return rows


def test_fixture_actually_builds_the_divergent_population():
    """Vacuity guard: the shared read must be multi-compatible AND clear the old gate.

    If EM does not concentrate its mass, the old threshold would have excluded it too,
    the two definitions would agree here, and the assertions below would prove nothing --
    which is the blind spot any corpus has where the two criteria coincide.
    """
    sg, dominant, scarce, shared = _fixture()
    rows = _decisions([dominant, scarce])
    d = rows["t_dominant"]

    assert shared.get_read_count() == 1
    assert int(d["near_unique"]) == N_EXCLUSIVE + 1, (
        "the shared read must clear unique_read_filter_min_frac for the dominant "
        "isoform, or this fixture does not exercise the changed behaviour: with "
        "equal fractional assignment its share is 0.5 against a gate of 0.4"
    )
    assert int(d["effectively_unique"]) == 1, (
        "the shared read must be recorded as clearing the old gate while not being "
        "exclusively assignable -- that is the whole divergent population"
    )


def test_a_multi_compatible_read_contributes_to_neither_magnitude():
    """The behaviour under test: exclusivity, not allocation.

    The shared read clears the allocation gate for the dominant isoform, so the old
    rule counted it in BOTH the literal tally and the weighted support. It is compatible with two
    isoforms, so exclusivity counts it in neither.
    """
    sg, dominant, scarce, shared = _fixture()
    rows = _decisions([dominant, scarce])
    d = rows["t_dominant"]

    assert int(d["unique_reads"]) == N_EXCLUSIVE, (
        "the literal unique count included a read compatible with two isoforms: "
        "get_isoform_unique_support is gating on fractional assignment, not on "
        "compatible-map membership"
    )

    gene_reads = float(d["gene_reads"])
    weighted = float(d["frac_unique"]) * gene_reads
    assert weighted == pytest.approx(float(N_EXCLUSIVE), rel=1e-6), (
        "the WEIGHTED unique support included the shared read's weight of "
        f"{SHARED_READ_WEIGHT}: the membership test did not reach the weighted "
        "accumulator, so frac_gene_unique_reads is still on the old definition"
    )


def test_a_sole_compatible_read_does_contribute():
    """Converse, so the assertions above are not passing because nothing counts.

    Same graph, but the extra read is reachable only by the scarce isoform. It is
    exclusively assignable, so it must appear in both of that isoform's magnitudes.
    """
    sg, node_ids = _graph()
    dominant = _transcript(sg, node_ids, [0, 1, 2], "t_dominant")
    scarce = _transcript(sg, node_ids, [0, 1, 3], "t_scarce")

    dominant.add_multipaths_evidence_assigned(
        [_mp(sg, node_ids, [0, 1, 2], f"ex{i}") for i in range(N_EXCLUSIVE)]
    )
    sole = _mp(sg, node_ids, [0, 1, 3], "sole", weight=SHARED_READ_WEIGHT)
    scarce.add_multipaths_evidence_assigned([sole])

    rows = _decisions([dominant, scarce])
    s = rows["t_scarce"]

    assert int(s["unique_reads"]) == 1, (
        "a read only one isoform can be assigned was not counted as unique support"
    )
    weighted = float(s["frac_unique"]) * float(s["gene_reads"])
    assert weighted == pytest.approx(SHARED_READ_WEIGHT, rel=1e-6), (
        "the exclusively-assignable read's weight did not reach the weighted support"
    )
