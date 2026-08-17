#!/usr/bin/env python3

"""The order a cluster's members are handed to gene clustering must be intrinsic.

`recluster_transcripts_to_genes` indexes transcripts by their position in a
contig-wide structural sort, so the index of any given transcript depends on how
many transcripts precede it -- i.e. on how much of the contig the run processed.
Materialising a cluster by iterating the raw `set` that
`networkx.connected_components` yields exposed that: CPython iterates a set of
ints in ascending `value mod table_size`, so shifting every index by a constant
rotates the member order.  That order becomes igraph vertex ids and
`leidenalg.find_partition` is order-sensitive, so the same transcripts were
called as different genes depending on the processed extent, and every isoform
filter downstream is gene-scoped.

Measured consequence before the fix, chr21 ref-guided HiFi: 4 models lost and 2
gained between a 46.7 Mb whole-contig run and a 9.56 Mb chunk over the same
extent, from 2 of the 121 clusters the two arms shared.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import GeneCommunityCluster
import LRAA_Globals
from Transcript import Transcript

# The real chr21 geometry: the divergent 30-member cluster sat at contig-wide
# index base 616 in the chunked arm and 2528 in the whole-contig arm.
CLUSTER_SIZE = 30
INDEX_BASES = (0, 40, 616, 2528)


def _transcript(exons, simple_path):
    transcript = Transcript("chr21", exons, "+")
    transcript.set_simple_path(simple_path)
    return transcript


def _structure_token(transcript):
    """Identity of a transcript independent of any code under test."""
    simple_path = transcript.get_simple_path()
    return (
        transcript.get_exons_string(),
        ",".join(str(node_id) for node_id in simple_path or []),
    )


def _fixture(index_base):
    """One overlap cluster of CLUSTER_SIZE, preceded by `index_base` singletons.

    The padding transcripts are upstream, single-exon and share no simple-path
    node with anything, so each is its own initial cluster and their only effect
    is to push the cluster of interest to contig-wide indices
    `index_base .. index_base + CLUSTER_SIZE - 1`.
    """
    padding = [
        _transcript([[100 + 10 * k, 104 + 10 * k]], [f"E:pad-{k}"])
        for k in range(index_base)
    ]
    # All members share node "E:hub", so the first-pass graph makes them one
    # connected component; distinct right ends give them a distinct intrinsic
    # order under structural_sort_key.
    cluster = [
        _transcript(
            [[1000000, 1000100], [1000000 + 1000 * m, 1000200 + 1000 * m]],
            ["E:hub", f"I:{m}", f"E:{m}"],
        )
        for m in range(1, CLUSTER_SIZE + 1)
    ]
    return padding + cluster


def _recorded_member_order(index_base, monkeypatch):
    """The member order `partition_with_leiden` is actually handed."""
    recorded = []

    def _recording_partition(transcripts, contig_acc, contig_strand, **kwargs):
        recorded.append([_structure_token(t) for t in transcripts])
        return list(range(len(transcripts)))

    monkeypatch.setattr(
        GeneCommunityCluster, "partition_with_leiden", _recording_partition
    )
    Transcript.recluster_transcripts_to_genes(_fixture(index_base), "chr21", "+")

    assert len(recorded) == 1, (
        f"expected exactly one multi-member cluster, got {len(recorded)}"
    )
    assert len(recorded[0]) == CLUSTER_SIZE
    return recorded[0]


def test_cluster_member_order_is_extent_invariant(monkeypatch):
    monkeypatch.setitem(LRAA_Globals.config, "debug_write_init_clusters", False)
    monkeypatch.setitem(LRAA_Globals.config, "use_community_clustering", True)

    orders = {
        base: _recorded_member_order(base, monkeypatch) for base in INDEX_BASES
    }

    reference_base = INDEX_BASES[0]
    reference = orders[reference_base]
    for base, order in orders.items():
        assert order == reference, (
            f"cluster member order at index base {base} differs from index base "
            f"{reference_base}: the order handed to gene clustering depends on "
            f"processed extent"
        )


def test_cluster_member_order_is_the_structural_sort_order(monkeypatch):
    """Invariance comes from the members' own key, not from a happy rotation."""
    monkeypatch.setitem(LRAA_Globals.config, "debug_write_init_clusters", False)
    monkeypatch.setitem(LRAA_Globals.config, "use_community_clustering", True)

    for base in INDEX_BASES:
        order = _recorded_member_order(base, monkeypatch)
        expected = [
            _structure_token(t)
            for t in sorted(
                (t for t in _fixture(base) if "E:hub" in t.get_simple_path()),
                key=Transcript.structural_sort_key,
            )
        ]
        assert order == expected, f"index base {base}"
