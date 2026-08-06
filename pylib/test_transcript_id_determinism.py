#!/usr/bin/env python3

import random
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA_Globals
from Transcript import Transcript


def _transcript(exons, simple_path):
    transcript = Transcript("chr1", exons, "+")
    transcript.set_simple_path(simple_path)
    return transcript


def _fixture_transcripts():
    """Several loci, including two clusters sharing a leftmost coordinate."""
    return [
        # locus A: two isoforms sharing a node, so one cluster
        _transcript([[1000, 1100], [1300, 1400]], ["E:1", "I:1", "E:2"]),
        _transcript([[1000, 1100], [1300, 1450]], ["E:1", "I:1", "E:3"]),
        # locus B: same leftmost coordinate as locus C, different right end
        _transcript([[2000, 2100], [2300, 2400]], ["E:4", "I:2", "E:5"]),
        # locus C: shares the leftmost coordinate of locus B but is shorter
        _transcript([[2000, 2050]], ["E:6"]),
        # loci D and E: identical spans, distinct structures
        _transcript([[3000, 3050]], ["E:7"]),
        _transcript([[3000, 3050]], ["E:8"]),
    ]


def _structure_token(transcript):
    """Identity of a transcript independent of any code under test."""
    simple_path = transcript.get_simple_path()
    return (
        transcript.get_exons_string(),
        ",".join(str(node_id) for node_id in simple_path or []),
    )


def _assign_ids(transcripts):
    revised = Transcript.recluster_transcripts_to_genes(transcripts, "chr1", "+")
    return {
        _structure_token(t): (t.get_gene_id(), t.get_transcript_id()) for t in revised
    }


def test_identifiers_do_not_depend_on_input_order(monkeypatch):
    monkeypatch.setitem(LRAA_Globals.config, "debug_write_init_clusters", False)

    baseline = _assign_ids(_fixture_transcripts())

    rng = random.Random(17)
    for _ in range(5):
        shuffled = _fixture_transcripts()
        rng.shuffle(shuffled)
        assert _assign_ids(shuffled) == baseline


def test_component_numbering_follows_genomic_order(monkeypatch):
    monkeypatch.setitem(LRAA_Globals.config, "debug_write_init_clusters", False)

    transcripts = _fixture_transcripts()
    leftmost_token = _structure_token(min(transcripts, key=lambda t: t.get_coords()))

    assigned = _assign_ids(transcripts)

    assert assigned[leftmost_token][0] == "g:chr1:+:comp-1"
    # every structure keeps its own identifier, including the same-span pair
    assert len(set(assigned.values())) == len(assigned)
