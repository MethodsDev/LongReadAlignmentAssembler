#!/usr/bin/env python3

import sys
from pathlib import Path

import networkx as nx

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA_Globals
import TranscriptFiltering
from MultiPathGraph import MultiPathGraph
from Transcript import Transcript


class _FakeMultiPath:
    def __init__(self, read_types):
        self._read_types = set(read_types)

    def get_read_types(self):
        return set(self._read_types)


class _FakeMPGN:
    def __init__(self, seq_length, read_types):
        self._seq_length = seq_length
        self._multipath = _FakeMultiPath(read_types)

    def get_seq_length(self):
        return self._seq_length

    def get_multiPathObj(self):
        return self._multipath


def _graph_shell():
    """A MultiPathGraph with only the node container remove_small_components touches."""
    mpg = MultiPathGraph.__new__(MultiPathGraph)
    mpg._mp_graph = nx.DiGraph()
    return mpg


def test_short_read_only_component_is_removed():
    mpg = _graph_shell()
    component = [_FakeMPGN(150, {"PacBio"})]

    surviving = mpg.remove_small_components([component], 200)

    assert surviving == []


def test_short_component_with_input_transcript_path_is_retained():
    mpg = _graph_shell()
    reference_component = [_FakeMPGN(161, {"reftranscript"})]
    merge_component = [_FakeMPGN(161, {"fake_for_merge"})]
    read_only_component = [_FakeMPGN(161, {"PacBio"})]

    surviving = mpg.remove_small_components(
        [reference_component, merge_component, read_only_component], 200
    )

    assert surviving == [reference_component, merge_component]


def test_long_read_only_component_is_retained():
    mpg = _graph_shell()
    component = [_FakeMPGN(150, {"PacBio"}), _FakeMPGN(400, {"PacBio"})]

    surviving = mpg.remove_small_components([component], 200)

    assert surviving == [component]


def _transcript(transcript_id, exons, reference_ids=(), read_counts=10.0):
    transcript = Transcript("chr1", exons, "+")
    transcript.set_gene_id("g1")
    transcript.set_transcript_id(transcript_id)
    transcript.set_read_counts_assigned(read_counts)
    if reference_ids:
        transcript.set_source_reference_transcript_ids(set(reference_ids))
    return transcript


def test_min_length_filter_exempts_expressed_annotation_contained_model(monkeypatch):
    monkeypatch.setitem(
        LRAA_Globals.config, "ref_trans_filter_mode", "retain_expressed"
    )
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 1000)

    short_with_reference = _transcript(
        "short.ref", [[100, 260]], reference_ids=["ref1"]
    )
    short_without_reference = _transcript("short.novel", [[400, 560]])
    long_novel = _transcript("long.novel", [[1000, 1400]])

    retained = TranscriptFiltering.filter_transcripts_by_min_length(
        [short_with_reference, short_without_reference, long_novel], 200
    )

    retained_ids = {transcript.get_transcript_id() for transcript in retained}
    assert retained_ids == {"short.ref", "long.novel"}


def test_min_length_filter_drops_unexpressed_annotation_contained_model(monkeypatch):
    monkeypatch.setitem(
        LRAA_Globals.config, "ref_trans_filter_mode", "retain_expressed"
    )
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 1000)

    short_with_reference = _transcript(
        "short.ref", [[100, 260]], reference_ids=["ref1"], read_counts=0.0
    )

    retained = TranscriptFiltering.filter_transcripts_by_min_length(
        [short_with_reference], 200
    )

    assert retained == []
