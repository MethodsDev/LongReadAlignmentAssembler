#!/usr/bin/env python3

import TranscriptFiltering
from Transcript import Transcript


class _DummySpliceGraph:
    def get_contig_acc(self):
        return "chr1"

    def get_contig_strand(self):
        return "-"


class _DummyMP:
    def __init__(self, read_count):
        self._read_count = read_count

    def get_read_count(self):
        return self._read_count


def _build_transcript(transcript_id, exons, simple_path, read_count):
    transcript = Transcript("chr1", exons, "-")
    transcript.set_gene_id("g1")
    transcript.set_transcript_id(transcript_id)
    transcript.set_simple_path(simple_path)
    transcript.set_read_counts_assigned(read_count)
    return transcript


def test_prune_degradation_products_handles_terminal_split_exon_paths():
    transcript_with_boundaries = _build_transcript(
        "t.keep",
        [[100, 150], [200, 250]],
        ["POLYA:1", "E:1", "E:2", "I:1", "E:3", "TSS:1"],
        10.0,
    )
    transcript_terminal_trunc = _build_transcript(
        "t.prune",
        [[101, 150], [200, 260]],
        ["E:2", "I:1", "E:3", "E:4"],
        3.0,
    )

    frac_read_assignments = {
        "t.keep": {_DummyMP(10): 1.0},
        "t.prune": {_DummyMP(3): 1.0},
    }

    retained = TranscriptFiltering.prune_likely_degradation_products(
        [transcript_with_boundaries, transcript_terminal_trunc],
        _DummySpliceGraph(),
        frac_read_assignments,
    )

    retained_ids = {transcript.get_transcript_id() for transcript in retained}
    assert retained_ids == {"t.keep"}
