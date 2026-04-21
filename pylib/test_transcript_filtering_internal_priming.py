#!/usr/bin/env python3

import TranscriptFiltering
from Transcript import Transcript


def _build_transcript(transcript_id, exons, strand="+"):
    transcript = Transcript("chr1", exons, strand)
    transcript.set_gene_id("g1")
    transcript.set_transcript_id(transcript_id)
    return transcript


def test_restrict_to_monoexonic_does_not_force_multiexonic_internal_priming_true():
    multiexonic = _build_transcript("t.multi", [[10, 20], [40, 50]], "+")

    retained = TranscriptFiltering.filter_internally_primed_transcripts(
        [multiexonic],
        contig_seq_str="C" * 100,
        contig_strand="+",
        known_transcripts=None,
        restrict_filter_to_monoexonic=True,
    )

    assert [transcript.get_transcript_id() for transcript in retained] == ["t.multi"]
    assert retained[0]._likely_internal_primed is False


def test_restrict_to_monoexonic_preserves_true_annotation_for_multiexonic_hits():
    multiexonic = _build_transcript("t.multi", [[10, 20], [40, 50]], "+")

    retained = TranscriptFiltering.filter_internally_primed_transcripts(
        [multiexonic],
        contig_seq_str="C" * 50 + "AAAAAAAAAAAA" + "C" * 38,
        contig_strand="+",
        known_transcripts=None,
        restrict_filter_to_monoexonic=True,
    )

    assert [transcript.get_transcript_id() for transcript in retained] == ["t.multi"]
    assert retained[0]._likely_internal_primed is True
