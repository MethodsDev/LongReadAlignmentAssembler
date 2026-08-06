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


def test_internal_priming_ignores_21st_downstream_base_on_forward_strand():
    downstream_twenty = "ACACACAC" + "AAAAAAACCCCC"
    contig_seq = "C" * 10 + downstream_twenty + "A" + "C" * 10

    assert TranscriptFiltering._looks_internally_primed(1, 10, "+", contig_seq) is False


def test_internal_priming_ignores_21st_downstream_base_on_reverse_strand():
    downstream_twenty = "TTTTTTTCCCCC" + "TCTCTCTC"
    contig_seq = "C" * 9 + "T" + downstream_twenty + "C" * 10

    assert (
        TranscriptFiltering._looks_internally_primed(31, 40, "-", contig_seq) is False
    )
# With restrict_filter_to_monoexonic disabled, multi-exonic transcripts reach the
# known-3'-end check instead of being retained unconditionally. That branch is the only
# consumer of the known-3'-end interval tree, and it is unreachable under the default
# config (restrict_internal_priming_filter_to_monoexonic is True), so it is covered here.
# The pair brackets the branch: same transcript and sequence, differing only in whether a
# known 3' end vouches for the end position.

_PRIMED_CONTIG = "C" * 50 + "A" * 12 + "C" * 38


def test_multiexonic_internal_priming_filtered_without_a_known_3prime_end():
    multiexonic = _build_transcript("t.multi", [[10, 20], [40, 50]], "+")

    retained = TranscriptFiltering.filter_internally_primed_transcripts(
        [multiexonic],
        contig_seq_str=_PRIMED_CONTIG,
        contig_strand="+",
        known_transcripts=None,
        restrict_filter_to_monoexonic=False,
    )

    assert retained == []
    assert multiexonic._likely_internal_primed is True


def test_multiexonic_internal_priming_spared_by_a_known_3prime_end():
    multiexonic = _build_transcript("t.multi", [[10, 20], [40, 50]], "+")
    # a known model ending at the same 3' coordinate vouches for that end
    known = _build_transcript("t.known", [[30, 50]], "+")

    retained = TranscriptFiltering.filter_internally_primed_transcripts(
        [multiexonic],
        contig_seq_str=_PRIMED_CONTIG,
        contig_strand="+",
        known_transcripts=[known],
        restrict_filter_to_monoexonic=False,
    )

    assert [transcript.get_transcript_id() for transcript in retained] == ["t.multi"]
    # the transcript is still judged internally primed; the known end is what spares it
    assert retained[0]._likely_internal_primed is True


def test_monoexonic_internal_priming_filtered_despite_known_3prime_by_default():
    """Default: a known 3' end does not spare a monoexonic model."""
    monoexonic = _build_transcript("t.mono", [[10, 50]], "+")
    known = _build_transcript("t.known", [[30, 50]], "+")

    retained = TranscriptFiltering.filter_internally_primed_transcripts(
        [monoexonic],
        contig_seq_str=_PRIMED_CONTIG,
        contig_strand="+",
        known_transcripts=[known],
        restrict_filter_to_monoexonic=True,
    )

    assert retained == []
    assert monoexonic._likely_internal_primed is True


def test_monoexonic_internal_priming_spared_when_option_enabled():
    monoexonic = _build_transcript("t.mono", [[10, 50]], "+")
    known = _build_transcript("t.known", [[30, 50]], "+")

    retained = TranscriptFiltering.filter_internally_primed_transcripts(
        [monoexonic],
        contig_seq_str=_PRIMED_CONTIG,
        contig_strand="+",
        known_transcripts=[known],
        restrict_filter_to_monoexonic=True,
        spare_monoexonic_with_known_3prime=True,
    )

    assert [transcript.get_transcript_id() for transcript in retained] == ["t.mono"]
    # spared, not re-evaluated: it still looks internally primed
    assert retained[0]._likely_internal_primed is True


def test_monoexonic_internal_priming_still_filtered_without_known_3prime():
    """The option spares only models whose 3' end a reference vouches for."""
    monoexonic = _build_transcript("t.mono", [[10, 50]], "+")

    retained = TranscriptFiltering.filter_internally_primed_transcripts(
        [monoexonic],
        contig_seq_str=_PRIMED_CONTIG,
        contig_strand="+",
        known_transcripts=None,
        restrict_filter_to_monoexonic=True,
        spare_monoexonic_with_known_3prime=True,
    )

    assert retained == []
