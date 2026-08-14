#!/usr/bin/env python3

"""The internal-priming veto waives itself where the reference calls a 3' end.

A read-derived PolyA candidate in A-rich genomic context is normally rejected before it
can become a graph vertex. But if the supplied reference annotation also calls a 3' end
at that position, the reference is evidence about cleavage that is independent of the
genomic context the veto inspects, so the veto is waived.

This is the DCT case from chr13: the veto discarded a site GENCODE calls to the base,
because 15 of the 20 bases downstream in transcript sense were T.

"At that position" means within max_dist_between_alt_polyA_sites / 2, i.e. +/-25 nt, not
an exact coordinate match; see Splice_graph._reference_endorses_polyA_site.

Inert without a reference, so ref-free behaviour is unchanged.
"""

import logging

import pytest

import LRAA_Globals
import Splice_graph
from Transcript import Transcript


def _guide(transcript_id, exons, strand):
    transcript = Transcript("chr1", exons, strand)
    transcript.set_gene_id("g")
    transcript.set_transcript_id(transcript_id)
    return transcript


def _polyA_sites(contig_seq_str, counter, strand="+", guides=None, spare=True):
    sg = Splice_graph.Splice_graph()
    sg._contig_seq_str = contig_seq_str
    LRAA_Globals.config["spare_polyA_veto_at_known_3prime"] = spare
    sg._reference_three_prime_ends = sg._collect_reference_three_prime_ends(
        guides, strand
    )
    sg._incorporate_PolyA_objects("chr1", strand, counter, from_reads=True)
    return sorted(obj.get_coords()[0] for obj in sg._PolyA_objs)


A_RICH_AT_1000 = "C" * 1000 + "A" * 20 + "C" * 1980


@pytest.fixture(autouse=True)
def _restore_config():
    original = LRAA_Globals.config.get("spare_polyA_veto_at_known_3prime")
    yield
    LRAA_Globals.config["spare_polyA_veto_at_known_3prime"] = original


def test_without_a_reference_the_veto_still_fires():
    """Ref-free behaviour must be untouched: no reference, no exemption."""
    assert _polyA_sites(A_RICH_AT_1000, {1000: 40}, guides=None) == []


def test_a_reference_3prime_end_at_the_same_position_waives_the_veto():
    guides = [_guide("known.1", [[500, 1000]], "+")]
    assert _polyA_sites(A_RICH_AT_1000, {1000: 40}, guides=guides) == [1000]


def test_a_reference_end_elsewhere_does_not_waive_it():
    guides = [_guide("known.far", [[500, 2500]], "+")]
    assert _polyA_sites(A_RICH_AT_1000, {1000: 40}, guides=guides) == []


def test_the_exemption_can_be_switched_off():
    guides = [_guide("known.1", [[500, 1000]], "+")]
    assert _polyA_sites(A_RICH_AT_1000, {1000: 40}, guides=guides, spare=False) == []


def test_only_same_strand_reference_ends_count():
    """A minus-strand annotation end must not waive a plus-strand candidate's veto."""
    guides = [_guide("known.minus", [[1000, 1500]], "-")]
    assert _polyA_sites(A_RICH_AT_1000, {1000: 40}, guides=guides) == []


def test_the_tolerance_is_the_alt_polyA_window_half():
    """+/-25 (max_dist_between_alt_polyA_sites / 2), inclusive, on both sides.

    Not an exact coordinate match: `position` reaches the veto from
    aggregate_sites_within_window, which names the single most-supported read end in a
    50 nt window, so it is not base-precise the way the annotation is.
    """
    tolerance = int(LRAA_Globals.config["max_dist_between_alt_polyA_sites"] / 2)
    for offset in (tolerance, -tolerance):
        inside = [_guide("k", [[500, 1000 + offset]], "+")]
        assert _polyA_sites(A_RICH_AT_1000, {1000: 40}, guides=inside) == [1000], offset
    for offset in (tolerance + 1, -(tolerance + 1)):
        outside = [_guide("k", [[500, 1000 + offset]], "+")]
        assert _polyA_sites(A_RICH_AT_1000, {1000: 40}, guides=outside) == [], offset


def test_the_tolerance_tracks_the_configured_window(monkeypatch):
    """The window is derived, not a literal, so the two cannot drift apart."""
    monkeypatch.setitem(LRAA_Globals.config, "max_dist_between_alt_polyA_sites", 4)
    assert _polyA_sites(A_RICH_AT_1000, {1000: 40},
                        guides=[_guide("k", [[500, 1002]], "+")]) == [1000]
    assert _polyA_sites(A_RICH_AT_1000, {1000: 40},
                        guides=[_guide("k", [[500, 1003]], "+")]) == []


def test_a_clean_context_candidate_is_unaffected_either_way():
    """The exemption must only ever loosen the veto, never tighten anything."""
    clean = "C" * 3000
    assert _polyA_sites(clean, {1000: 40}, guides=None) == [1000]
    guides = [_guide("k", [[500, 1000]], "+")]
    assert _polyA_sites(clean, {1000: 40}, guides=guides) == [1000]


def test_reference_ends_are_collected_on_the_transcript_three_prime_side():
    """On '-' a transcript's 3' end is its lend, not its rend."""
    sg = Splice_graph.Splice_graph()
    LRAA_Globals.config["spare_polyA_veto_at_known_3prime"] = True
    minus = [_guide("k", [[400, 900]], "-")]
    assert sg._collect_reference_three_prime_ends(minus, "-") == [400]
    plus = [_guide("k", [[400, 900]], "+")]
    assert sg._collect_reference_three_prime_ends(plus, "+") == [900]


class _Unusable:
    """A reference transcript whose 3' end cannot be read."""

    def get_transcript_id(self):
        return "broken.1"

    def get_strand(self):
        return "+"

    def get_coords(self):
        raise ValueError("no exon segments")


def test_an_unusable_reference_transcript_is_reported_not_dropped_silently(caplog):
    """Losing one shrinks the exemption: sites the annotation endorses get vetoed as
    internal priming instead, and that is indistinguishable from the annotation not
    calling an end there unless it is said out loud."""
    sg = Splice_graph.Splice_graph()
    LRAA_Globals.config["spare_polyA_veto_at_known_3prime"] = True
    with caplog.at_level(logging.WARNING, logger="Splice_graph"):
        ends = sg._collect_reference_three_prime_ends(
            [_guide("good", [[400, 900]], "+"), _Unusable()], "+"
        )
    assert ends == [900], "the usable transcript must still contribute its 3' end"
    warnings = [r for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1
    message = warnings[0].getMessage()
    assert "1 of 2 reference transcripts" in message
    assert "broken.1" in message
    assert "no exon segments" in message


def test_a_wrong_strand_reference_transcript_is_not_reported_as_unusable(caplog):
    """Nothing is wrong with it; it just has nothing to say about this strand."""
    sg = Splice_graph.Splice_graph()
    LRAA_Globals.config["spare_polyA_veto_at_known_3prime"] = True
    with caplog.at_level(logging.WARNING, logger="Splice_graph"):
        ends = sg._collect_reference_three_prime_ends(
            [_guide("minus", [[400, 900]], "-")], "+"
        )
    assert ends == []
    assert [r for r in caplog.records if r.levelno == logging.WARNING] == []
