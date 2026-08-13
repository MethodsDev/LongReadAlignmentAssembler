#!/usr/bin/env python3

"""Internal-priming coverage, across both stages that apply the rule.

`Splice_graph._incorporate_PolyA_objects` rejects an A-rich READ-DERIVED candidate
before it becomes a vertex. `TranscriptFiltering.filter_internally_primed_transcripts`
then annotates every emitted model and deletes by policy: monoexonic always,
multi-exonic only when `restrict_internal_priming_filter_to_monoexonic` is False.

Both stages share one predicate in `Util_funcs`, and a test here pins that they cannot
drift apart. The gate tests cover the property everything rests on: a read-derived
candidate is vetoed while a guide-derived terminus at the same coordinate is not.
"""

import LRAA_Globals
import Splice_graph
import TranscriptFiltering
import Util_funcs
from Transcript import Transcript


def _build_transcript(transcript_id, exons, strand="+"):
    transcript = Transcript("chr1", exons, strand)
    transcript.set_gene_id("g1")
    transcript.set_transcript_id(transcript_id)
    return transcript


def _filter(
    transcripts,
    contig_seq_str,
    contig_strand="+",
    known_transcripts=None,
    restrict_to_monoexonic=True,
    spare_monoexonic_with_known_3prime=False,
    delete=True,
):
    return TranscriptFiltering.filter_internally_primed_transcripts(
        transcripts,
        contig_seq_str,
        contig_strand,
        known_transcripts,
        restrict_filter_to_monoexonic=restrict_to_monoexonic,
        spare_monoexonic_with_known_3prime=spare_monoexonic_with_known_3prime,
        delete=delete,
    )


def _ids(transcripts):
    return [t.get_transcript_id() for t in transcripts]


# ---------------------------------------------------------------- the shared predicate

def test_internal_priming_ignores_21st_downstream_base_on_forward_strand():
    downstream_twenty = "ACACACAC" + "AAAAAAACCCCC"
    contig_seq = "C" * 10 + downstream_twenty + "A"
    assert Util_funcs.looks_internally_primed(contig_seq, 10, "+") is False


def test_internal_priming_ignores_21st_downstream_base_on_reverse_strand():
    downstream_twenty = "CCCCCTTTTTTT" + "TGTGTGTG"
    contig_seq = "T" + downstream_twenty + "C" * 10
    assert Util_funcs.looks_internally_primed(contig_seq, 22, "-") is False


def test_predicate_fires_on_the_count_arm_and_the_hexamer_arm():
    # 12 of 20 downstream A's, with no hexamer inside the first 8
    count_arm = "C" * 10 + ("AC" * 4) + "A" * 12
    assert count_arm[10:30].count("A") >= 12
    assert "AAAAAA" not in count_arm[10:18]
    assert Util_funcs.looks_internally_primed(count_arm, 10, "+") is True
    # a hexamer inside the first 8, but well under 12 A's overall
    hexamer_arm = "C" * 10 + "AAAAAA" + "C" * 14
    assert hexamer_arm[10:30].count("A") < 12
    assert Util_funcs.looks_internally_primed(hexamer_arm, 10, "+") is True
    # and neither
    assert Util_funcs.looks_internally_primed("C" * 40, 10, "+") is False


def test_transcript_stage_uses_the_same_predicate_as_the_gate():
    """The two stages must never drift apart; this pins them to one definition."""
    contig_seq = "C" * 50 + "A" * 12 + "C" * 38
    primed = _build_transcript("t.mono", [[10, 50]], "+")
    _filter([primed], contig_seq, delete=False)
    assert primed._likely_internal_primed is Util_funcs.looks_internally_primed(
        contig_seq, 50, "+"
    )
    assert primed._likely_internal_primed is True


# ------------------------------------------------------------- deletion policy

PRIMED_CONTIG = "C" * 50 + "A" * 12 + "C" * 38


def test_monoexonic_internal_priming_is_deleted():
    monoexonic = _build_transcript("t.mono", [[10, 50]], "+")
    assert _ids(_filter([monoexonic], PRIMED_CONTIG)) == []


def test_multiexonic_internal_priming_is_retained_by_default():
    """Default policy restricts deletion to monoexonic models; a flagged multi-exonic
    model is annotated and kept, its intron chain being corroborating evidence."""
    multiexonic = _build_transcript("t.multi", [[10, 20], [40, 50]], "+")
    retained = _filter([multiexonic], PRIMED_CONTIG)
    assert _ids(retained) == ["t.multi"]
    assert retained[0]._likely_internal_primed is True


def test_multiexonic_internal_priming_is_deleted_when_not_restricted():
    multiexonic = _build_transcript("t.multi", [[10, 20], [40, 50]], "+")
    assert _ids(_filter([multiexonic], PRIMED_CONTIG, restrict_to_monoexonic=False)) == []


def test_multiexonic_clean_context_is_retained_and_annotated_false():
    multiexonic = _build_transcript("t.multi", [[10, 20], [40, 50]], "+")
    retained = _filter([multiexonic], "C" * 100, restrict_to_monoexonic=False)
    assert _ids(retained) == ["t.multi"]
    assert retained[0]._likely_internal_primed is False


def test_multiexonic_clean_context_is_annotated_false_under_either_policy():
    """Restricting deletion to monoexonic models must not make the annotation lie about
    a multi-exonic model that was never in A-rich context to begin with."""
    multiexonic = _build_transcript("t.multi", [[10, 20], [40, 50]], "+")
    retained = _filter([multiexonic], "C" * 100, restrict_to_monoexonic=True)
    assert _ids(retained) == ["t.multi"]
    assert retained[0]._likely_internal_primed is False


def test_reverse_strand_transcript_is_judged_from_its_own_three_prime_end():
    # on '-' the 3' end is the transcript's lend, and the window runs 5'-ward of it
    monoexonic = _build_transcript("t.rev", [[51, 90]], "-")
    contig_seq = "C" * 38 + "T" * 12 + "C" * 50
    assert _ids(_filter([monoexonic], contig_seq, contig_strand="-")) == []


def test_only_the_flagged_model_is_deleted():
    models = [
        _build_transcript("t.mono.primed", [[10, 50]], "+"),
        _build_transcript("t.mono.clean", [[10, 20]], "+"),
    ]
    retained = _filter(models, PRIMED_CONTIG)
    assert _ids(retained) == ["t.mono.clean"]


# ------------------------------------------------- annotation survives without deletion


def test_annotation_is_applied_even_when_deletion_is_disabled():
    """`--no_filter_internal_priming`. The attribute is load-bearing beyond logging --
    Transcript.py reads InternalPriming back off an input GTF to re-export it -- so it
    must be set whether or not anything is removed."""
    models = [
        _build_transcript("t.mono", [[10, 50]], "+"),
        _build_transcript("t.multi", [[10, 20], [40, 50]], "+"),
    ]
    retained = _filter(models, PRIMED_CONTIG, delete=False)
    assert _ids(retained) == ["t.mono", "t.multi"]
    assert all(t._likely_internal_primed for t in retained)


# ------------------------------------------- the known-3'-end reprieve (not an atlas)


def _spare(model_end, known_end):
    """Does a monoexonic model ending at `model_end` survive, given one known
    transcript ending at `known_end`?"""
    monoexonic = _build_transcript("t.mono", [[10, model_end]], "+")
    known = [_build_transcript("known.1", [[5, known_end]], "+")]
    return _ids(
        _filter(
            [monoexonic],
            PRIMED_CONTIG,
            known_transcripts=known,
            spare_monoexonic_with_known_3prime=True,
        )
    ) == ["t.mono"]


def test_monoexonic_spared_by_a_coincident_known_transcript_end():
    """The reprieve keys on proximity to a `known_transcripts` 3' end -- the reference
    annotation handed to the filter -- and not on any measured cleavage atlas."""
    assert _spare(50, 50) is True


def test_the_reprieve_still_annotates_what_it_spares():
    monoexonic = _build_transcript("t.mono", [[10, 50]], "+")
    known = [_build_transcript("known.1", [[5, 50]], "+")]
    retained = _filter(
        [monoexonic],
        PRIMED_CONTIG,
        known_transcripts=known,
        spare_monoexonic_with_known_3prime=True,
    )
    assert _ids(retained) == ["t.mono"]
    assert retained[0]._likely_internal_primed is True


def test_the_reprieve_window_extends_beyond_a_coincident_end():
    """Proximity, not equality -- otherwise the coincident-end test above would pass on
    an implementation that only ever compared for equality.

    Deliberately probed well inside the window. The exact edge is NOT asserted here:
    intervals go in as [K - half, K + half) but are probed at [E, E + 1), which makes
    the window asymmetric by one base ([E - half + 1, E + half] for a half of 25). That
    looks like an off-by-one rather than a decision, it predates this change, and it is
    unreachable at shipped defaults, so it is reported rather than frozen into a test.
    """
    half = LRAA_Globals.config["max_dist_between_alt_polyA_sites"] // 2
    E = 50
    # not vacuous: with no known end nearby this model is deleted
    assert _ids(_filter([_build_transcript("t.mono", [[10, E]], "+")], PRIMED_CONTIG)) == []
    assert _spare(E, E + half // 2) is True
    assert _spare(E, E - half // 2) is True
    assert _spare(E, E + 4 * half) is False


def test_monoexonic_not_spared_when_no_known_end_is_nearby():
    assert _spare(50, 5000) is False


def test_monoexonic_not_spared_when_the_option_is_off():
    monoexonic = _build_transcript("t.mono", [[10, 50]], "+")
    known = [_build_transcript("known.1", [[5, 50]], "+")]
    assert _ids(_filter([monoexonic], PRIMED_CONTIG, known_transcripts=known)) == []


def test_multiexonic_spared_by_a_known_end_when_deletion_is_unrestricted():
    multiexonic = _build_transcript("t.multi", [[10, 20], [40, 50]], "+")
    known = [_build_transcript("known.1", [[5, 50]], "+")]
    retained = _filter(
        [multiexonic],
        PRIMED_CONTIG,
        known_transcripts=known,
        restrict_to_monoexonic=False,
    )
    assert _ids(retained) == ["t.multi"]


# ------------------------------------------------------------------------------ the gate

def _polyA_sites(contig_seq_str, counter, strand="+", from_reads=True):
    sg = Splice_graph.Splice_graph()
    sg._contig_seq_str = contig_seq_str
    sg._incorporate_PolyA_objects("chr1", strand, counter, from_reads=from_reads)
    return sorted(obj.get_coords()[0] for obj in sg._PolyA_objs)


def test_read_derived_A_rich_candidate_is_vetoed_at_identification():
    # 1000 is followed by 20 A's; 2000 is not. min_alignments_define_polyA_site is 5.
    contig_seq = "C" * 1000 + "A" * 20 + "C" * 980 + "C" * 20 + "C" * 980
    kept = _polyA_sites(contig_seq, {1000: 40, 2000: 40})
    assert 1000 not in kept, "A-rich read-derived candidate should be rejected"
    assert 2000 in kept, "clean-context candidate should survive"


def test_guide_derived_terminus_is_never_vetoed_at_the_same_coordinate():
    """The property everything rests on. A reference annotation's 3' end is an assertion
    about a transcript, not an inference from a soft clip, so it must survive the gate
    even where a read-derived candidate at the identical coordinate is rejected."""
    contig_seq = "C" * 1000 + "A" * 20 + "C" * 1980

    from_reads = _polyA_sites(contig_seq, {1000: 40}, from_reads=True)
    from_guide = _polyA_sites(contig_seq, {1000: 40}, from_reads=False)

    assert from_reads == [], "read-derived A-rich candidate must be rejected"
    assert from_guide == [1000], "guide-derived terminus must be retained"


def test_gate_is_strand_aware():
    # on '-' the window runs 5'-ward, so the A-run must be genomic T upstream of the site
    contig_seq = "C" * 980 + "T" * 20 + "C" * 2000
    kept_minus = _polyA_sites(contig_seq, {1001: 40}, strand="-")
    assert kept_minus == [], "A-rich in transcript sense on '-' should be rejected"
    # the same coordinate on '+' looks downstream instead, which is clean here
    kept_plus = _polyA_sites(contig_seq, {1001: 40}, strand="+")
    assert kept_plus == [1001]
