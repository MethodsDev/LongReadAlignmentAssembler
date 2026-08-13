#!/usr/bin/env python3

"""Internal-priming coverage.

The internal-priming DECISION now happens during PolyA site identification
(Splice_graph._incorporate_PolyA_objects); the transcript stage only annotates. So the
seven tests here that used to assert deletion now assert the annotation, over the same
fixtures, and two new tests cover the gate itself -- including the property everything
rests on, that a read-derived candidate is vetoed while a guide-derived terminus at the
same coordinate is not.
"""

import Splice_graph
import TranscriptFiltering
import Util_funcs
from Transcript import Transcript


def _build_transcript(transcript_id, exons, strand="+"):
    transcript = Transcript("chr1", exons, strand)
    transcript.set_gene_id("g1")
    transcript.set_transcript_id(transcript_id)
    return transcript


def _annotate(transcripts, contig_seq_str, contig_strand="+"):
    return TranscriptFiltering.annotate_internally_primed_transcripts(
        transcripts, contig_seq_str, contig_strand
    )


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
    contig_seq = "C" * 50 + "A" * 12 + "C" * 38
    primed = _build_transcript("t.mono", [[10, 50]], "+")
    _annotate([primed], contig_seq)
    assert primed._likely_internal_primed is Util_funcs.looks_internally_primed(
        contig_seq, 50, "+"
    )
    assert primed._likely_internal_primed is True


# ------------------------------------------------------- annotation, nothing is deleted

def test_multiexonic_clean_context_annotated_false_and_retained():
    multiexonic = _build_transcript("t.multi", [[10, 20], [40, 50]], "+")
    retained = _annotate([multiexonic], "C" * 100)
    assert [t.get_transcript_id() for t in retained] == ["t.multi"]
    assert retained[0]._likely_internal_primed is False


def test_multiexonic_internal_priming_annotated_true_and_retained():
    multiexonic = _build_transcript("t.multi", [[10, 20], [40, 50]], "+")
    retained = _annotate([multiexonic], "C" * 50 + "A" * 12 + "C" * 38)
    assert [t.get_transcript_id() for t in retained] == ["t.multi"]
    assert retained[0]._likely_internal_primed is True


def test_monoexonic_internal_priming_annotated_true_and_retained():
    # This used to be deleted. It is now annotated and kept: the decision is made at
    # PolyA site identification, and re-applying the rule to the emitted terminus would
    # delete models for which no spurious vertex was ever created.
    monoexonic = _build_transcript("t.mono", [[10, 50]], "+")
    retained = _annotate([monoexonic], "C" * 50 + "A" * 12 + "C" * 38)
    assert [t.get_transcript_id() for t in retained] == ["t.mono"]
    assert retained[0]._likely_internal_primed is True


def test_reverse_strand_transcript_annotated_from_its_own_three_prime_end():
    # on '-' the 3' end is the transcript's lend, and the window runs 5'-ward of it
    monoexonic = _build_transcript("t.rev", [[51, 90]], "-")
    retained = _annotate([monoexonic], "C" * 38 + "T" * 12 + "C" * 50, contig_strand="-")
    assert [t.get_transcript_id() for t in retained] == ["t.rev"]
    assert retained[0]._likely_internal_primed is True


def test_annotation_pass_deletes_nothing_even_when_every_model_is_flagged():
    contig_seq = "C" * 50 + "A" * 12 + "C" * 38
    models = [
        _build_transcript("t.mono", [[10, 50]], "+"),
        _build_transcript("t.multi", [[10, 20], [40, 50]], "+"),
    ]
    retained = _annotate(models, contig_seq)
    assert len(retained) == 2
    assert all(t._likely_internal_primed for t in retained)


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
