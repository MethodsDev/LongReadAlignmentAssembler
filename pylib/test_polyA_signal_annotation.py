#!/usr/bin/env python3

"""The canonical PAS upstream of a 3' terminus, recorded as a transcript attribute.

Independent of internal priming by design: that predicate reads the genome DOWNSTREAM
of the terminus looking for an oligo-dT template, this one reads UPSTREAM looking for
the signal a genuine cleavage site requires. A terminus can carry both, so the two are
separate attributes rather than one verdict.
"""

import pytest

import TranscriptFiltering
import Util_funcs
from Transcript import Transcript


def _build_transcript(transcript_id, exons, strand="+"):
    transcript = Transcript("chr1", exons, strand)
    transcript.set_gene_id("g1")
    transcript.set_transcript_id(transcript_id)
    # to_GTF_format() asserts a simple path exists; these tests are about attributes,
    # not graph structure, so a placeholder is enough
    transcript.set_simple_path(["n1"])
    return transcript


def _contig_with_signal(site, offset, hexamer="AATAAA", strand="+", length=400):
    """Place `hexamer` so its transcript-sense first base sits `offset` nt from `site`."""
    seq = ["C"] * length
    if strand == "+":
        start = site + offset  # 1-based genomic position of the first base
        for i, base in enumerate(hexamer):
            seq[start - 1 + i] = base
    else:
        # on '-' the transcript-sense first base is the RIGHTMOST genomic base
        first_base_pos = site - offset
        genomic_start = first_base_pos - (len(hexamer) - 1)
        spelling = "".join(
            reversed([{"A": "T", "T": "A", "C": "G", "G": "C"}[b] for b in hexamer])
        )
        for i, base in enumerate(spelling):
            seq[genomic_start - 1 + i] = base
    return "".join(seq)


# ------------------------------------------------------------------- the predicate


def test_signal_found_upstream_on_the_forward_strand():
    contig = _contig_with_signal(100, -21)
    assert Util_funcs.find_polyA_signal(contig, 100, "+") == ("AATAAA", -21)


def test_signal_found_upstream_on_the_reverse_strand():
    """Reported in transcript sense: the genome spells TTTATT, the attribute says
    AATAAA. Upstream on '-' means HIGHER genomic coordinates."""
    contig = _contig_with_signal(100, -21, strand="-")
    assert contig[115:121] == "TTTATT"
    assert Util_funcs.find_polyA_signal(contig, 100, "-") == ("AATAAA", -21)


@pytest.mark.parametrize(
    "offset,expected",
    [(-40, "AATAAA"), (-41, None), (-15, "AATAAA"), (-14, None), (-10, None)],
)
def test_only_fully_contained_motifs_count(offset, expected):
    """The window bounds the extracted region and the motif must fit entirely inside it,
    so the default -40..-10 admits first-base offsets of -40..-15.

    Worth stating plainly because it is NOT the same set as "starts 10 to 30 nt
    upstream": a 6-mer starting at -10 would run through -5, outside the window it is
    supposed to occupy.
    """
    contig = _contig_with_signal(100, offset)
    hexamer, _ = Util_funcs.find_polyA_signal(contig, 100, "+")
    assert hexamer == expected


def test_the_closest_signal_to_the_site_wins():
    seq = list(_contig_with_signal(100, -41))  # outside, must not be chosen
    for i, base in enumerate("AATAAA"):
        seq[100 - 21 - 1 + i] = base
    for i, base in enumerate("ATTAAA"):
        seq[100 - 15 - 1 + i] = base
    assert Util_funcs.find_polyA_signal("".join(seq), 100, "+") == ("ATTAAA", -15)


def test_both_canonical_hexamers_are_recognised_and_variants_are_not():
    assert Util_funcs.find_polyA_signal(_contig_with_signal(100, -20, "AATAAA"), 100, "+")[0] == "AATAAA"
    assert Util_funcs.find_polyA_signal(_contig_with_signal(100, -20, "ATTAAA"), 100, "+")[0] == "ATTAAA"
    # a single-substitution variant is deliberately NOT a canonical call
    assert Util_funcs.find_polyA_signal(_contig_with_signal(100, -20, "AATACA"), 100, "+") == (None, None)


def test_no_signal_returns_a_pair_of_nones():
    assert Util_funcs.find_polyA_signal("C" * 300, 100, "+") == (None, None)


def test_a_site_near_the_contig_edge_does_not_raise():
    assert Util_funcs.find_polyA_signal("AATAAACCCC", 3, "+") == (None, None)
    assert Util_funcs.find_polyA_signal("CCCC", 2, "-") == (None, None)


def test_rejects_an_unknown_strand():
    with pytest.raises(ValueError):
        Util_funcs.find_polyA_signal("C" * 100, 50, "?")


# ------------------------------------------------------------------ the annotation


def test_annotation_sets_the_attribute_and_deletes_nothing():
    contig = _contig_with_signal(50, -20)
    models = [
        _build_transcript("t.with", [[10, 50]], "+"),
        _build_transcript("t.without", [[10, 300]], "+"),
    ]
    retained = TranscriptFiltering.annotate_polyA_signal(models, contig, "+")
    assert [t.get_transcript_id() for t in retained] == ["t.with", "t.without"]
    assert retained[0]._polyA_signal == "AATAAA"
    assert retained[0]._polyA_signal_offset == -20
    assert retained[1]._polyA_signal is None


def test_absence_is_recorded_explicitly_in_the_gtf():
    """`PAS "none"` rather than an omitted key: a missing attribute cannot be told
    apart from a GTF written before this existed."""
    model = _build_transcript("t.none", [[10, 300]], "+")
    TranscriptFiltering.annotate_polyA_signal([model], "C" * 400, "+")
    gtf = model.to_GTF_format()
    assert 'PAS "none"' in gtf
    assert 'PAS_offset "' not in gtf


def test_a_found_signal_is_emitted_with_its_offset():
    contig = _contig_with_signal(50, -20)
    model = _build_transcript("t.with", [[10, 50]], "+")
    TranscriptFiltering.annotate_polyA_signal([model], contig, "+")
    gtf = model.to_GTF_format()
    assert 'PAS "AATAAA"' in gtf
    assert 'PAS_offset "-20"' in gtf


def test_the_attribute_is_absent_when_the_pass_never_ran():
    """No annotation must mean no key, so a GTF from a run without this pass is not
    misread as evidence of absence."""
    model = _build_transcript("t.unevaluated", [[10, 50]], "+")
    # `PAS "` and not bare "PAS": the GTF carries the tool's legacy PASA_SALRAA naming
    assert 'PAS "' not in model.to_GTF_format()


def test_pas_and_internal_priming_are_independent_on_the_same_terminus():
    """The case that motivates two attributes: a canonical signal upstream AND an
    A-rich stretch downstream. Both must be reported."""
    seq = ["C"] * 400
    for i, base in enumerate("AATAAA"):
        seq[50 - 20 - 1 + i] = base
    for i in range(50, 70):
        seq[i] = "A"
    contig = "".join(seq)

    model = _build_transcript("t.both", [[10, 50]], "+")
    TranscriptFiltering.annotate_polyA_signal([model], contig, "+")
    model.set_likely_internal_primed(
        TranscriptFiltering._looks_internally_primed(10, 50, "+", contig)
    )

    gtf = model.to_GTF_format()
    assert 'PAS "AATAAA"' in gtf
    assert 'InternalPriming "True"' in gtf


def test_pas_survives_a_gtf_round_trip(tmp_path):
    """Transcript.py re-exports imported attributes, so a value written by an assembly
    run must survive a later quant-only pass over the same GTF."""
    from Transcript import GTF_contig_to_transcripts

    gtf_path = tmp_path / "pas.gtf"
    gtf_path.write_text(
        'chr1\tLRAA\ttranscript\t100\t210\t.\t+\t.\t'
        'gene_id "g1"; transcript_id "t1"; PAS "AATAAA"; PAS_offset "-21";\n'
        'chr1\tLRAA\texon\t100\t210\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    )
    by_contig = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(str(gtf_path))
    transcript = by_contig["chr1"][0]

    assert transcript._polyA_signal == "AATAAA"
    assert transcript._polyA_signal_offset == -21

    reexported = transcript.to_GTF_format()
    assert 'PAS "AATAAA"' in reexported
    assert 'PAS_offset "-21"' in reexported


def _reference_implementation(contig_seq_str, three_prime_pos, strand, lo=-40, hi=-10):
    """How investigations/PIP_work/pip_profile.py does it: extract the transcript-sense
    window, then substring-search for either hexamer. Deliberately a separate, simpler
    formulation -- if the optimised scan in Util_funcs agrees with this over random
    sequence, the conventions are the same."""
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    n = len(contig_seq_str)
    if strand == "+":
        start, end = three_prime_pos + lo, three_prime_pos + hi
    else:
        start, end = three_prime_pos - hi, three_prime_pos - lo
    if start < 1 or end > n:
        return False
    region = contig_seq_str[start - 1 : end].upper()
    if strand == "-":
        region = "".join(reversed([complement[b] for b in region]))
    return any(hexamer in region for hexamer in ("AATAAA", "ATTAAA"))


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("window", [(-40, -10), (-30, -10), (-25, -12)])
def test_agrees_with_extract_then_substring_over_random_sequence(strand, window):
    """Equivalence with the straightforward formulation, at the default window and at
    the -40..-10 the offline analysis (investigations/PIP_work/pip_profile.py) used, so
    numbers from either convention remain reproducible from this helper."""
    import random

    rng = random.Random(20260813)
    site = 500
    disagreements = []
    for _ in range(2000):
        contig = "".join(rng.choice("AATCG") for _ in range(1000))
        mine = Util_funcs.find_polyA_signal(contig, site, strand, window=window)[0] is not None
        theirs = _reference_implementation(contig, site, strand, *window)
        if mine != theirs:
            disagreements.append((mine, theirs, contig[site - 60 : site + 10]))
    assert not disagreements, f"{len(disagreements)} disagreements, first: {disagreements[0]}"


# ------------------------------------------------------- the motif set is a real knob


@pytest.mark.parametrize("strand", ["+", "-"])
def test_a_custom_motif_set_works_on_both_strands(strand):
    """`hexamers` is documented as adjustable, so a non-canonical set must actually
    work -- including its minus-strand spelling, which is derived rather than looked up
    from a canonical table."""
    motif = "TGTGTG"
    contig = _contig_with_signal(200, -20, hexamer=motif, strand=strand, length=600)
    assert Util_funcs.find_polyA_signal(
        contig, 200, strand, hexamers=(motif,)
    ) == (motif, -20)
    # and the canonical set must NOT find it
    assert Util_funcs.find_polyA_signal(contig, 200, strand)[0] is None


def test_a_longer_motif_respects_containment():
    """Containment is derived from the motif width, so a 7-mer in -30..-10 has valid
    first-base offsets of -30..-16, one narrower than a 6-mer."""
    motif = "AATAAAG"
    inside = _contig_with_signal(200, -16, hexamer=motif, length=600)
    outside = _contig_with_signal(200, -15, hexamer=motif, length=600)
    assert Util_funcs.find_polyA_signal(inside, 200, "+", hexamers=(motif,)) == (motif, -16)
    assert Util_funcs.find_polyA_signal(outside, 200, "+", hexamers=(motif,)) == (None, None)


@pytest.mark.parametrize(
    "motifs,reason",
    [
        ((), "empty"),
        (("AATAAA", "ATTAA"), "mixed widths"),
        (("aataaa",), "lowercase"),
        (("AATAAN",), "non-ACGT"),
    ],
)
def test_an_unusable_motif_set_is_rejected(motifs, reason):
    """Mixed widths in particular must not be accepted: the scan derives its
    containment bound from a single width, so a shorter member would be silently
    mis-placed rather than mis-matched."""
    with pytest.raises(ValueError):
        Util_funcs.find_polyA_signal("C" * 400, 200, "+", hexamers=motifs)


# ------------------------------------------------- configurable for other organisms


def test_defaults_are_the_human_canonical_settings():
    motifs, window = Util_funcs.resolve_polyA_signal_settings(
        {"polyA_signal_motifs": ["AATAAA", "ATTAAA"], "polyA_signal_window": [-40, -10]}
    )
    assert motifs == ("AATAAA", "ATTAAA")
    assert window == (-40, -10)


def test_config_accepts_a_comma_separated_string_as_the_cli_supplies_it():
    motifs, window = Util_funcs.resolve_polyA_signal_settings(
        {"polyA_signal_motifs": "AATAAA, ATTAAA ,AATGAA", "polyA_signal_window": "-35,-12"}
    )
    assert motifs == ("AATAAA", "ATTAAA", "AATGAA")
    assert window == (-35, -12)


def test_motifs_are_upcased_so_a_lowercase_config_still_works():
    motifs, _ = Util_funcs.resolve_polyA_signal_settings(
        {"polyA_signal_motifs": ["aataaa"], "polyA_signal_window": [-40, -10]}
    )
    assert motifs == ("AATAAA",)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_a_non_human_configuration_annotates_its_own_motif(strand):
    """The point of making these settable: an organism whose signal is neither
    AATAAA nor ATTAAA must be annotatable without touching the code."""
    motif = "AATGAA"
    contig = _contig_with_signal(200, -22, hexamer=motif, strand=strand, length=600)
    motifs, window = Util_funcs.resolve_polyA_signal_settings(
        {"polyA_signal_motifs": [motif], "polyA_signal_window": [-40, -10]}
    )
    assert Util_funcs.find_polyA_signal(
        contig, 200, strand, window=window, hexamers=motifs
    ) == (motif, -22)
    # the human default must not find it
    assert Util_funcs.find_polyA_signal(contig, 200, strand)[0] is None


def test_a_wider_organism_window_is_honoured():
    contig = _contig_with_signal(200, -55, length=600)
    assert Util_funcs.find_polyA_signal(contig, 200, "+")[0] is None
    motifs, window = Util_funcs.resolve_polyA_signal_settings(
        {"polyA_signal_motifs": ["AATAAA"], "polyA_signal_window": [-60, -10]}
    )
    assert Util_funcs.find_polyA_signal(
        contig, 200, "+", window=window, hexamers=motifs
    ) == ("AATAAA", -55)


@pytest.mark.parametrize(
    "config,fragment",
    [
        ({"polyA_signal_motifs": []}, "at least one motif"),
        ({"polyA_signal_motifs": ["AATAAA", "ATTAA"]}, "share one length"),
        ({"polyA_signal_motifs": ["AATAAN"]}, "only A, C, G or T"),
        ({"polyA_signal_window": [-40]}, "exactly two integers"),
        ({"polyA_signal_window": ["a", "b"]}, "two integers"),
        ({"polyA_signal_window": [-10, -40]}, "low <= high"),
        ({"polyA_signal_window": [-13, -10]}, "too narrow"),
    ],
)
def test_a_bad_configuration_is_rejected_with_a_message_naming_the_key(config, fragment):
    """Startup validation: LRAA resolves these before any contig runs, so an organism
    misconfiguration must fail immediately and say which key is wrong."""
    full = {"polyA_signal_motifs": ["AATAAA"], "polyA_signal_window": [-40, -10]}
    full.update(config)
    with pytest.raises(ValueError, match=fragment):
        Util_funcs.resolve_polyA_signal_settings(full)


def test_the_annotation_pass_honours_configured_settings(monkeypatch):
    """End of the wire: the pass must read config rather than the module defaults."""
    import LRAA_Globals

    motif = "AATGAA"
    contig = _contig_with_signal(50, -22, hexamer=motif, length=400)
    monkeypatch.setitem(LRAA_Globals.config, "polyA_signal_motifs", [motif])
    monkeypatch.setitem(LRAA_Globals.config, "polyA_signal_window", [-40, -10])

    model = _build_transcript("t.plant", [[10, 50]], "+")
    TranscriptFiltering.annotate_polyA_signal([model], contig, "+")
    gtf = model.to_GTF_format()
    assert f'PAS "{motif}"' in gtf
    assert 'PAS_offset "-22"' in gtf
