#!/usr/bin/env python3

"""What retain_expressed demands of a reference model before exempting it.

The reprieve used to ask `get_TPM() > 0`. get_TPM() is
read_counts_assigned / num_total_reads, so on any real library that is "did EM put
ANY mass here", and a 52M-read PBMC run admitted 14,193 multi-exon reference chains
whose entire assigned mass was under 0.05 of one read. write_expr prints all_reads to
one decimal, so they surfaced as "0.0".

The gate is now a COUNT of assigned reads, min_reads_retain_reference, default 1.0.
Pinned here because the failure it prevents is invisible in the output it produces:
too low and unsupported annotation is reported as recovered isoforms, and the only
symptom is a column of zeros.
"""

import pytest

import LRAA_Globals
import TranscriptFiltering
from Transcript import Transcript


def _reference_model(read_counts, num_exons=2):
    """A transcript that contains a reference model and has been quantified."""
    if num_exons == 1:
        exons = [[100, 400]]
    else:
        exons = [[100, 150], [200, 250]]
    transcript = Transcript("chr1", exons, "+")
    transcript.set_transcript_id("draft")
    transcript.set_source_reference_transcript_ids({"reference"})
    transcript.set_read_counts_assigned(read_counts)
    return transcript


@pytest.fixture(autouse=True)
def _retain_expressed(monkeypatch):
    monkeypatch.setitem(
        LRAA_Globals.config, "ref_trans_filter_mode", "retain_expressed"
    )
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 52_000_000)
    monkeypatch.setitem(LRAA_Globals.config, "min_reads_retain_reference", 1.0)


@pytest.mark.parametrize(
    ("read_counts", "reprieved"),
    [
        (0.0, False),
        (0.000579, False),  # the measured median of the PBMC zero-read population
        (0.041831, False),  # the measured maximum; still prints as all_reads 0.0
        (0.999999, False),
        (1.0, True),  # a whole read's worth of assignment, boundary is inclusive
        (7.5, True),
    ],
)
def test_reprieve_requires_a_whole_assigned_read(read_counts, reprieved):
    transcript = _reference_model(read_counts)
    assert TranscriptFiltering.reference_model_reprieved(transcript) is reprieved


def test_reprieve_does_not_require_unique_reads(monkeypatch):
    """Total assigned mass qualifies; uniqueness is not asked for.

    An isoform no single read can distinguish from its neighbours never accumulates a
    uniquely assigned read, and must still earn an output row. The unique-read
    thresholds that gate novel models are set punitively here: if the reprieve ever
    starts consulting a unique-read quantity, this fails.
    """
    monkeypatch.setitem(LRAA_Globals.config, "min_unique_reads_novel_isoform", 1000)
    monkeypatch.setitem(LRAA_Globals.config, "min_frac_gene_unique_reads", 1.0)

    transcript = _reference_model(1.0)
    assert TranscriptFiltering.reference_model_reprieved(transcript) is True


def test_non_reference_model_is_never_reprieved():
    transcript = Transcript("chr1", [[100, 150], [200, 250]], "+")
    transcript.set_transcript_id("novel")
    transcript.set_read_counts_assigned(1000.0)
    assert TranscriptFiltering.reference_model_reprieved(transcript) is False


def test_retain_filtered_mode_disables_the_reprieve(monkeypatch):
    monkeypatch.setitem(
        LRAA_Globals.config, "ref_trans_filter_mode", "retain_filtered"
    )
    transcript = _reference_model(1000.0)
    assert TranscriptFiltering.reference_model_reprieved(transcript) is False


def test_threshold_of_zero_restores_any_nonzero_mass_behaviour(monkeypatch):
    """The old semantics remain reachable, so a run can be reproduced."""
    monkeypatch.setitem(LRAA_Globals.config, "min_reads_retain_reference", 0.0)
    assert TranscriptFiltering.reference_model_reprieved(_reference_model(0.0)) is True
    assert (
        TranscriptFiltering.reference_model_reprieved(_reference_model(0.000579))
        is True
    )


def test_min_length_filter_honours_the_floor():
    """The reprieve is what keeps a short reference model; below the floor it goes.

    filter_transcripts_by_min_length is the site where the reprieve exists to stop a
    reference model's reads being discarded without an output row. A model carrying no
    read has no reads to discard.
    """
    supported = _reference_model(1.0)
    starved = _reference_model(0.02)

    retained = TranscriptFiltering.filter_transcripts_by_min_length(
        [supported, starved], min_transcript_length=10_000
    )

    assert retained == [supported]


def _novel_model(read_counts, monoexonic=False):
    exons = [[100, 400]] if monoexonic else [[100, 150], [200, 250]]
    transcript = Transcript("chr1", exons, "+")
    transcript.set_transcript_id("novel")
    transcript.set_read_counts_assigned(read_counts)
    return transcript


# -- the absolute floor: every model, no exemptions ---------------------------


@pytest.mark.parametrize(
    ("read_counts", "retained"),
    [(0.0, False), (0.041831, False), (0.999999, False), (1.0, True), (12.0, True)],
)
def test_the_floor_applies_to_novel_models(read_counts, retained):
    kept = TranscriptFiltering.filter_isoforms_by_min_assigned_reads(
        [_novel_model(read_counts)], 1.0
    )
    assert bool(kept) is retained


@pytest.mark.parametrize(
    ("read_counts", "retained"),
    [(0.0, False), (0.041831, False), (0.999999, False), (1.0, True)],
)
def test_the_floor_applies_to_reference_models_with_no_reprieve(read_counts, retained):
    """The reference reprieve does not reach here, and must not.

    Every other threshold asks something the annotation can answer -- is this long
    enough, is it a big enough share of its gene, is it in enough cells. "Did a read
    support this" is not such a question, so a reference model gets no exemption from
    it. This is the assertion that keeps the 14,193 out.
    """
    kept = TranscriptFiltering.filter_isoforms_by_min_assigned_reads(
        [_reference_model(read_counts)], 1.0
    )
    assert bool(kept) is retained


def test_the_floor_applies_to_monoexonic_models():
    """Monoexonic models are exempted from the multi-exonic TPM gate; not from this."""
    kept = TranscriptFiltering.filter_isoforms_by_min_assigned_reads(
        [_novel_model(0.5, monoexonic=True), _novel_model(2.0, monoexonic=True)], 1.0
    )
    assert [t.get_read_counts_assigned() for t in kept] == [2.0]


@pytest.mark.parametrize("disabled", [0, 0.0, None, -1])
def test_a_non_positive_floor_disables_the_filter(disabled):
    models = [_novel_model(0.0), _reference_model(0.0)]
    assert (
        TranscriptFiltering.filter_isoforms_by_min_assigned_reads(models, disabled)
        == models
    )


def test_the_floor_preserves_order():
    """Filtering runs weakest-first downstream, so order carries meaning."""
    a, b, c = _novel_model(5.0), _novel_model(0.1), _novel_model(3.0)
    kept = TranscriptFiltering.filter_isoforms_by_min_assigned_reads([a, b, c], 1.0)
    assert kept == [a, c]


# -- an input GTF's TPM attribute is not a read count -------------------------


def _model_with_imported_tpm(tpm, read_counts, reference=False):
    """A model parsed from a GTF carrying a TPM attribute, then quantified.

    Transcript._imported_TPM_val is set by the GTF parser (Transcript.py:1054) and
    get_read_counts_assigned() returns it in preference to the quantified count, so
    these two numbers disagree on purpose.
    """
    transcript = Transcript("chr1", [[100, 150], [200, 250]], "+")
    transcript.set_transcript_id("imported")
    if reference:
        transcript.set_source_reference_transcript_ids({"reference"})
    transcript._imported_TPM_val = tpm
    transcript.set_read_counts_assigned(read_counts)
    return transcript


def test_a_high_imported_tpm_does_not_satisfy_the_read_floor():
    """A TPM is a rate against a library this run never measured.

    get_read_counts_assigned() would answer 5000.0 here -- the GTF's TPM -- and let a
    model through a floor it never cleared. The floor asks get_assigned_read_count().
    """
    starved = _model_with_imported_tpm(tpm=5000.0, read_counts=0.02)
    assert starved.get_read_counts_assigned() == 5000.0  # the trap
    assert starved.get_assigned_read_count() == 0.02  # what the floor sees

    kept = TranscriptFiltering.filter_isoforms_by_min_assigned_reads([starved], 1.0)
    assert kept == []


def test_a_zero_imported_tpm_does_not_veto_a_supported_model():
    """The converse: a GTF asserting TPM 0 must not delete a model reads support."""
    supported = _model_with_imported_tpm(tpm=0.0, read_counts=42.0)
    assert supported.get_read_counts_assigned() == 0.0  # the trap
    assert supported.get_assigned_read_count() == 42.0

    kept = TranscriptFiltering.filter_isoforms_by_min_assigned_reads([supported], 1.0)
    assert kept == [supported]


def test_the_reference_reprieve_also_ignores_an_imported_tpm():
    """Same substitution, same defect, so the reprieve reads the same accessor."""
    starved = _model_with_imported_tpm(tpm=5000.0, read_counts=0.02, reference=True)
    supported = _model_with_imported_tpm(tpm=0.0, read_counts=42.0, reference=True)

    assert TranscriptFiltering.reference_model_reprieved(starved) is False
    assert TranscriptFiltering.reference_model_reprieved(supported) is True


def test_an_unquantified_model_has_no_assigned_reads():
    """A floor must act on a model quant never reached, not raise or guess.

    get_read_counts_assigned() asserts here; the floor cannot, because it is applied
    to every model and a model that failed to map onto the graph is exactly the kind
    the floor exists to remove.
    """
    transcript = Transcript("chr1", [[100, 150], [200, 250]], "+")
    transcript.set_transcript_id("never_quantified")

    assert transcript.get_assigned_read_count() == 0.0
    assert (
        TranscriptFiltering.filter_isoforms_by_min_assigned_reads([transcript], 1.0)
        == []
    )
