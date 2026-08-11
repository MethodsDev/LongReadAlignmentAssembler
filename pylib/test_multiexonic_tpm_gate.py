#!/usr/bin/env python3

"""The post-EM gate that decides whether a selected multi-exonic model is reported.

Supplied models reach path selection on a synthetic template read, so selection can no
longer be read as "something supports this". Quantification can, and this gate is where
that judgement is made: a model EM gave no expression to is dropped, whatever put it in
the candidate set.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import TranscriptFiltering
from Transcript import Transcript


def _tx(exon_segments, tpm, tid):
    t = Transcript("chr1", exon_segments, "+")
    t.set_transcript_id(tid)
    # the field get_TPM() honours as an override; Transcript.py sets it the same way
    # when reading a TPM off an input GTF
    t._imported_TPM_val = tpm
    return t


def _multiexonic(tpm, tid):
    return _tx([[100, 200], [300, 400]], tpm, tid)


def _monoexonic(tpm, tid):
    return _tx([[100, 200]], tpm, tid)


def _ids(transcripts):
    return [t.get_transcript_id() for t in transcripts]


def test_multiexonic_tpm_gate_drops_unexpressed_and_keeps_expressed():
    """The default threshold of 0 means 'any expression at all'."""
    kept = _multiexonic(0.001, "expressed")
    dropped = _multiexonic(0.0, "unexpressed")

    retained = TranscriptFiltering.filter_multiexonic_isoforms_by_TPM_threshold(
        [kept, dropped], 0
    )

    assert _ids(retained) == ["expressed"]


def test_multiexonic_tpm_gate_leaves_monoexonic_models_alone():
    """Monoexonic models have their own thresholds; this gate must not double-filter."""
    mono_zero = _monoexonic(0.0, "mono_zero")
    multi_zero = _multiexonic(0.0, "multi_zero")

    retained = TranscriptFiltering.filter_multiexonic_isoforms_by_TPM_threshold(
        [mono_zero, multi_zero], 0
    )

    assert _ids(retained) == ["mono_zero"]


def test_multiexonic_tpm_gate_honours_a_raised_threshold():
    """A trace of expression can be required to be more than a trace."""
    trace = _multiexonic(0.4, "trace")
    real = _multiexonic(2.5, "real")

    retained = TranscriptFiltering.filter_multiexonic_isoforms_by_TPM_threshold(
        [trace, real], 1.0
    )

    assert _ids(retained) == ["real"]


def test_multiexonic_tpm_gate_is_disabled_by_a_negative_threshold():
    """Negative disables the check, reporting every selected multi-exonic model."""
    unexpressed = _multiexonic(0.0, "unexpressed")

    retained = TranscriptFiltering.filter_multiexonic_isoforms_by_TPM_threshold(
        [unexpressed], -1
    )

    assert _ids(retained) == ["unexpressed"]


def test_multiexonic_tpm_gate_is_strictly_greater_than():
    """A model exactly at the threshold does not clear it.

    This matters at the default: TPM == 0 is the unexpressed case the gate exists to
    remove, so the comparison cannot be >=.
    """
    at_threshold = _multiexonic(1.0, "at")
    above = _multiexonic(1.0001, "above")

    retained = TranscriptFiltering.filter_multiexonic_isoforms_by_TPM_threshold(
        [at_threshold, above], 1.0
    )

    assert _ids(retained) == ["above"]
