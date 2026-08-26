#!/usr/bin/env python3

"""`reference_transcripts` was added to public entry points, so it goes LAST.

It was originally inserted between `input_transcripts` and `quant_mode` in
Splice_graph.build_splice_graph_for_contig. Every in-tree caller passes `quant_mode` by
keyword, so nothing failed -- but any caller passing it positionally silently got
`quant_mode=False` AND had its transcript list adopted as the reference set, which is the
one thing the exemption must never do (see the circularity note in
build_splice_graph_for_contig). The failure is invisible: a discovery-mode build that
also waives the internal-priming veto on whatever it was handed.

These pin the argument ORDER, which is the part of a public signature a positional caller
depends on and no behavioural test in this tree exercises.
"""

import importlib.machinery
import importlib.util
import inspect
import os

import Splice_graph


def _load_lraa():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    loader = importlib.machinery.SourceFileLoader(
        "lraa_signature_test", os.path.join(here, "LRAA")
    )
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


lraa = _load_lraa()


def _last_parameter(func):
    return list(inspect.signature(func).parameters)[-1]


def test_a_positional_quant_mode_still_reaches_quant_mode():
    """The exhibit. Eight positional arguments, the eighth meaning quant_mode."""
    bound = inspect.signature(
        Splice_graph.Splice_graph.build_splice_graph_for_contig
    ).bind(None, "chr1", "+", "A" * 100, "reads.bam", None, None, ["guide"], True)
    bound.apply_defaults()
    assert bound.arguments["quant_mode"] is True
    assert bound.arguments["reference_transcripts"] is None


# The ORDERED tail of each public signature, ending at the last parameter.
# Pinned as a sequence, not a set: order is the thing a positional caller depends
# on, and only an ordered assertion catches BOTH failure directions --
# a parameter INSERTED before reference_transcripts (which shifts every
# positional argument onto the wrong name) and one APPENDED after it without a
# deliberate edit here.
_EXPECTED_TAIL = {
    "build_splice_graph_for_contig": (
        "restrict_splice_type",
        "SE_read_encapsulation_mask",
        "reference_transcripts",
    ),
    "_compute_splice_graph_cache_entry": (
        "transcripts",
        "quant_mode",
        "reference_transcripts",
    ),
    # bam_file_for_priors is the pass-1 theta bam: in cluster-guided single-cell
    # the splice graph comes from the shared sg bam while priors must come from
    # this cluster's own reads. Appended, never inserted.
    "run_quant_only": (
        "rescue_summary_path",
        "reference_transcripts",
        "bam_file_for_priors",
    ),
}


def test_the_parameter_tail_is_exactly_as_declared():
    """Pins the ordered tail of every signature reference_transcripts was added to.

    Two earlier forms of this test were weaker and both are worth recording,
    because each looked adequate:

    "is the last parameter" -- right property, but it forbids ever adding an
    argument, so the only way to add one and stay green is to insert BEFORE
    reference_transcripts, which is the exact mistake being guarded.

    "everything after it has a default" -- VACUOUS. Python forbids a required
    parameter after a defaulted one, so it can never fail; three negative
    controls, including the real hazard, all passed against it.

    An ordered tail fails on both an insertion and an undeclared append, and
    names what it expected.
    """

    for func in (
        Splice_graph.Splice_graph.build_splice_graph_for_contig,
        lraa._compute_splice_graph_cache_entry,
        lraa.run_quant_only,
    ):
        key = func.__name__
        expected = _EXPECTED_TAIL[key]
        names = list(inspect.signature(func).parameters)
        assert tuple(names[-len(expected) :]) == expected, (
            "{}: parameter tail is {} but {} was declared. Appending is safe -- "
            "update _EXPECTED_TAIL. Inserting before reference_transcripts is "
            "not, and shifts positional callers onto the reference "
            "set.".format(key, tuple(names[-len(expected) :]), expected)
        )

def test_reference_transcripts_defaults_to_no_reference():
    """Ref-free is the default at every entry point, so a caller cannot acquire a
    reference set by omission."""
    for func in (
        Splice_graph.Splice_graph.build_splice_graph_for_contig,
        lraa._compute_splice_graph_cache_entry,
        lraa.run_quant_only,
    ):
        assert (
            inspect.signature(func).parameters["reference_transcripts"].default is None
        )
