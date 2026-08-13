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


def test_reference_transcripts_is_the_last_parameter_everywhere_it_was_added():
    assert (
        _last_parameter(Splice_graph.Splice_graph.build_splice_graph_for_contig)
        == "reference_transcripts"
    )
    assert (
        _last_parameter(lraa._compute_splice_graph_cache_entry)
        == "reference_transcripts"
    )
    assert _last_parameter(lraa.run_quant_only) == "reference_transcripts"


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
