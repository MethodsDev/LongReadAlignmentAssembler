#!/usr/bin/env python3

"""The merge must not prune a TSS/PolyA node for being asserted by few inputs.

Splice_graph's boundary-support filters are written for read pileups, where a minority 5'
end is plausibly a degradation product. During a merge the counter holds input-model
assertions, so "minority support" means "rebuilt in few of the inputs" -- a cell-type
restricted model, precisely what the merge exists to carry forward. Once its boundary node
is purged its path is indistinguishable from an unannotated contained path, and
LRAA._validate_pairwise_incompatibilities absorbs it into any longer model spanning it.

Measured on the v0.35.0 PBMC cluster-guided run: of the FSM splice chains present in the
pooled init GTF and absent from the merged cluster GTF, 979 were built at cluster depth and
removed by the merge, and re-running one contig's merge with these thresholds zeroed brought
back chains that the stock settings dropped.

The zeroing is asserted here rather than left to the config defaults because it is a merge
policy: LRAA_Globals carries the read-pass values, and a reader of that file has no way to
know the merge overrides them.
"""

import importlib.machinery
import importlib.util
import os
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if os.path.join(REPO_ROOT, "pylib") not in sys.path:
    sys.path.insert(0, os.path.join(REPO_ROOT, "pylib"))

import LRAA_Globals


def _load_merge_script():
    path = os.path.join(REPO_ROOT, "util", "merge_LRAA_GTFs.py")
    loader = importlib.machinery.SourceFileLoader("merge_script_under_test", path)
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


@pytest.fixture
def merge_module():
    module = _load_merge_script()
    saved = {key: LRAA_Globals.config[key] for key in module.BOUNDARY_SUPPORT_FILTERS}
    yield module
    LRAA_Globals.config.update(saved)


def test_thresholds_are_zeroed_by_default(merge_module):
    configured = merge_module.configure_boundary_support_filters
    configured(apply_filters=False)

    for key in merge_module.BOUNDARY_SUPPORT_FILTERS:
        assert LRAA_Globals.config[key] == 0.0, key


def test_zeroed_thresholds_purge_nothing(merge_module):
    """Zero has to mean "off" for the comparisons that read these keys.

    _eliminate_low_support_TSS drops a site when its share is `< min_TSS_iso_fraction`, and
    the degradation walk drops a neighbour at `<= max_frac_alt_TSS_from_degradation` of a
    dominant site. A share is never negative, so zero disables the first; the second is an
    inclusive comparison, so zero still purges a site holding no support at all, which is
    the only case it may drop.
    """
    merge_module.configure_boundary_support_filters(apply_filters=False)

    min_fraction = LRAA_Globals.config["min_TSS_iso_fraction"]
    degradation = LRAA_Globals.config["max_frac_alt_TSS_from_degradation"]

    # a site asserted by 1 input against a dominant site asserted by 500
    lone_share = 1 / (1 + 500)
    assert not (lone_share < min_fraction)
    assert not (1 / 500 <= degradation)

    # and a site with no support behind it is still eligible for removal
    assert 0 / 500 <= degradation


def test_filters_restorable_for_reproducing_older_runs(merge_module):
    defaults = {
        "min_TSS_iso_fraction": 0.05,
        "min_PolyA_iso_fraction": 0.05,
        "max_frac_alt_TSS_from_degradation": 0.20,
    }
    merge_module.configure_boundary_support_filters(apply_filters=False)
    LRAA_Globals.config.update(defaults)

    merge_module.configure_boundary_support_filters(apply_filters=True)

    for key, value in defaults.items():
        assert LRAA_Globals.config[key] == value, key


def test_flag_defaults_to_disabling_the_filters(merge_module):
    """The CLI default is the policy: absent flag means filters off."""
    parser = None
    for action in _merge_parser_actions(merge_module):
        if action.dest == "apply_boundary_support_filters":
            parser = action
            break

    assert parser is not None, "--apply_boundary_support_filters not registered"
    assert parser.default is False


def _merge_parser_actions(module):
    """Build the merge parser without running a merge."""
    import argparse

    captured = {}
    real_parse = argparse.ArgumentParser.parse_args

    def capture(self, *args, **kwargs):
        captured["parser"] = self
        raise SystemExit(0)

    argparse.ArgumentParser.parse_args = capture
    try:
        module.main()
    except SystemExit:
        pass
    finally:
        argparse.ArgumentParser.parse_args = real_parse

    return captured["parser"]._actions
