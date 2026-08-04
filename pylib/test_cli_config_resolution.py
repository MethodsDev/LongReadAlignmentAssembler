#!/usr/bin/env python3

import argparse
import importlib.util
import json
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path
from types import SimpleNamespace

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))


def _load_lraa_cli_module():
    loader = SourceFileLoader("lraa_cli_config_test_module", str(REPO_ROOT / "LRAA"))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


LRAA_CLI = _load_lraa_cli_module()


def _base_config():
    return {
        "normalize_max_cov_level": 1000,
        "min_alt_splice_freq": 0.03,
        "min_alt_unspliced_freq": 0.01,
        "min_isoform_fraction": 0.01,
        "run_EM": True,
        "EM_alpha": 0.01,
    }


def _args(**updates):
    values = {
        "normalize_max_cov_level": 1000,
        "min_isoform_fraction": 0.01,
        "no_EM": False,
        "EM_alpha": 0.01,
        "HiFi": False,
    }
    values.update(updates)
    return SimpleNamespace(**values)


def _parse_splice_args(argv):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    LRAA_CLI._add_splice_frequency_arguments(parser)
    return parser, parser.parse_args(argv)


def test_splice_options_track_omission_without_changing_documented_defaults():
    parser, omitted = _parse_splice_args([])
    assert not hasattr(omitted, "min_alt_splice_freq")
    assert not hasattr(omitted, "min_alt_unspliced_freq")

    help_text = parser.format_help()
    assert "(default: 0.03)" in help_text
    assert "(default: 0.01)" in help_text

    _, explicit = _parse_splice_args(
        ["--min_alt_splice_freq", "0.03", "--min_alt_unspliced_freq", "0.02"]
    )
    assert explicit.min_alt_splice_freq == 0.03
    assert explicit.min_alt_unspliced_freq == 0.02


def test_ont_and_hifi_splice_defaults_preserve_explicit_hifi_value():
    ont_config = _base_config()
    ont_args = _args()
    LRAA_CLI._seed_authoritative_config_from_args(ont_config, ont_args)
    LRAA_CLI._apply_hifi_splice_default(ont_config, ont_args)
    assert ont_config["min_alt_splice_freq"] == 0.03

    hifi_config = _base_config()
    hifi_args = _args(HiFi=True)
    LRAA_CLI._seed_authoritative_config_from_args(hifi_config, hifi_args)
    LRAA_CLI._apply_hifi_splice_default(hifi_config, hifi_args)
    assert hifi_config["min_alt_splice_freq"] == 0.01

    explicit_config = _base_config()
    explicit_hifi_args = _args(HiFi=True, min_alt_splice_freq=0.03)
    LRAA_CLI._seed_authoritative_config_from_args(explicit_config, explicit_hifi_args)
    LRAA_CLI._apply_hifi_splice_default(explicit_config, explicit_hifi_args)
    assert explicit_config["min_alt_splice_freq"] == 0.03


def test_named_values_seed_config_and_json_wins_for_all_consumers():
    config = _base_config()
    named_args = _args(
        normalize_max_cov_level=321,
        min_alt_splice_freq=0.07,
        min_alt_unspliced_freq=0.08,
        min_isoform_fraction=0.09,
        no_EM=True,
        EM_alpha=0.10,
    )
    LRAA_CLI._seed_authoritative_config_from_args(config, named_args)
    assert config == {
        "normalize_max_cov_level": 321,
        "min_alt_splice_freq": 0.07,
        "min_alt_unspliced_freq": 0.08,
        "min_isoform_fraction": 0.09,
        "run_EM": False,
        "EM_alpha": 0.10,
    }

    json_update = json.loads(
        """{
            "normalize_max_cov_level": 654,
            "min_alt_splice_freq": 0.02,
            "min_alt_unspliced_freq": 0.04,
            "min_isoform_fraction": 0.05,
            "run_EM": true,
            "EM_alpha": 0.06
        }"""
    )
    applied, skipped = LRAA_CLI._apply_config_overrides(config, json_update)
    assert set(applied) == set(json_update)
    assert skipped == []

    runtime = LRAA_CLI._snapshot_runtime_config(config)
    assert runtime == {
        "normalize_max_cov_level": 654,
        "min_isoform_fraction": 0.05,
        "run_EM": True,
        "EM_alpha": 0.06,
    }
    splice_graph_params = LRAA_CLI._build_splice_graph_params(config)
    assert splice_graph_params["min_alt_splice_freq"] == 0.02
    assert splice_graph_params["min_alt_unspliced_freq"] == 0.04
