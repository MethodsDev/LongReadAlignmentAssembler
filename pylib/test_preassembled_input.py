#!/usr/bin/env python3

"""--preassembled_input: read minimums of 1 for preassembled transcripts.

Long-read runs keep the defaults (2); an explicit --min_reads_novel or
--config_update still wins; chunked mode forwards the resolved values.
"""

import importlib.util
import sys
from argparse import Namespace
from importlib.machinery import SourceFileLoader
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA_Globals


def _load_lraa():
    loader = SourceFileLoader("lraa_preassembled_input_test", str(REPO_ROOT / "LRAA"))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


LRAA = _load_lraa()
KEYS = ("min_intron_support", "min_reads_novel_isoform", "min_unique_reads_novel_isoform")


def _config(min_reads_novel=2):
    config = dict(LRAA_Globals.config)
    config["min_reads_novel_isoform"] = min_reads_novel
    config["min_unique_reads_novel_isoform"] = min_reads_novel
    return config


def _args(specified=(), **kw):
    return Namespace(user_specified_dests=set(specified), config_update=None, **kw)


def test_long_read_defaults_stay_at_two():
    assert all(LRAA_Globals.config[k] == 2 for k in KEYS)
    config = _config()
    LRAA._apply_preassembled_input_config_overrides(config, _args(preassembled_input=False))
    assert all(config[k] == 2 for k in KEYS)
    assert config["preassembled_input"] is False


def test_preassembled_input_lowers_minimums_to_one():
    config = _config()
    args = _args({"preassembled_input"}, preassembled_input=True)
    LRAA._apply_preassembled_input_config_overrides(config, args)
    assert all(config[k] == 1 for k in KEYS)
    assert config["preassembled_input"] is True


def test_explicit_min_reads_novel_wins():
    config = _config(min_reads_novel=3)
    args = _args({"preassembled_input", "min_reads_novel"}, preassembled_input=True,
                 min_reads_novel=3)
    LRAA._apply_preassembled_input_config_overrides(config, args)
    assert config["min_intron_support"] == 1
    assert config["min_reads_novel_isoform"] == 3
    assert config["min_unique_reads_novel_isoform"] == 3


def test_chunk_workers_receive_resolved_values():
    config = _config()
    args = _args({"preassembled_input", "pasa_mode"}, preassembled_input=True,
                 pasa_mode=True)
    LRAA._apply_pasa_mode_config_overrides(config, args)
    LRAA._apply_preassembled_input_config_overrides(config, args)
    forwarded = LRAA._collect_explicit_config_overrides(args, None, config)
    assert forwarded["preassembled_input"] is True and forwarded["pasa_mode"] is True
    assert all(forwarded[k] == 1 for k in KEYS)
    assert all(forwarded[k] == v for k, v in LRAA_Globals.PASA_MODE_CONFIG.items())
