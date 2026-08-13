#!/usr/bin/env python3

"""`Util_funcs.available_cpus` must report CPUs this process may RUN on, not CPUs the
machine happens to have. Sizing a pool from the latter oversubscribes by the ratio
between them under any cpuset -- containers, Slurm `--cpus-per-task`, taskset.
"""

import os

import pytest

import Util_funcs


def test_prefers_scheduler_affinity_over_machine_cpu_count(monkeypatch):
    """The whole point: a 4-of-16 cpuset must yield 4, not 16."""
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: {0, 1, 2, 3}, raising=False)
    monkeypatch.setattr(os, "cpu_count", lambda: 16)
    assert Util_funcs.available_cpus() == 4


def test_falls_back_to_cpu_count_where_affinity_is_unavailable(monkeypatch):
    """sched_getaffinity is Linux-only; on macOS and Windows it does not exist, and
    referencing it unguarded would raise AttributeError rather than degrade."""
    monkeypatch.delattr(os, "sched_getaffinity", raising=False)
    monkeypatch.setattr(os, "cpu_count", lambda: 8)
    assert Util_funcs.available_cpus() == 8


def test_falls_back_when_affinity_raises(monkeypatch):
    def boom(pid):
        raise OSError("no affinity for you")

    monkeypatch.setattr(os, "sched_getaffinity", boom, raising=False)
    monkeypatch.setattr(os, "cpu_count", lambda: 6)
    assert Util_funcs.available_cpus() == 6


def test_never_returns_less_than_one(monkeypatch):
    """os.cpu_count() returns None when the count is indeterminable; a thread pool sized
    at None or 0 is worse than one sized conservatively."""
    monkeypatch.delattr(os, "sched_getaffinity", raising=False)
    monkeypatch.setattr(os, "cpu_count", lambda: None)
    assert Util_funcs.available_cpus() == 1

    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: set(), raising=False)
    assert Util_funcs.available_cpus() == 1


def test_reports_something_sane_unmocked():
    n = Util_funcs.available_cpus()
    assert isinstance(n, int)
    assert n >= 1
    if hasattr(os, "sched_getaffinity"):
        assert n == len(os.sched_getaffinity(0))


def _load_lraa_script():
    """The clamp under test lives in the `LRAA` driver script, which has no .py suffix."""
    import importlib.machinery
    import importlib.util

    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(here, "LRAA")
    loader = importlib.machinery.SourceFileLoader("lraa_script_under_test", path)
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def test_singleton_stage_threads_respect_the_cpuset(monkeypatch):
    """The defect this fixes, at its call site.

    `_resolve_singleton_stage_threads` clamps a requested thread count to the hardware.
    Reading that ceiling from `os.cpu_count()` ignores the cpuset, so a request for 64
    threads inside a 4-of-16 container was clamped to 16 -- a 4x oversubscription of the
    CPUs actually granted. On `devel` this returns 16; it must now return 4.
    """
    lraa = _load_lraa_script()
    monkeypatch.setattr(lraa, "TOTAL_THREADS_FOR_SINGLETON_STAGE", 64, raising=False)
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: {0, 1, 2, 3}, raising=False)
    monkeypatch.setattr(os, "cpu_count", lambda: 16)

    assert lraa._resolve_singleton_stage_threads() == 4
