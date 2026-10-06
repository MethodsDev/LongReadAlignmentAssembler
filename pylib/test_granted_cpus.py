#!/usr/bin/env python3

"""Util_funcs.granted_cpus: the smaller of the affinity count and the cgroup quota.

The quota is what docker --cpus and Terra enforce, and it is invisible to
sched_getaffinity, so a pool sized from affinity alone oversubscribes there.
"""

import os
import Util_funcs


def _write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def test_cgroup_v2_quota(tmp_path):
    _write(str(tmp_path / "cpu.max"), "800000 100000\n")
    assert Util_funcs.cgroup_cpu_quota(str(tmp_path)) == 8


def test_cgroup_v2_unlimited(tmp_path):
    _write(str(tmp_path / "cpu.max"), "max 100000\n")
    assert Util_funcs.cgroup_cpu_quota(str(tmp_path)) is None


def test_cgroup_v1_quota_rounds_down(tmp_path):
    _write(str(tmp_path / "cpu" / "cpu.cfs_quota_us"), "250000\n")
    _write(str(tmp_path / "cpu" / "cpu.cfs_period_us"), "100000\n")
    assert Util_funcs.cgroup_cpu_quota(str(tmp_path)) == 2


def test_fractional_grant_floors_at_one(tmp_path):
    _write(str(tmp_path / "cpu.max"), "50000 100000\n")
    assert Util_funcs.cgroup_cpu_quota(str(tmp_path)) == 1


def test_no_cgroup_files(tmp_path):
    assert Util_funcs.cgroup_cpu_quota(str(tmp_path)) is None


def test_quota_below_affinity_wins(monkeypatch):
    monkeypatch.setattr(Util_funcs, "available_cpus", lambda: 16)
    monkeypatch.setattr(Util_funcs, "cgroup_cpu_quota", lambda root=None: 4)
    assert Util_funcs.granted_cpus() == 4


def test_affinity_below_quota_wins(monkeypatch):
    monkeypatch.setattr(Util_funcs, "available_cpus", lambda: 3)
    monkeypatch.setattr(Util_funcs, "cgroup_cpu_quota", lambda root=None: 8)
    assert Util_funcs.granted_cpus() == 3


def test_no_quota_uses_affinity(monkeypatch):
    monkeypatch.setattr(Util_funcs, "available_cpus", lambda: 12)
    monkeypatch.setattr(Util_funcs, "cgroup_cpu_quota", lambda root=None: None)
    assert Util_funcs.granted_cpus() == 12
