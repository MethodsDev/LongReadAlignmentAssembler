#!/usr/bin/env python3

"""Every WDL task must name a container image.

This exists because a missing one is invisible to this repo's entire local test
surface. miniwdl does not require `docker` in a runtime block -- it substitutes a
default image and runs the task -- while Cromwell on GCP Batch rejects the job
outright:

    No container image found in either 'container' or 'docker' runtime attributes.

So `require_annot_gtf` in WDL/LRAA-cell_cluster_guided.wdl shipped with cpu,
memory and disks but no docker, passed `miniwdl check`, passed a full
`make test_wdls`, passed three complete by_chunk single-cell smoke runs, and then
failed a real Terra submission. It was reachable only in quant_only mode, which
is why nothing local that ran it happened to be a Cromwell run.

`womtool validate` does not catch it either: an absent runtime attribute is legal
WDL. The requirement is a BACKEND one, so it needs a repo-level check, and this
is it -- cheap, static, and it fails on the exact shape that reached production.
"""

import re
from pathlib import Path

import pytest

WDL_DIR = Path(__file__).resolve().parents[1] / "WDL"

# `task NAME {` at column 0 -- every task in this repo is declared top-level.
TASK_RE = re.compile(r"^task\s+(\w+)\s*\{", re.M)
BOUNDARY_RE = re.compile(r"^(task|workflow)\s+\w+\s*\{", re.M)
RUNTIME_RE = re.compile(r"^\s*runtime\s*\{(.*?)^\s*\}", re.M | re.S)
# Cromwell accepts either spelling; WDL 1.1 renamed docker to container.
IMAGE_RE = re.compile(r"^\s*(docker|container)\s*:", re.M)


def _tasks(path):
    """(task name, its runtime body or None) for every task in ``path``."""

    src = path.read_text()
    out = []
    for match in TASK_RE.finditer(src):
        after = src[match.end() :]
        nxt = BOUNDARY_RE.search(after)
        body = src[match.start() : match.end() + (nxt.start() if nxt else len(after))]
        runtime = RUNTIME_RE.search(body)
        out.append((match.group(1), runtime.group(1) if runtime else None))
    return out


def _wdl_files():
    files = sorted(WDL_DIR.rglob("*.wdl"))
    assert files, "no WDL files found under {}".format(WDL_DIR)
    return files


@pytest.mark.parametrize("path", _wdl_files(), ids=lambda p: p.name)
def test_every_task_names_a_container_image(path):
    """A task with no image is a job Cromwell refuses to schedule."""

    offenders = []
    for name, runtime in _tasks(path):
        if runtime is None:
            offenders.append("{}: no runtime block at all".format(name))
        elif not IMAGE_RE.search(runtime):
            offenders.append("{}: runtime block with no docker/container".format(name))

    assert not offenders, (
        "{} has task(s) Cromwell cannot schedule -- GCP Batch reports 'No "
        "container image found in either container or docker runtime "
        "attributes' and the job fails. miniwdl runs these fine on a default "
        "image, so local green proves nothing here:\n  {}".format(
            path.name, "\n  ".join(offenders)
        )
    )


def test_the_check_can_actually_fail(tmp_path):
    """Guard the guard: a task with no image must be detected.

    Without this, a refactor of the parsing above could silently stop finding
    anything and every file would pass vacuously -- which is precisely how the
    original defect survived a green suite.
    """

    bad = tmp_path / "bad.wdl"
    bad.write_text(
        "version 1.0\n\n"
        "task no_image {\n"
        "    command <<<\n        echo hi\n    >>>\n"
        "    runtime {\n        cpu: 1\n        memory: \"1 GiB\"\n    }\n"
        "}\n"
    )

    names = [n for n, rt in _tasks(bad)]
    assert names == ["no_image"], names

    offenders = [
        n for n, rt in _tasks(bad) if rt is None or not IMAGE_RE.search(rt)
    ]
    assert offenders == ["no_image"], offenders
