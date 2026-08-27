#!/usr/bin/env python3

"""`diskSizeGB` must survive every hop from the cluster-guided entrypoint to `disks:`.

A declared-but-unbound resource input is the worst failure mode this repo has:
the caller sets a number, the run accepts it, and then a task dies on a full
disk while its inputs report that the request was honoured. That is not
hypothetical -- LRAA-cell_cluster_guided.wdl declared `Int diskSizeGB = 256` and
bound it nowhere, so raising it bought nothing. miniwdl reported it the whole
time ("UnusedDeclaration, nothing references Int diskSizeGB") and it stayed.

`UnusedDeclaration` only catches the FIRST hop, though. Once a value is bound
once, miniwdl is satisfied, and a mid-chain workflow that forwards it to some
calls but not others -- LRAA_quant_by_cluster.wdl's final quant was exactly
that -- is silent again. So the assertion here is the whole path, hop by hop,
ending at the runtime attribute that actually spends the number:

    LRAA_cell_cluster_guided
      -> LRAA_quant_by_cluster           (final quant)
        -> LRAA_wf                       (per cluster)
          -> LRAA_runner
            -> LRAA_runner_task          disks: "local-disk ~{diskSizeGB} HDD"

Each hop must bind `diskSizeGB` to a plain reference to the caller's own
same-named input. Anything else -- dropped, renamed, or hard-coded to a literal
-- means a caller's value stops there, which is the defect.

Static parse only: no docker, no fixtures, no execution.
"""

import pytest

WDL = pytest.importorskip("WDL", reason="miniwdl not installed")

from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
WDL_DIR = REPO / "WDL"

# (importing document, the call to inspect). Each entry's callee is the next
# entry's document, so the list IS the chain.
HOPS = [
    ("LRAA-cell_cluster_guided.wdl", "LRAA_quant_final_bamlist"),
    ("LRAA_quant_by_cluster.wdl", "LRAA_quant_cluster"),
    # LRAA.wdl reaches the runner two ways -- chromosome-sharded and direct --
    # and a run takes exactly one of them, so both have to carry the value.
    ("LRAA.wdl", "LRAA_scatter"),
    ("LRAA.wdl", "LRAA_direct"),
    ("subwdls/LRAA_runner.wdl", "LRAA_runner_task"),
]

RESOURCE = "diskSizeGB"


def _calls(node):
    """Every Call anywhere under ``node``, including inside if/scatter bodies."""

    for child in getattr(node, "body", None) or []:
        if isinstance(child, WDL.Tree.Call):
            yield child
        else:
            yield from _calls(child)


def _load(rel):
    return WDL.load(str(WDL_DIR / rel), path=[str(WDL_DIR)])


@pytest.mark.parametrize("rel,call_name", HOPS, ids=[f"{r}:{c}" for r, c in HOPS])
def test_hop_forwards_disksizegb(rel, call_name):
    workflow = _load(rel).workflow
    call = next((c for c in _calls(workflow) if c.name == call_name), None)
    assert call is not None, f"{rel} has no call named {call_name}"

    declared = {d.name for d in workflow.inputs or []}
    assert RESOURCE in declared, f"{workflow.name} in {rel} declares no {RESOURCE}"

    binding = call.inputs.get(RESOURCE)
    assert binding is not None, (
        f"{rel}: call {call_name} does not pass {RESOURCE}, so a caller's value "
        f"stops at {workflow.name} and everything below sizes its own disks"
    )
    assert str(binding) == RESOURCE, (
        f"{rel}: call {call_name} binds {RESOURCE} = {binding}, not the "
        f"workflow's own {RESOURCE} input"
    )

    accepts = {d.name for d in call.callee.inputs or []}
    assert RESOURCE in accepts, f"{rel}: callee of {call_name} accepts no {RESOURCE}"


def test_runner_task_spends_it_on_disks():
    """The end of the chain: the value must land in the `disks:` attribute."""

    doc = _load("subwdls/LRAA_runner.wdl")
    task = next(t for t in doc.tasks if t.name == "LRAA_runner_task")
    assert RESOURCE in {d.name for d in task.inputs or []}

    disks = task.runtime.get("disks")
    assert disks is not None, "LRAA_runner_task declares no disks: attribute"
    assert RESOURCE in str(disks), f"disks: does not interpolate {RESOURCE}: {disks}"


def test_chain_defaults_agree():
    """One number, quoted identically at every level.

    Divergent defaults would mean the unset behaviour of a nested call silently
    differs from its caller's, which is how a "threaded" input turns into a
    number nobody can predict from the top-level inputs alone.
    """

    seen = {}
    for rel in ["LRAA.wdl", "LRAA-singlecell.wdl", "LRAA-cell_cluster_guided.wdl",
                "LRAA_quant_by_cluster.wdl"]:
        workflow = _load(rel).workflow
        decl = next((d for d in workflow.inputs or [] if d.name == RESOURCE), None)
        assert decl is not None, f"{rel} declares no {RESOURCE}"
        assert str(decl.type) == "Int", f"{rel}: {RESOURCE} is {decl.type}, not Int"
        seen[rel] = str(decl.expr)

    assert set(seen.values()) == {"256"}, f"defaults disagree: {seen}"
