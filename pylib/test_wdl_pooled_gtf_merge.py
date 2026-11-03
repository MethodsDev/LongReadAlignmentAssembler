#!/usr/bin/env python3

"""The single-cell workflow merges the pooled catalog into the cluster-guided one.

LRAA-singlecell.wdl passes the initial pass's gtf to LRAA-cell_cluster_guided.wdl as
pooled_gtf, and the merge task hands it to merge_LRAA_GTFs.py --pooled_gtf. That is what
keeps isoforms and TSS/PolyA sites the per-cluster runs dilute, and why the workflow no
longer integrates the initial run's site beds separately. Dropping either link would
silently fall back to a clusters-only merge: every WDL still validates, nothing errors.
"""

from pathlib import Path

import WDL

WDL_DIR = Path(__file__).resolve().parents[1] / "WDL"


def _call(workflow, name):
    def walk(nodes):
        for node in nodes:
            if isinstance(node, WDL.Tree.Call) and node.name == name:
                return node
            found = walk(getattr(node, "body", []))
            if found is not None:
                return found
        return None

    return walk(workflow.body)


def test_singlecell_passes_init_gtf_as_pooled_gtf():
    doc = WDL.load(str(WDL_DIR / "LRAA-singlecell.wdl"))
    call = _call(doc.workflow, "cluster_guided")
    assert call is not None
    assert str(call.inputs["pooled_gtf"]) == "init_gtf_file"


def test_cluster_guided_merge_forwards_pooled_gtf():
    doc = WDL.load(str(WDL_DIR / "LRAA-cell_cluster_guided.wdl"))
    call = _call(doc.workflow, "lraa_merge_gtf_task")
    assert call is not None
    assert str(call.inputs["pooled_gtf"]) == "pooled_gtf"
    task = next(t for t in doc.tasks if t.name == "lraa_merge_gtf_task")
    assert '"--pooled_gtf " + pooled_gtf' in str(task.command)


def test_singlecell_no_longer_integrates_site_beds():
    doc = WDL.load(str(WDL_DIR / "LRAA-singlecell.wdl"))
    assert _call(doc.workflow, "integrate_sites") is None
    outputs = {o.name for o in doc.workflow.outputs}
    assert not outputs & {"integrated_TSS_bed", "integrated_PolyA_bed", "integrated_sites_summary"}


# ---- merge config (v0.50.0): the merge task takes --config_update ---------------------


def render_merge_task_command(inputs):
    """The lraa_merge_gtf_task command text for these task inputs, as miniwdl would run it."""
    doc = WDL.load(str(WDL_DIR / "LRAA-cell_cluster_guided.wdl"))
    task = next(t for t in doc.tasks if t.name == "lraa_merge_gtf_task")
    env = WDL.values_from_json(inputs, task.available_inputs, task.required_inputs)
    stdlib = WDL.StdLib.Base(wdl_version="1.0")
    for decl in (task.inputs or []) + (task.postinputs or []):
        if not env.has_binding(decl.name):
            value = decl.expr.eval(env, stdlib) if decl.expr is not None else WDL.Value.Null()
            env = env.bind(decl.name, value)
    return task.command.eval(env, stdlib).value


_MERGE_INPUTS = {
    "sample_id": "S",
    "LRAA_cell_cluster_gtfs": ["c1.gtf"],
    "pooled_gtf": "pooled.gtf",
    "referenceGenome": "g.fa",
    "docker": "d",
    "memoryGB": 8,
}


def test_merge_task_passes_config_update_only_when_set():
    unset = render_merge_task_command(_MERGE_INPUTS)
    assert "--config_update" not in unset
    cfg = '{"allow_unflagged_LRAA_guide": true}'
    command = render_merge_task_command(dict(_MERGE_INPUTS, config_update_json=cfg))
    assert '--config_update "$config_update_file"' in command
    # the JSON reaches the file the flag names, verbatim, through a quoted heredoc
    assert "<<'LRAA_CONFIG_UPDATE_EOF'\n" + cfg + "\nLRAA_CONFIG_UPDATE_EOF" in command


def test_merge_config_is_plumbed_from_both_workflows():
    sc = WDL.load(str(WDL_DIR / "LRAA-singlecell.wdl"))
    assert str(_call(sc.workflow, "cluster_guided").inputs["merge_config_update_json"]) == "merge_config_update_json"
    assert "merge_config_update_json" in {d.name for d in sc.workflow.available_inputs}

    guided = WDL.load(str(WDL_DIR / "LRAA-cell_cluster_guided.wdl"))
    call = _call(guided.workflow, "lraa_merge_gtf_task")
    assert str(call.inputs["config_update_json"]) == "merge_config_update_json"

    sample_set = WDL.load(str(WDL_DIR / "LRAA-merge_gtfs_from_sample_set.wdl"))
    call = _call(sample_set.workflow, "lraa_merge_gtf_task")
    assert str(call.inputs["config_update_json"]) == "config_update_json"
