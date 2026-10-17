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
