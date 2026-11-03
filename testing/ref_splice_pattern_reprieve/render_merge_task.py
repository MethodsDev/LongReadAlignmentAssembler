#!/usr/bin/env python3
"""Print the Terra merge task's command (WDL/LRAA-cell_cluster_guided.wdl lraa_merge_gtf_task)
for these inputs, as miniwdl renders it, so the Makefile can run it outside a container.

usage: render_merge_task.py <genome.fa> <pooled.gtf> <cluster.gtf> [config_update_json]
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../pylib"))
from test_wdl_pooled_gtf_merge import render_merge_task_command  # noqa: E402

genome, pooled, cluster = sys.argv[1:4]
inputs = {
    "sample_id": "S",
    "LRAA_cell_cluster_gtfs": [cluster],
    "pooled_gtf": pooled,
    "referenceGenome": genome,
    "docker": "unused",
    "memoryGB": 8,
}
if len(sys.argv) > 4:
    inputs["config_update_json"] = sys.argv[4]
print(render_merge_task_command(inputs))
