#!/usr/bin/env python3
"""Scripts inlined into WDL tasks must equal their tested util sources.

A WDL task that calls a newly added helper by name breaks whenever its pinned image
predates the helper, so some helpers are inlined into the task command as a quoted
heredoc instead. The util script stays the source and carries the unit tests; this test
keeps the inlined copy from drifting. To fix a failure, paste the util script, unchanged,
between the heredoc markers in the WDL task.
"""

from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
MARKER = "LRAA_INLINED_SCRIPT_EOF"

# (wdl, task, util script)
INLINED = [
    ("WDL/LRAA-cell_cluster_guided.wdl", "split_pooled_quant_by_cluster",
     "util/sc/tracking_to_cluster_pseudobulk_exprs.py"),
]


def _inlined(wdl, task):
    src = (REPO / wdl).read_text()
    start = src.index("task {} {{".format(task))
    opener = "<<'{}'\n".format(MARKER)
    i = src.index(opener, start) + len(opener)
    j = src.index("\n{}\n".format(MARKER), i) + 1
    return src[i:j]


@pytest.mark.parametrize("wdl,task,script", INLINED, ids=[t for _, t, _ in INLINED])
def test_inlined_copy_matches_util_script(wdl, task, script):
    assert _inlined(wdl, task) == (REPO / script).read_text(), (
        "{} task {} inlines a copy of {} that differs from it; paste the script into "
        "the heredoc unchanged".format(wdl, task, script))
