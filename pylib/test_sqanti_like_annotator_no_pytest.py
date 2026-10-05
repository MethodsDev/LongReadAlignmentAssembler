#!/usr/bin/env python3

"""SQANTI_like_annotator must import without pytest.

util/SQANTI-like_cats_for_reads_or_isoforms.py runs in lraa-core, which
deliberately ships without pytest. The annotator and MockTranscript once
imported pytest at module level, so the sqanti-like WDL failed there with
ModuleNotFoundError before doing any work.
"""

import os
import subprocess
import sys

PYLIB = os.path.dirname(os.path.abspath(__file__))


def test_annotator_imports_with_pytest_unavailable():
    # sys.modules[name] = None makes any later `import name` raise ImportError.
    code = (
        "import sys; sys.modules['pytest'] = None; "
        "sys.path.insert(0, {!r}); "
        "import SQANTI_like_annotator, MockTranscript".format(PYLIB)
    )
    proc = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True
    )
    assert proc.returncode == 0, proc.stderr
