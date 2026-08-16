#!/usr/bin/env python3

"""Command-line front end for the chunked quant pipeline and its control arm.

The pipeline itself lives in ``pylib/ChunkedRun.py``, where ``LRAA --chunk``
also reaches it. This script stays because it is what the parity evaluation and
the runtime-performance harness invoke, and because ``--arm baseline`` is a
measurement device rather than a user-facing mode: LRAA has no reason to offer
"run the whole-contig control", but the comparison that keeps chunking honest
has every reason to.
"""

import os
import sys

sys.path.insert(
    0,
    os.path.sep.join(
        [os.path.dirname(os.path.realpath(__file__)), "..", "..", "pylib"]
    ),
)

import ChunkedRun  # noqa: E402  (path insert must precede the import)


if __name__ == "__main__":
    try:
        ChunkedRun.main()
        sys.exit(0)
    except ChunkedRun.PipelineError as err:
        print("\nPIPELINE FAILED\n{}".format(err), file=sys.stderr)
        sys.exit(1)
