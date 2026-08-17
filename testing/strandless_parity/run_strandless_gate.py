#!/usr/bin/env python3

"""Command-line front end for the strandless-chunking equivalence gate.

The gate itself lives in ``pylib/StrandlessParity.py``, beside the pipeline it
measures, for the same reason ``util/misc/run_chunked_quant_pipeline.py`` is a
front end onto ``pylib/ChunkedRun.py``: the comparison has to import the
pipeline's own severed-read accounting, table reader and merge arithmetic rather
than reimplement them, and a second implementation of any of those would be a
second parity notion that could agree with itself and with nothing else.

Subcommands, cheapest and most falsifiable first:

    strand-invariant   per-read agreement between the whole-contig strand split
                       and a split run inside the chunk holding that read
    ordering-cost      what a strand-filtered extraction from a raw BAM would
                       misassign, executed on purpose, named read by read
    compare-arms       two finished pipeline runs, diffed
    gate               all of the above over one substrate, with the strand-first
                       arm's own figures pinned so a strandless figure means
                       something

Run ``--help`` on any of them.  ``make test`` in this directory is the tiny-corpus
form; the real-chromosome form is the same ``gate`` invocation pointed at a real
chromosome, and is not part of routine testing for the reason
``docs/chunked_quant_parity_evaluation.md`` gives: it needs a real chromosome with
real annotation to be worth running.
"""

import os
import sys

sys.path.insert(
    0,
    os.path.sep.join(
        [os.path.dirname(os.path.realpath(__file__)), "..", "..", "pylib"]
    ),
)

import StrandlessParity  # noqa: E402  (path insert must precede the import)


if __name__ == "__main__":
    try:
        sys.exit(StrandlessParity.main())
    except StrandlessParity.ParityError as err:
        print("\nGATE FAILED\n{}".format(err), file=sys.stderr)
        sys.exit(2)
