#!/usr/bin/env python3

"""Calibrate a library's TSS cap signature before running LRAA.

The same calibration LRAA runs itself at the start of a discovery run
(pylib/TSS_cap_calibration.py), as a separate step: a workflow runs it once on the full
bam and hands <prefix>.TSS_cap_calibration.config.json to every LRAA run of that
library (--config_update), so shards and cell clusters share one decision. Also writes
the report: <prefix>.TSS_cap_calibration.json (decision and statistics), .hist.tsv (the
per-peak capped-fraction histogram) and .svg (its plot).
"""

import sys, os

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../pylib"])
)

import TSS_cap_calibration

if __name__ == "__main__":
    TSS_cap_calibration.main()
