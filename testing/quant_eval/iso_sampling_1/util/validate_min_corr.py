#!/usr/bin/env python3

import sys

fname = sys.argv[1]
min_val = sys.argv[2]

min_val = float(min_val)

with open(fname, "rt") as fh:
    header = next(fh)
    val = next(fh)
    val = float(val)

    if val < min_val:
        raise RuntimeError(
            "Error, {} corr for {} < min val of {}".format(val, fname, min_val)
        )

    print("Corr of {} in {} is acceptable".format(val, fname))
    sys.exit(0)
