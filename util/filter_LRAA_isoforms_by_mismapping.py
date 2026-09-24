#!/usr/bin/env python3
# encoding: utf-8

"""Standalone whole-genome alignment-mismapping filter.

Removes isoforms that are alignment/strand-mismapping artifacts of a
much-higher-expressed transcript (see pylib/AlignmentMismappingFilter.py). Runs
on a MERGED, genome-wide LRAA gtf + quant.expr (never a per-chunk slice) and
writes a filtered gtf, a filtered+renormalized quant.expr, and a removal log.

This is the standalone entry point used by the WDL; the same pylib module is
invoked in-process by the LRAA driver as a post-merge step.
"""

import sys
import os
import argparse
import logging
import tempfile

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../pylib"])
)

import LRAA_Globals
import AlignmentMismappingFilter

FORMAT = (
    "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s:\n\t%(message)s\n"
)
logger = logging.getLogger()
logging.basicConfig(format=FORMAT, level=logging.INFO)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--gtf", required=True, help="merged whole-genome LRAA gtf")
    parser.add_argument(
        "--quant_expr", required=True, help="merged whole-genome LRAA quant.expr"
    )
    parser.add_argument(
        "--genome", required=True, help="reference genome FASTA (for cDNA extraction)"
    )
    parser.add_argument(
        "--output_gtf", required=True, help="filtered gtf output filename"
    )
    parser.add_argument(
        "--output_quant_expr",
        required=True,
        help="filtered + renormalized quant.expr output filename",
    )
    parser.add_argument(
        "--output_log",
        default=None,
        help="removal log (default: <output_gtf>.mismapping_filter.log)",
    )
    parser.add_argument("--work_dir", default=None, help="scratch dir for cDNA fasta")
    parser.add_argument("--threads", type=int, default=4, help="minimap2 threads")

    # Threshold overrides (default from LRAA_Globals.config).
    parser.add_argument("--min_seq_identity", type=float, default=None)
    parser.add_argument("--min_seq_coverage", type=float, default=None)
    parser.add_argument("--max_expr_fraction", type=float, default=None)
    parser.add_argument("--junction_tolerance", type=int, default=None)
    parser.add_argument("--min_base_overlap", type=float, default=None)
    parser.add_argument("--debug", "-d", action="store_true")

    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)
        LRAA_Globals.DEBUG = True

    overrides = {
        "mismap_min_seq_identity": args.min_seq_identity,
        "mismap_min_seq_coverage": args.min_seq_coverage,
        "mismap_max_expr_fraction": args.max_expr_fraction,
        "mismap_junction_tolerance": args.junction_tolerance,
        "mismap_min_base_overlap": args.min_base_overlap,
    }
    for k, v in overrides.items():
        if v is not None:
            LRAA_Globals.config[k] = v

    log_out = args.output_log or (args.output_gtf + ".mismapping_filter.log")

    workdir = args.work_dir
    tmp = None
    if workdir is None:
        tmp = tempfile.TemporaryDirectory(prefix="lraa_mismap_")
        workdir = tmp.name

    try:
        drop_set = AlignmentMismappingFilter.run_mismapping_filter(
            gtf_in=args.gtf,
            quant_in=args.quant_expr,
            genome_fasta=args.genome,
            gtf_out=args.output_gtf,
            quant_out=args.output_quant_expr,
            log_out=log_out,
            workdir=workdir,
            threads=args.threads,
        )
    finally:
        if tmp is not None:
            tmp.cleanup()

    logger.info(
        "Done. Removed %d models. Filtered gtf: %s ; quant: %s ; log: %s",
        len(drop_set),
        args.output_gtf,
        args.output_quant_expr,
        log_out,
    )


if __name__ == "__main__":
    main()
