#!/usr/bin/env python3

import sys, os

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../pylib"])
)

import SplicePatternCollapse
from SplicePatternCollapse import GeneConflictError, SiteSupportConflictError
import LRAA_Globals
import logging
import argparse

FORMAT = (
    "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s:\n\t%(message)s\n"
)

logger = logging.getLogger()
logging.basicConfig(format=FORMAT, level=logging.INFO)


def main():

    parser = argparse.ArgumentParser(
        description="Collapse an LRAA gtf by splice pattern, and report the TSS/PolyA "
        "sites its models describe. Isoforms are merged when they share BOTH the same "
        "gene_id AND the same intron (splice) pattern, so collapsed identifiers stay "
        "consistent with the expression-matrix collapse (build_LRAA_expr_matrices.py "
        "keys on gene_id^splice_hash). gene_symbol^ prefixes (if present) are retained. "
        "An intron pattern carried by more than one gene_id is REFUSED before anything "
        "is written, because identical splice patterns are one gene by definition and "
        "the collapsed transcript_id is minted from the pattern, so two such gene_ids "
        "would emit duplicate transcript_ids. Four files accompany the collapsed gtf: a "
        "merge report (each collapsed isoform and the isoforms merged into it), a "
        "gene-conflicts report (what the refusal names), and a TSS and a PolyA bed of "
        "the distinct boundary sites with their support -- and, for PolyA, the "
        "polyadenylation signal and internal-priming annotation. The beds describe the "
        "INPUT models: a merged model spans several terminal variants, so no single "
        "boundary count describes it. Merged models instead carry the same evidence as "
        "index-aligned per-site lists.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--gtf",
        type=str,
        required=True,
        help="LRAA gtf to collapse",
    )

    parser.add_argument(
        "--output_gtf",
        type=str,
        required=True,
        help="output collapsed gtf filename",
    )

    parser.add_argument(
        "--merge_report",
        type=str,
        default=None,
        help="TSV listing each collapsed isoform and the isoforms merged into it "
        "(default: <output_gtf>.isoform_merge_report.tsv)",
    )

    parser.add_argument(
        "--gene_conflicts_report",
        type=str,
        default=None,
        help="TSV listing intron patterns carried by more than one gene_id. Written, "
        "and the run then refused, when any exist (default: "
        "<output_gtf>.gene_conflicts.tsv)",
    )

    parser.add_argument(
        "--TSS_bed",
        type=str,
        default=None,
        help="bed of distinct TSS sites with their support "
        "(default: <output_gtf>.TSS.bed)",
    )

    parser.add_argument(
        "--PolyA_bed",
        type=str,
        default=None,
        help="bed of distinct PolyA sites with their support, polyadenylation signal, "
        "and internal-priming annotation (default: <output_gtf>.PolyA.bed)",
    )

    parser.add_argument(
        "--nonfatal_gene_conflicts",
        action="store_true",
        default=False,
        help="on a cross-gene splice pattern, report it and continue rather than "
        "refusing: the conflicts report and both bed files are still written, the "
        "collapsed gtf and merge report are not, and the exit status is 0. For callers "
        "that cannot afford a nonzero exit -- a failed WDL task publishes nothing, so "
        "refusing there would discard a completed run's outputs along with the report "
        "needed to diagnose the conflict",
    )

    parser.add_argument(
        "--debug",
        "-d",
        action="store_true",
        default=False,
        help="debug mode, more verbose",
    )

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
        LRAA_Globals.DEBUG = True

    logger.info(f"-capturing input transcripts from gtf {args.gtf}")

    try:
        written = SplicePatternCollapse.collapse_gtf(
            args.gtf,
            args.output_gtf,
            merge_report=args.merge_report,
            gene_conflicts_report=args.gene_conflicts_report,
            TSS_bed=args.TSS_bed,
            PolyA_bed=args.PolyA_bed,
            nonfatal_gene_conflicts=args.nonfatal_gene_conflicts,
        )
    except GeneConflictError as conflict:
        exit("Error, {}".format(conflict))
    except SiteSupportConflictError as conflict:
        exit("Error, {}".format(conflict))

    if "collapsed_gtf" not in written:
        logger.error(
            "-cross-gene splice patterns found; see %s. No collapsed gtf was written.",
            written["gene_conflicts_report"],
        )

    for label, path in written.items():
        logger.info(f"-wrote {label}: {path}")

    logger.info("Done.")

    sys.exit(0)


if __name__ == "__main__":
    main()
