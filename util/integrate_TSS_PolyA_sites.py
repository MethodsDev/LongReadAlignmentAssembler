#!/usr/bin/env python3

import sys, os

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../pylib"])
)

import SiteIntegration
import logging
import argparse

FORMAT = (
    "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s:\n\t%(message)s\n"
)

logger = logging.getLogger()
logging.basicConfig(format=FORMAT, level=logging.INFO)


def main():

    parser = argparse.ArgumentParser(
        description="Integrate the cluster-guided TSS/PolyA site beds with the basic "
        "(initial-catalog) ones. Every cluster-guided site is kept; a basic site is "
        "added only when it lies farther than the window from every "
        "cluster-guided site of the same type, contig and strand -- the distance over "
        "which LRAA's site definition absorbs read ends into one site, so no two sites "
        "of one run are closer. Inputs are the site beds "
        "collapse_LRAA_GTF_by_splice_pattern.py writes; outputs keep their columns and "
        "append 'source' (cluster_guided | basic). Writes "
        "<prefix>.integrated.TSS.bed, <prefix>.integrated.PolyA.bed and "
        "<prefix>.integrated_sites.summary.tsv.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--primary_TSS_bed", required=True,
                        help="cluster-guided TSS bed (all sites kept)")
    parser.add_argument("--primary_PolyA_bed", required=True,
                        help="cluster-guided PolyA bed (all sites kept)")
    parser.add_argument("--supplement_TSS_bed", required=True,
                        help="basic TSS bed (sites outside the window are added)")
    parser.add_argument("--supplement_PolyA_bed", required=True,
                        help="basic PolyA bed (sites outside the window are added)")
    parser.add_argument("--output_prefix", required=True, help="output file prefix")
    parser.add_argument(
        "--TSS_window", type=int, default=None,
        help="site-aggregation distance for TSS; default: LRAA config "
        "max_dist_between_alt_TSS_sites ({})".format(SiteIntegration.default_window("TSS")),
    )
    parser.add_argument(
        "--PolyA_window", type=int, default=None,
        help="site-aggregation distance for PolyA; default: LRAA config "
        "max_dist_between_alt_polyA_sites ({})".format(SiteIntegration.default_window("PolyA")),
    )

    args = parser.parse_args()

    all_counts = []
    for site_type, primary, supplement, window in (
        ("TSS", args.primary_TSS_bed, args.supplement_TSS_bed, args.TSS_window),
        ("PolyA", args.primary_PolyA_bed, args.supplement_PolyA_bed, args.PolyA_window),
    ):
        out = "{}.integrated.{}.bed".format(args.output_prefix, site_type)
        counts = SiteIntegration.integrate_site_beds(
            primary, supplement, site_type, out, window=window
        )
        logger.info(
            "%s: %d cluster_guided + %d basic supplement (%d basic dropped within %d nt of a "
            "cluster_guided site, %d collapsed as near-duplicates) -> %s",
            site_type, counts["primary"], counts["supplement_kept"],
            counts["supplement_dropped"], counts["window"], counts["supplement_collapsed"], out,
        )
        all_counts.append(counts)

    SiteIntegration.write_summary(
        "{}.integrated_sites.summary.tsv".format(args.output_prefix), all_counts
    )


if __name__ == "__main__":
    main()
