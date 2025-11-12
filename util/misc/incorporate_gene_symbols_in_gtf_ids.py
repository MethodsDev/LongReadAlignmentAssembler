#!/usr/bin/env python3

import sys, os, re
import logging
import argparse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = (
    os.path.dirname(os.path.dirname(SCRIPT_DIR))
    if os.path.basename(SCRIPT_DIR) == "misc"
    else os.path.dirname(SCRIPT_DIR)
)
PYLIB_DIR = os.path.join(ROOT_DIR, "pylib")
if PYLIB_DIR not in sys.path:
    sys.path.insert(0, PYLIB_DIR)

from gene_symbol_utils import (
    get_ref_gene_names,
    parse_gffcompare_mappings,
    resolve_gene_symbol,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description="incorporate gene symbols into single cell feature names as gene_name^identifier",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--target_gtf",
        type=str,
        required=True,
        help="target gtf file",
    )

    parser.add_argument(
        "--ref_gtf", type=str, required=True, help="reference annotation GTF file"
    )

    parser.add_argument(
        "--gffcompare_tracking",
        type=str,
        required=True,
        help="provide the tracking file from running: gffcompare -r ref_gtf target_gtf",
    )

    args = parser.parse_args()

    ref_gtf_filename = args.ref_gtf
    target_gtf_filename = args.target_gtf
    gffcompare_tracking_filename = args.gffcompare_tracking

    # begin

    ref_id_to_gene_name = get_ref_gene_names(ref_gtf_filename)

    gff_compare_target_id_to_REF_id_mapping = parse_gffcompare_mappings(
        gffcompare_tracking_filename
    )

    update_gtf_feature_ids(
        target_gtf_filename,
        ref_id_to_gene_name,
        gff_compare_target_id_to_REF_id_mapping,
    )

    sys.exit(0)


def update_gtf_feature_ids(
    target_gtf_filename, ref_id_to_gene_name, gff_compare_target_id_to_REF_id_mapping
):

    num_fields_updated = 0

    with open(target_gtf_filename, "rt") as fh:

        for line in fh:
            vals = line.split("\t")
            if len(vals) < 8:
                continue
            info = vals[8]
            m = re.search(
                'gene_id \\"([^\\"]+)\\"; transcript_id \\"([^\\"]+)\\";', info
            )
            if m:
                gene_id = m.group(1)
                transcript_id = m.group(2)

                gene_name = resolve_gene_symbol(
                    transcript_id,
                    ref_id_to_gene_name,
                    gff_compare_target_id_to_REF_id_mapping,
                    prefer_transcript_first=True,
                )

                if gene_name is None:
                    gene_name = resolve_gene_symbol(
                        gene_id,
                        ref_id_to_gene_name,
                        gff_compare_target_id_to_REF_id_mapping,
                        prefer_transcript_first=False,
                    )

                if gene_name is not None:
                    new_gene_id = gene_name + "^" + gene_id
                    line = line.replace(gene_id, new_gene_id)

                    new_transcript_id = gene_name + "^" + transcript_id
                    line = line.replace(transcript_id, new_transcript_id)

                    num_fields_updated += 1

            print(line, end="")

    if num_fields_updated > 0:
        logger.info(
            "-updated {} feature ids in gtf file {}".format(
                num_fields_updated, target_gtf_filename
            )
        )

    else:
        logger.error(
            "-no feature ids were updated from gtf file {}  - something wrong here...".format(
                target_gtf_filename
            )
        )
        sys.exit(1)

    return


if __name__ == "__main__":
    main()
