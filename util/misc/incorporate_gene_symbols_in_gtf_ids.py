#!/usr/bin/env python3

import sys, os, re
import logging
import argparse
import gzip
import subprocess

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

    gff_compare_target_id_to_REF_id_mapping = parse_GFFcompare_mappings(
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

                gene_name = None
                if gene_id in gff_compare_target_id_to_REF_id_mapping:
                    ref_gene_id = gff_compare_target_id_to_REF_id_mapping[gene_id][0]
                    if ref_gene_id in ref_id_to_gene_name:
                        gene_name = ref_id_to_gene_name[ref_gene_id]

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


def parse_GFFcompare_mappings(gffcompare_tracking_filename):

    logger.info("-parsing gffcompare output: {}".format(gffcompare_tracking_filename))

    gff_compare_mappings = dict()

    with open(gffcompare_tracking_filename, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            tcons, xloc, ref_info, compare_code, target_info = line.split("\t")

            if ref_info == "-":
                continue

            ensg_id, enst_id = ref_info.split("|")

            target_info = ":".join(
                target_info.split(":")[1:]
            )  # get rid of first q1 token

            target_vals = target_info.split("|")
            target_gene_id = target_vals[0]
            target_trans_id = target_vals[1]

            gff_compare_mappings[target_gene_id] = [ensg_id, enst_id]
            gff_compare_mappings[target_trans_id] = [ensg_id, enst_id]

    return gff_compare_mappings


def get_ref_gene_names(ref_gtf):

    logger.info(
        "-extracting gene_names and identifiers from reference gtf: {}".format(ref_gtf)
    )

    ref_id_to_gene_name = dict()

    with open(ref_gtf, "rt") as fh:
        for line in fh:
            vals = line.split("\t")
            if len(vals) < 8:
                continue
            info = vals[8]

            if vals[2] != "transcript":
                continue

            m = re.search(
                'gene_id \\"([^\\"]+)\\";.*transcript_id \\"([^\\"]+)\\";', info
            )
            if m:
                gene_id = m.group(1)
                transcript_id = m.group(2)
                gene_name = gene_id

                m2 = re.search(' gene_name "([^\\"]+)\\";', info)
                if m2:
                    gene_name = m2.group(1)

                    ref_id_to_gene_name[transcript_id] = gene_name
                    ref_id_to_gene_name[gene_id] = gene_name

    return ref_id_to_gene_name


if __name__ == "__main__":
    main()
