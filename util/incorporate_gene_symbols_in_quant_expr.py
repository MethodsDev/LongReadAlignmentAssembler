#!/usr/bin/env python3

import sys, os, re
import logging
import argparse
import gzip
import subprocess
import csv

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(
        description="""
        
        #############################################################################################################
        #
        # Incorporate gene symbols into gene_id and transcript_id as gene_name^identifier"
        # 
        # If LRAA was run in quant-only mode, then add gene symbols like so:
        # 
        #   incorporate_gene_symbols_in_sc_features.py --ref_gtf GRCh38.gencode.annotation.gtf --quant_expr LRAA.quant.A.expr [LRAA.quant.B.expr ...]
        #
        # If LRAA was run in isoform-discovery mode, then first use gffcompare to assign LRAA isoforms to reference annotation isoforms like so:
        #
        #     gffcompare -r  GRCh38.gencode.annotation.gtf  LRAA_gtf
        #
        #  then incorporate gene names like so:
        #
        #     incorporate_gene_symbols_in_sc_features.py --ref_gtf GRCh38.gencode.annotation.gtf \\
        #                                                --LRAA_gtf LRAA.gtf \\
        #                                                --gffcompare_tracking gffcmp.tracking \\
        #                                                 --quant_expr LRAA.quant.A.expr [LRAA.quant.B.expr ...] 
        #
        ##############################################################################################################


        A LRAA.quant.expr.wAnnot file will be generated for every input LRAA.quant.expr file.

        """,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # required - for ref-quant-only or invovling novel isoform ID discovery (ref-guided or de novo LRAA)
    parser.add_argument(
        "--ref_gtf", type=str, required=True, help="reference annotation GTF file"
    )

    parser.add_argument(
        "--quant_expr",
        type=str,
        required=True,
        help="LRAA quant expr file(s)y",
        nargs="+",
    )

    # optional - for ref-guided
    parser.add_argument(
        "--LRAA_gtf",
        type=str,
        required=False,
        default=None,
        help="LRAA gtf file",
    )

    parser.add_argument(
        "--gffcompare_tracking",
        type=str,
        required=False,
        default=None,
        help="provide the tracking file from running: gffcompare -r ref_gtf LRAA_gtf",
    )

    args = parser.parse_args()

    ref_gtf_filename = args.ref_gtf
    LRAA_gtf_filename = args.LRAA_gtf
    quant_expr_filenames = args.quant_expr
    gffcompare_tracking_filename = args.gffcompare_tracking

    # begin

    ref_id_to_gene_name = get_ref_gene_names(ref_gtf_filename)

    # print(str(ref_id_to_gene_name))

    gff_compare_LRAA_id_to_REF_id_mapping = None

    if gffcompare_tracking_filename is not None:
        gff_compare_LRAA_id_to_REF_id_mapping = parse_GFFcompare_mappings(
            gffcompare_tracking_filename
        )

    # print(str(gff_compare_LRAA_id_to_REF_id_mapping))

    feature_ids_updated = dict()

    for quant_expr in quant_expr_filenames:
        logger.info("-updating {}".format(quant_expr))
        update_quant_expr_feature_names(
            quant_expr,
            gff_compare_LRAA_id_to_REF_id_mapping,
            ref_id_to_gene_name,
            feature_ids_updated,
        )

    if LRAA_gtf_filename is not None:
        logger.info(
            "-writing " + LRAA_gtf_filename + ".updated.gtf including gene names"
        )
        update_LRAA_gff_feature_ids(
            LRAA_gtf_filename, LRAA_gtf_filename + ".wAnnot.gtf", feature_ids_updated
        )

    sys.exit(0)


def update_LRAA_gff_feature_ids(
    LRAA_gtf_filename, new_LRAA_gtf_filename, feature_ids_updated
):

    num_fields_updated = 0

    with open(LRAA_gtf_filename, "rt") as fh:
        with open(new_LRAA_gtf_filename, "wt") as ofh:

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

                    if gene_id in feature_ids_updated:
                        new_gene_id = feature_ids_updated[gene_id]
                        line = line.replace(gene_id, new_gene_id)
                        num_fields_updated += 1

                    if transcript_id in feature_ids_updated:
                        new_transcript_id = feature_ids_updated[transcript_id]
                        line = line.replace(transcript_id, new_transcript_id)
                        num_fields_updated += 1

                print(line, file=ofh, end="")

    if num_fields_updated > 0:
        logger.info(
            "-updated {} feature ids in gtf file {}, generated output file {}".format(
                num_fields_updated, LRAA_gtf_filename, new_LRAA_gtf_filename
            )
        )

    else:
        logger.error(
            "-no feature ids were updated from gtf file {}  - something wrong here...".format(
                LRAA_gtf_filename
            )
        )
        sys.exit(1)

    return


def update_quant_expr_feature_names(
    quant_expr_filename,
    gff_compare_LRAA_id_to_REF_id_mapping,
    ref_id_to_gene_name,
    feature_ids_updated,
):

    logger.info("-reassigning feature ids  {}".format(quant_expr_filename))

    annot_quant_expr_file = quant_expr_filename + ".wAnnot"

    num_gene_names_added = 0

    feature_ids_udpated = dict()

    with open(quant_expr_filename, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames)
        with open(annot_quant_expr_file, "wt") as ofh:
            writer = csv.DictWriter(
                ofh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            for row in reader:
                gene_name = None
                feature_id = None
                if "gene_id" in fieldnames:
                    feature_id = row["gene_id"]
                elif "transcript_id" in fieldnames:
                    feature_id = row["transcript_id"]

                if feature_id in ref_id_to_gene_name:
                    gene_name = ref_id_to_gene_name[feature_id]

                elif (
                    gff_compare_LRAA_id_to_REF_id_mapping is not None
                    and feature_id in gff_compare_LRAA_id_to_REF_id_mapping
                ):
                    ensg_id, enst_id = gff_compare_LRAA_id_to_REF_id_mapping[feature_id]

                    if ensg_id in ref_id_to_gene_name:
                        gene_name = ref_id_to_gene_name[ensg_id]
                    elif enst_id in ref_id_to_gene_name:
                        gene_name = ref_id_to_gene_name[enst_id]

                if gene_name is not None:
                    for feature_id_type in ("gene_id", "transcript_id"):
                        if feature_id_type in fieldnames:
                            new_feature_id = "^".join([gene_name, row[feature_id_type]])
                            feature_ids_updated[row[feature_id_type]] = new_feature_id
                            row[feature_id_type] = new_feature_id
                            num_gene_names_added += 1

                writer.writerow(row)

    if num_gene_names_added > 0:
        logger.info("- added {} gene names to feature ids".format(num_gene_names_added))
    else:
        logger.error(
            "-no gene names were assigned to feature ids... suggests a problem"
        )
        sys.exit(1)

    return


def parse_GFFcompare_mappings(gffcompare_tracking_filename):

    logger.info("-parsing gffcompare output: {}".format(gffcompare_tracking_filename))

    gff_compare_mappings = dict()

    with open(gffcompare_tracking_filename, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            tcons, xloc, ref_info, compare_code, lraa_info = line.split("\t")

            if ref_info == "-":
                continue

            ensg_id, enst_id = ref_info.split("|")

            lraa_info = ":".join(lraa_info.split(":")[1:])  # get rid of first q1 token

            lraa_vals = lraa_info.split("|")
            lraa_gene_id = lraa_vals[0]
            lraa_trans_id = lraa_vals[1]

            gff_compare_mappings[lraa_gene_id] = [ensg_id, enst_id]
            gff_compare_mappings[lraa_trans_id] = [ensg_id, enst_id]

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
                'gene_id \\"([^\\"]+)\\"; transcript_id \\"([^\\"]+)\\";', info
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
