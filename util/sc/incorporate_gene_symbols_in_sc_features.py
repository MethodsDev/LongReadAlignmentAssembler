#!/usr/bin/env python3

import sys, os, re
import logging
import argparse
import gzip
import subprocess
import csv
from collections import defaultdict

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
        # Incorporate gene symbols into single cell feature names as gene_name^identifier"
        # 
        # If LRAA was run in quant-only mode, then add gene symbols like so:
        # 
        #   incorporate_gene_symbols_in_sc_features.py --ref_gtf GRCh38.gencode.annotation.gtf \\
        #                                              --sparseM_dirs dataset^gene-sparseM dataset^isoform-sparseM
        #
        # If LRAA was run in isoform-discovery mode, then first use gffcompare to assign LRAA isoforms to reference annotation isoforms like so:
        #
        #     gffcompare -r  GRCh38.gencode.annotation.gtf  LRAA_gtf
        #
        #  then incorporate gene names like so:
        #
        #     incorporate_gene_symbols_in_sc_features.py --ref_gtf GRCh38.gencode.annotation.gtf \\
        #                                                --sparseM_dirs dataset^gene-sparseM dataset^isoform-sparseM \\
        #                                                --LRAA_gtf LRAA.gtf \\
        #                                                --gffcompare_tracking gffcmp.tracking 
        #
        ##############################################################################################################

        """,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # required - for ref-quant-only or invovling novel isoform ID discovery (ref-guided or de novo LRAA)
    parser.add_argument(
        "--ref_gtf", type=str, required=True, help="reference annotation GTF file"
    )

    parser.add_argument(
        "--id_mappings",
        type=str,
        required=True,
        help="id mappings file: *.gene_transcript_splicehashcode.tsv file",
    )

    parser.add_argument(
        "--sparseM_dirs",
        type=str,
        required=True,
        help="sparse matrix directory",
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
    id_mappings_file = args.id_mappings
    LRAA_gtf_filename = args.LRAA_gtf
    sparseM_dirnames = args.sparseM_dirs
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

    id_mappings = parse_id_mappings(id_mappings_file)

    feature_ids_updated = dict()

    for sparseM_dirname in sparseM_dirnames:
        logger.info("-updating {}".format(sparseM_dirname))
        update_sparseM_feature_names(
            sparseM_dirname,
            gff_compare_LRAA_id_to_REF_id_mapping,
            ref_id_to_gene_name,
            id_mappings,
            feature_ids_updated,
        )

    revised_id_mappings = write_annotated_id_mappings_file(
        id_mappings_file, feature_ids_updated, ref_id_to_gene_name
    )

    if LRAA_gtf_filename is not None:
        logger.info(
            "-writing " + LRAA_gtf_filename + ".updated.gtf including gene names"
        )
        update_LRAA_gff_feature_ids(
            LRAA_gtf_filename,
            LRAA_gtf_filename + ".updated.gtf",
            feature_ids_updated,
            revised_id_mappings,
        )

    sys.exit(0)


def write_annotated_id_mappings_file(
    id_mappings_file, feature_ids_updated, ref_id_to_gene_name
):

    revised_id_mappings = dict()

    with open(id_mappings_file, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = [
            "gene_id",
            "transcript_id",
            "transcript_splice_hash_code",
            "new_gene_id",
            "new_transcript_id",
            "new_transcript_splice_hash_code",
        ]
        new_id_mappings_file = id_mappings_file + ".wAnnotIDs"
        with open(new_id_mappings_file, "wt") as ofh:
            writer = csv.DictWriter(
                ofh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()

            for row in reader:
                row["new_gene_id"] = feature_ids_updated.get(
                    row["gene_id"], row["gene_id"]
                )
                row["new_transcript_id"] = feature_ids_updated.get(
                    row["transcript_id"], row["transcript_id"]
                )
                row["new_transcript_splice_hash_code"] = feature_ids_updated.get(
                    row["transcript_splice_hash_code"],
                    row["transcript_splice_hash_code"],
                )
                writer.writerow(row)

                revised_id_mappings[row["transcript_id"]] = row

    return revised_id_mappings


def parse_id_mappings(id_mappings_file):

    id_mappings = dict()

    with open(id_mappings_file, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            transcript_id = row["transcript_id"]
            transcript_splice_hash_code = row["transcript_splice_hash_code"]
            id_mappings[transcript_id] = row
            id_mappings[transcript_splice_hash_code] = row

    return id_mappings


def update_LRAA_gff_feature_ids(
    LRAA_gtf_filename, new_LRAA_gtf_filename, feature_ids_updated, revised_id_mappings
):

    num_fields_updated = 0

    with open(LRAA_gtf_filename, "rt") as fh:
        with open(new_LRAA_gtf_filename, "wt") as ofh:

            for line in fh:
                line = line.rstrip()
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

                    if transcript_id in revised_id_mappings:
                        new_transcript_splice_hash_code = revised_id_mappings[
                            transcript_id
                        ]["new_transcript_splice_hash_code"]
                        line += f' splice_pattern "{new_transcript_splice_hash_code}";'

                print(line, file=ofh)

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


def update_sparseM_feature_names(
    sparseM_dirname,
    gff_compare_LRAA_id_to_REF_id_mapping,
    ref_id_to_gene_name,
    id_mappings,
    feature_ids_updated,
):

    features_file = os.path.join(sparseM_dirname, "features.tsv.gz")

    logger.info("-reassigning feature ids  {}".format(features_file))

    revised_features_file = features_file + ".revised"

    num_gene_names_added = 0

    with gzip.open(features_file, "rt") as fh:
        with open(revised_features_file, "wt") as ofh:
            for feature_id in fh:
                feature_id = feature_id.rstrip()
                gene_name = None

                # check for splice hash code
                transcript_splice_hash_code = None
                if (
                    feature_id in id_mappings
                    and feature_id
                    == id_mappings[feature_id]["transcript_splice_hash_code"]
                ):
                    transcript_splice_hash_code = feature_id
                    feature_id = id_mappings[feature_id]["transcript_id"]

                if feature_id in ref_id_to_gene_name:
                    gene_name = ref_id_to_gene_name[feature_id]

                # check gffcompare results
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
                    if transcript_splice_hash_code is not None:
                        # restore splice hash code
                        feature_id = transcript_splice_hash_code

                    new_feature_id = "^".join([gene_name, feature_id])
                    print(new_feature_id, file=ofh)
                    num_gene_names_added += 1
                    feature_ids_updated[feature_id] = new_feature_id
                else:
                    print(feature_id, file=ofh)  # no change

    if num_gene_names_added > 0:
        logger.info("- added {} gene names to feature ids".format(num_gene_names_added))
    else:
        logger.error(
            "-no gene names were assigned to feature ids... suggests a problem"
        )
        sys.exit(1)

    logger.info("-gzipping new features file: {}".format(revised_features_file))
    # gzip new features file.
    subprocess.check_call("gzip -f {}".format(revised_features_file), shell=True)
    revised_features_file += ".gz"

    os.rename(features_file, features_file + ".orig")
    os.rename(revised_features_file, features_file)

    return


def parse_GFFcompare_mappings(gffcompare_tracking_filename):

    logger.info("-parsing gffcompare output: {}".format(gffcompare_tracking_filename))

    gff_compare_mappings = dict()

    LRAA_gene_id_eq_assigned = set()
    LRAA_trans_id_eq_assigned = set()

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

            if lraa_gene_id not in LRAA_gene_id_eq_assigned:
                gff_compare_mappings[lraa_gene_id] = [ensg_id, enst_id]

            if lraa_trans_id not in LRAA_trans_id_eq_assigned:
                gff_compare_mappings[lraa_trans_id] = [ensg_id, enst_id]

            if compare_code == "=":
                LRAA_gene_id_eq_assigned.add(lraa_gene_id)
                LRAA_trans_id_eq_assigned.add(lraa_trans_id)

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
                'gene_id \\"([^\\"]+)\\";.* transcript_id \\"([^\\"]+)\\";', info
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
