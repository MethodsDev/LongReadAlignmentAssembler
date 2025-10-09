#!/usr/bin/env python3
import os, sys
import argparse
import logging

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../pylib"])
)

from Pipeliner import Pipeliner, Command

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="Run ORF prediction pipeline with Pipeliner"
    )
    parser.add_argument("--genome", required=True, help="Genome fasta file")
    parser.add_argument("--gtf", required=True, help="Input GTF file")
    parser.add_argument("--output_prefix", required=True, help="Prefix for outputs")
    parser.add_argument(
        "--td_dir", required=True, help="Path to TransDecoder executables"
    )
    parser.add_argument(
        "--complete_orfs_only",
        action="store_true",
        default=False,
        help="require complete ORFs starting with Met",
    )
    parser.add_argument(
        "--single_best_only",
        action="store_true",
        default=False,
        help="single best orf per transcript only",
    )

    parser.add_argument(
        "-m",
        "--min_prot_length",
        default=100,
        help="minimum protein length to be considered",
    )

    args = parser.parse_args()

    utildir = args.td_dir + "/util"

    pipeliner = Pipeliner("__chckpts")

    # Step 1: generate alignment gff3 formatted output
    cmd = f"{utildir}/gtf_to_alignment_gff3.pl {args.gtf} > {args.output_prefix}.input_converted.gff3"
    pipeliner.add_commands([Command(cmd, "gtf_to_alignment_gff3.ok")])

    # Step 2: generate cdna fasta
    cmd = f"{utildir}/gtf_genome_to_cdna_fasta.pl {args.gtf} {args.genome} > {args.output_prefix}.cdna.fasta"
    pipeliner.add_commands([Command(cmd, "gtf_genome_to_cdna_fasta.ok")])

    # Step 3: extract the long ORFs
    cmd = f"{args.td_dir}/TransDecoder.LongOrfs -t {args.output_prefix}.cdna.fasta -S -m {args.min_prot_length}"
    if args.complete_orfs_only:
        cmd += " --complete_orfs_only "
    pipeliner.add_commands([Command(cmd, "transdecoder_longorfs.ok")])

    # Step 4: predict likely ORFs
    cmd = f"{args.td_dir}/TransDecoder.Predict -t {args.output_prefix}.cdna.fasta --no_refine_starts "
    if args.single_best_only:
        cmd += " --single_best_only"

    pipeliner.add_commands([Command(cmd, "transdecoder_predict.ok")])

    # Step 5: convert to genome coordinates
    cmd = (
        f"{utildir}/cdna_alignment_orf_to_genome_orf.pl "
        f"{args.output_prefix}.cdna.fasta.transdecoder.gff3 "
        f"{args.output_prefix}.input_converted.gff3 "
        f"{args.output_prefix}.cdna.fasta "
        f"> {args.output_prefix}.transcripts.fasta.transdecoder.genome.gff3"
    )
    pipeliner.add_commands([Command(cmd, "cdna_orf_to_genome_orf.ok")])

    # Step 6a: gtf to bed
    cmd = (
        f"{utildir}/gtf_to_bed.pl {args.gtf} "
        f"| sort -k1,1 -k2,2g -k3,3g > {args.output_prefix}.input_converted.bed"
    )
    pipeliner.add_commands([Command(cmd, "gtf_to_bed.ok")])

    # Step 6b: bgzip + tabix gtf bed
    cmd = (
        f"bgzip -f {args.output_prefix}.input_converted.bed && "
        f"tabix -f {args.output_prefix}.input_converted.bed.gz"
    )
    pipeliner.add_commands([Command(cmd, "index_input_bed.ok")])

    # Step 6c: genome-based gff3 to bed
    cmd = (
        f"{utildir}/gff3_file_to_bed.pl "
        f"{args.output_prefix}.transcripts.fasta.transdecoder.genome.gff3 "
        f"| sort -k1,1 -k2,2g -k3,3g "
        f"> {args.output_prefix}.transcripts.fasta.transdecoder.genome.bed"
    )
    pipeliner.add_commands([Command(cmd, "gff3_to_bed.ok")])

    # Step 6d: bgzip + tabix genome bed
    cmd = (
        f"bgzip -f {args.output_prefix}.transcripts.fasta.transdecoder.genome.bed && "
        f"tabix -f {args.output_prefix}.transcripts.fasta.transdecoder.genome.bed.gz"
    )
    pipeliner.add_commands([Command(cmd, "index_genome_bed.ok")])

    # Run all steps
    pipeliner.run()


if __name__ == "__main__":
    sys.exit(main())
