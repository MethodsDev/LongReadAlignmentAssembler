#!/usr/bin/env python3
import os, sys
import argparse
import logging
import shlex

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
        description="Run ORF prediction pipeline with Pipeliner. "
        "Uses the TransDecoder v6 pipeline wrapper and optionally performs "
        "homology search during ORF prediction."
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
        "--no_refine_starts",
        action="store_true",
        default=False,
        help="Disable TransDecoder v6 start-codon refinement.",
    )

    parser.add_argument(
        "-m",
        "--min_prot_length",
        default=100,
        help="minimum protein length to be considered",
    )

    parser.add_argument(
        "--blast_search_pep",
        "--blastp_db",
        dest="blast_search_pep",
        type=str,
        default=None,
        help="Protein FASTA file for TransDecoder v6 homology search. "
        "TransDecoder now builds the search database automatically.",
    )

    parser.add_argument(
        "--blast_evalue",
        "--blastp_evalue",
        dest="blast_evalue",
        type=float,
        default=1e-5,
        help="E-value threshold for blastp/diamond search (default: 1e-5)",
    )

    parser.add_argument(
        "--blastp_max_target_seqs",
        type=int,
        default=1,
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--blast_threads",
        "--blastp_num_threads",
        dest="blast_threads",
        type=int,
        default=1,
        help="Number of threads for blastp/diamond (default: 1)",
    )

    parser.add_argument(
        "--blast_tool",
        "--search_method",
        dest="blast_tool",
        type=str,
        choices=["diamond", "blastp"],
        default="diamond",
        help="Homology search method to use: diamond (default, faster) or blastp (slower, NCBI blast+)",
    )

    parser.add_argument(
        "--pfam_search_db",
        "--pfam-search-db",
        dest="pfam_search_db",
        type=str,
        default=None,
        help="Pfam HMM database for TransDecoder v6 hmmsearch-based ORF retention.",
    )

    args = parser.parse_args()

    utildir = args.td_dir + "/util"
    transdecoder_exe = os.path.join(args.td_dir, "TransDecoder")

    pipeliner = Pipeliner("__chckpts")

    if not os.path.exists(transdecoder_exe):
        raise FileNotFoundError(
            f"TransDecoder v6 wrapper not found at expected path: {transdecoder_exe}"
        )

    if args.blastp_max_target_seqs != 1:
        logger.warning(
            "--blastp_max_target_seqs=%s requested, but TransDecoder v6 uses the "
            "single top target internally; proceeding with 1",
            args.blastp_max_target_seqs,
        )

    # Step 1: run the full TransDecoder v6 pipeline wrapper
    cmd_parts = [
        shlex.quote(transdecoder_exe),
        "--genome",
        shlex.quote(args.genome),
        "--gtf",
        shlex.quote(args.gtf),
        "-t",
        shlex.quote(f"{args.output_prefix}.cdna.fasta"),
        "-O",
        shlex.quote("."),
        "-S",
        "-m",
        shlex.quote(str(args.min_prot_length)),
    ]

    if args.no_refine_starts:
        cmd_parts.append("--no_refine_starts")

    if args.complete_orfs_only:
        cmd_parts.append("--complete_orfs_only")
    if args.single_best_only:
        cmd_parts.append("--single_best_only")

    blastp_outfile = None
    if args.blast_search_pep:
        if not os.path.exists(args.blast_search_pep):
            raise FileNotFoundError(
                f"Protein FASTA for homology search not found: {args.blast_search_pep}"
            )
        if args.blast_search_pep.endswith(".dmnd"):
            raise ValueError(
                "TransDecoder v6 expects --blast_search_pep to reference a protein FASTA, "
                "not a prebuilt Diamond database (*.dmnd)."
            )

        blastp_outfile = f"{args.output_prefix}.blastp.outfmt6"
        cmd_parts.extend(
            [
                "--blast_search_pep",
                shlex.quote(args.blast_search_pep),
                "--blast_tool",
                shlex.quote(args.blast_tool),
                "--blast_evalue",
                shlex.quote(str(args.blast_evalue)),
                "--blast_threads",
                shlex.quote(str(args.blast_threads)),
            ]
        )

    pfam_outfile = None
    if args.pfam_search_db:
        if not os.path.exists(args.pfam_search_db):
            raise FileNotFoundError(
                f"Pfam HMM database for ORF retention not found: {args.pfam_search_db}"
            )

        pfam_outfile = f"{args.output_prefix}.pfam.domtblout"
        cmd_parts.extend(
            [
                "--pfam-search-db",
                shlex.quote(args.pfam_search_db),
            ]
        )

    cmd = " ".join(cmd_parts)
    pipeliner.add_commands([Command(cmd, "transdecoder_pipeline.ok")])

    # Step 2: keep blast output filename stable for downstream consumers
    if blastp_outfile:
        td_blastp_outfile = (
            f"{args.output_prefix}.cdna.fasta.transdecoder_dir/blastp.outfmt6"
        )
        cmd = f"cp {shlex.quote(td_blastp_outfile)} {shlex.quote(blastp_outfile)}"
        pipeliner.add_commands([Command(cmd, "copy_blastp_output.ok")])

    if pfam_outfile:
        td_pfam_outfile = f"{args.output_prefix}.cdna.fasta.transdecoder_dir/pfam.domtblout"
        cmd = f"cp {shlex.quote(td_pfam_outfile)} {shlex.quote(pfam_outfile)}"
        pipeliner.add_commands([Command(cmd, "copy_pfam_output.ok")])

    # Step 3a: gtf to bed
    cmd = (
        f"{utildir}/gtf_to_bed.pl {args.gtf} "
        f"| sort -k1,1 -k2,2g -k3,3g > {args.output_prefix}.input_converted.bed"
    )
    pipeliner.add_commands([Command(cmd, "gtf_to_bed.ok")])

    # Step 3b: bgzip + tabix gtf bed
    cmd = (
        f"bgzip -f {args.output_prefix}.input_converted.bed && "
        f"tabix -f {args.output_prefix}.input_converted.bed.gz"
    )
    pipeliner.add_commands([Command(cmd, "index_input_bed.ok")])

    # Step 3c: transcript-coordinate gff3 to bed
    cmd = (
        f"{utildir}/gff3_file_to_bed.pl "
        f"{args.output_prefix}.cdna.fasta.transdecoder.gff3 "
        f"| sort -k1,1 -k2,2g -k3,3g "
        f"> {args.output_prefix}.cdna.fasta.transdecoder.bed"
    )
    pipeliner.add_commands([Command(cmd, "transdecoder_gff3_to_bed.ok")])

    # Step 3d: genome-based gff3 to bed
    cmd = (
        f"{utildir}/gff3_file_to_bed.pl "
        f"{args.output_prefix}.cdna.fasta.transdecoder.genome.gff3 "
        f"| sort -k1,1 -k2,2g -k3,3g "
        f"> {args.output_prefix}.cdna.fasta.transdecoder.genome.bed"
    )
    pipeliner.add_commands([Command(cmd, "gff3_to_bed.ok")])

    # Step 3e: bgzip + tabix genome bed
    cmd = (
        f"bgzip -f {args.output_prefix}.cdna.fasta.transdecoder.genome.bed && "
        f"tabix -f {args.output_prefix}.cdna.fasta.transdecoder.genome.bed.gz"
    )
    pipeliner.add_commands([Command(cmd, "index_genome_bed.ok")])

    # Run all steps
    pipeliner.run()


if __name__ == "__main__":
    sys.exit(main())
