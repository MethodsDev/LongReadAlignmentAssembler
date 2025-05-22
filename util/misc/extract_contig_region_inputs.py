#!/usr/bin/env python

import sys, os, re
import argparse
import subprocess

def main():

    parser = argparse.ArgumentParser(description="extract genome contig and bam for target contig of interest", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--genome_fa", type=str, required=True, help="genome fasta file")

    parser.add_argument("--bam", type=str, required=True, help="bam aligned reads")

    parser.add_argument("--contig", type=str, required=True, help="contig of interest")
    
    args = parser.parse_args()

    genome_fa_filename = args.genome_fa
    bam_filename = args.bam
    contig = args.contig


    # extract contig
    cmd = f"samtools faidx {genome_fa_filename} {contig} > {contig}.genome.fa"
    subprocess.check_call(cmd, shell=True)

    cmd = f"samtools faidx {contig}.genome.fa"
    subprocess.check_call(cmd, shell=True)
    
    # extract bam
    cmd = f"samtools view {bam_filename} {contig} -bo {contig}.bam"
    subprocess.check_call(cmd, shell=True)

    cmd = f"samtools index {contig}.bam"
    subprocess.check_call(cmd, shell=True)



    sys.exit(0)


if __name__=='__main__':
    main()



    
