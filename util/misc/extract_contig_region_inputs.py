#!/usr/bin/env python

import sys, os, re
import argparse
import subprocess
import pysam


def main():

    parser = argparse.ArgumentParser(
        description="extract genome contig and bam for target contig of interest",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--genome_fa", type=str, required=True, help="genome fasta file"
    )

    parser.add_argument("--bam", type=str, required=True, help="bam aligned reads")

    parser.add_argument("--gtf", type=str, required=False, help="gtf annotation")

    parser.add_argument(
        "--region",
        type=str,
        required=True,
        help="contig region of interest. Formatted chr\\d[+-]?:\\d+-\\d+ ",
    )

    parser.add_argument(
        "--output_prefix",
        type=str,
        required=False,
        default=None,
        help="prefix for output files. Default: region string.",
    )

    args = parser.parse_args()

    genome_fa_filename = args.genome_fa
    bam_filename = args.bam
    gtf_filename = args.gtf
    region = args.region
    output_prefix = args.output_prefix

    if output_prefix is None:
        output_prefix = region

    m = re.match("([^\\:\\+\\-]+)([\\+\\-]?):(\\d+)-(\\d+)", region)
    if m is None:
        raise RuntimeError("Cannot parse region string: {}".format(region))

    chrom = m.group(1)
    strand = m.group(2)
    lend = int(m.group(3))
    rend = int(m.group(4))

    # extract contig
    genome_region_fa_filename = f"{output_prefix}.fa"
    cmd = f"samtools faidx {genome_fa_filename} {chrom}:{lend}-{rend}"
    fa_seq = subprocess.check_output(cmd, shell=True).decode()
    # remove header line
    fa_seq = "\n".join([f">{chrom}"] + fa_seq.split("\n")[1:])
    with open(genome_region_fa_filename, "wt") as ofh:
        print(fa_seq, file=ofh)

    cmd = f"samtools faidx {genome_region_fa_filename}"
    subprocess.check_call(cmd, shell=True)

    # extract bam
    bamreader = pysam.AlignmentFile(bam_filename, "rb")

    output_bam_filename = output_prefix + ".bam"
    bamwriter = pysam.AlignmentFile(output_bam_filename, "wb", template=bamreader)

    for read in bamreader.fetch(chrom, lend, rend):
        if strand != "":
            if (strand == "+" and not read.is_forward) or (
                strand == "-" and not read.is_reverse
            ):
                continue

        if read.reference_start >= lend:
            read.reference_start -= lend - 1
            bamwriter.write(read)

    bamwriter.close()

    cmd = f"samtools index {output_bam_filename}"
    subprocess.check_call(cmd, shell=True)

    # extract gtf
    if gtf_filename is not None:
        gtf_output_filename = output_prefix + ".gtf"
        with open(gtf_filename, "rt") as fh, open(gtf_output_filename, "wt") as ofh:

            for line in fh:
                vals = line.split("\t")
                if len(vals) < 8:
                    continue
                gtf_chrom, gtf_lend, gtf_rend, gtf_strand = (
                    vals[0],
                    int(vals[3]),
                    int(vals[4]),
                    vals[6],
                )
                if gtf_chrom != chrom:
                    continue
                if strand != "" and strand != gtf_strand:
                    continue

                if gtf_lend >= lend and gtf_rend <= gtf_rend:
                    gtf_lend -= lend - 1
                    gtf_rend -= lend - 1
                    vals[3] = str(gtf_lend)
                    vals[4] = str(gtf_rend)
                    print("\t".join(vals), file=ofh, end="")

    sys.exit(0)


if __name__ == "__main__":
    main()
