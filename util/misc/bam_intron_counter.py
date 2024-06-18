#!/usr/bin/env python3

import sys, os, re
import pysam
from collections import defaultdict

MIN_MAPPING_QUALITY = 60

OK_SPLICES = ("GT--AG:+", "GC--AG:+", "AT--AC:+", # forward strand
              "CT--AC:-", "CT--GC:-", "GT--AT:-" # reverse strand
              )

def main():

    usage = "\n\nusage: {} alignments.bam genome.fasta\n\n".format(sys.argv[0])

    if len(sys.argv) < 3:
        exit(usage)

    bam_filename = sys.argv[1]
    genome_fasta_filename = sys.argv[2]
    
    intron_counter = defaultdict( lambda: defaultdict(int) )
    
    fasta_reader = pysam.FastaFile(genome_fasta_filename)
    
    bam_reader = pysam.AlignmentFile(bam_filename, "rb")
    for read in bam_reader.fetch():
        if read.mapping_quality < MIN_MAPPING_QUALITY:
            continue

        if read.is_secondary:
            continue


        chrom = bam_reader.get_reference_name(read.reference_id)
        
        strand = '+' if read.is_forward else '-'

        introns = bam_reader.find_introns([read])

        for coordpair in introns.keys():
            lend, rend = coordpair
            lend += 1
            
            intron_key = f"{lend}\t{rend}\t{strand}"

            intron_counter[chrom][intron_key] += 1


    for chrom, chrom_icounter in intron_counter.items():
        
        chrom_seq = fasta_reader.fetch(chrom)
        
        for intron, count in chrom_icounter.items():

            lend, rend, strand = intron.split("\t")
            lend = int(lend)
            rend = int(rend)

            left_dinuc = chrom_seq[lend-1:lend+1]
            right_dinuc = chrom_seq[rend-1-1:rend]

            splice_tok = f"{left_dinuc}--{right_dinuc}"

            splice_tok_adj = splice_tok + f":{strand}"

            splice_flag = "OK" if splice_tok_adj in OK_SPLICES else "NON"
            
            print("\t".join([chrom, str(lend), str(rend), str(strand), splice_tok, str(count), splice_flag]))
            
        
    






if __name__=='__main__':
    main()
