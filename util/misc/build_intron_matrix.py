#!/usr/bin/env python3

import sys, os, re
import glob
from collections import defaultdict


def main():


    samplename_to_intron_to_counts = dict()

    samplenames = list()
    
    for file in glob.glob("*.introns"):

        sample_name = file
        sample_name = sample_name.replace(".aligned.sorted.bam.introns", "")

        samplenames.append(sample_name)
        
        intron_counts = parse_intron_counts(file)
        samplename_to_intron_to_counts[sample_name] = intron_counts


    print("\t" + "\t".join(samplenames))

    all_introns = set()
    for sample_name, intron_counts in samplename_to_intron_to_counts.items():
        intron_toks = intron_counts.keys()
        all_introns.update(intron_toks)

    for intron in intron_toks:
        vals = list()
        vals.append(intron)
        for sample_name in samplenames:
            count = str(samplename_to_intron_to_counts[sample_name][intron])
            vals.append(count)

        print("\t".join(vals))
        
    sys.exit(0)

        

def parse_intron_counts(introns_filename):



    intron_counts = defaultdict(int)
    
    with open(introns_filename) as fh:

        for line in fh:
            """
            0       chr1
            1       14830
            2       14969
            3       -
            4       CT--AC
            5       52
            6       OK
            """

            line = line.rstrip()
            
            vals = line.split("\t")
            count = vals.pop(5)
            tok = "^".join(vals)
            intron_counts[tok] = count

    return intron_counts






if __name__=='__main__':
    main()
