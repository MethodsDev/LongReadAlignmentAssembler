#!/usr/bin/env python3

import sys, os, re
import gzip

def main():

    usage = "\n\n\tusage: {} isoquant.gtf.gz\n\n".format(sys.argv[0])
    if len(sys.argv) != 2:
        exit(usage)

    gtf_file = sys.argv[1]

    with gzip.open(gtf_file, "rt") as fh:
        for line in fh:
            vals = line.split("\t")
            if len(vals) < 8:
                continue
            info = vals[8]

            if vals[2] != 'transcript':
                continue
            
            m = re.search('gene_id \\"([^\\"]+)\\"; transcript_id \\"([^\\"]+)\\";', info)
            if m:
                gene_id = m.group(1)
                transcript_id = m.group(2)
                gene_name = gene_id

                m2 = re.search(' gene_name \"([^\\"]+)\\";', info)
                if m2:
                    gene_name = m2.group(1)
                
                print("\t".join([gene_id, transcript_id, gene_name]))


    sys.exit(0)


if __name__=='__main__':
    main()

    
