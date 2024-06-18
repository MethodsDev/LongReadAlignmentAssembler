#!/usr/bin/env python3


import sys, os, re
import random


usage = "\n\n\tusage: {} input_transcripts.fasta.tab lowrange highrange\n\n".format(sys.argv[0])
if len(sys.argv) < 4:
    exit(usage)


input_trans_filename = sys.argv[1]
lowrange = int(sys.argv[2])
highrange = int(sys.argv[3])

if not re.search("\\.tab$", input_trans_filename):
    exit(usage)

assert lowrange < highrange, "Error, lowrange {} must be < highrange {}".format(lowrange, highrange)


ofh_fasta = open("sim.fasta", "wt")
ofh_vals = open("sim.fasta.abundance_vals.tsv", "wt")

with open(input_trans_filename) as fh:
    for line in fh:
        line = line.rstrip()
        acc, sequence = line.split("\t")
        num_seqs = random.randint(lowrange, highrange)
        print("\t".join([acc, str(num_seqs)]), file=ofh_vals)
        for i in range(num_seqs):
            read_acc = f"{i}-{acc}"
            print(f">{read_acc}\n{sequence}", file=ofh_fasta)

ofh_fasta.close()
ofh_vals.close()


