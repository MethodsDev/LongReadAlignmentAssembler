#!/usr/bin/env python3

import sys, os, re
import logging
import argparse
from collections import defaultdict
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)


def main():

    parser = argparse.ArgumentParser(description="extract GTF records according to benchmarking categories",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--reco_gtf", type=str, required=True, help="reco gtf file")
    parser.add_argument("--truth_gtf", type=str, required=True, help="reco gtf file")

    parser.add_argument("--refmap_full", type=str, required=True, help="gffcompare.LRAA.gtf.refmap.full.tsv")
    
    args = parser.parse_args()
    
    reco_gtf_file = args.reco_gtf
    truth_gtf_file = args.truth_gtf
    
    refmap_full = args.refmap_full
    

    reco_transcript_to_gtf_txt = parse_transcripts_from_gtf(reco_gtf_file)
    
    truth_transcript_to_gtf_txt = parse_transcripts_from_gtf(truth_gtf_file)
    
    base_gtf_name = os.path.basename(reco_gtf_file)

    refmap_df = pd.read_csv(refmap_full, sep="\t")

    
    
    types = list(refmap_df.class_code.unique())

    fhs = dict()
    for ctype in types:
        fhs[ctype] = open(base_gtf_name + "." + ctype + ".gtf", "wt")

        
    for _, row in refmap_df.iterrows():

        ctype = row['class_code']
        qry_id_list = row['qry_id_list']
        gtf_txt = None
        if pd.isna(qry_id_list):
            transcript_id = row['transcript_id']
            gtf_txt = truth_transcript_to_gtf_txt[transcript_id]
        else:
            transcript_id = qry_id_list.split("|")[1]
            gtf_txt = reco_transcript_to_gtf_txt[transcript_id]

        print(gtf_txt, file=fhs[ctype])


    sys.exit(0)







def parse_transcripts_from_gtf(gtf_file):

    transcript_to_gtf_txt = defaultdict(str)
    
    with open(gtf_file) as fh:
        for line in fh:
            m = re.search("transcript_id \"([^\"]+)\"", line)
            if m:
                transcript_id = m.group(1)
                transcript_to_gtf_txt[transcript_id] += line
                
    return transcript_to_gtf_txt

    

if __name__=='__main__':
    main()
