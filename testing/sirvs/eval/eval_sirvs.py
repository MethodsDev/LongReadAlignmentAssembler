#!/usr/bin/env python3

import sys, os
import pandas as pd

usage = "\n\n\tusage: {} gene_trans_map gffcompare.refmap\n\n".format(sys.argv[0])

if len(sys.argv) < 3:
    exit(usage)

ref_trans_map_file = sys.argv[1]
gffcompare_refmap_file = sys.argv[2]

    
ref_trans_map = pd.read_csv(ref_trans_map_file, header=None, names=['gene_id', 'transcript_id'], sep="\t")

gffcompare_refmap = pd.read_csv(gffcompare_refmap_file, sep="\t")

df = ref_trans_map.merge(gffcompare_refmap, how='left', left_on=['gene_id', 'transcript_id'], right_on=['ref_gene', 'ref_id'] )

df['class_code'] = df['class_code'].apply(lambda x: "!" if pd.isna(x) else x)

df_summary = df[ ['gene_id', 'transcript_id', 'class_code'] ].groupby(['gene_id','class_code']).count()

df.to_csv(gffcompare_refmap_file + ".full.tsv", sep="\t", index=False)

df_summary.to_csv(gffcompare_refmap_file + ".full.summary.tsv", sep="\t", index=True)  

sys.exit(0)

