#!/bin/bash

set -ex

gffcompare  -r SIRV2.annot.gtf LRAA.gtf -o gffcompare

../eval/eval_sirvs.py SIRV2.gene_trans_map gffcompare.LRAA.gtf.refmap

../eval/eval_sirvs.extract_gtfs_by_category.py  --reco_gtf LRAA.gtf --truth_gtf SIRV2.annot.gtf --refmap_full gffcompare.LRAA.gtf.refmap.full.tsv

