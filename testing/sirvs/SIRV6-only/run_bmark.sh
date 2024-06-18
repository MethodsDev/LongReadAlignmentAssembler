#!/bin/bash

set -ex

gffcompare  -r SIRV6.annot.gtf LRAA.gtf -o gffcompare

../eval/eval_sirvs.py SIRV6.gene_trans_map gffcompare.LRAA.gtf.refmap

../eval/eval_sirvs.extract_gtfs_by_category.py  --reco_gtf LRAA.gtf --truth_gtf SIRV6.annot.gtf --refmap_full gffcompare.LRAA.gtf.refmap.full.tsv

