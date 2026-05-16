#!/bin/bash

set -ex

gffcompare  -r SIRV4.annot.gtf LRAA.LRAA.ref-free.gtf -o gffcompare

../eval/eval_sirvs.py SIRV4.gene_trans_map gffcompare.LRAA.LRAA.ref-free.gtf.refmap

../eval/eval_sirvs.extract_gtfs_by_category.py  --reco_gtf LRAA.LRAA.ref-free.gtf --truth_gtf SIRV4.annot.gtf --refmap_full gffcompare.LRAA.LRAA.ref-free.gtf.refmap.full.tsv
