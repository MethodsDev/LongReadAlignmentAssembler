for file in *.bam; do ../../../../util/SQANTI-like_cats_for_reads_or_isoforms.py --input_bam $file  --ref_gtf ../SIRVs1-7.annot.gtf   --output_prefix $file; done
