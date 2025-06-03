test:
	../../LRAA --genome chr19.genome.fa --bam chr19.bam
	../../util/sc/singlecell_tracking_to_sparse_matrix.R --tracking LRAA.quant.tracking --output_prefix test


clean:
	rm -f ./chr19.bam.SS* ./LRAA* ./chr19.bam.norm*
	rm -rf ./test* ./__chckpts ./__alignment_cache
