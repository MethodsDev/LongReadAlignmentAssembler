test:
	../../../LRAA --genome ../data/chr19.genome.fa --bam ../data/chr19.bam --HiFi --verbose --no_cleanup --gtf ../data/chr19.gtf  --quant_only
	../../../util/sc/singlecell_tracking_to_sparse_matrix.py --tracking LRAA.quant.tracking --output_prefix test
	../../../util/sc/incorporate_gene_symbols_in_sc_features.py  --ref_gtf ../data/chr19.gtf --id_mappings test.gene_transcript_splicehashcode.tsv   --sparseM_dirs 'test^gene-sparseM' 'test^isoform-sparseM' 'test^splice_pattern-sparseM'


test_wdls:


clean:
	rm -rf ./chr19.bam.SS* ./LRAA* ./chr19.bam.norm*
	rm -rf ./test* ./__chckpts ./__alignment_cache
	rm -rf ./__*  ./*.count



