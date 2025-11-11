test:
	../../LRAA --genome chr19.genome.fa --bam chr19.bam --HiFi --verbose --no_cleanup
	../../util/sc/singlecell_tracking_to_sparse_matrix.R --tracking LRAA.quant.tracking --output_prefix test
	gffcompare -r chr19.gtf LRAA.gtf
	../../util/sc/incorporate_gene_symbols_in_sc_features.py  --ref_gtf chr19.gtf --LRAA_gtf LRAA.gtf --gffcompare_tracking gffcmp.tracking  --id_mappings test.gene_transcript_splicehashcode.tsv   --sparseM_dirs 'test^gene-sparseM' 'test^isoform-sparseM' test^splice_pattern-sparseM


test_full_sc_pipe_wdl:
	 miniwdl run ../../WDL/LRAA-singlecell.wdl -i singlecell-test.inputs.json -d sc_pipe_outdir

clean:
	rm -rf ./chr19.bam.SS* ./LRAA* ./chr19.bam.norm*
	rm -rf ./test* ./__chckpts ./__alignment_cache
	rm -rf ./sc_pipe_outdir
	rm -rf ./__* ./gffcomp*
	cd sc_cluster_guided_pipe && $(MAKE) clean


