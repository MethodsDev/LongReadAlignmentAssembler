LRAA = ../../../LRAA
SAME_GENE_GTF = ablim1_region.all_same_gene_id.gtf
PRE_EXPR = LRAA.LRAA.quant-only.pre-cross-gene-EM.quant.expr
FINAL_EXPR = LRAA.LRAA.quant-only.quant.expr

test:
	python rewrite_gtf_gene_id.py data/ablim1_region.gtf $(SAME_GENE_GTF) ABLIM1_ALL_ISOFORMS
	$(LRAA) --bam data/ablim1_region.bam --genome data/ablim1_region.fa --gtf $(SAME_GENE_GTF) --quant_only --HiFi --num_total_reads 21509677 --no_shuffle_parallel_jobs > run.log 2>&1
	python validate_same_gene_pre_post.py run.log $(PRE_EXPR) $(FINAL_EXPR)

clean:
	rm -rf $(SAME_GENE_GTF) run.log LRAA.LRAA.quant-only* __LRAA.LRAA.quant-only* __alignment_cache tx_rescue.*
