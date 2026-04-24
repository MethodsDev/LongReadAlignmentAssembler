test:
	../../../LRAA --genome data/minigenome.fa --gtf data/minigenome.annot.gtf --bam data/sim.fasta.mm2.cs.hq.bam --quant_only --tag_bam --HiFi --verbose --no_cleanup --no_shuffle_parallel_jobs
	./validate_quant_expr.py LRAA.quant.expr data/perfect_quant.tsv


clean:
	rm -rf ./LRAA.quant.* ./LRAA.tagged.bam ./LRAA.contigtmp ./__LRAA.contigtmp ./__LRAA.sgcache ./__LRAA.norm_cache ./__alignment_cache
