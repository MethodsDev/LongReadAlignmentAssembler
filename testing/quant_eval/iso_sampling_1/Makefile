test: no_EM with_EM

no_EM:
	../../../LRAA --genome data/minigenome.fa --gtf data/annot.gtf --bam data/IsoSeqFaked.bam --quant_only --no_EM --output_prefix LRAA.no_EM --HiFi
	./util/get_corr.Rscript LRAA.no_EM.quant.expr
	./util/validate_min_corr.py LRAA.no_EM.quant.expr.corr_val 0.96

with_EM:
	../../../LRAA --genome data/minigenome.fa --gtf data/annot.gtf --bam data/IsoSeqFaked.bam --quant_only --output_prefix LRAA.with_EM --HiFi
	./util/get_corr.Rscript LRAA.with_EM.quant.expr 
	./util/validate_min_corr.py LRAA.with_EM.quant.expr.corr_val 0.97



with_EM_no_read_weighting:
	../../../LRAA --genome data/minigenome.fa --gtf data/annot.gtf --bam data/IsoSeqFaked.bam --quant_only --output_prefix LRAA.with_EM --no_weight_reads_by_3prime_agreement --HiFi
	./util/get_corr.Rscript LRAA.with_EM.quant.expr 
	./util/validate_min_corr.py LRAA.with_EM.quant.expr.corr_val 0.97
clean:
	rm -rf ./LRAA.* ./__alignment_cache
