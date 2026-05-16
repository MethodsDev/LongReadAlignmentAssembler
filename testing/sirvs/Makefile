

test: test_ref_free_ID test_quant test_ref_guided_ID


test_ref_free_ID:
	make clean
	../../LRAA --genome data/SIRVs1-7.genome.fa --bam data/sim.fasta.mm2.bam --norm 0 --HiFi --verbose --no_cleanup --no_shuffle_parallel_jobs --config_update SIRVs.LRAA_config_update.json
	gffcompare -r data/SIRVs1-7.annot.gtf LRAA.LRAA.ref-free.gtf -o gffcompare
	./eval/eval_sirvs.py data/SIRVs1-7.gene_trans_map gffcompare.LRAA.LRAA.ref-free.gtf.refmap
	diff gffcompare.LRAA.LRAA.ref-free.gtf.refmap.full.summary.tsv data/latest_summary.tsv

test_quant:
	make clean
	../../LRAA --quant_only --tag_bam --genome data/SIRVs1-7.genome.fa --bam data/sim.fasta.mm2.bam --gtf data/SIRVs1-7.annot.gtf --HiFi --verbose --no_cleanup --no_shuffle_parallel_jobs
	./eval/validate_quant_expr.py LRAA.LRAA.quant-only.quant.expr data/latest_quant_info.tsv


test_ref_guided_ID:
	make clean
	../../LRAA --genome data/SIRVs1-7.genome.fa --bam data/sim.fasta.mm2.bam --norm 0 --gtf data/SIRVs1-7.annot.reduced.gtf --HiFi --verbose --no_cleanup --no_shuffle_parallel_jobs


MINIWDL_CFG = ../../miniwdl.test.cfg
MINIWDL_RUN = miniwdl run --cfg $(MINIWDL_CFG)

test_docker:
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i data/input_refFree_asm.json
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i data/input_refGuided_asm.json
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i data/input_quant_only.json


clean:
	rm -rf ./gffcompare* ./LRAA.* ./__* ./*_lraa_wf ./_LAST ./__miniwdl_call_cache ./__miniwdl_download_cache
