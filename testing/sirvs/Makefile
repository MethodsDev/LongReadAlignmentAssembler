

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


# Guided by data/SIRVs1-7.annot.half.gtf, which supplies 37 of the 69 SIRV
# isoforms, including the two under 200 bases (SIRV311, SIRV618). Recovery is
# scored against the complete annotation. The baseline pins per-gene class-code
# counts, currently 67 exact matches -- all 37 supplied models plus 30 of the 32
# withheld ones -- leaving SIRV107 and SIRV206 undiscovered. It is a snapshot of
# current behavior, not a target: regenerate it when a change moves those counts
# intentionally.
#
# Comparing against test_ref_free_ID needs care: that target loosens the TSS and
# PolyA site thresholds (SIRVs.LRAA_config_update.json) and this one does not.
# Under matched settings guided recovery is a superset of annotation-free -- 67 vs
# 64 at the defaults, 68 vs 66 with both loosened. As the two targets actually
# run, neither set contains the other: annotation-free reaches SIRV206 on the
# looser thresholds, while this target reaches the two short models. SIRV107 is
# recovered by no mode in either configuration, despite carrying 11 reads.
test_ref_guided_ID:
	make clean
	../../LRAA --genome data/SIRVs1-7.genome.fa --bam data/sim.fasta.mm2.bam --norm 0 --gtf data/SIRVs1-7.annot.half.gtf --HiFi --verbose --no_cleanup --no_shuffle_parallel_jobs
	gffcompare -r data/SIRVs1-7.annot.gtf LRAA.LRAA.ref-guided.gtf -o gffcompare
	./eval/eval_sirvs.py data/SIRVs1-7.gene_trans_map gffcompare.LRAA.LRAA.ref-guided.gtf.refmap
	diff gffcompare.LRAA.LRAA.ref-guided.gtf.refmap.full.summary.tsv data/latest_refguided_summary.tsv


MINIWDL_CFG = ../../miniwdl.test.cfg
MINIWDL_RUN = miniwdl run --cfg $(MINIWDL_CFG)

test_docker:
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i data/input_refFree_asm.json
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i data/input_refGuided_asm.json
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i data/input_quant_only.json


clean:
	rm -rf ./gffcompare* ./LRAA.* ./__* ./*_lraa_wf ./_LAST ./__miniwdl_call_cache ./__miniwdl_download_cache
