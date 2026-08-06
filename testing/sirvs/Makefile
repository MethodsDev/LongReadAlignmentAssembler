

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
# counts, currently 60 exact matches -- all 37 supplied models plus 23 of the 32
# withheld ones -- with the remaining 9 withheld models undiscovered. It is a
# snapshot of current behavior, not a target: regenerate it when a change moves
# those counts intentionally.
# Of those 9, exactly 7 (SIRV202, SIRV304, SIRV306, SIRV404, SIRV506, SIRV606,
# SIRV608) are recovered by the annotation-free target, under both default and
# the loosened TSS/PolyA settings, so guided discovery is currently the weaker
# of the two on this dataset. That asymmetry is unexplained; this baseline
# records it rather than endorsing it.
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
