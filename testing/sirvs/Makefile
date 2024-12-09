

test: test_ref_free_ID test_quant test_ref_guided_ID


test_ref_free_ID:
	../../LRAA --genome data/SIRVs1-7.genome.fa --bam data/sim.fasta.mm2.bam --norm 0 
	gffcompare -r data/SIRVs1-7.annot.gtf LRAA.gtf -o gffcompare
	./eval/eval_sirvs.py data/SIRVs1-7.gene_trans_map gffcompare.LRAA.gtf.refmap
	diff gffcompare.LRAA.gtf.refmap.full.summary.tsv data/latest_summary.tsv

test_quant:
	../../LRAA --quant_only --tag_bam --genome data/SIRVs1-7.genome.fa --bam data/sim.fasta.mm2.bam --gtf data/SIRVs1-7.annot.gtf
	diff LRAA.quant.expr data/latest_quant_info.tsv


test_ref_guided_ID:
	../../LRAA --genome data/SIRVs1-7.genome.fa --bam data/sim.fasta.mm2.bam --norm 0 --gtf data/SIRVs1-7.annot.reduced.gtf 



test_docker:
	miniwdl run ../../WDL/LRAA.wdl -i data/input_refFree_asm.json --no-cache
	miniwdl run ../../WDL/LRAA.wdl -i data/input_refGuided_asm.json --no-cache
	miniwdl run ../../WDL/LRAA.wdl -i data/input_quant_only.json --no-cache


clean:
	rm -rf ./gffcompare* ./LRAA.gtf ./__* ./*_lraa_wf ./_LAST ./LRAA.*.lmdb ./LRAA.tagged.bam


