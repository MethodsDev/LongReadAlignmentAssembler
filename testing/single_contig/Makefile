
test_quick:
	../../LRAA --bam pacbio.PBLR.small.bam --genome ref_genome.fa --CPU 4 --output_prefix LRAA.small  --no_norm --tag_bam --HiFi --verbose --no_cleanup


test: test_quick test_reads_only test_include_ref_annot_guide test_quant test_oversimplify


test_reads_only:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --CPU 4 --output_prefix LRAA  --no_norm --HiFi --verbose --no_cleanup


test_include_ref_annot_guide:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --CPU 4 --output_prefix LRAA.guided  --no_norm --gtf ref_annot.gtf --HiFi --verbose --no_cleanup


test_quant:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --CPU 4 --output_prefix  LRAA.ref --gtf ref_annot.gtf --quant_only --tag_bam --HiFi --verbose --no_cleanup





test_oversimplify: test_oversimplify_quantifyonly test_oversimplify_denovo test_oversimplify_refguided


test_oversimplify_quantifyonly:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --CPU 4 --output_prefix  LRAA.quantonly.oversimplify --gtf ref_annot.gtf --quant_only --HiFi --oversimplify minigenome --verbose --no_cleanup

test_oversimplify_denovo:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --CPU 4 --output_prefix  LRAA.denovo.oversimplify --HiFi --oversimplify minigenome --verbose --no_cleanup

test_oversimplify_refguided:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --CPU 4 --output_prefix  LRAA.refguided.oversimplify  --gtf ref_annot.gtf  --HiFi --oversimplify minigenome --verbose --no_cleanup




## WDLs

test_wdls: test_wdl_ID_ref_free test_wdl_ID_ref_guided test_wdl_Quant_only

test_wdl_ID_ref_free: clean_miniwdl
	miniwdl run ../../WDL/LRAA.wdl -i wdl.ID_ref_free_Quant_mode.json --verbose --no-cache

test_wdl_ID_ref_guided: clean_miniwdl
	miniwdl run ../../WDL/LRAA.wdl -i wdl.ID_ref_guided_Quant_mode.json --verbose --no-cache

test_wdl_Quant_only: clean_miniwdl
	miniwdl run ../../WDL/LRAA.wdl -i wdl.Quant_only_mode.json --verbose --no-cache



## Cleaning

clean_miniwdl: 
	rm -rf ./_LAST ./*_LRAA_wf 


clean: clean_miniwdl
	rm -rf ./__* ./LRAA*


