

test: test_quick

all: test_reads_only test_include_ref_annot_guide test_quant

test_quick:
	../../LRAA --bam pacbio.PBLR.small.bam --genome ref_genome.fa --CPU 4 --output_prefix LRAA.small  --no_norm


test_reads_only:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --CPU 4 --output_prefix LRAA  --no_norm


test_include_ref_annot_guide:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --CPU 4 --output_prefix LRAA.guided  --no_norm --gtf ref_annot.gtf


test_quant:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --CPU 4 --output_prefix  LRAA.ref --gtf ref_annot.gtf --quant_only --tag_bam


## WDLs

test_wdls: test_wdl_ID_ref_free test_wdl_ID_ref_guided test_wdl_Quant_only

test_wdl_ID_ref_free: clean_miniwdl
	miniwdl run ../../WDL/LRAA.wdl -i wdl.ID_ref_free_Quant_mode.json

test_wdl_ID_ref_guided: clean_miniwdl
	miniwdl run ../../WDL/LRAA.wdl -i wdl.ID_ref_guided_Quant_mode.json

test_wdl_Quant_only: clean_miniwdl
#	miniwdl run ../../WDL/LRAA.scatter.wdl -i wdl.Quant_only_mode.json
# - note miniwdl errors out... but cromwell works fine
	java -jar ~/BIN/cromwell-85.jar run ../../WDL/LRAA.wdl -i wdl.Quant_only_mode.json 



## Cleaning

clean_miniwdl: 
	rm -rf ./*CombinedWorkflow ./_LAST

clean_cromwell: clean_miniwdl
	sudo rm -rf ./cromwell-*

clean: clean_cromwell
	rm -rf ./__* ./LRAA*.gtf ./LRAA*.quant ./LRAA*.lmdb ./LRAA.ref.tagged.bam


