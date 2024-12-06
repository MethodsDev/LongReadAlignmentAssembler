

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


test_wdl_ID_ref_free:
	miniwdl run ../../WDL/LRAA.scatter.wdl -i wdl.ID_ref_free_Quant_mode.json



clean:
	rm -rf ./__* ./LRAA*.gtf ./LRAA*.quant ./LRAA*.lmdb ./LRAA.ref.tagged.bam


