
test_quick:
	../../LRAA --bam pacbio.PBLR.small.bam --genome ref_genome.fa --cpu_budget 4 --output_prefix LRAA.small  --no_norm --tag_bam --HiFi --verbose --no_cleanup


test: test_quick test_reads_only test_include_ref_annot_guide test_quant test_chunked test_oversimplify


test_reads_only:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --cpu_budget 4 --output_prefix LRAA  --no_norm --HiFi --verbose --no_cleanup


test_include_ref_annot_guide:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --cpu_budget 4 --output_prefix LRAA.guided  --no_norm --gtf ref_annot.gtf --HiFi --verbose --no_cleanup


test_quant:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --cpu_budget 4 --output_prefix  LRAA.ref --gtf ref_annot.gtf --quant_only --tag_bam --HiFi --verbose --no_cleanup


# Chunked mode, on both routes into it.
#
# Nothing here exercised chunking until v0.21.0, and the cost of that was two
# releases in which every chunked run aborted: v0.20.0 gzipped quant tracking,
# the driver still looked for the uncompressed name, and so each chunk that
# succeeded was reported FAILED. A test that merely runs the thing would have
# caught it on the day.
#
# 3.4 Mb yields no cuts at the default spacing, so the cut size is forced down
# to make this a genuinely multi-chunk run rather than a one-chunk no-op.
test_chunked: test_chunked_entry_point test_chunked_route_parity

test_chunked_entry_point:
	../../LRAA --chunk --quant_only --bam pacbio.PBLR.bam --genome ref_genome.fa --gtf ref_annot.gtf \
	    --cpu_budget 4 --HiFi --approx_MB_per_cut 0.5 --approx_MB_per_cut_wiggle_window 0.1 \
	    --output_prefix LRAA.chunked
	test -s LRAA.chunked.LRAA.quant-only.quant.expr
	test -s LRAA.chunked.LRAA.quant-only.quant.tracking.gz

# The two routes must not drift: LRAA --chunk and the parity harness call the
# same orchestration, so their merged output is compared rather than pinned to
# expected numbers -- a relationship survives a legitimate change to the corpus
# or the algorithm, and a hardcoded row count does not.
test_chunked_route_parity: test_chunked_entry_point
	python3 ../../util/misc/run_chunked_quant_pipeline.py \
	    --bam pacbio.PBLR.bam --genome_fa ref_genome.fa --gtf ref_annot.gtf \
	    --output_dir __chunked_pipeline --arm chunked --cpu_budget 4 --HiFi \
	    --approx_MB_per_cut 0.5 --approx_MB_per_cut_wiggle_window 0.1
	diff __chunked_pipeline/merged/chunked.quant.expr LRAA.chunked.LRAA.quant-only.quant.expr
	gunzip -c __chunked_pipeline/merged/chunked.quant.tracking.gz > __chunked_pipeline.tracking
	gunzip -c LRAA.chunked.LRAA.quant-only.quant.tracking.gz > __chunked_lraa.tracking
	diff __chunked_pipeline.tracking __chunked_lraa.tracking





test_oversimplify: test_oversimplify_quantifyonly test_oversimplify_denovo test_oversimplify_refguided


test_oversimplify_quantifyonly:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --cpu_budget 4 --output_prefix  LRAA.quantonly.oversimplify --gtf ref_annot.gtf --quant_only --HiFi --oversimplify minigenome --verbose --no_cleanup

test_oversimplify_denovo:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --cpu_budget 4 --output_prefix  LRAA.denovo.oversimplify --HiFi --oversimplify minigenome --verbose --no_cleanup

test_oversimplify_refguided:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --cpu_budget 4 --output_prefix  LRAA.refguided.oversimplify  --gtf ref_annot.gtf  --HiFi --oversimplify minigenome --verbose --no_cleanup




## WDLs

include ../lraa_test_docker.mk

MINIWDL_CFG = ../../miniwdl.test.cfg
MINIWDL_RUN = miniwdl run --cfg $(MINIWDL_CFG)

test_wdls: check_test_image_revision test_wdl_ID_ref_free test_wdl_ID_ref_guided test_wdl_Quant_only test_wdl_ORF_prediction

test_wdl_ID_ref_free:
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i wdl.ID_ref_free_Quant_mode.json $(LRAA_WDL_DOCKER) --verbose

test_wdl_ID_ref_guided:
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i wdl.ID_ref_guided_Quant_mode.json $(LRAA_WDL_DOCKER) --verbose

test_wdl_Quant_only:
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i wdl.Quant_only_mode.json $(LRAA_WDL_DOCKER) --verbose

test_wdl_ORF_prediction:
	$(MINIWDL_RUN) ../../WDL/LRAA-ORF-prediction.wdl -i wdl.ORF_prediction.json $(LRAA_WDL_DOCKER_ORF) --verbose



## Cleaning

clean_miniwdl: 
	rm -rf ./_LAST ./*_LRAA_wf ./__miniwdl_call_cache ./__miniwdl_download_cache


clean: clean_miniwdl
	rm -rf ./__* ./LRAA*
