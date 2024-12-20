


test:
	../../LRAA --bam pacbioccsminimap.minig.PBLR.bam --genome minig.fa --no_norm


test_wdls: test_wdl_quant_only test_wdl_ref_guided test_wdl_ref_free


##########
# wdls

test_wdl_quant_only: wdl_clean
	 miniwdl run ../../WDL/LRAA.wdl -i wf.quant-only.json  --verbose --no-cache

test_wdl_ref_guided: wdl_clean
	miniwdl run ../../WDL/LRAA.wdl -i wf.ref-guided.json  --verbose --no-cache

test_wdl_ref_free: wdl_clean
	miniwdl run ../../WDL/LRAA.wdl -i wf.ref-free.json  --verbose --no-cache



######
# cleaning

wdl_clean:
	rm -rf  ./*_LRAA_wf ./_LAST

clean: wdl_clean
	rm -rf ./__*  ./LRAA.gtf


