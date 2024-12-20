


test:
	../../LRAA --bam pacbioccsminimap.minig.PBLR.bam --genome minig.fa --no_norm



test_wdl_quant_only:
	 miniwdl run ../../WDL/LRAA.wdl -i wf.quant-only.json

test_wdl_ref_guided:
	miniwdl run ../../WDL/LRAA.wdl -i wf.ref-guided.json

clean:
	rm -rf ./__*  ./LRAA.gtf *_LRAA_wf _LAST


