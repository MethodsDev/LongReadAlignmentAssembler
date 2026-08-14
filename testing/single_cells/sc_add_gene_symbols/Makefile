test:


include ../../lraa_test_docker.mk

MINIWDL_CFG = ../../../miniwdl.test.cfg
MINIWDL_RUN = miniwdl run --cfg $(MINIWDL_CFG)

test_wdls:
	$(MINIWDL_RUN) ../../../WDL/subwdls/Incorporate_gene_symbols.wdl  -i inputs.json -d test_incorporate_symbols_outdir $(LRAA_WDL_DOCKER) --verbose

clean:
	rm -rf ./test_incorporate_symbols_outdir ./__miniwdl_call_cache ./__miniwdl_download_cache
