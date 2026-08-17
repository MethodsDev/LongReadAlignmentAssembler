test:


include ../../lraa_test_docker.mk

# MINIWDL_CFG and MINIWDL_RUN come from the include above.

test_wdls: check_test_image_revision
	$(MINIWDL_RUN) ../../../WDL/subwdls/Incorporate_gene_symbols.wdl  -i inputs.json -d test_incorporate_symbols_outdir $(LRAA_WDL_DOCKER) --verbose

clean: reclaim_root_owned
	rm -rf ./test_incorporate_symbols_outdir ./__miniwdl_call_cache ./__miniwdl_download_cache
