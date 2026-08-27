
test:


include ../../lraa_test_docker.mk

# MINIWDL_CFG and MINIWDL_RUN come from the include above.

test_wdls: check_test_image_revision denovo_wdl ref_guided_wdl quant_only_wdl

quant_only_wdl:
	$(MINIWDL_RUN) ../../../WDL/LRAA-cell_cluster_guided.wdl -i sc_pipe_inputs-quant_only.json -d sc_pipe_quant_only_outdir $(LRAA_WDL_DOCKER_CORE_SC) $(LRAA_WDL_EXTRA)


ref_guided_wdl:
	$(MINIWDL_RUN) ../../../WDL/LRAA-cell_cluster_guided.wdl -i sc_pipe_inputs-ref_guided.json -d sc_pipe_ref_guided_outdir $(LRAA_WDL_DOCKER_CORE_SC) $(LRAA_WDL_EXTRA)

denovo_wdl:
	$(MINIWDL_RUN) ../../../WDL/LRAA-cell_cluster_guided.wdl -i sc_pipe_inputs-denovo.json -d sc_pipe_denovo_outdir $(LRAA_WDL_DOCKER_CORE_SC) $(LRAA_WDL_EXTRA)

clean: reclaim_root_owned
	rm -rf ./sc_pipe_ref_guided_outdir ./sc_pipe_denovo_outdir ./sc_pipe_quant_only_outdir ./__miniwdl_call_cache ./__miniwdl_download_cache
