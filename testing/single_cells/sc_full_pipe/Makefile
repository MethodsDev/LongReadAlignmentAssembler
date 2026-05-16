
test:

MINIWDL_CFG = ../../../miniwdl.test.cfg
MINIWDL_RUN = miniwdl run --cfg $(MINIWDL_CFG)

test_wdls:  test_full_sc_pipe_basic_wdl  test_full_sc_pipe_cluster_guided_wdl test_full_sc_pipe_basic_wdl_QOnly test_full_sc_pipe_cluster_guided_wdl_QOnly


test_full_sc_pipe_basic_wdl:
	$(MINIWDL_RUN) ../../../WDL/LRAA-singlecell.wdl -i singlecell-test-mode-basic.inputs.json -d sc_pipe_outdir-basic --verbose


test_full_sc_pipe_basic_wdl_QOnly:
	$(MINIWDL_RUN) ../../../WDL/LRAA-singlecell.wdl -i singlecell-test-mode-basic.QUANTONLY.inputs.json  -d sc_pipe_outdir-basic-QOnly  --verbose


test_full_sc_pipe_cluster_guided_wdl:
	$(MINIWDL_RUN) ../../../WDL/LRAA-singlecell.wdl -i singlecell-test-mode-cluster-guided.inputs.json -d sc_pipe_outdir-cluster-guided --verbose

test_full_sc_pipe_cluster_guided_wdl_QOnly:
	$(MINIWDL_RUN) ../../../WDL/LRAA-singlecell.wdl -i singlecell-test-mode-cluster-guided.QUANTONLY.inputs.json -d sc_pipe_outdir-cluster-guided-QOnly --verbose




clean:
	rm -rf ./sc_pipe_outdir* ./__miniwdl_call_cache ./__miniwdl_download_cache
