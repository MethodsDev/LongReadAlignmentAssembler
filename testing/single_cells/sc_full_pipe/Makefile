
test:

test_wdls:  test_full_sc_pipe_basic_wdl  test_full_sc_pipe_cluster_guided_wdl test_full_sc_pipe_basic_wdl_QOnly


test_full_sc_pipe_basic_wdl:
	 miniwdl run ../../../WDL/LRAA-singlecell.wdl -i singlecell-test-mode-basic.inputs.json -d sc_pipe_outdir-basic --verbose --no-cache


test_full_sc_pipe_basic_wdl_QOnly:
	miniwdl run ../../../WDL/LRAA-singlecell.wdl -i singlecell-test-mode-basic.QUANTONLY.inputs.json  -d sc_pipe_outdir-basic-QOnly  --verbose --no-cache


test_full_sc_pipe_cluster_guided_wdl:
	 miniwdl run ../../../WDL/LRAA-singlecell.wdl -i singlecell-test-mode-cluster-guided.inputs.json -d sc_pipe_outdir-cluster-guided --verbose --no-cache




clean:
	 rm -rf ./sc_pipe_outdir*

