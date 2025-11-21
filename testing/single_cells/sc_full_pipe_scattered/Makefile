
test:

test_wdls: test_full_sc_pipe_cluster_guided_scattered_wdl

test_full_sc_pipe_cluster_guided_scattered_wdl:
	 miniwdl run ../../../WDL/LRAA-singlecell.wdl -i singlecell-test-mode-cluster-guided.scattered.inputs.json -d sc_pipe_outdir-cluster-guided-scattered --verbose --no-cache


clean:
	 rm -rf ./sc_pipe_outdir*

