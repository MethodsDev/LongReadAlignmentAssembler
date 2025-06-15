
ref_guided:
	miniwdl run ../../../WDL/LRAA-cell_cluster_guided.wdl -i sc_pipe_inputs-ref_guided.json

quant_only:
	miniwdl run ../../../WDL/LRAA-cell_cluster_guided.wdl -i sc_pipe_inputs-quant_only.json

clean:
	rm -rf ./202*
	rm -f ./_LAST
