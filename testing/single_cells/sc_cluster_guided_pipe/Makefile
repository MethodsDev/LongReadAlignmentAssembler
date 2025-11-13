
test:


test_wdls: denovo_wdl ref_guided_wdl quant_only_wdl

quant_only_wdl:
	miniwdl run ../../../WDL/LRAA-cell_cluster_guided.wdl -i sc_pipe_inputs-quant_only.json -d sc_pipe_quant_only_outdir


ref_guided_wdl:
	miniwdl run ../../../WDL/LRAA-cell_cluster_guided.wdl -i sc_pipe_inputs-ref_guided.json -d sc_pipe_ref_guided_outdir

denovo_wdl:
	miniwdl run ../../../WDL/LRAA-cell_cluster_guided.wdl -i sc_pipe_inputs-denovo.json -d sc_pipe_denovo_outdir

clean:
	rm -rf ./sc_pipe_ref_guided_outdir ./sc_pipe_denovo_outdir ./sc_pipe_quant_only_outdir

