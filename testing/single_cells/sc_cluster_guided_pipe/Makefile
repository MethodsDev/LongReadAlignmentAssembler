
all: ref_guided quant_only

ref_guided:
	miniwdl run ../../../WDL/LRAA-cell_cluster_guided.wdl -i sc_pipe_inputs-ref_guided.json -d sc_pipe_ref_guided_outdir

denovo:
	miniwdl run ../../../WDL/LRAA-cell_cluster_guided.wdl -i sc_pipe_inputs-denovo.json -d sc_pipe_denovo_outdir

clean:
	rm -rf ./sc_pipe_ref_guided_outdir ./sc_pipe_denovo_outdir

