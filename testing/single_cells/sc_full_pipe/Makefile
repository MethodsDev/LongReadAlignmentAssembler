
test:

test_wdls: test_full_sc_pipe_wdl


test_full_sc_pipe_wdl:
	 miniwdl run ../../../WDL/LRAA-singlecell.wdl -i singlecell-test.inputs.json -d sc_pipe_outdir --verbose --no-cache


clean:
	 rm -rf ./sc_pipe_outdir

