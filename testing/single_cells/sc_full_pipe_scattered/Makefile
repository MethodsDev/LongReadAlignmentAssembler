
test:

include ../../lraa_test_docker.mk

MINIWDL_CFG = ../../../miniwdl.test.cfg
MINIWDL_RUN = miniwdl run --cfg $(MINIWDL_CFG)

test_wdls: test_full_sc_pipe_cluster_guided_scattered_wdl

test_full_sc_pipe_cluster_guided_scattered_wdl:
	$(MINIWDL_RUN) ../../../WDL/LRAA-singlecell.wdl -i singlecell-test-mode-cluster-guided.scattered.inputs.json -d sc_pipe_outdir-cluster-guided-scattered $(LRAA_WDL_DOCKER_CORE_SC) --verbose


clean:
	rm -rf ./sc_pipe_outdir* ./__miniwdl_call_cache ./__miniwdl_download_cache
