
test:

include ../../lraa_test_docker.mk

# MINIWDL_CFG and MINIWDL_RUN come from the include above.

test_wdls: check_test_image_revision test_full_sc_pipe_basic_wdl  test_full_sc_pipe_cluster_guided_wdl test_full_sc_pipe_basic_wdl_QOnly test_full_sc_pipe_cluster_guided_wdl_QOnly


test_full_sc_pipe_basic_wdl:
	$(MINIWDL_RUN) ../../../WDL/LRAA-singlecell.wdl -i singlecell-test-mode-basic.inputs.json -d sc_pipe_outdir-basic $(LRAA_WDL_DOCKER_CORE_SC) $(LRAA_WDL_EXTRA) --verbose


test_full_sc_pipe_basic_wdl_QOnly:
	$(MINIWDL_RUN) ../../../WDL/LRAA-singlecell.wdl -i singlecell-test-mode-basic.QUANTONLY.inputs.json  -d sc_pipe_outdir-basic-QOnly $(LRAA_WDL_DOCKER_CORE_SC) $(LRAA_WDL_EXTRA) --verbose


test_full_sc_pipe_cluster_guided_wdl:
	$(MINIWDL_RUN) ../../../WDL/LRAA-singlecell.wdl -i singlecell-test-mode-cluster-guided.inputs.json -d sc_pipe_outdir-cluster-guided $(LRAA_WDL_DOCKER_CORE_SC) $(LRAA_WDL_EXTRA) --verbose

test_full_sc_pipe_cluster_guided_wdl_QOnly:
	$(MINIWDL_RUN) ../../../WDL/LRAA-singlecell.wdl -i singlecell-test-mode-cluster-guided.QUANTONLY.inputs.json -d sc_pipe_outdir-cluster-guided-QOnly $(LRAA_WDL_DOCKER_CORE_SC) $(LRAA_WDL_EXTRA) --verbose




clean: reclaim_root_owned
	rm -rf ./sc_pipe_outdir* ./__miniwdl_call_cache ./__miniwdl_download_cache

# NOTE: the whole-genome PBMC measurement target that lived here has been
# removed. Its inputs JSON names absolute paths on one workstation
# (/home/unix/bhaas/projects/LRAA_PAPER_Analyses/...), so a tracked target
# pointing at it is dead in any other clone -- which is the same defect as a
# target referencing an ignored file. The measurement it produced is recorded in
# the chunk resource defaults and their comments; to repeat it, write a local
# inputs JSON against your own library and invoke miniwdl directly.