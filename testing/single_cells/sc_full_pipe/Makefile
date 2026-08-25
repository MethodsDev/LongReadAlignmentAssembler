
test:

include ../../lraa_test_docker.mk

# MINIWDL_CFG and MINIWDL_RUN come from the include above.

test_wdls: check_test_image_revision test_full_sc_pipe_basic_wdl  test_full_sc_pipe_cluster_guided_wdl test_full_sc_pipe_basic_wdl_QOnly test_full_sc_pipe_cluster_guided_wdl_QOnly


test_full_sc_pipe_basic_wdl:
	$(MINIWDL_RUN) ../../../WDL/LRAA-singlecell.wdl -i singlecell-test-mode-basic.inputs.json -d sc_pipe_outdir-basic $(LRAA_WDL_DOCKER_CORE_SC) --verbose


test_full_sc_pipe_basic_wdl_QOnly:
	$(MINIWDL_RUN) ../../../WDL/LRAA-singlecell.wdl -i singlecell-test-mode-basic.QUANTONLY.inputs.json  -d sc_pipe_outdir-basic-QOnly $(LRAA_WDL_DOCKER_CORE_SC) --verbose


test_full_sc_pipe_cluster_guided_wdl:
	$(MINIWDL_RUN) ../../../WDL/LRAA-singlecell.wdl -i singlecell-test-mode-cluster-guided.inputs.json -d sc_pipe_outdir-cluster-guided $(LRAA_WDL_DOCKER_CORE_SC) --verbose

test_full_sc_pipe_cluster_guided_wdl_QOnly:
	$(MINIWDL_RUN) ../../../WDL/LRAA-singlecell.wdl -i singlecell-test-mode-cluster-guided.QUANTONLY.inputs.json -d sc_pipe_outdir-cluster-guided-QOnly $(LRAA_WDL_DOCKER_CORE_SC) --verbose




clean: reclaim_root_owned
	rm -rf ./sc_pipe_outdir* ./__miniwdl_call_cache ./__miniwdl_download_cache

# NOTE: the whole-genome PBMC measurement target that lived here has been
# removed. Its inputs JSON names absolute paths on one workstation
# (/home/unix/bhaas/projects/LRAA_PAPER_Analyses/...), so a tracked target
# pointing at it is dead in any other clone -- which is the same defect as a
# target referencing an ignored file. The measurement it produced is recorded in
# the chunk resource defaults and their comments; to repeat it, write a local
# inputs JSON against your own library and invoke miniwdl directly.

# Cluster-guided fixture run, REF-GUIDED discovery (quant_only=false with an
# initial_annot_gtf), with all three scatter knobs set to by_chunk -- init,
# per-cluster discovery, and final quant.  The fixture contig is 2 Mb, so the
# general approx_MB_per_cut must also be dropped to 0.2 (not just the _init
# variant): the per-cluster and final-quant phases read the general one, and at
# the 10 Mb default every phase would produce exactly one chunk and exercise
# the scatter plumbing trivially.
test_sc_cg_by_chunk_refguided:
	$(MINIWDL_RUN) ../../../WDL/LRAA-singlecell.wdl -i singlecell-cluster-guided.by_chunk.REFGUIDED.inputs.json -d sc_cg_outdir-by_chunk-refguided $(LRAA_WDL_DOCKER_CORE_SC) --verbose

# Cluster-guided fixture run with ALL THREE scatter phases set to by_chunk,
# QUANT-ONLY.  The fixture contig is 2 Mb and the default cut spacing is 10 Mb,
# so approx_MB_per_cut{,_init} are dropped to 0.2 -- at defaults every phase
# would produce exactly one chunk and the fan-out would be untested.
test_sc_cg_by_chunk_qonly:
	$(MINIWDL_RUN) ../../../WDL/LRAA-singlecell.wdl -i singlecell-cluster-guided.by_chunk.QUANTONLY.inputs.json -d sc_cg_outdir-by_chunk-qonly $(LRAA_WDL_DOCKER_CORE_SC) --verbose

# Cluster-guided fixture run, DE NOVO discovery, with all three scatter knobs
# set to by_chunk -- init, per-cluster discovery, and final quant.
#
# De novo means NO annotation anywhere: initial_annot_gtf is omitted entirely,
# so the init phase runs annotation-free (ChunkedRun gets no --gtf, the chunk
# extractor writes no mini GTF, and has_gtf is false in every chunk manifest).
# ref_annot_gtf_source_gene_symbols is omitted too -- it is File? and gated on
# defined() in LRAA-singlecell.wdl, so nothing requires it; supplying the
# reference GTF there would put a reference annotation into a run whose point
# is not having one.  The consequence is that the gene-symbol incorporation
# step does not run, which is correct for this test and is covered by the
# ref-guided target above.
#
# The final-quant phase still selects its shared cut plan against an
# annotation: LRAA_quant_by_cluster.emit_shared_chunk_plan takes a REQUIRED
# annot_gtf, and in cluster-guided mode that is lraa_merge_gtf_task's
# consolidated GTF -- here a purely de novo one.  So annotation_present is true
# on both the plan and its consumers, which is what keeps ChunkedRun's
# annotation-identity refusal from firing.
#
# The fixture contig is 2 Mb, so the general approx_MB_per_cut is dropped to
# 0.2 alongside the _init variant: the per-cluster and final-quant phases read
# the general one, and at the 10 Mb default every phase would produce exactly
# one chunk and exercise the scatter plumbing trivially.
test_sc_cg_by_chunk_denovo:
	$(MINIWDL_RUN) ../../../WDL/LRAA-singlecell.wdl -i singlecell-cluster-guided.by_chunk.DENOVO.inputs.json -d sc_cg_outdir-by_chunk-denovo $(LRAA_WDL_DOCKER_CORE_SC) --verbose