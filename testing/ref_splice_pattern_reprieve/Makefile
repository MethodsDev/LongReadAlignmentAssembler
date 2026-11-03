# Reference identity in ref-guided runs (v0.50.0): the splice_hash / ref_splice_pattern /
# ref_FSM_reads / ref_supported GTF attributes, guiding a run with an LRAA GTF (as the
# single-cell cluster runs are guided by the init GTF), the merge carrying the flags, and
# the refusal of an LRAA guide without them. SIRV simulated reads (testing/sirvs/data);
# each run takes seconds. check_ref_identity.py asserts properties, not a golden GTF.

LRAA = ../../LRAA
MERGE = ../../util/merge_LRAA_GTFs.py
D = ../sirvs/data
RUN = --genome ../$(D)/SIRVs1-7.genome.fa --bam ../$(D)/sim.fasta.mm2.bam --HiFi
CHECK = ./check_ref_identity.py

.PHONY: test clean

test: run_rg/done run_rg_chunked/done run_guided/done run_unflagged/done run_unflagged_chunked/done \
      run_unflagged_allowed/done run_unflagged_quant/done run_merge/done run_merge_unflagged/done \
      run_merge_task_unflagged/done run_merge_task_config/done
	$(CHECK) $(D)/sim.fasta.mm2.bam

# 1. guided by half of the SIRV annotation (GTF column 2 is not LRAA: a reference)
run_rg/done:
	rm -rf run_rg; mkdir -p run_rg
	cd run_rg && ../$(LRAA) $(RUN) --gtf ../$(D)/SIRVs1-7.annot.half.gtf --no_chunk --output_prefix rg > run.log 2>&1
	touch $@

# 1b. the same through the chunked pipeline (chunk workers scan their own chunk bam)
run_rg_chunked/done:
	rm -rf run_rg_chunked; mkdir -p run_rg_chunked
	cd run_rg_chunked && ../$(LRAA) $(RUN) --gtf ../$(D)/SIRVs1-7.annot.half.gtf --chunk --output_prefix rg > run.log 2>&1
	touch $@

# 2. guided by run 1's GTF, as a cluster run is guided by the init GTF: the pooled
#    ref_FSM_reads come from the guide
run_guided/done: run_rg/done
	rm -rf run_guided; mkdir -p run_guided
	cd run_guided && ../$(LRAA) $(RUN) --gtf ../run_rg/rg.LRAA.ref-guided.gtf --no_chunk --output_prefix guided > run.log 2>&1
	touch $@

# run 1's GTF as an LRAA GTF from before v0.50.0 (attributes stripped)
unflagged.gtf: run_rg/done
	sed -E 's/ (splice_hash|ref_splice_pattern|ref_FSM_reads|ref_supported) "[^"]*";//g' run_rg/rg.LRAA.ref-guided.gtf > $@

# 3. discovery guided by it is refused (the run must fail)
run_unflagged/done: unflagged.gtf
	rm -rf run_unflagged; mkdir -p run_unflagged
	cd run_unflagged && ! ../$(LRAA) $(RUN) --gtf ../unflagged.gtf --no_chunk --output_prefix unflagged > run.log 2>&1
	touch $@

# 3b. the chunked pipeline refuses it at startup, before any chunk work: the run dir
#     must hold nothing but the log
run_unflagged_chunked/done: unflagged.gtf
	rm -rf run_unflagged_chunked; mkdir -p run_unflagged_chunked
	cd run_unflagged_chunked && ! ../$(LRAA) $(RUN) --gtf ../unflagged.gtf --chunk --output_prefix unflagged > run.log 2>&1
	test "$$(ls -A run_unflagged_chunked)" = "run.log"
	touch $@

# 4. ... unless allow_unflagged_LRAA_guide
run_unflagged_allowed/done: unflagged.gtf
	rm -rf run_unflagged_allowed; mkdir -p run_unflagged_allowed
	cd run_unflagged_allowed && ../$(LRAA) $(RUN) --gtf ../unflagged.gtf --no_chunk --output_prefix allowed \
	    --config_update ../allow_unflagged.json > run.log 2>&1
	touch $@

# 5. quant-only filters nothing, so it takes any GTF
run_unflagged_quant/done: unflagged.gtf
	rm -rf run_unflagged_quant; mkdir -p run_unflagged_quant
	cd run_unflagged_quant && ../$(LRAA) $(RUN) --gtf ../unflagged.gtf --quant_only --no_chunk --output_prefix quant > run.log 2>&1
	touch $@

# 6. the cross-cluster merge: run 2's GTF as a cluster GTF, run 1's as --pooled_gtf
run_merge/done: run_rg/done run_guided/done
	rm -rf run_merge; mkdir -p run_merge
	cd run_merge && ../$(MERGE) --genome ../$(D)/SIRVs1-7.genome.fa --gtf ../run_guided/guided.LRAA.ref-guided.gtf \
	    --pooled_gtf ../run_rg/rg.LRAA.ref-guided.gtf --output_gtf merged.gtf > merge.log 2>&1
	touch $@

# 7. the merge refuses unflagged LRAA inputs (must fail)
run_merge_unflagged/done: unflagged.gtf run_guided/done
	rm -rf run_merge_unflagged; mkdir -p run_merge_unflagged
	cd run_merge_unflagged && ! ../$(MERGE) --genome ../$(D)/SIRVs1-7.genome.fa --gtf ../run_guided/guided.LRAA.ref-guided.gtf \
	    --pooled_gtf ../unflagged.gtf --output_gtf merged.gtf > merge.log 2>&1
	touch $@

# 8. the Terra merge task's own command (rendered from the WDL, run here without a
#    container): refuses the unflagged pooled GTF, and merges it with config_update_json
TASK_ENV = PATH=$(abspath ../..):$(abspath ../../util):$$PATH
TASK_ARGS = ../$(D)/SIRVs1-7.genome.fa ../unflagged.gtf ../run_guided/guided.LRAA.ref-guided.gtf

run_merge_task_unflagged/done: unflagged.gtf run_guided/done
	rm -rf run_merge_task_unflagged; mkdir -p run_merge_task_unflagged
	cd run_merge_task_unflagged && ../render_merge_task.py $(TASK_ARGS) > task.sh && \
	    ! $(TASK_ENV) bash task.sh > task.log 2>&1
	grep -q allow_unflagged_LRAA_guide run_merge_task_unflagged/command_output.log
	touch $@

run_merge_task_config/done: unflagged.gtf run_guided/done
	rm -rf run_merge_task_config; mkdir -p run_merge_task_config
	cd run_merge_task_config && ../render_merge_task.py $(TASK_ARGS) '{"allow_unflagged_LRAA_guide": true}' > task.sh && \
	    $(TASK_ENV) bash task.sh > task.log 2>&1
	test -s run_merge_task_config/S.LRAA.sc_merged.gtf
	touch $@

clean:
	rm -rf run_* unflagged.gtf
