

test: test_ref_free_ID test_quant test_ref_guided_ID test_unsupported_ref_model


test_ref_free_ID:
	make clean
	../../LRAA --genome data/SIRVs1-7.genome.fa --bam data/sim.fasta.mm2.bam --norm 0 --HiFi --verbose --no_cleanup --no_shuffle_parallel_jobs --config_update SIRVs.LRAA_config_update.json
	gffcompare -r data/SIRVs1-7.annot.gtf LRAA.LRAA.ref-free.gtf -o gffcompare
	./eval/eval_sirvs.py data/SIRVs1-7.gene_trans_map gffcompare.LRAA.LRAA.ref-free.gtf.refmap
	diff gffcompare.LRAA.LRAA.ref-free.gtf.refmap.full.summary.tsv data/latest_summary.tsv

test_quant:
	make clean
	../../LRAA --quant_only --tag_bam --genome data/SIRVs1-7.genome.fa --bam data/sim.fasta.mm2.bam --gtf data/SIRVs1-7.annot.gtf --HiFi --verbose --no_cleanup --no_shuffle_parallel_jobs
	./eval/validate_quant_expr.py LRAA.LRAA.quant-only.quant.expr data/latest_quant_info.tsv


# Guided by data/SIRVs1-7.annot.half.gtf, which supplies 37 of the 69 SIRV
# isoforms, including the two under 200 bases (SIRV311, SIRV618). Recovery is
# scored against the complete annotation. The baseline pins per-gene class-code
# counts, currently 67 exact matches -- all 37 supplied models plus 30 of the 32
# withheld ones -- leaving SIRV107 and SIRV206 undiscovered. It is a snapshot of
# current behavior, not a target: regenerate it when a change moves those counts
# intentionally.
#
# Comparing against test_ref_free_ID needs care: that target loosens the TSS and
# PolyA site thresholds (SIRVs.LRAA_config_update.json) and this one does not.
# Under matched settings guided recovery is a superset of annotation-free: 67 exact
# against 63 at the defaults, 68 against 65 with both loosened, and annotation-free
# also leaves one model merely contained rather than exact. As the two targets
# actually run, neither set contains the other: annotation-free reaches SIRV206 on
# the looser thresholds, while this target reaches the two short models.
#
# SIRV107 is a settled case, not an open one. It is undiscoverable rather than
# under-supported: its 11 reads are pristine full-length alignments, but it is
# annotated on the minus strand while its junction at 10792-10882 carries a GT-AG
# motif implying the plus strand (that junction is shared with plus-strand SIRV109),
# and its junction at 11243-11403 is CT-AG, canonical in neither orientation. A
# read-derived junction is admitted to the splice graph only when its motif implies
# the strand being built, so those reads have no graph path and the locus cannot be
# reconstructed from reads in any mode. Junctions supplied by an annotation skip that
# check, so SIRV107 is recovered whenever it is provided: test_quant quantifies it at
# its full 11 reads, and a guided run against the complete annotation rebuilds its
# exact structure. Its absence here follows from the reference sequence and is not a
# regression to chase.
#
# The class-code summary above keys on intron chains, so it cannot see terminal
# coordinates: a boundary shift leaves every count unchanged. The three checks after
# the diff pin terminals directly. SIRV311 and SIRV618 are asserted exactly, matching
# their annotated spans.
#
# SIRV7 pins the harder case. SIRV701 and SIRV705 are terminal twins: identical intron
# chains, identical cDNA length, differing by 2 bases at each end, and both supplied
# here. Both are reconstructed, deterministically and with stable read counts of 2 and
# 3, and degradation pruning then removes one. Each twin is a valid degradation
# product of the other, so this used to be decided by whichever was examined first,
# which was set iteration over Transcript objects -- address order, since Transcript
# defines no __hash__. The run-to-run flip that caused is fixed: pruning now breaks a
# cdna-length tie on assigned reads, so SIRV705 wins on its 3 reads against 2 and the
# terminals follow the read evidence.
#
# The last check pins total assigned reads at 537 of the 551 in the BAM. The 14
# unassigned are SIRV107's 11 and SIRV206's 3: both are withheld from the half
# annotation and neither is reconstructed, so neither is ever a rescue target.
# SIRV206's 3 were already unassigned when this baseline read 548. SIRV107's 11
# joined them once rescue gained an aligned-length floor: transcriptome rescue had
# been crediting those reads fractionally to unrelated SIRV1 isoforms on alignments
# covering under 30% of the read, which is why any nonzero
# rescue_unassigned_min_aligned_read_frac now rejects them and only 0 restores 548.
# A read this run cannot place is left unassigned rather than attached to a model
# that is not its own. test_quant, which supplies the complete annotation, still
# accounts for all 551.
test_ref_guided_ID:
	make clean
	../../LRAA --genome data/SIRVs1-7.genome.fa --bam data/sim.fasta.mm2.bam --norm 0 --gtf data/SIRVs1-7.annot.half.gtf --HiFi --verbose --no_cleanup --no_shuffle_parallel_jobs
	gffcompare -r data/SIRVs1-7.annot.gtf LRAA.LRAA.ref-guided.gtf -o gffcompare
	./eval/eval_sirvs.py data/SIRVs1-7.gene_trans_map gffcompare.LRAA.LRAA.ref-guided.gtf.refmap
	diff gffcompare.LRAA.LRAA.ref-guided.gtf.refmap.full.summary.tsv data/latest_refguided_summary.tsv
	@echo "checking exact terminal coordinates, which the class-code summary cannot see"
	awk '$$1=="SIRV3" && $$3=="transcript" && $$4==4602 && $$5==4762 {n++} END{exit !(n==1)}' LRAA.LRAA.ref-guided.gtf
	awk '$$1=="SIRV6" && $$3=="transcript" && $$4==2359 && $$5==2547 {n++} END{exit !(n==1)}' LRAA.LRAA.ref-guided.gtf
	awk '$$1=="SIRV7" && $$3=="transcript" && $$7=="-" && $$4==1006 && $$5==147925 {n++} END{exit !(n==1)}' LRAA.LRAA.ref-guided.gtf
	@echo "checking total assigned reads, which the class-code summary also cannot see"
	awk 'BEGIN{FS="\t"} !/^#/ && NR>1 {s+=$$4} END{printf "  assigned reads: %.0f\n", s; exit !(s>536.5 && s<537.5)}' LRAA.LRAA.ref-guided.quant.expr


# A supplied annotation model with no read evidence must never become a candidate.
# data/sim.fasta.mm2.no_SIRV311.bam is the standard BAM with the 13 reads overlapping
# SIRV311 (SIRV3:4602-4762, minus strand, 161 bases) removed. SIRV311 is the only
# minus-strand model with exons in that window and nothing spans it across an intron,
# so those are the only reads dropped and every other locus is untouched. SIRV311 stays
# in the guide annotation, so a path can be built for it only from the synthetic read
# standing in for the annotation -- which is not evidence.
#
# The check must read the pre-filtering transcripts: downstream filtering drops a
# zero-support model from the final GTF either way, so asserting on the final output
# would pass even when the model is wrongly reconstructed.
test_unsupported_ref_model:
	make clean
	../../LRAA --genome data/SIRVs1-7.genome.fa --bam data/sim.fasta.mm2.no_SIRV311.bam --norm 0 --gtf data/SIRVs1-7.annot.half.gtf --HiFi --verbose --no_cleanup --no_shuffle_parallel_jobs --include_prelim_transcripts
	@echo "checking that unsupported SIRV311 was never reconstructed"
	! awk '$$1=="SIRV3" && $$3=="transcript" && $$4>=4590 && $$5<=4780 {found=1} END{exit !found}' __PRE-Tx-Filtering.transcripts.gtf
	@echo "checking that the neighboring supported SIRV3 locus still reconstructs"
	awk '$$1=="SIRV3" && $$3=="transcript" && $$4>=8700 && $$5<=10000 {n++} END{exit !(n>=2)}' __PRE-Tx-Filtering.transcripts.gtf


MINIWDL_CFG = ../../miniwdl.test.cfg
MINIWDL_RUN = miniwdl run --cfg $(MINIWDL_CFG)
# The gating test below asserts on task behavior, so it must never be served from the
# call cache: a cache hit skips the task and the assertion checks nothing.
MINIWDL_RUN_UNCACHED = miniwdl run

test_docker: test_docker_direct_cross_gene_EM_gating
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i data/input_refFree_asm.json
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i data/input_refGuided_asm.json
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i data/input_quant_only.json


# Both runs reuse data/input_quant_only.json, which leaves main_chromosomes unset, so
# LRAA.wdl takes the direct (non-scattered) branch; the flag under test is supplied as a
# miniwdl CLI override rather than a new input file, since testing/** is gitignored.
# LRAA writes *.pre-cross-gene-EM.quant.* only when secondary alignments are on, EM is on,
# AND the cross-gene pass is enabled -- see
# _maybe_run_cross_gene_em_for_secondary_alignments in ../../LRAA. LRAA.wdl used to hardcode
# run_final_cross_gene_EM=true on the direct branch while LRAA_runner asserted the file's
# existence off that flag alone, so --no_EM and --no_allow_secondary_alignments each killed
# an otherwise successful run in the task's post-processing block, after all compute was
# done and with every output discarded. The scattered path was never affected: it passes
# run_final_cross_gene_EM=false and gates the merge on allow_secondary_alignments && !no_EM.
# Each run must therefore succeed and report a null preCrossGeneEMQuantExpr.
test_docker_direct_cross_gene_EM_gating:
	$(MINIWDL_RUN_UNCACHED) ../../WDL/LRAA.wdl -i data/input_quant_only.json no_EM=true
	python3 -c 'import json; o=json.load(open("_LAST/outputs.json")); assert o.get("LRAA_wf.mergedQuantExpr"), "no_EM: quant expr missing"; assert o.get("LRAA_wf.preCrossGeneEMQuantExpr") is None, "no_EM: pre-cross-gene-EM output should be absent"; print("  no_EM gating ok")'
	$(MINIWDL_RUN_UNCACHED) ../../WDL/LRAA.wdl -i data/input_quant_only.json allow_secondary_alignments=false
	python3 -c 'import json; o=json.load(open("_LAST/outputs.json")); assert o.get("LRAA_wf.mergedQuantExpr"), "no_secondary: quant expr missing"; assert o.get("LRAA_wf.preCrossGeneEMQuantExpr") is None, "no_secondary: pre-cross-gene-EM output should be absent"; print("  no_secondary gating ok")'




clean:
	rm -rf ./gffcompare* ./LRAA.* ./__* ./*_lraa_wf ./_LAST ./__miniwdl_call_cache ./__miniwdl_download_cache \
		./tx_rescue.* ./reassign_multigene_tracking.* ./*.bam.count
