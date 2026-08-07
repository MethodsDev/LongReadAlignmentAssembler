

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
test_ref_guided_ID:
	make clean
	../../LRAA --genome data/SIRVs1-7.genome.fa --bam data/sim.fasta.mm2.bam --norm 0 --gtf data/SIRVs1-7.annot.half.gtf --HiFi --verbose --no_cleanup --no_shuffle_parallel_jobs
	gffcompare -r data/SIRVs1-7.annot.gtf LRAA.LRAA.ref-guided.gtf -o gffcompare
	./eval/eval_sirvs.py data/SIRVs1-7.gene_trans_map gffcompare.LRAA.LRAA.ref-guided.gtf.refmap
	diff gffcompare.LRAA.LRAA.ref-guided.gtf.refmap.full.summary.tsv data/latest_refguided_summary.tsv


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

test_docker:
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i data/input_refFree_asm.json
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i data/input_refGuided_asm.json
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i data/input_quant_only.json


clean:
	rm -rf ./gffcompare* ./LRAA.* ./__* ./*_lraa_wf ./_LAST ./__miniwdl_call_cache ./__miniwdl_download_cache
