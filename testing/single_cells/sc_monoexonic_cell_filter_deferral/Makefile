# Regression test for min_monoexonic_supporting_cells's cell-count evidence source.
#
# The gate counts DISTINCT CELLS, not reads: filter_novel_monoexonic_isoforms_by_min_cells
# reads transcript.get_multipaths_evidence_assigned() and pulls a barcode off each read
# name, so its evidence is whichever bam the draft quantification read. That bam now
# mirrors the final quant's estimating pass -- bam_file_for_sg under --stream_reads,
# bam_file_for_quant without it -- so THE GATE'S EVIDENCE IS MODE-DEPENDENT, and this
# directory pins both arms rather than asserting one answer.
#
# Cells are the one quantity XW weighting cannot repair. A read retained at p=1/15 stands
# for fifteen reads but for exactly one cell, and you cannot recover which distinct
# barcodes existed from reads that were discarded. So under streaming the count is of
# cells surviving in the thinned evidence.
#
# That is an ACCEPTED TRADEOFF, not a proof of safety. The reasoning is that 5 distinct
# cells is a low bar at a locus retaining a thousand reads or more, and that it scales --
# sourcing cells from the full bam would cost discovery the full-bam read that reading the
# normalized one saves. But retained DEPTH does not guarantee per-cell retention: a barcode
# contributing few reads to a locus can lose all of them to thinning however deep the locus
# is overall, and sparse barcodes are exactly the ones a min-cells bar is deciding about.
# The exposure named at the bottom of this comment is the accurate statement of the risk.
#
# What this fixture engineers is therefore a WORST CASE, not a typical one, and does so
# deterministically via --bam_for_sg/--no_norm rather than the real normalizer's
# probabilistic window sampling (see data/build_fixture.py):
#
#   locus A (comp-1)  8 distinct cells on the full bam, 2 on the stand-in sg bam
#                     -> kept without streaming (8 >= 5), filtered with it (2 < 5)
#   locus B (comp-2)  genuinely 2 distinct cells on both -> filtered in both arms
#   locus C (comp-3)  8 distinct cells on both -> kept in both arms, the control that
#                     isolates thinning as the only variable
#
# The residual exposure comp-1 stands for is narrow but real: a lowly-expressed monoexonic
# novel isoform inside a deeply-sequenced locus, whose handful of reads come from distinct
# cells. Coverage normalization's scarce-evidence exemption cannot protect it, because that
# exemption is keyed on junctions and a monoexonic model has none. Raising
# --normalize_max_cov_level retains more reads and is the lever if this ever bites.
#
# --bam_for_sg is stamped XW:f:1.0 because LRAA refuses a supplied splice-graph bam that
# carries no XW, that flag naming already-normalized evidence. 1.0 is what normalization
# writes for a read it kept whole, so the stamp declares "nothing thinned" and changes no
# number here; see util/misc/stamp_unit_xw.py.
#
# --HiFi is deliberately NOT passed: it gates monoexonic retention on inferred TSS/PolyA
# annotation, which this fixture's exact-CIGAR, unclipped synthetic reads never produce,
# and that gate is unrelated to what this test checks.

LRAA = $(abspath ../../../LRAA)
CHECK = python3 check_fixture.py
STAMP = $(abspath ../../../util/misc/stamp_unit_xw.py)

test: test_nonstreaming test_streaming

test_nonstreaming:
	rm -rf run_nonstreaming
	mkdir run_nonstreaming
	$(STAMP) data/chrTest.sg.bam run_nonstreaming/sg.tagged.bam
	cd run_nonstreaming && $(LRAA) \
	    --genome ../data/chrTest.genome.fa --bam ../data/chrTest.full.bam \
	    --bam_for_sg sg.tagged.bam --no_norm \
	    --no_chunk --no_stream_reads \
	    --debug --no_cleanup --no_parallelize_contigs \
	    --output_prefix nonstreaming > run.log 2>&1
	$(CHECK) --label non-streaming --locus_a retained \
	    --quant_expr run_nonstreaming/nonstreaming.LRAA.ref-free.quant.expr \
	    --tracking run_nonstreaming/nonstreaming.LRAA.ref-free.quant.tracking.gz
	grep -q 'FILTERING novel monoexonic transcript t:chrTest:+:comp-2:iso-1: supported by 2 cells' run_nonstreaming/run.log
	! grep -q 'FILTERING novel monoexonic transcript t:chrTest:+:comp-1:iso-1' run_nonstreaming/run.log
	! grep -q 'FILTERING novel monoexonic transcript t:chrTest:+:comp-3:iso-1' run_nonstreaming/run.log

# Streaming's draft quant reads the splice-graph bam, so comp-1's cells are the 2 present
# there rather than the 8 in the library, and it is filtered. That difference from the arm
# above is the whole point of pinning both.
test_streaming:
	rm -rf run_streaming
	mkdir run_streaming
	$(STAMP) data/chrTest.sg.bam run_streaming/sg.tagged.bam
	cd run_streaming && $(LRAA) \
	    --genome ../data/chrTest.genome.fa --bam ../data/chrTest.full.bam \
	    --bam_for_sg sg.tagged.bam --no_norm --stream_reads \
	    --no_rescue_unassigned_reads_via_transcriptome_alignment \
	    --no_chunk \
	    --debug --no_cleanup --no_parallelize_contigs \
	    --output_prefix streaming > run.log 2>&1
	$(CHECK) --label streaming --locus_a filtered \
	    --quant_expr run_streaming/streaming.LRAA.ref-free.quant.expr \
	    --tracking run_streaming/streaming.LRAA.ref-free.quant.tracking.gz
	grep -q 'FILTERING novel monoexonic transcript t:chrTest:+:comp-2:iso-1: supported by 2 cells' run_streaming/run.log
	grep -q 'FILTERING novel monoexonic transcript t:chrTest:+:comp-1:iso-1: supported by 2 cells' run_streaming/run.log
	! grep -q 'FILTERING novel monoexonic transcript t:chrTest:+:comp-3:iso-1' run_streaming/run.log

test_wdls:

clean:
	rm -rf run_nonstreaming run_streaming __*
