# Regression test for min_monoexonic_supporting_cells's cell-count evidence source.
#
# Investigated as the single_cell_monoexonic_cell_filter_deferral_plan.md scoping
# document's proposed fix: deferring the filter (LRAA:6008, inside
# run_transcript_assembly) to run after final quant, on the theory that its cell
# counts came from the coverage-normalized splice-graph bam rather than the full one.
# That diagnosis does not hold: filter_novel_monoexonic_isoforms_by_min_cells reads
# transcript.get_multipaths_evidence_assigned(), populated by the draft-mode
# quantify() calls (LRAA:5742/5835/5867) against mp_all_counter -- built at
# LRAA:5726-5735 from bam_file_for_quant, the full un-normalized bam, unconditionally
# in both --stream_reads and default mode (confirmed via `git blame`: that argument
# has been bam_file_for_quant since 2025-02-19, not a recent change). No file
# patching was applied; this directory exists to pin the correct behavior down as a
# regression test rather than to exercise a fix.
#
# data/chrTest.{full,sg}.bam (see data/build_fixture.py) engineer the exact scenario
# deterministically via --bam_for_sg/--no_norm rather than the real normalizer's
# probabilistic window sampling: locus A has 8 distinct cells on the full bam but
# only 2 on the stand-in splice-graph bam (default min_monoexonic_supporting_cells is
# 5), locus B has genuinely 2 distinct cells on both (must stay filtered), locus C
# has 8 distinct cells on both (always-kept, unaffected control).
#
# --HiFi is deliberately NOT passed: it gates monoexonic retention on inferred
# TSS/PolyA annotation, which this fixture's exact-CIGAR, unclipped synthetic reads
# never produce, and that gate is unrelated to what this test checks.

LRAA = $(abspath ../../../LRAA)
CHECK = python3 check_fixture.py

test: test_nonstreaming test_streaming

test_nonstreaming:
	rm -rf run_nonstreaming
	mkdir run_nonstreaming
	cd run_nonstreaming && $(LRAA) \
	    --genome ../data/chrTest.genome.fa --bam ../data/chrTest.full.bam \
	    --bam_for_sg ../data/chrTest.sg.bam --no_norm \
	    --no_chunk --no_stream_reads \
	    --debug --no_cleanup --no_parallelize_contigs \
	    --output_prefix nonstreaming > run.log 2>&1
	$(CHECK) --label non-streaming \
	    --quant_expr run_nonstreaming/nonstreaming.LRAA.ref-free.quant.expr \
	    --tracking run_nonstreaming/nonstreaming.LRAA.ref-free.quant.tracking.gz
	grep -q 'FILTERING novel monoexonic transcript t:chrTest:+:comp-2:iso-1: supported by 2 cells' run_nonstreaming/run.log
	! grep -q 'FILTERING novel monoexonic transcript t:chrTest:+:comp-1:iso-1' run_nonstreaming/run.log
	! grep -q 'FILTERING novel monoexonic transcript t:chrTest:+:comp-3:iso-1' run_nonstreaming/run.log

test_streaming:
	rm -rf run_streaming
	mkdir run_streaming
	cd run_streaming && $(LRAA) \
	    --genome ../data/chrTest.genome.fa --bam ../data/chrTest.full.bam \
	    --bam_for_sg ../data/chrTest.sg.bam --no_norm --stream_reads \
	    --no_rescue_unassigned_reads_via_transcriptome_alignment \
	    --no_chunk \
	    --debug --no_cleanup --no_parallelize_contigs \
	    --output_prefix streaming > run.log 2>&1
	$(CHECK) --label streaming \
	    --quant_expr run_streaming/streaming.LRAA.ref-free.quant.expr \
	    --tracking run_streaming/streaming.LRAA.ref-free.quant.tracking.gz
	grep -q 'FILTERING novel monoexonic transcript t:chrTest:+:comp-2:iso-1: supported by 2 cells' run_streaming/run.log
	! grep -q 'FILTERING novel monoexonic transcript t:chrTest:+:comp-1:iso-1' run_streaming/run.log
	! grep -q 'FILTERING novel monoexonic transcript t:chrTest:+:comp-3:iso-1' run_streaming/run.log

test_wdls:

clean:
	rm -rf run_nonstreaming run_streaming __*
