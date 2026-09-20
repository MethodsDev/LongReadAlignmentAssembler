.PHONY: test validate clean

PREFIX = LRAA.AT3G28550_repetitive
OUT_PREFIX = $(PREFIX).LRAA.ref-guided

test: clean
	../../../../LRAA \
		--genome data/AT3G28550.fa \
		--gtf data/AT3G28550.gtf \
		--bam data/AT3G28550.2.rescued_original_alignments.bam \
		--no_norm \
		--no_parallelize_contigs \
		--cpu_budget 1 \
		--no_chunk --no_stream_reads \
		--output_prefix $(PREFIX)
	$(MAKE) validate

# The failed_genome / unassigned_quant split was re-pinned in v0.37.0, 0/141 ->
# 55/86.  The OUTCOME is unchanged -- all 141 reads are still requested and all
# 141 still rescued, against the same expected read-name list -- so what moved is
# the route a read takes to rescue, not whether it gets there.
#
# Cause: v0.37.0 enables infer_PolyA for every platform, so PolyA vertices now
# exist in the graph.  At this locus 55 of the 141 reads acquire a SPACER in
# their graph path and are handed to transcriptome rescue as "no graph path"
# rather than reaching quant and being rescued as "unassigned"
# (run.log: `total=141 kept=86 discard_spacer=55 no_path=55`).
#
# This locus is deliberately repetitive and amplifies the effect.  Measured on
# SGNex MCF7 chr1:1-4Mb, real human ONT, the same change moves spacer discards
# from 77 to 85 of 2,129 reads (+0.4 points), so the 39-point swing here is a
# property of this fixture rather than of the change.  If this split moves again,
# check the OUTCOME line first: 141/141 rescued is the assertion that matters.
validate:
	python ../../validate_txalign_path_rescue.py \
		--case_dir . \
		--output_prefix $(OUT_PREFIX) \
		--expected_read_names data/AT3G28550.2.rescued_read_names.txt \
		--expected_total_reads 141 \
		--expected_rescue_requested 141 \
		--expected_rescued 141 \
		--expected_rescue_requested_failed_genome 55 \
		--expected_rescue_requested_unassigned_quant 86

clean:
	rm -rf __$(OUT_PREFIX).* \
		$(OUT_PREFIX).* \
		tx_rescue.* \
		__init_clusters.tsv \
		$(PREFIX).* \
		__$(PREFIX).*
