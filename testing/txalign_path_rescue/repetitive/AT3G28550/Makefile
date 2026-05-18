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
		--num_threads_per_worker 1 \
		--output_prefix $(PREFIX)
	$(MAKE) validate

validate:
	python ../../validate_txalign_path_rescue.py \
		--case_dir . \
		--output_prefix $(OUT_PREFIX) \
		--expected_read_names data/AT3G28550.2.rescued_read_names.txt \
		--expected_total_reads 141 \
		--expected_rescue_requested 141 \
		--expected_rescued 141 \
		--expected_rescue_requested_failed_genome 141 \
		--expected_rescue_requested_unassigned_quant 0

clean:
	rm -rf __$(OUT_PREFIX).* \
		$(OUT_PREFIX).* \
		tx_rescue.* \
		reassign_multigene_tracking.* \
		*.quant.secondary_rescue.bam.count \
		__init_clusters.tsv \
		$(PREFIX).* \
		__$(PREFIX).*
