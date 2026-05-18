.PHONY: test validate clean

PREFIX = LRAA.AT4G05320_repetitive
OUT_PREFIX = $(PREFIX).LRAA.ref-guided

test: clean
	../../../../LRAA \
		--genome data/AT4G05320.fa \
		--gtf data/AT4G05320.gtf \
		--bam data/AT4G05320.1.rescued_original_alignments.bam \
		--no_norm \
		--no_parallelize_contigs \
		--num_threads_per_worker 1 \
		--output_prefix $(PREFIX)
	$(MAKE) validate

validate:
	python ../../validate_txalign_path_rescue.py \
		--case_dir . \
		--output_prefix $(OUT_PREFIX) \
		--expected_read_names data/AT4G05320.1.rescued_read_names.txt \
		--expected_total_reads 70 \
		--expected_rescue_requested 70 \
		--expected_rescued 70 \
		--expected_rescue_requested_failed_genome 70 \
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
