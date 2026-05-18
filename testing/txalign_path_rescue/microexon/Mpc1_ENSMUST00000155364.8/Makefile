.PHONY: test validate clean

PREFIX = LRAA.Mpc1_microexon
OUT_PREFIX = $(PREFIX).LRAA.ref-guided

test: clean
	../../../../LRAA \
		--genome data/Mpc1.fa \
		--gtf data/Mpc1.gtf \
		--bam data/ENSMUST00000155364.8.rescued_original_alignments.bam \
		--no_norm \
		--no_parallelize_contigs \
		--num_threads_per_worker 1 \
		--output_prefix $(PREFIX)
	$(MAKE) validate

validate:
	python ../../validate_txalign_path_rescue.py \
		--case_dir . \
		--output_prefix $(OUT_PREFIX) \
		--expected_read_names data/ENSMUST00000155364.8.rescued_read_names.txt \
		--expected_total_reads 570 \
		--expected_rescue_requested 570 \
		--expected_rescued 570 \
		--expected_microexon chr17:10654-10657:+

clean:
	rm -rf __$(OUT_PREFIX).* \
		$(OUT_PREFIX).* \
		tx_rescue.* \
		reassign_multigene_tracking.* \
		*.quant.secondary_rescue.bam.count \
		__init_clusters.tsv \
		$(PREFIX).* \
		__$(PREFIX).*
