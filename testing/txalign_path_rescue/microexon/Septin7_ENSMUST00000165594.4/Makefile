.PHONY: test validate clean

PREFIX = LRAA.Septin7_microexon
OUT_PREFIX = $(PREFIX).LRAA.ref-guided

test: clean
	../../../../LRAA \
		--genome data/Septin7.fa \
		--gtf data/Septin7.gtf \
		--bam data/ENSMUST00000165594.4.rescued_original_alignments.bam \
		--no_norm \
		--no_parallelize_contigs \
		--num_threads_per_worker 1 \
		--output_prefix $(PREFIX)
	$(MAKE) validate

validate:
	python ../../validate_txalign_path_rescue.py \
		--case_dir . \
		--output_prefix $(OUT_PREFIX) \
		--expected_read_names data/ENSMUST00000165594.4.rescued_read_names.txt \
		--expected_total_reads 216 \
		--expected_rescue_requested 216 \
		--expected_rescued 216 \
		--expected_microexon chr9:12753-12757:+

clean:
	rm -rf __$(OUT_PREFIX).* \
		$(OUT_PREFIX).* \
		tx_rescue.* \
		reassign_multigene_tracking.* \
		*.quant.secondary_rescue.bam.count \
		__init_clusters.tsv \
		$(PREFIX).* \
		__$(PREFIX).*
