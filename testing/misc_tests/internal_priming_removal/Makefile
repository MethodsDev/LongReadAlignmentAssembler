

test: test_remove_iprime test_keep_iprime




test_remove_iprime:
	../../../LRAA --bam data/iprime_reads.bam --genome data/genome_region.fa --output_prefix test_remove_iprime --config_update data/config_adjust.json --HiFi --verbose --no_cleanup
	@if grep -v '^#' test_remove_iprime.LRAA.ref-free.gtf | grep -q .; then \
		echo "Error: test_remove_iprime.LRAA.ref-free.gtf has non-comment content!"; \
		exit 1; \
	else \
		echo "OK: test_remove_iprime.LRAA.ref-free.gtf has no non-comment content."; \
	fi


test_keep_iprime:
	../../../LRAA --bam data/iprime_reads.bam --genome data/genome_region.fa --gtf data/use_as_ref.gtf  --output_prefix test_keep_iprime --config_update data/config_adjust.json --HiFi --verbose --no_cleanup
	@if [ -s test_keep_iprime.LRAA.ref-guided.gtf ]; then \
		echo "OK: test_keep_iprime.LRAA.ref-guided.gtf is not empty."; \
	else \
		echo "Error: test_keep_iprime.LRAA.ref-guided.gtf is empty or missing!"; \
		exit 1; \
	fi



clean:
	rm -rf ./test_*  ./iprime* ./__*
