

test: test_remove_iprime test_keep_iprime




test_remove_iprime:
	../../../LRAA --bam data/iprime_reads.bam --genome data/genome_region.fa --output_prefix test_remove_iprime
	@if [ -s test_remove_iprime.gtf ]; then \
		echo "Error: test_remove_iprime.gtf is not empty!"; \
		exit 1; \
	else \
		echo "OK: test_remove_iprime.gtf is empty."; \
	fi


test_keep_iprime:
	../../../LRAA --bam data/iprime_reads.bam --genome data/genome_region.fa --gtf data/use_as_ref.gtf  --output_prefix test_keep_iprime
	@if [ -s test_keep_iprime.gtf ]; then \
		echo "OK: test_keep_iprime.gtf is not empty."; \
	else \
		echo "Error: test_keep_iprime.gtf is empty or missing!"; \
		exit 1; \
	fi



clean:
	rm -rf ./test_*  ./iprime* ./__*

