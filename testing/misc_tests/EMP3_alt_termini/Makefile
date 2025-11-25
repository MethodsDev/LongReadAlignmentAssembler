test:
	../../../util/merge_LRAA_GTFs.py --HiFi --genome data/EMP3.genome.fa --gtf data/EMP3.*gtf
	diff LRAA.merged.gtf data/expected_merged.GTF


clean:
	rm -f ./LRAA* 
	rm -rf ./__* 
