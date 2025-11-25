test:
	../../../util/merge_LRAA_GTFs.py --HiFi --genome data/EMP3.genome.fa --gtf data/EMP3.*gtf
	sort LRAA.merged.gtf > LRAA.merged.gtf.sorted
	diff LRAA.merged.gtf.sorted data/expected_merged.GTF.sorted


clean:
	rm -f ./LRAA* 
	rm -rf ./__* 


igv:
	igv.sh -g data/EMP3.genome.fa LRAA.merged.gtf
