test:
	../../../LRAA --genome minig-ENSG00000066135.genome.fa --bam minig-ENSG00000066135.mm2.bam -d 2>&1 | tee run.log
	gffcompare -r minig-ENSG00000066135.annot.gtf LRAA.gtf
	cat gffcmp.tracking | cut -f4 | sort > gffcmp.tracking.codes
#	diff gffcmp.tracking.codes data/gffcmp.tracking.codes.expected || \
#	(echo "Output does not match any expected result" && exit 1)
	cat gffcmp.tracking.codes

clean:
	rm -rf ./__*
	rm -rf ./LRAA.* ./gffcmp.*  ./run.log

