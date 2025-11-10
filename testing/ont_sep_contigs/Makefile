test:
	../../LRAA --genome minigenome.fa --bam ONT_sim_alignments.mm2.bam -d --verbose --no_cleanup 2>&1 | tee run.log
	gffcompare -r minigenome.gtf LRAA.gtf
	cat gffcmp.tracking | cut -f4 | sort > gffcmp.tracking.codes
	diff gffcmp.tracking.codes data/gffcmp.tracking.codes.expected || \
	(echo "Output does not match any expected result" && exit 1)

clean:
	rm -rf ./__*
	rm -rf ./LRAA.* ./gffcmp.* ./ONT_sim_alignments.mm2.bam.count ./ONT_sim_alignments.mm2.bam.SS* ./ONT_sim_alignments.mm2.bam.norm* ./run.log

