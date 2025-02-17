test:
	../../LRAA --genome minigenome.fa --bam ONT_sim_alignments.mm2.bam --LowFi
	gffcompare -r minigenome.gtf LRAA.gtf
	cat gffcmp.tracking | cut -f4 | sort > gffcmp.tracking.codes
	diff gffcmp.tracking.code __data/gffcmp.tracking.codes.expected


clean:
	rm -rf ./__chckpts
	rm -f ./LRAA.* ./gffcmp.* ./ONT_sim_alignments.mm2.bam.count ./ONT_sim_alignments.mm2.bam.SS* ./ONT_sim_alignments.mm2.bam.norm*

