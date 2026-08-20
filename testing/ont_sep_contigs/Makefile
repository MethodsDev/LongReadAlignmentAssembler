# `LRAA ... | tee run.log` reports tee's exit status, which is always 0, so a
# crashed run reads as a pass here.  What then decides the result is whether
# gffcompare happens to find a usable gtf -- present and complete from a run
# that died later, or left over from a previous one.  pipefail puts LRAA's own
# status back in charge.
SHELL := /bin/bash
.SHELLFLAGS := -o pipefail -c

DIRS = minig-ENSG00000066135

.PHONY: test test_minigenome test_subdirs igv clean

test: test_minigenome test_subdirs

test_subdirs:
	@for dir in $(DIRS); do \
		$(MAKE) -C $$dir test || exit $$?; \
	done

test_minigenome:
	../../LRAA --genome minigenome.fa --bam ONT_sim_alignments.mm2.bam -d --verbose --no_cleanup --no_chunk --no_stream_reads 2>&1 | tee run.log
	gffcompare -r minigenome.gtf LRAA.LRAA.ref-free.gtf
	cat gffcmp.tracking | cut -f4 | sort > gffcmp.tracking.codes
	diff gffcmp.tracking.codes data/gffcmp.tracking.codes.expected || \
	(echo "Output does not match any expected result" && exit 1)

igv:
	igv.sh -g minigenome.fa minigenome.gtf,LRAA.LRAA.ref-free.bed,ONT_sim_alignments.mm2.bam


clean:
	@for dir in $(DIRS); do \
		$(MAKE) -C $$dir clean || exit $$?; \
	done
	rm -rf ./__*
	rm -rf ./LRAA.* ./gffcmp.* ./ONT_sim_alignments.mm2.bam.count ./ONT_sim_alignments.mm2.bam.SS* ./ONT_sim_alignments.mm2.bam.norm* ./run.log
