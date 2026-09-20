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

# data/gffcmp.tracking.codes.expected was re-pinned in v0.37.0, 10 '=' -> 11 '='
# (the single 'u' is unchanged).  This corpus is ONT, and v0.37.0 turned on
# infer_PolyA together with strip_polyA_on_proximal_window: a read carrying a
# polyA tail can now define a polyA site, where before the unstripped tail left
# residual soft clipping and max_soft_clip_at_PolyA (0) barred it.  The result
# here is that all 11 reference transcripts are now recovered with exact
# structure rather than 10, which is the intended effect of that change rather
# than drift -- re-measure before re-pinning again if this moves.
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
