# One ONT locus, 1,145 reads over a 22-exon transcript, asserted to reconstruct
# as gffcompare class `=` -- structurally identical to ENST00000372396.  The
# parent's corpus is many genes shallowly; this is one gene with enough exons
# that a single missed junction changes the code to c, j or k.
#
# The directory arrived in "incorporate progress monitoring" as a reproducer,
# with this comparison commented out against a baseline that was never written
# and no target recursing into it, so it ran nothing for ten months while an
# output rename swept through it.  The baseline is now captured and ../Makefile
# recurses here.

# See ../Makefile: without pipefail, tee's exit status hides a crashed LRAA.
SHELL := /bin/bash
.SHELLFLAGS := -o pipefail -c

test:
	../../../LRAA --genome minig-ENSG00000066135.genome.fa --bam minig-ENSG00000066135.mm2.bam -d --no_chunk --no_stream_reads 2>&1 | tee run.log
	gffcompare -r minig-ENSG00000066135.annot.gtf LRAA.LRAA.ref-free.gtf
	cat gffcmp.tracking | cut -f4 | sort > gffcmp.tracking.codes
	diff gffcmp.tracking.codes data/gffcmp.tracking.codes.expected || \
	(echo "Output does not match the expected result" && exit 1)

clean:
	rm -rf ./__*
	rm -rf ./LRAA.* ./gffcmp.*  ./run.log
