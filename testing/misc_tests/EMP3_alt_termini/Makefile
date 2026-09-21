test: check_merged_coords check_merged_3prime_attributes

LRAA.merged.gtf:
	../../../util/merge_LRAA_GTFs.py --HiFi --genome data/EMP3.genome.fa --gtf data/EMP3.*gtf

check_merged_coords: LRAA.merged.gtf
	sort LRAA.merged.gtf | cut -f1-7 > LRAA.merged.gtf.sorted
	diff LRAA.merged.gtf.sorted data/expected_merged.GTF.sorted

# The coordinate diff above cuts to columns 1-7, so it says nothing about the
# attributes -- which is how the merge came to emit models with no PAS,
# PAS_offset or InternalPriming at all while this test stayed green.  These are
# 14 per-cluster gtfs of one alt-termini locus, so the merged models are exactly
# the case where boundary evidence has to survive reconstruction.
#
# Asserted as properties rather than as a golden column 9: a second expected
# file would have to be regenerated for every attribute ever added, and a
# regenerated golden file records whatever the code did that day.
check_merged_3prime_attributes: LRAA.merged.gtf
	@awk -F'\t' '$$3 == "transcript" { \
	    n++; \
	    if ($$9 !~ /PAS "/)             { no_pas++ } \
	    if ($$9 !~ /InternalPriming "/) { no_ip++ } \
	    if ($$9 ~ /PAS "[ACGT]/)        { motif++ } \
	    if ($$9 ~ /PolyA "True"/) { polyA++; if ($$9 !~ /PolyA_read_count "/) no_count++ } \
	  } END { \
	    if (n == 0)        { print "FAIL: no transcript records in merged gtf"; exit 1 } \
	    if (no_pas)        { print "FAIL: " no_pas " of " n " merged models carry no PAS attribute"; exit 1 } \
	    if (no_ip)         { print "FAIL: " no_ip " of " n " merged models carry no InternalPriming attribute"; exit 1 } \
	    if (motif == 0)    { print "FAIL: no merged model found a PAS hexamer; annotation ran but reported nothing"; exit 1 } \
	    if (polyA == 0)    { print "FAIL: no merged model retained a PolyA boundary"; exit 1 } \
	    if (no_count)      { print "FAIL: " no_count " of " polyA " PolyA models lost PolyA_read_count"; exit 1 } \
	    print "ok: " n " merged models annotated, " motif " with a PAS hexamer, " polyA " with a retained PolyA boundary" \
	  }' LRAA.merged.gtf


clean:
	rm -f ./LRAA* 
	rm -rf ./__* 


igv:
	igv.sh -g data/EMP3.genome.fa LRAA.merged.gtf
