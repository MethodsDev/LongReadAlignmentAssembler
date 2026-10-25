# LRAA --pasa_mode on synthetic staggered reads (fixture/, made by build_fixture.py):
# spliced reads sharing an intron chain into full-length models, reads overlapping
# only within an exon or joined only by an unspliced read stay separate, and short
# intron overhangs don't break the chain.

PREFIX = pasa_mode

test:
	../../LRAA --pasa_mode --bam fixture/reads.bam --genome fixture/genome.fa \
	    --cpu_budget 2 --output_prefix $(PREFIX)
	python3 check_pasa_mode.py fixture/expected.tsv $(PREFIX).LRAA.ref-free.gtf

clean:
	rm -rf $(PREFIX).LRAA.* __$(PREFIX).* __init_clusters.tsv __read_tracking __alignment_cache *.mapped_primary.count
