test:
	../../../LRAA --genome data/minigenome.fa --gtf data/minigenome.annot.gtf --bam data/sim.fasta.mm2.cs.hq.bam --quant_only --tag_bam --HiFi --verbose --no_cleanup
	diff LRAA.quant.expr data/perfect_quant.tsv


clean:
	rm -rf ./LRAA.quant.* ./LRAA.tagged.bam ./__alignment_cache
