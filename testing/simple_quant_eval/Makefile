test:
	../../LRAA --genome data/minigenome.fa --gtf data/minigenome.annot.gtf --bam data/sim.fasta.mm2.cs.hq.bam --quant_only
	diff LRAA.quant.expr data/perfect_quant.tsv


clean:
	rm -f ./LRAA.quant.*
