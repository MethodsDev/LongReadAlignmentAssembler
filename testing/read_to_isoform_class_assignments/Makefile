test:
	../../util/classify_read_alignments_to_SQANTI_cats.py  --gtf data/LRAA.gtf --bam data/sim.fasta.mm2.bam --output_prefix test


test_wdl:
	miniwdl run ../../WDL/LRAA-sqanti-like-reads-eval.wdl -i data/inputs.json --no-cache --verbose

clean:
	rm -f ./test.*
