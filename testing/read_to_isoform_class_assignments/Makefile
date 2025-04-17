test: test_single test_multi


test_single:
	../../util/classify_read_alignments_to_SQANTI_cats.py  --gtf data/LRAA.gtf --bam data/sim.fasta.mm2.bam --output_prefix test


test_multi:
	../../util/misc/plot_SQANTI_cats.summarize_mult_samples.Rscript --output_prefix test.multi --sample_stats data/example_summaries/SM-*


test_wdl:
	miniwdl run ../../WDL/LRAA-sqanti-like-reads-eval.wdl -i data/inputs.json --no-cache --verbose

clean:
	rm -f ./test.* ./_LAST
	rm -rf ./*sqanti_like_reads_eval_wf
