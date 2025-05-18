test: test_single test_multi


test_single:
	../../util/SQANTI-like_cats_for_reads_or_isoforms.py  --ref_gtf data/LRAA.gtf --input_bam data/sim.fasta.mm2.bam --output_prefix test.bam
	../../util/SQANTI-like_cats_for_reads_or_isoforms.py  --ref_gtf data/SIRVs1-7.annot.gtf --input_gtf data/LRAA.gtf --output_prefix test.gtf


test_multi:
	../../util/misc/plot_SQANTI_cats.summarize_mult_samples.Rscript --output_prefix test.multi --sample_stats data/example_summaries/*tsv.gz
	../../util/misc/plot_feature_lengths.summarize_mult_samples.Rscript --sample_iso_cats data/example_summaries/indiv_cats/*tsv.gz --output_prefix test.lengths

test_wdl: test_wdl_single test_wdl_multi


test_wdl_single:
	miniwdl run ../../WDL/LRAA-sqanti-like-reads-eval.wdl -i data/inputs.json --no-cache --verbose


test_wdl_multi:
	miniwdl run ../../WDL/LRAA-sqanti-like-reads-eval.summarize_mult_samples.wdl -i data/inputs.mult_samples_summary.json --no-cache --verbose


clean:
	rm -f ./test.* ./_LAST
	rm -rf ./*sqanti_like_reads_eval_wf ./202*


