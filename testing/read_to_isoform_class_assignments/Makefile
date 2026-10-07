test: test_single test_multi test_code test_parallel_parity


test_single:
	../../util/SQANTI-like_cats_for_reads_or_isoforms.py  --ref_gtf data/LRAA.gtf --input_bam data/sim.fasta.mm2.bam --output_prefix test.bam
	../../util/SQANTI-like_cats_for_reads_or_isoforms.py  --ref_gtf data/SIRVs1-7.annot.gtf --input_gtf data/LRAA.gtf --output_prefix test.gtf
	../../util/misc/plot_SQANTI_cats.Rscript data/misc_iso_cats/control_bc01.iso_cats.summary_counts.tsv.gz test.misc.iso_cats.summary_counts.pdf

test_multi:
	../../util/misc/plot_SQANTI_cats.summarize_mult_samples.Rscript --output_prefix test.multi --sample_stats data/example_summaries/*tsv.gz
	../../util/misc/plot_feature_lengths.summarize_mult_samples.Rscript --sample_iso_cats data/example_summaries/indiv_cats/*tsv.gz --output_prefix test.lengths

test_code:
	pytest ./sqanti-like_code_test.py 

# --CPU classifies one contig per worker and concatenates the parts; it must
# reproduce the serial run exactly: rows, summary order, and every bam record.
test_parallel_parity:
	python3 make_parallel_parity_fixture.py data/sim.fasta.mm2.bam test.par.fixture.bam data/SIRVs1-7.annot.gtf test.par.ref_noSIRV7.gtf
	../../util/SQANTI-like_cats_for_reads_or_isoforms.py --ref_gtf test.par.ref_noSIRV7.gtf --input_bam test.par.fixture.bam --output_prefix test.par.cpu1 --CPU 1
	../../util/SQANTI-like_cats_for_reads_or_isoforms.py --ref_gtf test.par.ref_noSIRV7.gtf --input_bam test.par.fixture.bam --output_prefix test.par.cpu4 --CPU 4
	cmp test.par.cpu1.iso_cats.tsv test.par.cpu4.iso_cats.tsv
	cmp test.par.cpu1.iso_cats.summary_counts.tsv test.par.cpu4.iso_cats.summary_counts.tsv
	samtools view -h --no-PG test.par.cpu1.iso_cats.bam > test.par.cpu1.sam
	samtools view -h --no-PG test.par.cpu4.iso_cats.bam > test.par.cpu4.sam
	cmp test.par.cpu1.sam test.par.cpu4.sam
	../../util/SQANTI-like_cats_for_reads_or_isoforms.py --ref_gtf test.par.ref_noSIRV7.gtf --input_bam test.par.fixture.bam --output_prefix test.par.auto
	cmp test.par.cpu1.iso_cats.tsv test.par.auto.iso_cats.tsv
	cp test.par.fixture.bam test.par.noindex.bam
	../../util/SQANTI-like_cats_for_reads_or_isoforms.py --ref_gtf test.par.ref_noSIRV7.gtf --input_bam test.par.noindex.bam --output_prefix test.par.auto_noindex
	cmp test.par.cpu1.iso_cats.tsv test.par.auto_noindex.iso_cats.tsv
	../../util/SQANTI-like_cats_for_reads_or_isoforms.py --ref_gtf test.par.ref_noSIRV7.gtf --input_bam test.par.fixture.bam --output_prefix test.par.gz4 --CPU 4 --gzip_tsv
	gzip -t test.par.gz4.iso_cats.tsv.gz
	zcat test.par.gz4.iso_cats.tsv.gz | cmp - test.par.cpu1.iso_cats.tsv
	../../util/SQANTI-like_cats_for_reads_or_isoforms.py --ref_gtf test.par.ref_noSIRV7.gtf --input_bam test.par.fixture.bam --output_prefix test.par.gz1 --CPU 1 --gzip_tsv
	zcat test.par.gz1.iso_cats.tsv.gz | cmp - test.par.cpu1.iso_cats.tsv
	@echo "ok: --CPU 4, --CPU auto (indexed and not) and --gzip_tsv (serial and parallel) match the serial run"

include ../lraa_test_docker.mk

# MINIWDL_CFG and MINIWDL_RUN come from the include above, and both runs pass the
# :testing lraa-sc image built from this commit. Without docker_sc the WDLs fall
# back to their default (lraa-sc:latest), which tests the last release instead.
# check_test_image_revision FIRST: make runs prerequisites in declaration order.
test_wdls: check_test_image_revision test_wdl_single test_wdl_multi

test_wdl_single:
	$(MINIWDL_RUN) ../../WDL/LRAA-sqanti-like-reads-eval.wdl -i data/inputs.json $(LRAA_WDL_DOCKER_SC) --verbose --dir wdl_single_out
	$(MAKE) check_wdl_single

test_wdl_multi:
	$(MINIWDL_RUN) ../../WDL/LRAA-sqanti-like-reads-eval.summarize_mult_samples.wdl -i data/inputs.mult_samples_summary.json $(LRAA_WDL_DOCKER_SC) --verbose --dir wdl_multi_out
	$(MAKE) check_wdl_multi

# A workflow can finish with an empty plot or no rows (the R step once failed only
# after every read had been classified), so assert on the outputs themselves.
check_wdl_single:
	@out=wdl_single_out/_LAST/out; \
	for f in iso_cats_tsv iso_cats_summary_counts_tsv iso_cats_summary_counts_pdf; do \
	  test -s $$out/$$f/* || { echo "FAIL: $$f missing or empty"; exit 1; }; \
	done; \
	n=`zcat $$out/iso_cats_tsv/*.gz | tail -n +2 | wc -l`; \
	test $$n -gt 0 || { echo "FAIL: iso_cats_tsv has no rows"; exit 1; }; \
	head -c 4 $$out/iso_cats_summary_counts_pdf/* | grep -q '%PDF' || { echo "FAIL: summary plot is not a PDF"; exit 1; }; \
	echo "ok: sqanti-like WDL classified $$n reads and wrote the summary PDF"

check_wdl_multi:
	@out=wdl_multi_out/_LAST/out; \
	n=`find $$out -name '*.pdf' -size +0 | wc -l`; \
	test $$n -gt 0 || { echo "FAIL: multi-sample summary wrote no PDF"; exit 1; }; \
	echo "ok: multi-sample summary wrote $$n PDF(s)"


clean:
	rm -f ./test.* ./_LAST
	rm -rf ./*sqanti_like_reads_eval_wf ./202* ./wdl_single_out ./wdl_multi_out ./__miniwdl_call_cache ./__miniwdl_download_cache
