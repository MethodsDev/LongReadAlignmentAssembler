
test: test_serial test_parallel


# Run serial across (contig,strand) using component-level multiprocessing only
test_serial:
	../../LRAA --bam pacbioccsminimap.minig.PBLR.bam --genome minig.fa --no_norm --output_prefix LRAA.serial


# Run with new contig/strand-level parallelization (plus component-level inside each worker)
test_parallel:
	../../LRAA --bam pacbioccsminimap.minig.PBLR.bam --genome minig.fa --no_norm --parallelize_contigs --output_prefix LRAA.par


# Quick timing harness to compare serial vs parallelized contig processing
# Produces: bench.summary.tsv with real/user/sys times for both modes
bench: bench_clean
	( time -p ../../LRAA --bam pacbioccsminimap.minig.PBLR.bam --genome minig.fa --no_norm --CPU 15 --output_prefix LRAA.serial ) 2> serial.time.txt
	( time -p ../../LRAA --bam pacbioccsminimap.minig.PBLR.bam --genome minig.fa --no_norm --CPU 15 --parallelize_contigs --output_prefix LRAA.par ) 2> parallel.time.txt
	@echo -e "mode\tCPU\tparallelize_contigs\treal\tuser\tsys" > bench.summary.tsv
	@SER_REAL=`grep '^real' serial.time.txt | awk '{print $$2}'`; \
	 SER_USER=`grep '^user' serial.time.txt | awk '{print $$2}'`; \
	 SER_SYS=`grep '^sys'  serial.time.txt | awk '{print $$2}'`; \
	 echo -e "serial\t15\tno\t$${SER_REAL}\t$${SER_USER}\t$${SER_SYS}" >> bench.summary.tsv
	@PAR_REAL=`grep '^real' parallel.time.txt | awk '{print $$2}'`; \
	 PAR_USER=`grep '^user' parallel.time.txt | awk '{print $$2}'`; \
	 PAR_SYS=`grep  '^sys' parallel.time.txt | awk '{print $$2}'`; \
	 echo -e "parallel\t15\tyes\t$${PAR_REAL}\t$${PAR_USER}\t$${PAR_SYS}" >> bench.summary.tsv
	@echo "Wrote bench.summary.tsv" && cat bench.summary.tsv


test_wdls: test_wdl_quant_only test_wdl_ref_guided test_wdl_ref_free


##########
# wdls

test_wdl_quant_only: wdl_clean
	 miniwdl run ../../WDL/LRAA.wdl -i wf.quant-only.json  --verbose --no-cache

test_wdl_ref_guided: wdl_clean
	miniwdl run ../../WDL/LRAA.wdl -i wf.ref-guided.json  --verbose --no-cache

test_wdl_ref_free: wdl_clean
	miniwdl run ../../WDL/LRAA.wdl -i wf.ref-free.json  --verbose --no-cache



######
# cleaning

wdl_clean:
	rm -rf  ./*_LRAA_wf ./_LAST

clean: wdl_clean
	rm -rf ./__*  ./LRAA.gtf ./LRAA.bed ./LRAA.quant.expr ./LRAA.quant.tracking \
		./LRAA.serial.* ./LRAA.par.* serial.time.txt parallel.time.txt bench.summary.tsv \
		./LRAA.contigtmp ./LRAA.serial.contigtmp ./LRAA.par.contigtmp

bench_clean:
	$(MAKE) clean


