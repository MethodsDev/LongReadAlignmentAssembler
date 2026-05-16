
test: test_serial test_parallel


# Run serial across (contig,strand) using component-level multiprocessing only
test_serial:
	../../LRAA --bam pacbioccsminimap.minig.PBLR.bam --genome minig.fa --no_norm --no_parallelize_contigs --output_prefix LRAA.serial --HiFi --verbose --no_cleanup


# Run with default contig/strand-level parallelization (1 core per worker)
test_parallel:
	../../LRAA --bam pacbioccsminimap.minig.PBLR.bam --genome minig.fa --no_norm --output_prefix LRAA.par --HiFi --verbose --no_cleanup


# Quick timing harness to compare serial vs parallelized contig processing
# Produces: bench.summary.tsv with real/user/sys times for both modes
bench: bench_clean
	( time -p ../../LRAA --bam pacbioccsminimap.minig.PBLR.bam --genome minig.fa --no_norm --num_threads_per_worker 15 --no_parallelize_contigs --output_prefix LRAA.serial --HiFi --verbose --no_cleanup ) 2> serial.time.txt
	( time -p ../../LRAA --bam pacbioccsminimap.minig.PBLR.bam --genome minig.fa --no_norm --num_threads_per_worker 15 --output_prefix LRAA.par --HiFi --verbose --no_cleanup ) 2> parallel.time.txt
	@echo -e "mode\tCPU\treal\tuser\tsys" > bench.summary.tsv
	@SER_REAL=`grep '^real' serial.time.txt | awk '{print $$2}'`; \
	 SER_USER=`grep '^user' serial.time.txt | awk '{print $$2}'`; \
	 SER_SYS=`grep '^sys'  serial.time.txt | awk '{print $$2}'`; \
	 echo -e "serial\t15\t$${SER_REAL}\t$${SER_USER}\t$${SER_SYS}" >> bench.summary.tsv
	@PAR_REAL=`grep '^real' parallel.time.txt | awk '{print $$2}'`; \
	 PAR_USER=`grep '^user' parallel.time.txt | awk '{print $$2}'`; \
	 PAR_SYS=`grep  '^sys' parallel.time.txt | awk '{print $$2}'`; \
	 echo -e "parallel\t15\t$${PAR_REAL}\t$${PAR_USER}\t$${PAR_SYS}" >> bench.summary.tsv
	@echo "Wrote bench.summary.tsv" && cat bench.summary.tsv


test_wdls: test_wdl_quant_only test_wdl_ref_guided test_wdl_ref_free


##########
# wdls

MINIWDL_CFG = ../../miniwdl.test.cfg
MINIWDL_RUN = miniwdl run --cfg $(MINIWDL_CFG)

test_wdls: test_wdl_quant_only test_wdl_ref_guided test_wdl_ref_free


test_wdl_quant_only:
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i wf.quant-only.json  --verbose

test_wdl_ref_guided:
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i wf.ref-guided.json  --verbose

test_wdl_ref_free:
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i wf.ref-free.json  --verbose



######
# cleaning

wdl_clean:
	rm -rf  ./*_LRAA_wf ./_LAST ./__miniwdl_call_cache ./__miniwdl_download_cache

clean: wdl_clean
	rm -rf ./__*  ./LRAA.gtf ./LRAA.bed ./LRAA.quant.expr ./LRAA.quant.tracking \
		./LRAA.serial.* ./LRAA.par.* serial.time.txt parallel.time.txt bench.summary.tsv \
		./LRAA.contigtmp ./LRAA.serial.contigtmp ./LRAA.par.contigtmp

bench_clean:
	$(MAKE) clean
