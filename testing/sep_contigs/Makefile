
test: test_serial test_parallel test_chunked_validation_recipe


# Run serial across (contig,strand) using component-level multiprocessing only
test_serial:
	../../LRAA --bam pacbioccsminimap.minig.PBLR.bam --genome minig.fa --no_norm --no_parallelize_contigs --output_prefix LRAA.serial --HiFi --verbose --no_cleanup --no_chunk --no_stream_reads


# Run with default contig/strand-level parallelization (1 core per worker)
test_parallel:
	../../LRAA --bam pacbioccsminimap.minig.PBLR.bam --genome minig.fa --no_norm --output_prefix LRAA.par --HiFi --verbose --no_cleanup --no_chunk --no_stream_reads


# C1: `--chunk` with a cut geometry no contig can exceed is a VALIDATION RECIPE
# (see CodeBacklog.chunking-single-cell.md, "Per-contig chunking is a
# validation recipe, not a feature" -- decided, not open), not a production
# mode: it runs the full chunked machinery -- extraction, per-chunk
# normalization, quant, merge -- while placing zero cuts, so a correct result
# here proves that machinery against the unchunked reference path with
# nothing to attribute a difference to.
#
# 300 Mb exceeds every contig in this fixture (each is one gene locus, the
# largest under 1 Mb) by construction, per select_contig_cut_points.py:452 --
# not a property this corpus happens to have. Separately, and not what this
# geometry is chosen for: every gene locus here turns out to be annotated
# edge-to-edge, so no legal cut position exists at ANY geometry down to 10 kb
# (checked). That is a property of this corpus, not what makes the recipe
# work in general.
#
# Does NOT assert extractions_reusing_source_bam == chunks. This bam has
# alignments overhanging their nominal contig length (worst on KLHL11: 377 of
# 399 records), and extract_contig_region_inputs.py declines a requested
# reuse whenever any alignment does, forcing a real extraction that then
# drops the overhang instead. That drop does not touch quant.expr (verified
# below), but it does mean reuse is not universal on this corpus, so the
# invariant would fail here for a reason unrelated to correctness.
test_chunked_validation_recipe:
	rm -rf __unchunked_ref* __chunked_validation*
	../../LRAA --bam pacbioccsminimap.minig.PBLR.bam --genome minig.fa --gtf minig.gtf \
	    --quant_only --cpu_budget 4 --HiFi --no_chunk --no_stream_reads --output_prefix __unchunked_ref
	../../LRAA --chunk --strandless_chunks --approx_MB_per_cut 300 \
	    --bam pacbioccsminimap.minig.PBLR.bam --genome minig.fa --gtf minig.gtf \
	    --quant_only --cpu_budget 4 --HiFi --no_stream_reads --output_prefix __chunked_validation
	python3 -c 'import json; d = json.load(open("__chunked_validation.chunked_work/timing.json")); cp = d["cut_placement"]; assert cp["cuts_placed"] == 0, "cuts were placed: %r" % (cp,); assert cp["alignments_severed"] == 0, "alignments were severed: %r" % (cp,); offsets = set(c["offset"] for c in d["chunk_manifests"]); assert offsets == {0}, "non-zero chunk offset(s): %r" % (offsets,); print("structural invariants hold: %d chunks, 0 cuts, 0 severed, all offsets 0" % cp["chunks"])'
	python3 -c 'a = {(f[0], f[1]): f[2:] for f in (ln.rstrip("\n").split("\t") for ln in open("__unchunked_ref.LRAA.quant-only.quant.expr") if not ln.startswith("#")) if f[0] != "gene_id"}; b = {(f[0], f[1]): f[2:] for f in (ln.rstrip("\n").split("\t") for ln in open("__chunked_validation.LRAA.quant-only.quant.expr") if not ln.startswith("#")) if f[0] != "gene_id"}; only_a = set(a) - set(b); only_b = set(b) - set(a); diffs = [k for k in a if k in b and a[k] != b[k]]; assert not only_a and not only_b and not diffs, "rows differ: only_a=%r only_b=%r differing=%d" % (only_a, only_b, len(diffs)); print("quant.expr identical: %d rows, 0 unique keys, 0 differing values" % len(a))'


# Quick timing harness to compare serial vs parallelized contig processing
# Produces: bench.summary.tsv with real/user/sys times for both modes
bench: bench_clean
	( time -p ../../LRAA --bam pacbioccsminimap.minig.PBLR.bam --genome minig.fa --no_norm --cpu_budget 15 --no_parallelize_contigs --output_prefix LRAA.serial --HiFi --verbose --no_cleanup --no_chunk --no_stream_reads ) 2> serial.time.txt
	( time -p ../../LRAA --bam pacbioccsminimap.minig.PBLR.bam --genome minig.fa --no_norm --cpu_budget 15 --output_prefix LRAA.par --HiFi --verbose --no_cleanup --no_chunk --no_stream_reads ) 2> parallel.time.txt
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


test_wdls: check_test_image_revision test_wdl_quant_only test_wdl_ref_guided test_wdl_ref_free


##########
# wdls

include ../lraa_test_docker.mk

# MINIWDL_CFG, MINIWDL_RUN and MINIWDL_RUN_UNCACHED come from the include above.

test_wdls: test_wdl_quant_only test_wdl_ref_guided test_wdl_ref_free


test_wdl_quant_only:
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i wf.quant-only.json $(LRAA_WDL_DOCKER) --verbose

test_wdl_ref_guided:
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i wf.ref-guided.json $(LRAA_WDL_DOCKER) --verbose

test_wdl_ref_free:
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i wf.ref-free.json $(LRAA_WDL_DOCKER) --verbose



######
# cleaning

wdl_clean: reclaim_root_owned
	rm -rf  ./*_LRAA_wf ./_LAST ./__miniwdl_call_cache ./__miniwdl_download_cache

clean: wdl_clean
	rm -rf ./__*  ./LRAA.gtf ./LRAA.bed ./LRAA.quant.expr ./LRAA.quant.tracking.gz \
		./LRAA.serial.* ./LRAA.par.* serial.time.txt parallel.time.txt bench.summary.tsv \
		./LRAA.contigtmp ./LRAA.serial.contigtmp ./LRAA.par.contigtmp

bench_clean:
	$(MAKE) clean
