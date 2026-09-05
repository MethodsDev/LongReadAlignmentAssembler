
test_quick:
	../../LRAA --bam pacbio.PBLR.small.bam --genome ref_genome.fa --cpu_budget 4 --output_prefix LRAA.small  --no_norm --tag_bam --HiFi --verbose --no_cleanup --no_chunk --no_stream_reads


test: test_quick test_reads_only test_include_ref_annot_guide test_quant test_chunked test_observability test_oversimplify


test_reads_only:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --cpu_budget 4 --output_prefix LRAA  --no_norm --HiFi --verbose --no_cleanup --no_chunk --no_stream_reads


test_include_ref_annot_guide:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --cpu_budget 4 --output_prefix LRAA.guided  --no_norm --gtf ref_annot.gtf --HiFi --verbose --no_cleanup --no_chunk --no_stream_reads


test_quant:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --cpu_budget 4 --output_prefix  LRAA.ref --gtf ref_annot.gtf --quant_only --tag_bam --HiFi --verbose --no_cleanup --no_chunk --no_stream_reads


# Chunked mode, on both routes into it.
#
# Nothing here exercised chunking until v0.21.0, and the cost of that was two
# releases in which every chunked run aborted: v0.20.0 gzipped quant tracking,
# the driver still looked for the uncompressed name, and so each chunk that
# succeeded was reported FAILED. A test that merely runs the thing would have
# caught it on the day.
#
# 3.4 Mb yields no cuts at the default spacing, so the cut size is forced down
# to make this a genuinely multi-chunk run rather than a one-chunk no-op.
test_chunked: test_chunked_entry_point test_chunked_route_parity test_chunked_strandless test_chunked_strandless_requires_N

test_chunked_entry_point:
	../../LRAA --chunk --quant_only --bam pacbio.PBLR.bam --genome ref_genome.fa --gtf ref_annot.gtf \
	    --cpu_budget 4 --HiFi --approx_MB_per_cut 0.5 --approx_MB_per_cut_wiggle_window 0.1 \
	    --no_stream_reads \
	    --output_prefix LRAA.chunked
	test -s LRAA.chunked.LRAA.quant-only.quant.expr
	test -s LRAA.chunked.LRAA.quant-only.quant.tracking.gz

# The two routes must not drift: LRAA --chunk and the parity harness call the
# same orchestration, so their merged output is compared rather than pinned to
# expected numbers -- a relationship survives a legitimate change to the corpus
# or the algorithm, and a hardcoded row count does not.
# --arm chunked skips stage 1 (strandless is the default since the
# chunk_by_strand flip), so ChunkedRun has no retained-record count to default
# -N to and requires it explicitly -- unlike `LRAA --chunk` above, which counts
# genome-mapped reads itself with the same -F 0x904 mask when -N is omitted.
# Passed explicitly here so both routes quantify against the same denominator;
# a mismatch would fail the diff below for a TPM reason, not a structural one.
test_chunked_route_parity: test_chunked_entry_point
	python3 ../../util/misc/run_chunked_quant_pipeline.py \
	    --bam pacbio.PBLR.bam --genome_fa ref_genome.fa --gtf ref_annot.gtf \
	    --num_total_reads $$(samtools view -c -F 0x904 pacbio.PBLR.bam) \
	    --output_dir __chunked_pipeline --arm chunked --cpu_budget 4 --HiFi \
	    --approx_MB_per_cut 0.5 --approx_MB_per_cut_wiggle_window 0.1 \
	    --no_stream_reads
	diff __chunked_pipeline/merged/chunked.quant.expr LRAA.chunked.LRAA.quant-only.quant.expr
	gunzip -c __chunked_pipeline/merged/chunked.quant.tracking.gz > __chunked_pipeline.tracking
	gunzip -c LRAA.chunked.LRAA.quant-only.quant.tracking.gz > __chunked_lraa.tracking
	diff __chunked_pipeline.tracking __chunked_lraa.tracking

# Strandless chunking, THROUGH LRAA.
#
# The mode shipped reachable only from util/misc/run_chunked_quant_pipeline.py:
# --strandless_chunks existed in ChunkedRun's own parser and LRAA never had the
# flag, so _run_chunked_mode could not pass it and no inputs JSON, CLI line or
# WDL could turn the mode on. Its parity gate ran the harness directly and so
# said nothing about that. This drives it through the entry point users have.
#
# Asserts the strandless SHAPE, not just a zero exit: cut selection is one
# strandless pass rather than one per orientation, stage 1 is skipped outright,
# and each interval carries the split it did internally. A --strandless_chunks
# that parsed and then quietly ran strand-first would satisfy neither.
#
# -N only sets the TPM denominator here and the assertions are structural, so it
# does not have to agree with any particular stage's count -- but it does have to
# be passed, which is the next test.
#
# Run from an UNINDEXED copy of the bam, deliberately. Only the strandless arm
# fetches the raw bam by region; strand-first reads the indexed stage-1 output, so
# an unindexed input was never a prerequisite of a chunked run, and every corpus
# here ships a committed .bai -- which is exactly why nobody noticed that the
# strandless arm needs one and that LRAA --chunk dispatches before LRAA's own
# index step. It surfaced as pysam's "fetch called on bamfile without index" the
# first time a caller had no .bai, which was the WDL layer. Testing against the
# indexed corpus copy would restore that blind spot.
test_chunked_strandless:
	cp pacbio.PBLR.bam __strandless_input.bam
	test ! -f __strandless_input.bam.bai
	../../LRAA --chunk --strandless_chunks --quant_only \
	    --num_total_reads $$(samtools view -c -F 0x904 __strandless_input.bam) \
	    --bam __strandless_input.bam --genome ref_genome.fa --gtf ref_annot.gtf \
	    --cpu_budget 4 --HiFi --approx_MB_per_cut 0.5 --approx_MB_per_cut_wiggle_window 0.1 \
	    --no_stream_reads \
	    --output_prefix __strandless
	test -f __strandless_input.bam.bai
	test -s __strandless.LRAA.quant-only.quant.expr
	python3 -c 'import json,collections; d=json.load(open("__strandless.chunked_work/timing.json")); assert "cuts_strandless" in d["stages"] and "cuts_plus" not in d["stages"], "cut selection was not strandless: %r" % (sorted(d["stages"]),); assert "skipped" in d["stages"]["strand_split"], "the whole-bam strand split still ran: %r" % (d["stages"]["strand_split"],); assert d["stage1_retained_records"] is None, "stage 1 still ran"; m=d["chunk_manifests"]; per=collections.Counter(c["region"].split(":")[0] for c in m); print("strandless intervals per contig: %r, quant units: %d" % (dict(per), d["arms"]["chunked"]["quant_units"])); assert max(per.values()) > 1, "no contig was cut: %r" % ([c["region"] for c in m],); assert all(c["split_counts"] for c in m), "intervals carry no per-chunk split accounting, so the split did not run inside them"'

# The refusal that goes with it -- but not through LRAA. LRAA's own wrapper
# auto-counts genome-mapped reads itself when -N is omitted (same -F 0x904
# helper, cached beside the bam), unconditionally of --strandless_chunks, so
# `LRAA --chunk --strandless_chunks` without -N now succeeds rather than
# refuses. The refusal above ChunkedRun's own parser doc still describes is
# real, just reachable only one level down: the --arm harness bypasses LRAA's
# auto-fill, and --arm chunked has no stage-1 count to default to either.
# Asserts the message says WHY, not merely that something is missing -- a bare
# argparse "required" would send the reader looking for a value the run cannot
# produce.
test_chunked_strandless_requires_N:
	! python3 ../../util/misc/run_chunked_quant_pipeline.py \
	    --bam pacbio.PBLR.bam --genome_fa ref_genome.fa --gtf ref_annot.gtf \
	    --output_dir __strandless_no_N --arm chunked --cpu_budget 4 --HiFi \
	    > __strandless_no_N.log 2>&1
	grep -q -- '-N is required with --strandless_chunks --arm chunked' __strandless_no_N.log
	grep -q 'does not run stage 1' __strandless_no_N.log

# Observability of a run that succeeds.
#
# Everything a worker logs goes to <tmp_prefix>.err.log inside the contigtmp
# tree, which is deleted on success -- so a successful run used to leave no
# per-contig record of anything, and streaming rescue's counters reached only
# `stats()`, which nothing outside the tests has ever called. A rescue rate of
# 6 reads in 17,897 was therefore indistinguishable from a rescue path that was
# silently disabled.
#
# Asserts the two summaries survive cleanup and carry real values: one rescue
# row per contig-strand, and a per-worker duration that is not zero. The
# duration is the specific regression -- the summary differenced a single
# timestamp against itself and reported every work unit shorter than the
# monitor's write interval as lasting 0.0 seconds.
test_observability:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --gtf ref_annot.gtf \
	    --quant_only --stream_reads \
	    --stream_reads_rescue_unassigned --cpu_budget 4 --HiFi \
	    --no_chunk \
	    --output_prefix LRAA.observability
	test -s LRAA.observability.LRAA.quant-only.rescue.summary.tsv
	test $$(grep -vc '^contig' LRAA.observability.LRAA.quant-only.rescue.summary.tsv) -ge 2
	awk -F'\t' 'NR>1 && $$3+0 > 0 {n++} END{exit !(n>=2)}' \
	    LRAA.observability.LRAA.quant-only.rescue.summary.tsv
	test -s LRAA.observability.LRAA.quant-only.resources.summary.tsv
	awk -F'\t' 'NR>1 && $$4+0 > 0 {n++} END{exit !(n>=1)}' \
	    LRAA.observability.LRAA.quant-only.resources.summary.tsv





test_oversimplify: test_oversimplify_quantifyonly test_oversimplify_denovo test_oversimplify_refguided


test_oversimplify_quantifyonly:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --cpu_budget 4 --output_prefix  LRAA.quantonly.oversimplify --gtf ref_annot.gtf --quant_only --HiFi --oversimplify minigenome --verbose --no_cleanup --no_chunk --no_stream_reads

test_oversimplify_denovo:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --cpu_budget 4 --output_prefix  LRAA.denovo.oversimplify --HiFi --oversimplify minigenome --verbose --no_cleanup --no_chunk --no_stream_reads

test_oversimplify_refguided:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --cpu_budget 4 --output_prefix  LRAA.refguided.oversimplify  --gtf ref_annot.gtf  --HiFi --oversimplify minigenome --verbose --no_cleanup --no_chunk --no_stream_reads




## WDLs

include ../lraa_test_docker.mk

# MINIWDL_CFG, MINIWDL_RUN and MINIWDL_RUN_UNCACHED come from the include above.

test_wdls: check_test_image_revision test_wdl_ID_ref_free test_wdl_ID_ref_guided test_wdl_Quant_only test_wdl_chunked test_wdl_chunked_strandless test_wdl_ORF_prediction

test_wdl_ID_ref_free:
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i wdl.ID_ref_free_Quant_mode.json $(LRAA_WDL_DOCKER) --verbose

test_wdl_ID_ref_guided:
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i wdl.ID_ref_guided_Quant_mode.json $(LRAA_WDL_DOCKER) --verbose

test_wdl_Quant_only:
	$(MINIWDL_RUN) ../../WDL/LRAA.wdl -i wdl.Quant_only_mode.json $(LRAA_WDL_DOCKER) --verbose

# Chunking, through the WDL layer.
#
# Chunking was unreachable from any WDL until these inputs existed: LRAA has always
# needed --chunk, no workflow passed it, and LRAA_runner.wdl enumerates every flag
# explicitly, so there was no generic passthrough to smuggle it through either.
# Reachability is the whole point of the wiring, so reachability is what this asserts,
# from the chunk manifests rather than from exit 0 -- a flag that plumbs through
# cleanly and never takes effect fails no test that only looks at the status.
#
# What it asserts is that some contig-strand was cut into more than one chunk, NOT
# that more than one chunk exists. Strand-first chunking emits at least one chunk per
# orientation no matter what, so an uncut run still reports two manifests and "more
# than one chunk" would pass on a run that placed no cuts at all. Grouping by the
# manifest's own region field is the difference between observing chunking and
# observing that LRAA has two strands.
#
# 0.5 Mb cuts over a 3.4 Mb contig, matching test_chunked_entry_point above: at the
# 10 Mb default this corpus is shorter than one chunk, so no cut is placed and chunked
# mode degenerates to the whole contig -- which is exactly the green-but-vacuous run
# the grouped assertion rejects.
#
# --dir __wdl_chunked/. pins the run directory instead of reading _LAST, so the
# assertion cannot end up inspecting whichever run finished most recently. UNCACHED
# because a cache hit skips the task, and a replayed chunk report would satisfy the
# assertion while proving nothing about this build.
test_wdl_chunked:
	rm -rf __wdl_chunked
	$(MINIWDL_RUN_UNCACHED) ../../WDL/LRAA.wdl -i wdl.Quant_only_chunked_mode.json $(LRAA_WDL_DOCKER) --dir __wdl_chunked/. --verbose
	python3 -c 'import json,collections; r=json.load(open("__wdl_chunked/outputs.json"))["LRAA_wf.chunkReports"]; assert len(r)==1, "expected one shard chunk report, got %r" % (r,); m=json.load(open(r[0]))["chunk_manifests"]; per=collections.Counter(c["region"].split(":")[0] for c in m); print("chunks per contig-strand: %r" % (dict(per),)); assert max(per.values()) > 1, "no contig-strand was cut, so chunking never happened: %r" % ([c["region"] for c in m],)'

# Strandless chunking, through the WDL layer.
#
# The mode was reachable from util/misc/run_chunked_quant_pipeline.py and from
# nothing else -- --strandless_chunks lived in ChunkedRun's parser alone, so LRAA
# could not pass it and no WDL could reach it. test_chunked_strandless above covers
# the LRAA half; this covers the WDL half, which is the half that was missing.
#
# No num_total_reads in the inputs, deliberately. The strandless arm runs no stage 1
# and so requires -N, and LRAA refuses without it (test_chunked_strandless_requires_N).
# The scattered path already counts the bam and passes the count to every shard, so
# the requirement is satisfied without the caller knowing it exists -- and if that
# ever stops being true this run fails on the refusal rather than quietly.
#
# Asserts the strandless SHAPE, not that a Boolean was accepted: cut selection is one
# strandless pass rather than one per orientation, stage 1 is skipped, some contig was
# cut, and each interval carries the split it did internally. A --strandless_chunks
# that plumbed through and then ran strand-first satisfies none of those.
test_wdl_chunked_strandless:
	rm -rf __wdl_chunked_strandless
	$(MINIWDL_RUN_UNCACHED) ../../WDL/LRAA.wdl -i wdl.Quant_only_chunked_strandless_mode.json $(LRAA_WDL_DOCKER) --dir __wdl_chunked_strandless/. --verbose
	python3 -c 'import json,collections; r=json.load(open("__wdl_chunked_strandless/outputs.json"))["LRAA_wf.chunkReports"]; assert len(r)==1, "expected one shard chunk report, got %r" % (r,); d=json.load(open(r[0])); assert "cuts_strandless" in d["stages"] and "cuts_plus" not in d["stages"], "cut selection was not strandless: %r" % (sorted(d["stages"]),); assert "skipped" in d["stages"]["strand_split"], "the whole-bam strand split still ran: %r" % (d["stages"]["strand_split"],); m=d["chunk_manifests"]; per=collections.Counter(c["region"].split(":")[0] for c in m); print("strandless intervals: %r, quant units: %d, -N from: %s" % (dict(per), d["arms"]["chunked"]["quant_units"], d["num_total_reads_source"])); assert max(per.values()) > 1, "no contig was cut: %r" % ([c["region"] for c in m],); assert all(c["split_counts"] for c in m), "intervals carry no per-chunk split accounting, so the split did not run inside them"'

test_wdl_ORF_prediction:
	$(MINIWDL_RUN) ../../WDL/LRAA-ORF-prediction.wdl -i wdl.ORF_prediction.json $(LRAA_WDL_DOCKER_ORF) --verbose



## Cleaning

clean_miniwdl: reclaim_root_owned
	rm -rf ./_LAST ./*_LRAA_wf ./__miniwdl_call_cache ./__miniwdl_download_cache


clean: clean_miniwdl
	rm -rf ./__* ./LRAA*
