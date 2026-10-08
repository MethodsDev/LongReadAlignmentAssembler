# Site-first differential TSS / PolyA usage across single-cell clusters: an execution
# test on seven PBMC genes (see README.md for the walk-through and the method).
#
#   make test     run the whole pipeline on data/ and check the expected switches
#   make figures  also draw the showcase events' read-track figures (R tidyverse, cowplot)
#   make clean
#
# Needs python3 with pysam, numpy, scipy, pandas; R with DEXSeq, stageR, data.table,
# argparse, BiocParallel, Matrix (the sc image). Part of testing/single_cells' `make test`.

SHELL      := /bin/bash
.SHELLFLAGS := -o pipefail -c
LRAA_HOME  ?= $(abspath ../../..)
UTIL       := $(LRAA_HOME)/util/sc/site_usage
CPU        ?= 2

D          := data
P          := test
DEX        := $(P).dexseq
KINDS      := TSS PolyA
EVENTS     := $(foreach k,$(KINDS),$(DEX).$(k).events.tsv)

.PHONY: test test_wdls figures clean
.DELETE_ON_ERROR:
.SECONDARY:

test: read_tracks/manifest.tsv
	python3 check_results.py --prefix $(P)

# 1. site table: LRAA's sites, assigned to genes
$(P).sites.tsv: $(D)/sites.TSS.bed $(D)/sites.PolyA.bed $(UTIL)/prep_site_table.py
	python3 $(UTIL)/prep_site_table.py --TSS_bed $(D)/sites.TSS.bed --PolyA_bed $(D)/sites.PolyA.bed \
		--window_TSS 25 --gtf $(D)/models.gtf --gene_trans_map $(D)/gene_trans_map.tsv \
		--output $@ --gene_spans_output $(P).gene_spans.tsv

# 2. per-cell read support with LRAA's site-support rules
support/$(P).TSS-sparseM/matrix.mtx.gz: $(D)/sites.TSS.bed $(D)/sites.PolyA.bed $(D)/reads.bam \
		$(LRAA_HOME)/util/sc/site_read_support_to_sparse_matrix.py
	mkdir -p support
	python3 $(LRAA_HOME)/util/sc/site_read_support_to_sparse_matrix.py --bam $(D)/reads.bam \
		--TSS_bed $(D)/sites.TSS.bed --PolyA_bed $(D)/sites.PolyA.bed --cell_clusters $(D)/cell_clusters.tsv \
		--HiFi --CPU $(CPU) --output_prefix support/$(P)

# 3. into the layout the test reads
$(P).site_counts.mtx.gz: $(P).sites.tsv support/$(P).TSS-sparseM/matrix.mtx.gz $(UTIL)/site_counts_from_lraa_support.py
	python3 $(UTIL)/site_counts_from_lraa_support.py --sites $< --support_prefix support/$(P) \
		--cell_clusters $(D)/cell_clusters.tsv --output_prefix $(P)

# 4. DEXSeq on pseudo-replicates, per site kind
$(DEX).%.genes.tsv: $(P).site_counts.mtx.gz $(UTIL)/site_usage_dexseq.R
	Rscript $(UTIL)/site_usage_dexseq.R --counts_prefix $(P) --sites $(P).sites.tsv --kind $* \
		--cores $(CPU) --output_prefix $(DEX) 2>&1 | grep -v "converting counts to integer mode" | tee dexseq.$*.log
	test -s $(DEX).$*.pairwise.tsv.gz

# 5. pairwise switch events, annotated
$(DEX).%.events.tsv: $(DEX).%.genes.tsv $(UTIL)/annotate_site_usage_events.py
	python3 $(UTIL)/annotate_site_usage_events.py --dexseq_prefix $(DEX) --kind $* --sites $(P).sites.tsv \
		--gtf $(D)/models.gtf --cluster_quant_tar $(D)/cluster_quant.tar.gz --output $@

# 6. alternative terminal usage vs switches coupled with alternative splicing, from the reads
$(DEX).site_pairs.splicing.tsv: $(EVENTS) $(UTIL)/classify_site_pairs_by_splicing.py
	python3 $(UTIL)/classify_site_pairs_by_splicing.py --events $(foreach k,$(KINDS),$(k)=$(DEX).$(k).events.tsv) \
		--sites $(P).sites.tsv --gene_spans $(P).gene_spans.tsv --bam $(D)/reads.bam \
		--cell_clusters $(D)/cell_clusters.tsv --CPU $(CPU) --output $@

# 7. showcase events and their read-track data
read_tracks/events.tsv: $(DEX).site_pairs.splicing.tsv $(UTIL)/select_showcase_events.py
	mkdir -p read_tracks
	python3 $(UTIL)/select_showcase_events.py --dexseq_prefix $(DEX) --splicing $< \
		--cell_clusters $(D)/cell_clusters.tsv --min_cluster_cells 50 --min_gene_reads 20 \
		--sites $(P).sites.tsv --gtf $(D)/models.gtf --cluster_quant_tar $(D)/cluster_quant.tar.gz --output $@

read_tracks/manifest.tsv: read_tracks/events.tsv $(UTIL)/build_site_event_read_tracks.py
	python3 $(UTIL)/build_site_event_read_tracks.py --events $< --sites $(P).sites.tsv --gtf $(D)/models.gtf \
		--cluster_quant_tar $(D)/cluster_quant.tar.gz --tracking $(D)/tracking.tsv.gz --bam $(D)/reads.bam \
		--cell_clusters $(D)/cell_clusters.tsv --max_reads 30 --outdir read_tracks

# no WDL in this test; the target exists for testing/single_cells' loop over test dirs
test_wdls:
	@echo "sc_site_usage: no WDL to test"

figures: read_tracks/manifest.tsv
	Rscript figures.R

clean:
	rm -rf support read_tracks read_tracks.*.pdf $(P).* dexseq.*.log make.log
