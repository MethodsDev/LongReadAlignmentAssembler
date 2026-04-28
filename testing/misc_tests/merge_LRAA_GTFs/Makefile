PYTHON ?= python
MERGE_TOOL = ../../../util/merge_LRAA_GTFs.py
PREP_DIR = prepared_inputs
GENOME = $(PREP_DIR)/chr19.genome.fa
GTFS = $(PREP_DIR)/*.gtf

test: test_nohifi test_hifi

prep_inputs:
	$(PYTHON) ./prepare_merge_test_inputs.py --src data --out $(PREP_DIR)

test_nohifi: prep_inputs
	$(PYTHON) $(MERGE_TOOL) --genome $(GENOME) --output_gtf lraa.merged.nohifi.gtf --gtf $(GTFS)

test_hifi: prep_inputs
	$(PYTHON) $(MERGE_TOOL) --HiFi --genome $(GENOME) --output_gtf lraa.merged.hifi.gtf --gtf $(GTFS)

clean:
	rm -rf ./__* lraa.merged.gtf lraa.merged.nohifi.gtf lraa.merged.hifi.gtf \
		lraa.merged.nohifi.gtf.tracking.tsv lraa.merged.hifi.gtf.tracking.tsv $(PREP_DIR)

