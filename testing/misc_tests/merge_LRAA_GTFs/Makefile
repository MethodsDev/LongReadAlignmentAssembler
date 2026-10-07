PYTHON ?= python
MERGE_TOOL = ../../../util/merge_LRAA_GTFs.py
DATA_DIR = data
GENOME = $(DATA_DIR)/chr19.genome.fa
GTFS = $(DATA_DIR)/*.gtf

# The merge respects input TSS/PolyA by default (the former --HiFi merge, v0.42.0);
# --ignore_TSS_POLYA restores the former non-HiFi merge that ignores them.
test: test_ignore_TSS_POLYA test_respect_TSS_POLYA test_parallel_parity

test_ignore_TSS_POLYA:
	$(PYTHON) $(MERGE_TOOL) --ignore_TSS_POLYA --genome $(GENOME) --output_gtf lraa.merged.ignore_TSS_POLYA.gtf --gtf $(GTFS)

test_respect_TSS_POLYA:
	$(PYTHON) $(MERGE_TOOL) --genome $(GENOME) --output_gtf lraa.merged.respect_TSS_POLYA.gtf --gtf $(GTFS)
	$(MAKE) check_TSS_PolyA_annotations

# --CPU N merges contig/strand units in forked workers; the output must match the
# serial (--CPU 1) run byte for byte, gtf and tracking tsv alike. Serial is run twice
# as well: the tracking rows were once ordered by string-hash set iteration, which
# changes with every process's hash seed, so a serial-vs-parallel match alone could
# be luck. Last, chr19+ merged on its own (--contig) must equal the chr19+ records of
# the full run, where it is merged AFTER chr19-: a unit's output may not depend on
# which units its process merged before it, or the workers could not reproduce it.
PARITY = lraa.merged.parity
test_parallel_parity:
	$(PYTHON) $(MERGE_TOOL) --CPU 1 --genome $(GENOME) --output_gtf $(PARITY).serial.gtf --gtf $(GTFS)
	$(PYTHON) $(MERGE_TOOL) --CPU 1 --genome $(GENOME) --output_gtf $(PARITY).serial2.gtf --gtf $(GTFS)
	$(PYTHON) $(MERGE_TOOL) --CPU 4 --genome $(GENOME) --output_gtf $(PARITY).cpu4.gtf --gtf $(GTFS)
	$(PYTHON) $(MERGE_TOOL) --CPU auto --genome $(GENOME) --output_gtf $(PARITY).auto.gtf --gtf $(GTFS)
	$(PYTHON) $(MERGE_TOOL) --CPU 1 --contig chr19+ --genome $(GENOME) --output_gtf $(PARITY).chr19plus.gtf --gtf $(GTFS)
	cmp $(PARITY).serial.gtf $(PARITY).serial2.gtf
	cmp $(PARITY).serial.gtf.tracking.tsv $(PARITY).serial2.gtf.tracking.tsv
	cmp $(PARITY).serial.gtf $(PARITY).cpu4.gtf
	cmp $(PARITY).serial.gtf.tracking.tsv $(PARITY).cpu4.gtf.tracking.tsv
	cmp $(PARITY).serial.gtf $(PARITY).auto.gtf
	cmp $(PARITY).serial.gtf.tracking.tsv $(PARITY).auto.gtf.tracking.tsv
	awk -F'\t' 'NF >= 9 && $$7 == "+"' $(PARITY).serial.gtf > $(PARITY).serial.chr19plus_records
	awk -F'\t' 'NF >= 9' $(PARITY).chr19plus.gtf > $(PARITY).chr19plus_records
	test -s $(PARITY).chr19plus_records
	cmp $(PARITY).serial.chr19plus_records $(PARITY).chr19plus_records
	awk -F'\t' 'NR == 1 || $$4 == "+"' $(PARITY).serial.gtf.tracking.tsv > $(PARITY).serial.chr19plus_tracking
	cmp $(PARITY).serial.chr19plus_tracking $(PARITY).chr19plus.gtf.tracking.tsv
	@echo "PASS: --CPU 4 and --CPU auto output identical to serial (gtf + tracking); serial reproducible; chr19+ alone == chr19+ in full run"

check_TSS_PolyA_annotations:
	@awk 'BEGIN{FS=OFS="\t"; tx=0;annotated=0;tss_true=0;polya_true=0} \
		$$3=="transcript"{ \
			tx++; \
			if (index($$9, "TSS \"")>0 || index($$9, "PolyA \"")>0) annotated++; \
			if (index($$9, "TSS \"True\"")>0) tss_true++; \
			if (index($$9, "PolyA \"True\"")>0) polya_true++; \
		} \
		END{ \
			if (tx==0) {print "ERROR: No transcript entries found in lraa.merged.respect_TSS_POLYA.gtf"; exit 1} \
			if (annotated==0) {print "ERROR: No TSS/PolyA annotations found in lraa.merged.respect_TSS_POLYA.gtf"; exit 1} \
			if (tss_true==0 || polya_true==0) {print "ERROR: No TSS/PolyA True annotations found in lraa.merged.respect_TSS_POLYA.gtf"; exit 1} \
			print "PASS: merged GTF retains TSS/PolyA attrs (transcripts=" tx ", annotated=" annotated ", TSS_true=" tss_true ", PolyA_true=" polya_true ")"; \
		}' lraa.merged.respect_TSS_POLYA.gtf
	@awk 'BEGIN{FS=OFS="\t"; tss_col=0; polya_col=0; bad=0; rows=0; tss_ones=0; polya_ones=0} \
		NR==1{ \
			for (i=1; i<=NF; i++) { \
				if ($$i=="source_has_TSS") tss_col=i; \
				if ($$i=="source_has_PolyA") polya_col=i; \
			} \
			if (!tss_col || !polya_col) {print "ERROR: tracking header missing source_has_TSS/source_has_PolyA columns"; exit 1} \
			next; \
		} \
		{ \
			rows++; \
			if (!(($$tss_col=="0" || $$tss_col=="1") && ($$polya_col=="0" || $$polya_col=="1"))) bad++; \
			if ($$tss_col=="1") tss_ones++; \
			if ($$polya_col=="1") polya_ones++; \
		} \
		END{ \
			if (rows==0) {print "ERROR: No rows found in lraa.merged.respect_TSS_POLYA.gtf.tracking.tsv"; exit 1} \
			if (bad>0) {print "ERROR: Non-binary source_has_TSS/source_has_PolyA values in " bad " tracking rows"; exit 1} \
			if (tss_ones==0 || polya_ones==0) {print "ERROR: No source rows with source_has_TSS/source_has_PolyA == 1 in tracking file"; exit 1} \
			print "PASS: merge tracking retains source TSS/PolyA status (rows=" rows ", source_has_TSS==1 rows=" tss_ones ", source_has_PolyA==1 rows=" polya_ones ")"; \
		}' lraa.merged.respect_TSS_POLYA.gtf.tracking.tsv

clean:
	rm -rf ./__* lraa.merged.gtf lraa.merged.ignore_TSS_POLYA.gtf lraa.merged.respect_TSS_POLYA.gtf \
		lraa.merged.ignore_TSS_POLYA.gtf.tracking.tsv lraa.merged.respect_TSS_POLYA.gtf.tracking.tsv \
		$(PARITY).*
