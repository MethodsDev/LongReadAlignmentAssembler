PYTHON ?= python
MERGE_TOOL = ../../../util/merge_LRAA_GTFs.py
DATA_DIR = data
GENOME = $(DATA_DIR)/chr19.genome.fa
GTFS = $(DATA_DIR)/*.gtf

test: test_nohifi test_hifi

test_nohifi:
	$(PYTHON) $(MERGE_TOOL) --genome $(GENOME) --output_gtf lraa.merged.nohifi.gtf --gtf $(GTFS)

test_hifi:
	$(PYTHON) $(MERGE_TOOL) --HiFi --genome $(GENOME) --output_gtf lraa.merged.hifi.gtf --gtf $(GTFS)
	$(MAKE) check_hifi_annotations

check_hifi_annotations:
	@awk 'BEGIN{FS=OFS="\t"; tx=0;annotated=0;tss_true=0;polya_true=0} \
		$$3=="transcript"{ \
			tx++; \
			if (index($$9, "TSS \"")>0 || index($$9, "PolyA \"")>0) annotated++; \
			if (index($$9, "TSS \"True\"")>0) tss_true++; \
			if (index($$9, "PolyA \"True\"")>0) polya_true++; \
		} \
		END{ \
			if (tx==0) {print "ERROR: No transcript entries found in lraa.merged.hifi.gtf"; exit 1} \
			if (annotated==0) {print "ERROR: No TSS/PolyA annotations found in lraa.merged.hifi.gtf"; exit 1} \
			if (tss_true==0 || polya_true==0) {print "ERROR: No TSS/PolyA True annotations found in lraa.merged.hifi.gtf"; exit 1} \
			print "PASS: HiFi merged GTF retains TSS/PolyA attrs (transcripts=" tx ", annotated=" annotated ", TSS_true=" tss_true ", PolyA_true=" polya_true ")"; \
		}' lraa.merged.hifi.gtf
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
			if (rows==0) {print "ERROR: No rows found in lraa.merged.hifi.gtf.tracking.tsv"; exit 1} \
			if (bad>0) {print "ERROR: Non-binary source_has_TSS/source_has_PolyA values in " bad " tracking rows"; exit 1} \
			if (tss_ones==0 || polya_ones==0) {print "ERROR: No source rows with source_has_TSS/source_has_PolyA == 1 in tracking file"; exit 1} \
			print "PASS: HiFi tracking retains source TSS/PolyA status (rows=" rows ", source_has_TSS==1 rows=" tss_ones ", source_has_PolyA==1 rows=" polya_ones ")"; \
		}' lraa.merged.hifi.gtf.tracking.tsv

clean:
	rm -rf ./__* lraa.merged.gtf lraa.merged.nohifi.gtf lraa.merged.hifi.gtf \
		lraa.merged.nohifi.gtf.tracking.tsv lraa.merged.hifi.gtf.tracking.tsv
