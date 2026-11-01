# Where isoform ends land: v0.48.0 5' end rule (own_vote + extension, SPACER revert), merge end
# logic, SE weak-junction guard, internal TSS of an ISM isoform. Five small BT474 minigenomes
# (README.md); each run takes seconds. check_end_placement.py asserts properties, not a golden GTF.

LRAA = ../../LRAA
MERGE = ../../util/merge_LRAA_GTFs.py
# The minigenomes are too small for the library cap calibration; this is BT474's calibration.
LRAA_ARGS = --HiFi --min_mapping_quality 0 --min_mapping_quality_for_final_quant 0 --cpu_budget 2 \
            --config_update ../data/cap_calib.json

LOCI = TIGD7 SPACER_revert SEPHS2 COMT

.PHONY: test clean $(addprefix check_,$(LOCI)) check_TIGD7_merge

test: $(addprefix check_,$(LOCI)) check_TIGD7_merge

run_%/done:
	rm -rf run_$*; mkdir -p run_$*
	cd run_$* && ../$(LRAA) --genome ../data/$*.fa --bam ../data/$*.bam --output_prefix $* $(LRAA_ARGS) > run.log 2>&1
	touch $@

$(addprefix check_,$(LOCI)): check_%: run_%/done
	./check_end_placement.py $* run_$*

run_TIGD7_merge/done:
	rm -rf run_TIGD7_merge; mkdir -p run_TIGD7_merge
	cd run_TIGD7_merge && ../$(MERGE) --genome ../data/TIGD7.fa \
	    --gtf ../data/TIGD7.merge_input.c1.gtf ../data/TIGD7.merge_input.c2.gtf \
	    --output_gtf merged.gtf > merge.log 2>&1
	touch $@

check_TIGD7_merge: run_TIGD7_merge/done
	./check_end_placement.py TIGD7_merge run_TIGD7_merge

clean:
	rm -rf run_*
