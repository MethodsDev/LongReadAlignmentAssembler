DIRS = sirvs quant_eval single_contig strandless_parity sep_contigs ont_sep_contigs read_to_isoform_class_assignments single_cells txalign_path_rescue call_cache_soundness

QC = read_to_isoform_class_assignments

MISC = misc_tests

.PHONY: test test_pytest test_QC test_docker test_wdls test_wdls_apptainer test_DTU clean

test: test_pytest
	@for i in $(DIRS) $(MISC); do \
    echo "Running example in $$i..."; \
    (cd $$i; $(MAKE) test) || exit $$?; done

test_pytest:
	cd .. && pytest -q pylib


all: test test_QC test_wdls test_DTU 

test_QC:
	@for i in $(QC); do \
    echo "Running example in $$i..."; \
    (cd $$i; $(MAKE) test) || exit $$?; done

test_docker:
	cd sirvs && ${MAKE} test_docker


# call_cache_soundness first, and deliberately.  It asserts that the runs below
# cannot be served stale outputs from the call cache after an image rebuild --
# the failure that made this suite report green in two minutes on 2026-08-16
# without starting a container.  A broken cache should be the first thing
# reported, not something discovered after forty minutes of runs that may all
# have been replays.
test_wdls:
	cd call_cache_soundness && $(MAKE) test_wdls
	cd sep_contigs && $(MAKE) test_wdls
	cd single_contig && $(MAKE) test_wdls
	cd single_cells && $(MAKE) test_wdls


# The whole wdl suite under Apptainer/Singularity rather than Docker.  One
# override, propagated by make into every sub-make below, so the recursion above
# is reused rather than duplicated.  Absolute, because each test directory sets
# MINIWDL_CFG relative to its own depth.
test_wdls_apptainer:
	$(MAKE) test_wdls MINIWDL_CFG=$(abspath $(CURDIR)/../miniwdl.apptainer.cfg)


test_DTU:
	cd DTU_eval && $(MAKE) clean test


clean:
	@for i in $(DIRS) $(QC) $(MISC); do \
    echo "Running example in $$i..."; \
    (cd $$i; $(MAKE) clean) || exit $$?; done
