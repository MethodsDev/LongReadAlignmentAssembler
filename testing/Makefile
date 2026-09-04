DIRS = sirvs quant_eval single_contig strandless_parity sep_contigs ont_sep_contigs read_to_isoform_class_assignments single_cells txalign_path_rescue call_cache_soundness

QC = read_to_isoform_class_assignments

MISC = misc_tests

.PHONY: test test_pytest test_QC test_docker test_wdls test_wdls_apptainer test_DTU clean

test: test_pytest
	@for i in $(DIRS) $(MISC); do \
    echo "Running example in $$i..."; \
    (cd $$i; $(MAKE) test) || exit $$?; done

# RUN THE SUITE IN lraa-sc, NOT lraa-core.  lraa-core is deliberately lean and
# carries neither pytest nor pandas/scipy, so the suite cannot even start there,
# and the tests that drive util/sc/singlecell_tracking_to_sparse_matrix.py as a
# subprocess report a missing module rather than anything about the code.  The sc
# image has all three; $(LRAA_TEST_DOCKER_SC) in lraa_test_docker.mk names it.
#
# Two further conditions, both of which produce failures that look like code
# defects and are not:
#   - Run from a WRITABLE checkout as a user who can write the working directory.
#     Several tests write outputs relative to cwd, and a read-only or root-owned
#     tree turns those into PermissionError.
#   - `import WDL` resolves to this repo's own WDL/ DIRECTORY when the repo root is
#     on sys.path, shadowing miniwdl's package.  test_disksizegb_reaches_final_quant.py
#     loads WDL through the miniwdl API and needs miniwdl installed AND a working
#     directory that does not shadow it.  The sc image does not ship miniwdl, and its
#     absence is the ONLY thing between that image and a clean suite: sc as published,
#     non-root, on a writable checkout leaves exactly that file's tests failing, and
#     they pass once miniwdl is installed.  So `pip install miniwdl` in the container,
#     or run the suite where it is already available.
#
# python3 -m pytest, not bare pytest: several modules import through the package path
# (e.g. test_alignment_per_id_mixed_cigar.py's `from pylib.Util_funcs import ...`), which
# needs the repo root on sys.path. `python -m` puts cwd there; the bare console script
# does not, so bare pytest fails at COLLECTION and takes the whole suite with it.
test_pytest:
	cd .. && python3 -m pytest -q pylib


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
