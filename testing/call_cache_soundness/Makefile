
# Nothing to run outside of the wdl suite.
test:


include ../lraa_test_docker.mk

# MINIWDL_CFG, MINIWDL_RUN and MINIWDL_RUN_UNCACHED come from the include above.

# No check_test_image_revision prerequisite: this directory never runs an LRAA
# image, it builds its own throwaway one.
#
# test_reclaim_root_owned comes from the include; it lives beside the target it
# tests and runs here because this is where harness self-checks run.  Both are
# seconds, and testing/Makefile runs this directory first so that a broken
# harness is reported before forty minutes of runs that depend on it.
test_wdls: test_reclaim_root_owned test_cache_soundness

test_cache_soundness:
	bash ./test_cache_soundness.sh '$(MINIWDL_CFG)' '$(LRAA_TEST_DOCKER)'

# Driven by test_cache_soundness.sh, which overrides LRAA_TEST_DOCKER to its
# throwaway probe image and reads the outputs file.  Goes through $(MINIWDL_RUN)
# on purpose: the thing under test is what that variable puts in front of
# miniwdl, so a copy of it here would test the copy.
probe_run:
	$(MINIWDL_RUN) probe.wdl docker=$(LRAA_TEST_DOCKER) \
		-d probe_outdir -o probe.outputs.json --verbose


clean: reclaim_root_owned
	rm -rf ./probe_outdir ./probe.outputs.json ./probe.run*.log \
		./__miniwdl_call_cache ./__miniwdl_download_cache
	-docker image rm -f lraa-call-cache-probe:testing >/dev/null 2>&1
