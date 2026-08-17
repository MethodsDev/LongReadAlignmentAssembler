# The images the WDL test targets run against.  Defined here and included by
# every Makefile under testing/ that invokes miniwdl, so the tag moves in one
# edit instead of once per directory.
#
# Why an override is needed at all: the WDLs default every image input to
# :latest, which the release script writes -- so it names whatever commit was
# last released, not the one in front of you.  A test run left on the defaults
# exercises released code, so it cannot fail on anything broken since, and it
# passes just as happily when the WDL and the code it drives have gone out of
# step.  :testing is built from the current commit by
# Docker/build_docker.testing.sh.
#
# Any of these can be overridden on the command line for a one-off image:
#
#     make test_docker LRAA_TEST_DOCKER=lraa-core:my-branch
#
# A locally built and tagged image works without a registry round trip; miniwdl
# only pulls when the tag is absent from the local daemon.

LRAA_TEST_REGISTRY ?= us-central1-docker.pkg.dev/methods-dev-lab/lraa
LRAA_TEST_TAG ?= testing

LRAA_TEST_DOCKER ?= $(LRAA_TEST_REGISTRY)/lraa-core:$(LRAA_TEST_TAG)
LRAA_TEST_DOCKER_SC ?= $(LRAA_TEST_REGISTRY)/lraa-sc:$(LRAA_TEST_TAG)
LRAA_TEST_DOCKER_ORF ?= $(LRAA_TEST_REGISTRY)/lraa-orf:$(LRAA_TEST_TAG)

# miniwdl input overrides, by the input names the workflows declare.  A workflow
# with more than one image input needs each one named: an unset one keeps its
# :latest default and half the run tests the last release.  Pass the set that
# matches the workflow being run --
#
#   LRAA.wdl, LRAA-merge_gtfs_from_sample_set.wdl,
#   subwdls/Incorporate_gene_symbols.wdl        $(LRAA_WDL_DOCKER)
#   LRAA-cell_cluster_guided.wdl,
#   LRAA-singlecell.wdl                         $(LRAA_WDL_DOCKER_CORE_SC)
#   LRAA-ORF-prediction.wdl                     $(LRAA_WDL_DOCKER_ORF)
#
# -- rather than pointing everything at core, which would run the single-cell
# tasks in an image without R and the ORF tasks in one without TransDecoder.
LRAA_WDL_DOCKER = docker=$(LRAA_TEST_DOCKER)
LRAA_WDL_DOCKER_SC = docker_sc=$(LRAA_TEST_DOCKER_SC)
LRAA_WDL_DOCKER_CORE_SC = $(LRAA_WDL_DOCKER) $(LRAA_WDL_DOCKER_SC)
LRAA_WDL_DOCKER_ORF = docker=$(LRAA_TEST_DOCKER_ORF)


# Repo root, derived from this file's own location rather than assumed, so it
# survives the directory being moved.  Every path built from it is absolute,
# because a test directory's depth varies -- testing/sep_contigs against
# testing/single_cells/sc_full_pipe -- and a relative path resolves differently
# in each.
LRAA_REPO_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))..)

# How every test directory invokes miniwdl.  Defined once here rather than
# copied into each Makefile, which is what it was: seven copies differing only
# in how many ../ they needed to reach the config, so seven places for one of
# them to drift.
#
# ?= so that `make test_wdls MINIWDL_CFG=...` and test_wdls_apptainer below
# still win: a command-line assignment beats a Makefile's `=` either way, but ?=
# also keeps an includer free to point one directory somewhere else.
MINIWDL_CFG ?= $(LRAA_REPO_ROOT)/miniwdl.test.cfg

# MINIWDL__CALL_CACHE__DIR overrides [call_cache] dir from the config file --
# miniwdl reads MINIWDL__<SECTION>__<KEY> ahead of any --cfg.  It names a cache
# directory namespaced by the CONTENT of the test images, because miniwdl's own
# cache key names the image only by the string it was given, and $(LRAA_TEST_TAG)
# is a moving tag.  Without it, `make test_wdls` after an image rebuild replays
# the previous build's task outputs, touches no container, and reports success
# in about two minutes.  Observed on 2026-08-16, green.
#
# testing/miniwdl_call_cache_dir.sh carries the full reasoning, including why
# this is not done by pinning the image by digest.
#
# Resolved inside the recipe, per invocation: no directory pays for it on `make
# clean` or a non-wdl target, and an image rebuilt between two targets of one
# suite run is still seen.  The `|| exit 1` is load-bearing -- the exit status of
# a command substitution used in a command's environment prefix is discarded, so
# a failure to resolve the namespace would otherwise reach miniwdl as an empty
# dir, which it resolves against [file_io] root.
LRAA_CALL_CACHE_DIR = $$('$(LRAA_REPO_ROOT)/testing/miniwdl_call_cache_dir.sh' $(LRAA_TEST_DOCKER) $(LRAA_TEST_DOCKER_SC) $(LRAA_TEST_DOCKER_ORF))

MINIWDL_RUN = _lraa_cc="$(LRAA_CALL_CACHE_DIR)" || exit 1; MINIWDL__CALL_CACHE__DIR="$$_lraa_cc" miniwdl run --cfg $(MINIWDL_CFG)

# For a target that asserts something about a task HAVING RUN: a chunk count, an
# output only a real execution produces, a log line.  A cache hit skips the task
# entirely, so such an assertion would be reading a replay and would pass with
# the task thoroughly broken.  Namespacing the cache does not help here -- it
# separates image builds, and this is a repeat run against one build, which is
# meant to hit.  Only call caching is disabled; downloads still cache.
#
# Not miniwdl's own --no-cache: that also turns off download-cache lookup, and
# it leaves call_cache put=true, so the run would still write entries -- into the
# un-namespaced default directory, where nothing will ever read them.  The two
# env overrides say exactly what is meant.
MINIWDL_RUN_UNCACHED = MINIWDL__CALL_CACHE__GET=false MINIWDL__CALL_CACHE__PUT=false miniwdl run --cfg $(MINIWDL_CFG)

# Container tasks run as root and bind-mount their work directory, so a miniwdl
# run directory routinely holds root-owned files -- one single-cell run here
# leaves upwards of a hundred.  rm -rf fails with Permission denied on the first
# root-owned directory it cannot descend into, so `make clean` aborts and the
# next suite run is blocked by the debris of the last one.  That has now cost
# two diagnoses.
#
# No sudo needed: a throwaway container is already root and can chown the tree
# back through the same bind mount.  Guarded on there being anything to reclaim,
# so the ordinary case costs one find and no container.
#
# Owner only, not owner:group.  These directories are group-shared and setgid, so
# the debris already carries the right group; chowning to `id -g` would replace
# it with the caller's private primary group and quietly take the tree away from
# everyone else in the project group.  rm only needs us to own the directory.
#
# Every wdl test directory lists this as the first prerequisite of its own clean
# target.  A new one should do the same rather than grow a second mechanism.
.PHONY: reclaim_root_owned
reclaim_root_owned:
	@foreign=`find . -xdev ! -user $$(id -u) -print -quit 2>/dev/null`; \
	if [ -z "$$foreign" ]; then exit 0; fi; \
	if ! command -v docker >/dev/null 2>&1; then \
	  echo "" >&2; \
	  echo "container-created files not owned by you are present under $$PWD" >&2; \
	  echo "  e.g. $$foreign" >&2; \
	  echo "and docker is not on PATH to reclaim them, so rm -rf may fail.  Under" >&2; \
	  echo "Apptainer this should not arise: --fakeroot maps container root back to" >&2; \
	  echo "the calling user outside the container." >&2; \
	  echo "" >&2; \
	  exit 0; \
	fi; \
	echo "reclaiming container-created files under $$PWD (e.g. $$foreign)"; \
	docker run --rm -v "$$PWD":/reclaim alpine \
	  chown -R "$$(id -u)" /reclaim \
	  || { echo "could not reclaim root-owned files under $$PWD; rm -rf will fail" >&2; exit 1; }

# The tag above is a pointer, and a miniwdl run never says which commit it
# resolved to.  A worktree that has moved since :testing was built therefore
# drives the current WDL against older code, and the resulting failure names
# whatever the drift happens to touch -- a script deleted two releases ago
# reads as `command not found`, which looks like a packaging bug rather than a
# stale checkout.  That cost a diagnosis once; the guard makes the skew itself
# the error message.
#
# Only lraa-core is inspected.  The build scripts build all four images from
# one checkout and verify each one's revision label before pushing any of them,
# so core's revision is the tag's revision.
#
# Escape hatch, for driving a hand-built or deliberately older image:
#
#     make test_wdls LRAA_TEST_SKIP_IMAGE_CHECK=1
#
.PHONY: check_test_image_revision
check_test_image_revision:
	@if [ -n "$(LRAA_TEST_SKIP_IMAGE_CHECK)" ]; then \
	  echo "image revision check skipped"; exit 0; \
	fi; \
	if ! command -v docker >/dev/null 2>&1; then \
	  echo "" >&2; \
	  echo "docker is not on PATH, so the image revision cannot be read and the" >&2; \
	  echo "worktree-vs-image check is being skipped.  This is expected under" >&2; \
	  echo "Apptainer, which can inspect labels only after converting an image to a" >&2; \
	  echo "SIF -- too expensive to do for a check.  The run continues, but nothing" >&2; \
	  echo "is verifying that $(LRAA_TEST_DOCKER) was built from this commit:" >&2; \
	  echo "  worktree at `git rev-parse HEAD`" >&2; \
	  echo "" >&2; \
	  exit 0; \
	fi; \
	docker image inspect $(LRAA_TEST_DOCKER) >/dev/null 2>&1 \
	  || docker pull $(LRAA_TEST_DOCKER) >/dev/null \
	  || { echo "cannot inspect or pull $(LRAA_TEST_DOCKER)" >&2; exit 1; }; \
	built=`docker inspect --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' $(LRAA_TEST_DOCKER)`; \
	head=`git rev-parse HEAD`; \
	if [ -z "$$built" ]; then \
	  echo "" >&2; \
	  echo "$(LRAA_TEST_DOCKER) carries no org.opencontainers.image.revision label," >&2; \
	  echo "so the commit it was built from cannot be established.  Images published" >&2; \
	  echo "before v0.19.1 predate the label -- :latest is one of them.  Rebuild the" >&2; \
	  echo "tag, or set LRAA_TEST_SKIP_IMAGE_CHECK=1 to run against it unverified." >&2; \
	  echo "" >&2; \
	  exit 1; \
	fi; \
	if [ "$$built" != "$$head" ]; then \
	  echo "" >&2; \
	  echo "$(LRAA_TEST_DOCKER) was built from $$built" >&2; \
	  echo "this worktree is at                    $$head" >&2; \
	  echo "" >&2; \
	  echo "The wdls come from the worktree and the code from the image, so a run" >&2; \
	  echo "now tests neither.  Rebuild the tag (bash Docker/build_docker.testing.sh)," >&2; \
	  echo "check out $$built, or set LRAA_TEST_SKIP_IMAGE_CHECK=1 to proceed anyway." >&2; \
	  echo "" >&2; \
	  exit 1; \
	fi; \
	dirty=`git status --porcelain --untracked-files=no | wc -l`; \
	if [ "$$dirty" -gt 0 ]; then \
	  echo "" >&2; \
	  echo "warning: $$dirty tracked file(s) modified since $$head." >&2; \
	  echo "The image was built from the commit, not the worktree, so none of those" >&2; \
	  echo "changes are under test -- including any fix you are trying to confirm." >&2; \
	  echo "" >&2; \
	fi; \
	echo "image and worktree agree at $$head"

# Each includer lists check_test_image_revision as the FIRST prerequisite of its
# own test_wdls.  Appending it from here instead would put it last -- make runs
# prerequisites in declaration order, so the check would report the skew after
# the runs it was meant to stop.  Declared here it would also look sufficient,
# which is worse than absent.
#
# .NOTPARALLEL because that ordering is only guaranteed for a serial make.
.NOTPARALLEL:


# Run the same wdl targets under Apptainer/Singularity instead of Docker.
#
#     make test_wdls_apptainer
#
# The only difference is which miniwdl config is used, and a command-line
# variable assignment beats a Makefile's own `=` AND propagates into sub-makes,
# so overriding MINIWDL_CFG once here reaches every directory the recursion
# touches.  Verified, because relying on unverified make semantics is how the
# prerequisite ordering above went wrong.
#
# The path is absolute, like MINIWDL_CFG above and for the same reason.
MINIWDL_APPTAINER_CFG ?= $(LRAA_REPO_ROOT)/miniwdl.apptainer.cfg

.PHONY: test_wdls_apptainer
test_wdls_apptainer:
	@test -f "$(MINIWDL_APPTAINER_CFG)" \
	  || { echo "no apptainer config at $(MINIWDL_APPTAINER_CFG)" >&2; exit 1; }
	$(MAKE) test_wdls MINIWDL_CFG=$(MINIWDL_APPTAINER_CFG)
