# The images the WDL test targets run against.  Defined here and included by
# every Makefile under testing/ that invokes miniwdl, so the tag moves in one
# edit instead of once per directory.
#
# Why an override is needed at all: the WDLs default every image input to
# :latest, which is built from Docker/LRAA_CO.txt -- the last release.  A test
# run left on the defaults exercises released code, so it cannot fail on
# anything broken since the release, and it passes just as happily when the WDL
# and the code it drives have gone out of step.  :testing is built from the
# current commit by Docker/build_docker.testing.sh.
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
