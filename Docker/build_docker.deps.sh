#!/bin/bash

set -ex

# The dependency images: everything that is NOT the LRAA checkout.
#
#   lraa-base      Dockerfile.base      samtools, htslib, minimap2, gffcompare, python
#   lraa-sc-base   Dockerfile.sc-base   FROM lraa-base, plus the R stack and scientific python
#
# Run this when the packages in those two files change.  DO NOT run it to cut a
# release: the release scripts pull these tags and build FROM them, which is the
# whole point of the split.  Nothing here depends on the commit being built, so
# a release that rebuilt them would be recompiling Seurat to ship an unchanged
# 2.3 GB of R.
#
# That is not hypothetical.  Before the split, the R layers lived in
# Dockerfile.sc above the checkout, so they were reused only while the LOCAL
# build cache held them -- and every release rebuilt lraa-base, which gave it a
# new image id and invalidated every layer below FROM in Dockerfile.sc.
# --cache-from could not save it either: BUILDKIT_INLINE_CACHE records only the
# final stage of a build, and Dockerfile.base is two-stage, so its `builder`
# stage always re-ran.  MEASURED on the v0.34.0 testing build: 3583 s in
# Dockerfile.sc, on a machine that had compiled the same packages an hour before.
#
# Tags: :latest is what the release scripts pull by default, and a dated tag is
# written beside it so a release can be pinned to the dependency set it was
# built against.  There is no version tag, deliberately -- these images do not
# track the LRAA version and naming them after one would imply they do.

# Newer Docker daemons reject old client API pins inherited from the shell.
unset DOCKER_API_VERSION

REGISTRY=us-central1-docker.pkg.dev/methods-dev-lab/lraa
DATE_TAG=`date +%Y%m%d`

cd "$(dirname "$0")"

# lraa-base first: lraa-sc-base is FROM it, and the local tag below is what that
# build resolves.  --cache-from seeds from the published image so an unchanged
# Dockerfile.base is layers off the registry rather than a fresh apt run.
docker build -f Dockerfile.base \
    --build-arg BUILDKIT_INLINE_CACHE=1 \
    --cache-from ${REGISTRY}/lraa-base:latest \
    -t lraa-base:latest \
    -t ${REGISTRY}/lraa-base:latest \
    -t ${REGISTRY}/lraa-base:${DATE_TAG} .

docker build -f Dockerfile.sc-base \
    --build-arg LRAA_BASE_IMAGE=lraa-base:latest \
    --build-arg BUILDKIT_INLINE_CACHE=1 \
    --cache-from ${REGISTRY}/lraa-sc-base:latest \
    -t lraa-sc-base:latest \
    -t ${REGISTRY}/lraa-sc-base:latest \
    -t ${REGISTRY}/lraa-sc-base:${DATE_TAG} .

# Both, or neither.  A release pulls both by tag, so publishing one of a pair
# that were built together leaves the next release building lraa-sc FROM an
# sc-base whose parent it never saw.
for name in lraa-base lraa-sc-base; do
    docker push ${REGISTRY}/${name}:latest
    docker push ${REGISTRY}/${name}:${DATE_TAG}
done

set +x
echo ""
echo "dependency images published:"
for name in lraa-base lraa-sc-base; do
    echo "  ${REGISTRY}/${name}:latest"
    echo "  ${REGISTRY}/${name}:${DATE_TAG}"
    echo "      `docker inspect --format '{{index .RepoDigests 0}}' ${REGISTRY}/${name}:latest`"
done
echo ""
echo "Releases pull these; they are not rebuilt by build_docker.{testing,latest,versioned}.sh."
echo ""
