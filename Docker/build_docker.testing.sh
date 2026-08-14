#!/bin/bash

set -ex

# Newer Docker daemons reject old client API pins inherited from the shell.
# Clear any stale override so the installed client can negotiate its default.
unset DOCKER_API_VERSION

VERSION=testing
LRAA_VERSION=`cat VERSION.txt`

# A versioned image built for testing carries the suffix so it cannot be mistaken
# for the published release of the same version.  Only the release scripts write
# a bare version tag.  Local test_wdls follow the moving :testing tag; a run on a
# VM should pin ${LRAA_VERSION}-testing instead, so the image under test is still
# identifiable after :testing has moved on.
VERSIONED_TAG=${LRAA_VERSION}-testing

# The checkout comes from the commit you are sitting on, which all three scripts
# now do.  They differ only in the tags they write.
LRAA_CO=`git rev-parse HEAD`

REGISTRY=us-central1-docker.pkg.dev/methods-dev-lab/lraa

# Must stay in step with the URL the Dockerfiles curl the checkout from; the
# reachability check below is worthless if it probes a different repository.
GITHUB_REPO=MethodsDev/LongReadAlignmentAssembler

BASE_IMAGE=lraa-base:${LRAA_VERSION}

# The Dockerfiles fetch the checkout by SHA from GitHub, so an unpushed commit
# does not build: the fetch 404s a few layers in, long after the expensive ones,
# and the error names a tarball rather than the mistake.  Worse, a SHA that
# happens to exist upstream but is not what you have locally builds something
# that is not your working tree.  Fail here instead.
if ! curl -sSf -o /dev/null -L https://github.com/${GITHUB_REPO}/archive/${LRAA_CO}.tar.gz; then
    set +x
    echo "" >&2
    echo "commit ${LRAA_CO} is not fetchable from ${GITHUB_REPO}." >&2
    echo "The Dockerfiles fetch the LRAA checkout by SHA from GitHub, so this" >&2
    echo "commit must be pushed before it can be built.  Push it and re-run:" >&2
    echo "" >&2
    echo "    git push origin HEAD" >&2
    echo "" >&2
    exit 1
fi

# lraa-base carries the dependencies the three published images share and is not
# pushed; its layers ship inside them.  It is built first because they are all
# FROM it.
#
#   image           dockerfile
#   lraa-base       Dockerfile.base   (build input, not published)
#   lraa-core       Dockerfile.core
#   lraa-sc         Dockerfile.sc
#   lraa-orf        Dockerfile.orf
#   lraa-combined   Dockerfile        (everything in one image)
#
# Every image takes the LRAA checkout as its final layer, so pointing at a
# different commit rebuilds one small layer per image rather than recompiling
# Seurat.  That is what makes a per-commit build cheap enough to be routine.
#
# Two tags are written, both unmistakably non-release.  :testing is the moving
# pointer test harnesses follow; ${LRAA_VERSION}-testing is version-bearing, for
# pinning a specific candidate or comparing two of them.  Both come out of the
# same build, so they cannot drift apart.
#
# The suffix is the point of the second one.  A bare version tag promises a
# published release, and an image built from an arbitrary commit must never be
# reachable under a name that makes that promise -- someone pulling 0.19.0 has
# no way to tell it was cut mid-development.  So a versioned image built for
# testing carries -testing, and :latest and the bare version tags stay the
# exclusive property of the release scripts.

docker build -f Dockerfile.base -t ${BASE_IMAGE} .

build_and_push() {
    local name=$1
    local dockerfile=$2

    docker build -f ${dockerfile} \
        --build-arg LRAA_BASE_IMAGE=${BASE_IMAGE} \
        --build-arg LRAA_VERSION=v${LRAA_VERSION} \
        --build-arg LRAA_CO=${LRAA_CO} \
        -t ${REGISTRY}/${name}:${VERSION} \
        -t ${REGISTRY}/${name}:${VERSIONED_TAG} .
    docker push ${REGISTRY}/${name}:${VERSION}
    docker push ${REGISTRY}/${name}:${VERSIONED_TAG}
}

build_and_push lraa-core     Dockerfile.core
build_and_push lraa-sc       Dockerfile.sc
build_and_push lraa-orf      Dockerfile.orf
build_and_push lraa-combined Dockerfile

# No 'lraa:testing' alias.  The plain name is what pipelines that have not named
# an image end up pulling, and its tags are release aliases of lraa-core; a
# pre-release commit does not belong under it.  Test workflows name lraa-core
# explicitly (see testing/lraa_test_docker.mk).

# verify: assert rather than print, because the failure this guards against is
# an image that built and runs but carries the wrong checkout -- a stale cached
# layer, or a commit that resolved to something other than what was intended.
# VERSION.txt is the host copy; the version reported is compiled into the
# checkout inside the image, so the two agreeing means the image really does
# hold the source this working tree is at.
EXPECTED="LRAA VERSION: v${LRAA_VERSION}"
REPORTED=`docker run --rm ${REGISTRY}/lraa-core:${VERSION} /usr/local/src/LRAA/LRAA --version`
if [ "${REPORTED}" != "${EXPECTED}" ]; then
    set +x
    echo "" >&2
    echo "lraa-core:${VERSION} reports '${REPORTED}', expected '${EXPECTED}'." >&2
    echo "Built from commit ${LRAA_CO}.  Either VERSION.txt disagrees with the" >&2
    echo "VERSION in the LRAA script at that commit, or the image did not pick" >&2
    echo "up the checkout it was told to." >&2
    exit 1
fi
echo "lraa-core:${VERSION} reports ${REPORTED}, built from ${LRAA_CO}"
