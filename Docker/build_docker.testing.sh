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

# NEITHER tag above is stable enough to identify a benchmark run, and that is what
# this third one is for.  ${VERSION} moves on every build; ${VERSIONED_TAG} moves
# on every build that does not bump VERSION.txt, which during a release cycle is
# most of them -- 0.22.0-testing named three different commits in one afternoon.
# A result attributed to either is unreproducible the moment someone rebuilds.
#
# So: one tag per commit, never reused, because the commit is in the name.  A
# benchmarking run pins THIS, records it beside its outputs, and stays
# reattributable years later.  It is still -testing-lineage rather than a release
# alias: the bare version tags remain the exclusive property of the release
# scripts, for the reason given further down.
#
# Short SHA rather than full, because it is a docker tag a human retypes into an
# inputs JSON; the immutability comes from never moving it, and the full SHA is
# recorded in the image's own org.opencontainers.image.revision label, which the
# loop below asserts against this build.  Resolve a tag to its digest with
#   docker inspect --format '{{index .RepoDigests 0}}' <image>:<tag>
# if a run needs a reference that cannot be re-tagged even by mistake.
COMMIT_TAG=${LRAA_VERSION}-`git rev-parse --short=7 HEAD`

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
# Retried, because the archive endpoint lags a push: a commit confirmed present
# by `git ls-remote` still 404s here for twenty to thirty seconds.  Failing on
# the first attempt turns "you pushed a moment ago" into "you did not push",
# which is the opposite of what this check is for.
# The commit must be pushed before it is built.  Not because the build fetches it
# -- it no longer does -- but because an image labelled with a commit nobody else
# can obtain is not reproducible.  Checked locally: a commit reachable from any
# remote-tracking branch is on the remote.  The previous check GET'd
# https://github.com/<repo>/archive/<sha>.tar.gz to /dev/null, which downloaded
# 88MB to answer a yes/no question, and retried it six times.
git -C .. fetch --quiet origin || true
if [ -z "`git -C .. branch -r --contains ${LRAA_CO} 2>/dev/null`" ]; then
    set +x
    echo "" >&2
    echo "commit ${LRAA_CO} is not on any remote-tracking branch of ${GITHUB_REPO}." >&2
    echo "Images are labelled with this SHA, so it must be pushed for the image to be" >&2
    echo "reproducible by anyone else.  Push it and re-run:" >&2
    echo "" >&2
    echo "    git push origin HEAD" >&2
    echo "" >&2
    exit 1
fi

# The checkout is staged into the build context by git archive rather than fetched
# from GitHub.  git archive reads content-addressed objects, so the tarball is
# byte-determined by the SHA -- the same guarantee fetching by SHA gave -- without
# the network, and without shipping testing/.  That matters more than it sounds:
# the GitHub archive of this repo is 88MB because it CONTAINS testing/, and every
# Dockerfile downloaded all of it and then ran `rm -rf LRAA/testing` to delete
# ~84MB of it.  Four images plus the old reachability GET made it ~440MB per build
# to install 1.2MB of content, and it put every build at the mercy of an endpoint
# that rate-limits: two builds were lost to a 429 body being written into the
# tarball and surfacing as "gzip: stdin: not in gzip format".
#
# The sidecar .sha is not decoration.  A stale tarball left in the context by an
# interrupted build would otherwise be reused under a new LRAA_CO, producing an
# image whose label names a commit it does not contain.  The Dockerfiles refuse
# unless the two agree.
rm -f lraa_checkout.tar.gz lraa_checkout.sha
git -C .. archive --format=tar.gz --prefix=LRAA/ ${LRAA_CO} \
    `git -C .. ls-tree --name-only ${LRAA_CO} | grep -vx testing` \
    > lraa_checkout.tar.gz
echo ${LRAA_CO} > lraa_checkout.sha
echo "staged checkout `du -h lraa_checkout.tar.gz | cut -f1` for ${LRAA_CO} (testing/ excluded)"

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

# --cache-from, and the inline cache that makes it work.
#
# The docker driver supports INLINE cache only: a pushed image carries its own
# cache metadata when built with BUILDKIT_INLINE_CACHE=1, and a later build seeds
# from it by naming that image in --cache-from.  `type=registry` needs a
# docker-container driver, which this script does not use.
#
# Why base is now PUSHED as well: the R stack in Dockerfile.sc is 2.3 GB and
# roughly an hour of compiling Seurat, cached across commits only because the
# checkout is the last layer.  That holds while the local cache holds.  Once it is
# pruned, base is rebuilt, and `apt-get update` resolves whatever package versions
# are current that day -- so base comes out with a DIFFERENT digest, every layer
# in Dockerfile.sc that descends from it misses, and Seurat recompiles.  Naming a
# published base as both the parent cache source and a --cache-from input is what
# breaks that chain: the layers come off the registry instead of the compiler.
#
# The release scripts deliberately do NOT do this.  A release should resolve its
# own apt and CRAN content rather than inherit a cached layer from whenever a devel
# build last ran; freshness is worth an hour there, and releases are rare.
#
# --cache-from on an image that does not exist is a warning, not an error, so a
# first run against an empty registry still builds.
docker build -f Dockerfile.base \
    --build-arg BUILDKIT_INLINE_CACHE=1 \
    --cache-from ${REGISTRY}/lraa-base:${VERSION} \
    -t ${BASE_IMAGE} \
    -t ${REGISTRY}/lraa-base:${VERSION} \
    -t ${REGISTRY}/lraa-base:${VERSIONED_TAG} \
    -t ${REGISTRY}/lraa-base:${COMMIT_TAG} .

build_image() {
    local name=$1
    local dockerfile=$2

    docker build -f ${dockerfile} \
        --build-arg LRAA_BASE_IMAGE=${BASE_IMAGE} \
        --build-arg LRAA_VERSION=v${LRAA_VERSION} \
        --build-arg LRAA_CO=${LRAA_CO} \
        --build-arg BUILDKIT_INLINE_CACHE=1 \
        --cache-from ${REGISTRY}/${name}:${VERSION} \
        -t ${REGISTRY}/${name}:${VERSION} \
        -t ${REGISTRY}/${name}:${VERSIONED_TAG} \
        -t ${REGISTRY}/${name}:${COMMIT_TAG} .
}

# Every image is built and checked before ANY of them is pushed.  Pushing inside
# the build loop publishes each image as it completes, so a check that fails on
# the last one leaves the registry holding a partial set: some tags moved to this
# commit, some still on the previous, and nothing recording which.  The tags a
# test harness follows are exactly the ones that must not be half updated.
IMAGES="lraa-core lraa-sc lraa-orf lraa-combined"

build_image lraa-core     Dockerfile.core
build_image lraa-sc       Dockerfile.sc
build_image lraa-orf      Dockerfile.orf
build_image lraa-combined Dockerfile

# Every image must carry the commit this build was told to use.  Checking one and
# assuming the rest is what let four images drift apart before: they share a base
# and a build arg but nothing structurally forces the checkout to match, so the
# assertion has to name each of them.
for name in ${IMAGES}; do
    STAMPED=`docker inspect ${REGISTRY}/${name}:${VERSION} \
        --format '{{index .Config.Labels "org.opencontainers.image.revision"}}'`
    if [ "${STAMPED}" != "${LRAA_CO}" ]; then
        set +x
        echo "" >&2
        echo "${name}:${VERSION} is stamped ${STAMPED}, expected ${LRAA_CO}." >&2
        echo "Nothing has been pushed.  A stale cached layer is the usual cause." >&2
        exit 1
    fi
done

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

# Only now.  Every image exists locally, every one is stamped with this commit,
# and lraa-core reports the expected version.
for name in ${IMAGES}; do
    docker push ${REGISTRY}/${name}:${VERSION}
    docker push ${REGISTRY}/${name}:${VERSIONED_TAG}
    docker push ${REGISTRY}/${name}:${COMMIT_TAG}
done

# lraa-base, pushed for one reason: to be the --cache-from source for the next
# build on a machine whose local cache is gone.  It is NOT in ${IMAGES} because
# that list drives the revision-label assertion above and base carries no
# LRAA_CO -- it holds no checkout, which is exactly why it is reusable across
# commits.  Pushed after the assertions for the same reason as the rest: a base
# published while a later image fails its check would seed future builds from a
# commit that never passed.
for tag in ${VERSION} ${VERSIONED_TAG} ${COMMIT_TAG}; do
    docker push ${REGISTRY}/lraa-base:${tag}
done

# The pin a benchmark run should record, printed rather than left to be
# reconstructed: the tag names the commit, and the digest is what cannot be moved
# even deliberately.  BENCHMARKING.md's WDLs default to :latest, so a run that
# wants this code has to pass `docker` explicitly -- these are the strings to pass.
echo ""
echo "immutable pins for this build (commit ${LRAA_CO}):"
for name in ${IMAGES}; do
    DIGEST=`docker inspect ${REGISTRY}/${name}:${COMMIT_TAG} \
        --format '{{if .RepoDigests}}{{index .RepoDigests 0}}{{else}}<unpushed>{{end}}'`
    echo "  ${REGISTRY}/${name}:${COMMIT_TAG}"
    echo "      ${DIGEST}"
done
