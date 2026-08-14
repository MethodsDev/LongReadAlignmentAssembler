#!/bin/bash

set -ex

# Newer Docker daemons reject old client API pins inherited from the shell.
# Clear any stale override so the installed client can negotiate its default.
unset DOCKER_API_VERSION


VERSION=latest
LRAA_VERSION=`cat VERSION.txt`
# Must stay in step with the URL the Dockerfiles curl the checkout from; the
# reachability check below is worthless if it probes a different repository.
GITHUB_REPO=MethodsDev/LongReadAlignmentAssembler

# The checkout is the commit you are sitting on.  There is no pinned SHA file:
# one existed, and it drifted -- it named a 0.18.3-era commit while VERSION.txt
# said 0.19.0, so this script would have published 0.19.0 images containing code
# from before the release.  A build describes the tree it was run from, and the
# SHA is stamped into every image as ENV LRAA_CO and as an OCI revision label.
LRAA_CO=`git rev-parse HEAD`

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

REGISTRY=us-central1-docker.pkg.dev/methods-dev-lab/lraa

BASE_IMAGE=lraa-base:${LRAA_VERSION}


#docker buildx create --name terra-builder --use

#docker buildx build --platform linux/amd64 \
#  -t us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:latest \
#  --push .


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
# Every image takes the LRAA checkout as its final layer, from the commit named
# by the LRAA_CO build arg, so bumping a version rebuilds one small layer rather
# than recompiling Seurat.

docker build -f Dockerfile.base -t ${BASE_IMAGE} .

build_and_push() {
    local name=$1
    local dockerfile=$2

    docker build -f ${dockerfile} \
        --build-arg LRAA_BASE_IMAGE=${BASE_IMAGE} \
        --build-arg LRAA_VERSION=v${LRAA_VERSION} \
        --build-arg LRAA_CO=${LRAA_CO} \
        -t ${REGISTRY}/${name}:${VERSION} .
    docker push ${REGISTRY}/${name}:${VERSION}
}

# 'lraa' is an alias for 'lraa-core', not a separate build.  Every pull off the
# plain name costs us egress, so the plain name has to be the smallest image
# that can run an assembly, not the one carrying R and TransDecoder.
alias_core() {
    local tag=$1

    docker tag ${REGISTRY}/lraa-core:${VERSION} ${REGISTRY}/lraa:${tag}
    docker push ${REGISTRY}/lraa:${tag}
}

build_and_push lraa-core     Dockerfile.core
build_and_push lraa-sc       Dockerfile.sc
build_and_push lraa-orf      Dockerfile.orf
build_and_push lraa-combined Dockerfile

alias_core ${VERSION}

# verify
docker run --rm ${REGISTRY}/lraa-core:${VERSION} /usr/local/src/LRAA/LRAA --version
