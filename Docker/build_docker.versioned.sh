#!/bin/bash

set -ex

# Newer Docker daemons reject old client API pins inherited from the shell.
# Clear any stale override so the installed client can negotiate its default.
unset DOCKER_API_VERSION

VERSION=`cat VERSION.txt`

REGISTRY=us-central1-docker.pkg.dev/methods-dev-lab/lraa

CORE_IMAGE=${REGISTRY}/lraa-core:${VERSION}

#docker buildx create --name terra-builder --use

#docker buildx build --platform linux/amd64 \
#  -t us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:${VERSION} \
#  --push .


# lraa-core must be built and pushed first: Dockerfile.sc and Dockerfile.orf
# are FROM ${LRAA_CORE_IMAGE} and are built against the core image just made
# here, not against whatever copy the registry happens to hold.
#
#   image           dockerfile
#   lraa-core       Dockerfile.core
#   lraa-sc         Dockerfile.sc
#   lraa-orf        Dockerfile.orf
#   lraa-combined   Dockerfile       (everything in one image)
#
# Each image gets both ${VERSION} and :testing.  The :testing tag is applied
# with 'docker tag' on the image just built; rebuilding it would cost about an
# hour of R layers for a byte-identical result.

build_and_push() {
    local name=$1
    local dockerfile=$2

    docker build -f ${dockerfile} --build-arg LRAA_CORE_IMAGE=${CORE_IMAGE} \
        -t ${REGISTRY}/${name}:${VERSION} .
    docker tag ${REGISTRY}/${name}:${VERSION} ${REGISTRY}/${name}:testing

    docker push ${REGISTRY}/${name}:${VERSION}
    docker push ${REGISTRY}/${name}:testing
}

# 'lraa' is an alias for 'lraa-core', not a separate build.  Every pull off the
# plain name costs us egress, so the plain name has to be the smallest image
# that can run an assembly, not the one carrying R and TransDecoder.
alias_core() {
    local tag=$1

    docker tag ${CORE_IMAGE} ${REGISTRY}/lraa:${tag}
    docker push ${REGISTRY}/lraa:${tag}
}

build_and_push lraa-core     Dockerfile.core
build_and_push lraa-sc       Dockerfile.sc
build_and_push lraa-orf      Dockerfile.orf
build_and_push lraa-combined Dockerfile

alias_core ${VERSION}
alias_core testing

# verify
docker run --rm ${CORE_IMAGE} /usr/local/src/LRAA/LRAA --version
