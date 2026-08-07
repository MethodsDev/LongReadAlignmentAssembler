#!/bin/bash

set -ex

# Newer Docker daemons reject old client API pins inherited from the shell.
# Clear any stale override so the installed client can negotiate its default.
unset DOCKER_API_VERSION


VERSION=latest

REGISTRY=us-central1-docker.pkg.dev/methods-dev-lab/lraa

CORE_IMAGE=${REGISTRY}/lraa-core:${VERSION}


#docker buildx create --name terra-builder --use

#docker buildx build --platform linux/amd64 \
#  -t us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest \
#  --push .


# lraa-core must be built and pushed first: Dockerfile.sc and Dockerfile.orf
# are FROM ${LRAA_CORE_IMAGE} and are built against the core image just made
# here, not against whatever copy the registry happens to hold.
#
#   image        dockerfile
#   lraa-core    Dockerfile.core
#   lraa-sc      Dockerfile.sc
#   lraa-orf     Dockerfile.orf
#   lraa         Dockerfile        (combined image, kept for anyone pinning it)

build_and_push() {
    local name=$1
    local dockerfile=$2

    docker build -f ${dockerfile} --build-arg LRAA_CORE_IMAGE=${CORE_IMAGE} \
        -t ${REGISTRY}/${name}:${VERSION} .
    docker push ${REGISTRY}/${name}:${VERSION}
}

build_and_push lraa-core Dockerfile.core
build_and_push lraa-sc   Dockerfile.sc
build_and_push lraa-orf  Dockerfile.orf
build_and_push lraa      Dockerfile

# verify
docker run --rm ${CORE_IMAGE} /usr/local/src/LRAA/LRAA --version
