#!/bin/bash

set -ex

# Newer Docker daemons reject old client API pins inherited from the shell.
# Clear any stale override so the installed client can negotiate its default.
unset DOCKER_API_VERSION

VERSION=`cat VERSION.txt`
LRAA_CO=`cat LRAA_CO.txt`

REGISTRY=us-central1-docker.pkg.dev/methods-dev-lab/lraa

BASE_IMAGE=lraa-base:${VERSION}

#docker buildx create --name terra-builder --use

#docker buildx build --platform linux/amd64 \
#  -t us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:${VERSION} \
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
# in LRAA_CO.txt, so bumping a version rebuilds one small layer per image rather
# than recompiling Seurat.  Passing it as a build arg keeps one copy of the SHA
# instead of four Dockerfiles that can drift apart.
#
# Each image gets both ${VERSION} and :testing.  The :testing tag is applied
# with 'docker tag' on the image just built; rebuilding it would cost about an
# hour of R layers for a byte-identical result.

docker build -f Dockerfile.base -t ${BASE_IMAGE} .

build_and_push() {
    local name=$1
    local dockerfile=$2

    docker build -f ${dockerfile} \
        --build-arg LRAA_BASE_IMAGE=${BASE_IMAGE} \
        --build-arg LRAA_VERSION=v${VERSION} \
        --build-arg LRAA_CO=${LRAA_CO} \
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

    docker tag ${REGISTRY}/lraa-core:${VERSION} ${REGISTRY}/lraa:${tag}
    docker push ${REGISTRY}/lraa:${tag}
}

build_and_push lraa-core     Dockerfile.core
build_and_push lraa-sc       Dockerfile.sc
build_and_push lraa-orf      Dockerfile.orf
build_and_push lraa-combined Dockerfile

alias_core ${VERSION}
alias_core testing

# verify
docker run --rm ${REGISTRY}/lraa-core:${VERSION} /usr/local/src/LRAA/LRAA --version
