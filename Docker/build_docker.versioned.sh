#!/bin/bash

set -ex

# Newer Docker daemons reject old client API pins inherited from the shell.
# Clear any stale override so the installed client can negotiate its default.
unset DOCKER_API_VERSION

VERSION=`cat VERSION.txt`
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
# Retried, because the archive endpoint lags a push: a commit confirmed present
# by `git ls-remote` still 404s here for twenty to thirty seconds.  Failing on
# the first attempt turns "you pushed a moment ago" into "you did not push",
# which is the opposite of what this check is for.
reachable=no
for attempt in 1 2 3 4 5 6; do
    if curl -sSf -o /dev/null -L https://github.com/${GITHUB_REPO}/archive/${LRAA_CO}.tar.gz; then
        reachable=yes
        break
    fi
    [ $attempt -lt 6 ] && sleep 10
done
if [ "${reachable}" != "yes" ]; then
    set +x
    echo "" >&2
    echo "commit ${LRAA_CO} is not fetchable from ${GITHUB_REPO} after 6 attempts" >&2
    echo "over ~50s.  The Dockerfiles fetch the LRAA checkout by SHA from GitHub, so" >&2
    echo "this commit must be pushed before it can be built.  Push it and re-run:" >&2
    echo "" >&2
    echo "    git push origin HEAD" >&2
    echo "" >&2
    exit 1
fi

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
# by the LRAA_CO build arg, so bumping a version rebuilds one small layer rather
# than recompiling Seurat.  Passing it as a build arg keeps one copy of the SHA
# instead of four Dockerfiles that can drift apart.
#
# Only ${VERSION} is written here.  :testing belongs to build_docker.testing.sh
# and means "built from the commit under test": a release build stamping it too
# would leave the WDL tests pointed at released code, which is the one thing
# they must not be.  A release is reachable as :latest and as its version.

docker build -f Dockerfile.base -t ${BASE_IMAGE} .

build_image() {
    local name=$1
    local dockerfile=$2

    docker build -f ${dockerfile} \
        --build-arg LRAA_BASE_IMAGE=${BASE_IMAGE} \
        --build-arg LRAA_VERSION=v${VERSION} \
        --build-arg LRAA_CO=${LRAA_CO} \
        -t ${REGISTRY}/${name}:${VERSION} .
}

# 'lraa' is an alias for 'lraa-core', not a separate build.  Every pull off the
# plain name costs us egress, so the plain name has to be the smallest image
# that can run an assembly, not the one carrying R and TransDecoder.
alias_core() {
    local tag=$1

    docker tag ${REGISTRY}/lraa-core:${VERSION} ${REGISTRY}/lraa:${tag}
    docker push ${REGISTRY}/lraa:${tag}
}

# Every image is built and stamped-checked before ANY of them is pushed.  Pushing
# inside the build loop publishes each as it completes, so a check failing on the
# last one leaves the registry holding a partial set: some tags on this commit,
# some on the previous, nothing recording which.  For release tags that is worse
# than for testing tags, because :latest is what every WDL default resolves to.
IMAGES="lraa-core lraa-sc lraa-orf lraa-combined"

build_image lraa-core     Dockerfile.core
build_image lraa-sc       Dockerfile.sc
build_image lraa-orf      Dockerfile.orf
build_image lraa-combined Dockerfile

# Each image separately: they share a base and a build arg, but nothing
# structurally forces the checkout to match, which is how four images drifted
# apart before.
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

# The version too, and also before the push: an image can build and run while
# carrying the wrong checkout, and the reported version is compiled into that
# checkout.  Printing it after publication tells you what you already shipped.
EXPECTED="LRAA VERSION: v${VERSION}"
REPORTED=`docker run --rm ${REGISTRY}/lraa-core:${VERSION} /usr/local/src/LRAA/LRAA --version`
if [ "${REPORTED}" != "${EXPECTED}" ]; then
    set +x
    echo "" >&2
    echo "lraa-core:${VERSION} reports '${REPORTED}', expected '${EXPECTED}'." >&2
    echo "Built from ${LRAA_CO}.  Nothing has been pushed." >&2
    exit 1
fi
echo "lraa-core:${VERSION} reports ${REPORTED}, built from ${LRAA_CO}"

for name in ${IMAGES}; do
    docker push ${REGISTRY}/${name}:${VERSION}
done

alias_core ${VERSION}

