#!/bin/bash

# Builds the CAS image and pushes two tags:
#   :<cellarium-cas version>-<git shortsha>   pinned, reproducible reference
#   :latest                                   what the WDL uses by default
# The CAS plugin image is exempt from the repo rule that reserves :latest for
# official LRAA releases; see CLAUDE.md.

set -ex

cd "$(dirname "$0")"

IMAGE=us-central1-docker.pkg.dev/methods-dev-lab/lraa/cas

if [ -n "$(git status --porcelain -- .)" ]; then
    echo "Error: uncommitted changes under Plugins/CAS; commit first so the tag's sha matches the image." >&2
    exit 1
fi

SHORTSHA=$(git rev-parse --short HEAD)

docker build -t ${IMAGE}:build-${SHORTSHA} .

CAS_VERSION=$(docker run --rm --entrypoint python ${IMAGE}:build-${SHORTSHA} \
    -c "from importlib.metadata import version; print(version('cellarium-cas'))")

VERSION_TAG=${CAS_VERSION}-${SHORTSHA}

docker tag ${IMAGE}:build-${SHORTSHA} ${IMAGE}:${VERSION_TAG}
docker tag ${IMAGE}:build-${SHORTSHA} ${IMAGE}:latest
docker rmi ${IMAGE}:build-${SHORTSHA}

docker push ${IMAGE}:${VERSION_TAG}
docker push ${IMAGE}:latest

echo "Pushed ${IMAGE}:${VERSION_TAG} and ${IMAGE}:latest"
