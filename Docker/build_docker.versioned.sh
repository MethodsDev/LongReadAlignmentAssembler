#!/bin/bash

set -ex

# Newer Docker daemons reject old client API pins inherited from the shell.
# Clear any stale override so the installed client can negotiate its default.
unset DOCKER_API_VERSION

VERSION=`cat VERSION.txt`

#docker buildx create --name terra-builder --use

#docker buildx build --platform linux/amd64 \
#  -t us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:${VERSION} \
#  --push .


docker build -t us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:${VERSION} .
docker push  us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:${VERSION}

docker build -t us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:testing .
docker push  us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:testing

# verify
docker run --rm us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:${VERSION} /usr/local/src/LRAA/LRAA --version
