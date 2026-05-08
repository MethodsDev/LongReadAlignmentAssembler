#!/bin/bash

set -ex

# Newer Docker daemons reject old client API pins inherited from the shell.
# Clear any stale override so the installed client can negotiate its default.
unset DOCKER_API_VERSION


VERSION=latest


#docker buildx create --name terra-builder --use

#docker buildx build --platform linux/amd64 \
#  -t us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest \
#  --push .


docker build -t us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:${VERSION} .
docker push  us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:${VERSION}

# verify
docker run --rm us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:${VERSION} /usr/local/src/LRAA/LRAA --version
