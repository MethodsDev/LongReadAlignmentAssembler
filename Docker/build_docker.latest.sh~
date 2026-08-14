#!/bin/bash

set -ex

VERSION=`cat VERSION.txt`

#docker buildx create --name terra-builder --use

docker buildx build --platform linux/amd64 \
  -t us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:${VERSION} \
  --push .


# verify
docker run --rm -it us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:${VERSION} /usr/local/src/LRAA/LRAA --version
