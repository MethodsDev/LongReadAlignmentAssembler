#!/bin/bash

set -ex


VERSION=latest


#docker buildx create --name terra-builder --use

#docker buildx build --platform linux/amd64 \
#  -t us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest \
#  --push .


docker build -t us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:${VERSION} .
docker push  us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:${VERSION}

# verify
docker run --rm -it us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:${VERSION} /usr/local/src/LRAA/LRAA --version


