#!/bin/bash

set -ex

#docker buildx create --name terra-builder --use

docker buildx build --platform linux/amd64 \
  -t us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest \
  --push .


# verify
docker pull us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest
docker run --rm -it us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest /usr/local/src/LRAA/LRAA --version
