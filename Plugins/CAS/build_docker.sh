#!/bin/bash

set -ex

docker build -t us-central1-docker.pkg.dev/methods-dev-lab/lraa/cas .
docker push us-central1-docker.pkg.dev/methods-dev-lab/lraa/cas 
