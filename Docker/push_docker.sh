#!/bin/bash

VERSION=`cat VERSION.txt`

docker push us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:${VERSION}
docker push us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest

