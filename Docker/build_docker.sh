#!/bin/bash

VERSION=`cat VERSION.txt`

docker build -t lraa/lraa:${VERSION} .
docker build -t lraa/lraa:latest .
docker build -t  us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:${VERSION} .
docker build -t  us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest .
