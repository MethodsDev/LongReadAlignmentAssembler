#!/bin/bash

VERSION=`cat VERSION.txt`

LRAA_CO=`cat LRAA_CO.txt`


if [ ! -d "LRAA" ]; then
    git clone --recursive git@github.com:MethodsDev/LongReadAlignmentAssembler.git LRAA
fi

cd LRAA && git checkout ${LRAA_CO} && cd ../

docker build -t lraa/lraa:${VERSION} .
docker build -t lraa/lraa:latest .
docker build -t  us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:${VERSION} .
docker build -t  us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest .
