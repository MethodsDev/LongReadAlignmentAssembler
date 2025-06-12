#!/bin/bash

VERSION=`cat VERSION.txt`

singularity build lraa.v${VERSION}.simg docker://us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:${VERSION} 

singularity exec -e lraa.v${VERSION}.simg LRAA --version


