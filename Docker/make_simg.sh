#!/bin/bash

VERSION=`cat VERSION.txt`

# lraa-core is what an LRAA run needs.  Swap in lraa-combined if you also want
# R, Seurat and TransDecoder inside the container.
singularity build lraa.v${VERSION}.simg docker://us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:${VERSION}

singularity exec -e lraa.v${VERSION}.simg LRAA --version


