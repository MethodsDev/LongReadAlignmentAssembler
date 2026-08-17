FROM ubuntu:24.04
MAINTAINER bhaas@broadinstitute.org


ENV DEBIAN_FRONTEND=noninteractive

# Every apt layer clears /var/lib/apt/lists afterwards; the lists are ~60MB and
# nothing reads them at runtime.
#
# Deliberately absent:
#   libghc-zlib-dev -- the Haskell zlib binding, which drags in the 667MB GHC
#     toolchain. The C builds here want zlib1g-dev, which is listed below.
#   default-jre -- nothing in the tree invokes java or a jar.
RUN apt-get -qq update && apt-get -qq -y install \
    automake \
    build-essential \
    bzip2 \
    cmake \
    curl \
    fort77 \
    ftp \
    g++ \
    gcc \
    gfortran \
    git \
    libblas-dev \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libdb-dev \
    libjpeg-dev \
    liblzma-dev \
    libncurses-dev \
    libncurses5-dev \
    libpcre3-dev \
    libpng-dev \
    libreadline-dev \
    libssl-dev \
    libtbb-dev \
    libx11-dev \
    libxml2-dev \
    libxt-dev \
    libzmq3-dev \
    make \
    nano \
    perl \
    pigz \
    pkg-config \
    python3 \
    python3-dev \
    python3-pip \
    python3-setuptools \
    rsync \
    texlive-latex-base \
    tzdata \
    unzip \
    wget \
    x11-common \
    zlib1g-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*


## Perl stuff
RUN curl -L https://cpanmin.us | perl - App::cpanminus && \
    cpanm DB_File URI::Escape && \
    rm -rf /root/.cpanm


## set up tool config and deployment area:

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

## Python 3 stuff
RUN ln -sf /usr/bin/python3 /usr/bin/python

WORKDIR $SRC

# One layer, and --no-cache-dir throughout: the wheel cache for the numeric
# stack alone is ~166MB and is never reused inside the image.
# igraph/leidenalg back the Leiden gene clustering in pylib/GeneCommunityCluster.py
# and install from manylinux wheels, so no igraph C headers are needed.
RUN pip install --no-cache-dir --break-system-packages \
    argparse \
    pysam \
    networkx \
    lmdb \
    psutil \
    tqdm \
    intervaltree \
    pytest \
    igraph==0.11.4 \
    leidenalg==0.10.2 \
    "numpy>=1.26" \
    "pandas>=2.2" \
    "pyarrow>=15" \
    "scipy>=1.11" \
    "matplotlib>=3.8"

## Samtools
ENV SAMTOOLS_VERSION=1.13
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure && make && make install && \
    cd $SRC && rm -rf samtools-${SAMTOOLS_VERSION} samtools-${SAMTOOLS_VERSION}.tar.bz2

## HTSlib (for tabix and bgzip)
ENV HTSLIB_VERSION=1.13
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    tar xf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure && make && make install && \
    cd $SRC && rm -rf htslib-${HTSLIB_VERSION} htslib-${HTSLIB_VERSION}.tar.bz2

## minimap2 -- transcriptome read rescue realigns against local transcript
## sequences with it (pylib/IsoformReadRescue.py), and rescue is on by default.
ENV MINIMAP2_VERSION=2.28
RUN wget https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 && \
    tar xf minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 && \
    cp minimap2-${MINIMAP2_VERSION}_x64-linux/minimap2 /usr/local/bin/ && \
    cd $SRC && rm -rf minimap2-${MINIMAP2_VERSION}_x64-linux minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2



## R-stuff
# Install R >= 4.5.0 from CRAN to meet ggrepel/Seurat requirements
RUN apt-get update && apt-get install -y --no-install-recommends \
    software-properties-common \
    dirmngr \
    gnupg \
    ca-certificates && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    apt-get update && \
    apt-get install -y r-base r-base-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# System libraries the CRAN/Bioconductor sources below compile against.
# The GDAL/GEOS/PROJ stack is deliberately absent: it costs ~400MB (libgraphblas
# and boost headers dominate) and no spatial R package (sf, terra, units, s2) is
# installed here.
RUN apt-get update && \
    apt-get install -y --fix-missing \
        libharfbuzz-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libfribidi-dev \
        libglpk-dev \
        libuv1-dev \
        pandoc \
        libtiff5-dev \
        libudunits2-dev \
        patch \
        python3-venv && \
    apt-get clean && rm -rf /var/lib/apt/lists/*


RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 1)' >> /usr/lib/R/etc/Rprofile.site

# One layer so the package caches can be dropped before it is committed, and
# Seurat installed with dependencies = NA (Depends/Imports/LinkingTo only).
# dependencies = TRUE pulled its Suggests as well -- BPCells, Rfast, Rfast2 and
# Rnanoflann alone are ~400MB, and nothing here loads them.
# devtools is not installed: remotes is what install_github needs, and devtools
# adds pak/gert/testthat for no caller.
#
# clustermole is not installed and never has been: it needs GSVA ->
# SpatialExperiment -> magick, magick needs libmagick++-dev, and without those
# headers BiocManager::install only warns, so the layer looked like it had
# succeeded while leaving nothing behind but its partial dependency tree.  Only
# util/sc/notebook_templates/Seurat_get_cell_clusters.Rmd loads it, and no
# workflow renders that notebook.  Add libmagick++-dev here to bring it back.
#
# The final check turns any such silent skip into a build failure.
RUN Rscript --vanilla -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")' && \
    Rscript --vanilla -e 'install.packages(c("tidyverse", "cowplot", "argparse", "remotes", "FNN", "SeuratObject"), repos="https://cloud.r-project.org")' && \
    Rscript --vanilla -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); BiocManager::install(c("DropletUtils", "seqLogo"), ask = FALSE, update = FALSE)' && \
    Rscript --vanilla -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); remotes::install_github("satijalab/seurat", ref = "seurat5", quiet = FALSE, upgrade = "never", dependencies = NA)' && \
    Rscript --vanilla -e 'required <- c("BiocManager", "tidyverse", "cowplot", "argparse", "remotes", "FNN", "SeuratObject", "DropletUtils", "seqLogo", "Seurat"); missing <- setdiff(required, rownames(installed.packages())); if (length(missing)) stop("R packages failed to install: ", paste(missing, collapse = ", ")); invisible(lapply(required, function(p) suppressPackageStartupMessages(library(p, character.only = TRUE)))); cat("Seurat", as.character(packageVersion("Seurat")), "\n")' && \
    rm -rf /root/.cache/R /tmp/Rtmp* /tmp/downloaded_packages


WORKDIR $SRC

## gffcompare
ENV GFFCOMPARE_VERSION=0.12.6
RUN cd $SRC && \
    git clone --depth 1 --branch v${GFFCOMPARE_VERSION} --recurse-submodules https://github.com/gpertea/gffcompare.git && \
    cd gffcompare && \
    make release && \
    cp gffcompare trmap /usr/local/bin/ && \
    cd $SRC && rm -rf gffcompare gclib


WORKDIR $SRC

## Diamond (for fast protein homology search)
ENV DIAMOND_VERSION=2.1.8
RUN wget https://github.com/bbuchfink/diamond/releases/download/v${DIAMOND_VERSION}/diamond-linux64.tar.gz && \
    tar xzf diamond-linux64.tar.gz && \
    mv diamond /usr/local/bin/ && \
    rm diamond-linux64.tar.gz

## NCBI BLAST+ (alternative to Diamond, via --blast_tool blastp)
# Only the two programs TransDecoder's blastp path runs are kept; the full
# release is 672MB of which the rest is never invoked.
ENV BLAST_VERSION=2.14.0
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz && \
    tar xzf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz && \
    cp ncbi-blast-${BLAST_VERSION}+/bin/blastp ncbi-blast-${BLAST_VERSION}+/bin/makeblastdb /usr/local/bin/ && \
    rm -rf ncbi-blast-${BLAST_VERSION}+ ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz

## TransDecoder (for ORF prediction)
ENV TRANSDECODER_VERSION=6.0.0
RUN wget https://github.com/TransDecoder/TransDecoder/archive/TransDecoder-v${TRANSDECODER_VERSION}.tar.gz && \
    tar xzf TransDecoder-v${TRANSDECODER_VERSION}.tar.gz && \
    mv TransDecoder-TransDecoder-v${TRANSDECODER_VERSION} /usr/local/src/TransDecoder && \
    cd /usr/local/src/TransDecoder && \
    make && \
    rm -f $SRC/TransDecoder-v${TRANSDECODER_VERSION}.tar.gz

ENV TRANSDECODER_DIR=/usr/local/src/TransDecoder
ENV PATH="${TRANSDECODER_DIR}:${PATH}"

WORKDIR $SRC


ARG LRAA_VERSION
ARG LRAA_CO
ENV LRAA_VERSION=${LRAA_VERSION}
ENV LRAA_CO=${LRAA_CO}

# Real provenance, and the only provenance: the SHA the checkout below was
# fetched with.  Readable without running the image, unlike an ENV.
LABEL org.opencontainers.image.revision=${LRAA_CO}


# The pinned commit as a source archive rather than a clone: no runtime code
# consults git, and a clone of this repository carries ~85MB of pack history.
# testing/ holds ~155MB of fixture data that only the repo's own Makefiles read,
# from a host checkout, so it is not shipped.
#
# The archive is now staged into the build context by the build script's
# `git archive`, not fetched from GitHub. Same reasoning, one step further: the
# GitHub archive of this repo is 88MB precisely BECAUSE it contains testing/, so
# fetching it and then deleting testing/ downloaded ~84MB in order to discard it.
# git archive reads content-addressed objects, so the tarball is still
# byte-determined by the SHA, and testing/ is excluded before it is ever sent.
COPY lraa_checkout.tar.gz lraa_checkout.sha /tmp/
RUN if [ -z "${LRAA_CO}" ]; then echo "build arg LRAA_CO is required; the build scripts pass git rev-parse HEAD" >&2; exit 1; fi; \
    if [ "`cat /tmp/lraa_checkout.sha`" != "${LRAA_CO}" ]; then \
        echo "the staged checkout is `cat /tmp/lraa_checkout.sha` but LRAA_CO is ${LRAA_CO}: the build context holds a stale tarball, so this image would be labelled with a commit it does not contain. Re-run the build script, which regenerates it." >&2; \
        exit 1; \
    fi; \
    cd ${SRC} && \
    tar xzf /tmp/lraa_checkout.tar.gz && \
    rm -f /tmp/lraa_checkout.tar.gz /tmp/lraa_checkout.sha



ENV PATH="$SRC/LRAA/:$SRC/LRAA/util/:$SRC/LRAA/util/sc/:$SRC/LRAA/util/misc:${PATH}"
