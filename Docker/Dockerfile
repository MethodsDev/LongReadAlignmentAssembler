FROM ubuntu:24.04
MAINTAINER bhaas@broadinstitute.org


ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update && apt-get -qq -y install \
    automake \
    build-essential \
    bzip2 \
    cmake \
    curl \
    default-jre \
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
    libghc-zlib-dev \
    libjpeg-dev \
    liblzma-dev \
    libncurses-dev \
    libncurses5-dev \
    libpcre3-dev \
    libpng-dev \
    libreadline-dev \
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
    zlib1g-dev



## Perl stuff
RUN curl -L https://cpanmin.us | perl - App::cpanminus
RUN cpanm install DB_File
RUN cpanm install URI::Escape


## set up tool config and deployment area:

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

## Python 3 stuff
RUN ln -sf /usr/bin/python3 /usr/bin/python

RUN pip install argparse pysam networkx lmdb psutil tqdm --break-system-packages

WORKDIR $SRC

RUN pip install intervaltree --break-system-packages
## Leiden community clustering dependencies (Python)
# igraph requires C core; install build and runtime deps first
RUN apt-get -qq update && apt-get -qq -y install \
    libxml2-dev \
    libigraph0-dev || true
RUN pip install --break-system-packages igraph==0.11.4 leidenalg==0.10.2

## Samtools
ENV SAMTOOLS_VERSION=1.13
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure && make && make install

## HTSlib (for tabix and bgzip)
ENV HTSLIB_VERSION=1.13
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    tar xvf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure && make && make install



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
    apt-get install -y r-base r-base-dev

RUN R -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")'

RUN apt-get update && \
    apt-get install -y --fix-missing \
        libharfbuzz-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libfribidi-dev \
        libgdal-dev \
        libgeos-dev \
        libglpk-dev \
        libuv1-dev \
        libproj-dev \
        pandoc \
        libtiff5-dev \
        libudunits2-dev \
        patch \
        python3-venv

    
RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 1)' >> /usr/lib/R/etc/Rprofile.site
RUN Rscript --vanilla -e 'install.packages("tidyverse", repos="https://cloud.r-project.org")'
RUN Rscript --vanilla -e 'install.packages("cowplot", repos="https://cloud.r-project.org")'
RUN Rscript --vanilla -e 'install.packages("argparse", repos="https://cloud.r-project.org")'
RUN Rscript --vanilla -e 'install.packages("devtools", repos="https://cloud.r-project.org")'
RUN Rscript --vanilla -e 'install.packages("remotes", repos="https://cloud.r-project.org")'
RUN Rscript --vanilla -e 'install.packages("FNN", repos="https://cloud.r-project.org")'
RUN Rscript --vanilla -e 'BiocManager::install("DropletUtils", ask = FALSE, update = FALSE)'
RUN Rscript --vanilla -e 'BiocManager::install("seqLogo", ask = FALSE, update = FALSE)'
RUN Rscript --vanilla -e 'BiocManager::install("clustermole", ask = FALSE, update = FALSE)'
RUN Rscript --vanilla -e 'install.packages("SeuratObject", repos="https://cloud.r-project.org")'
RUN Rscript --vanilla -e 'remotes::install_github("satijalab/seurat", ref = "seurat5", quiet = FALSE, upgrade = "never", dependencies = TRUE); stopifnot("Seurat" %in% rownames(installed.packages())); suppressPackageStartupMessages(library(Seurat)); cat(as.character(packageVersion("Seurat")), "\n")'

RUN pip install numpy --break-system-packages

RUN pip install pytest --break-system-packages

RUN python -m pip install \
    "numpy>=1.26" \
    "pandas>=2.2" \
    "pyarrow>=15" \
    "scipy>=1.11" \
    "matplotlib>=3.8" --break-system-packages

        
WORKDIR $SRC

## gffcompare
ENV GFFCOMPARE_VERSION=0.12.6
RUN cd $SRC && \
    git clone --depth 1 --branch v${GFFCOMPARE_VERSION} --recurse-submodules https://github.com/gpertea/gffcompare.git && \
    cd gffcompare && \
    make release && \
    cp gffcompare trmap /usr/local/bin/ && \
    cd $SRC && rm -rf gffcompare


WORKDIR $SRC

## Diamond (for fast protein homology search)
ENV DIAMOND_VERSION=2.1.8
RUN wget https://github.com/bbuchfink/diamond/releases/download/v${DIAMOND_VERSION}/diamond-linux64.tar.gz && \
    tar xzf diamond-linux64.tar.gz && \
    mv diamond /usr/local/bin/ && \
    rm diamond-linux64.tar.gz

## NCBI BLAST+ (alternative to Diamond)
ENV BLAST_VERSION=2.14.0
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz && \
    tar xzf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz && \
    cp ncbi-blast-${BLAST_VERSION}+/bin/* /usr/local/bin/ && \
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


ENV LRAA_VERSION="v0.16.3"
ENV LRAA_CO=52d305d428a4adc24323e5ba577bac2404e61bce


# update VERSION and LRAA_CO files with version info and checkout

RUN git clone --recursive https://github.com/MethodsDev/LongReadAlignmentAssembler.git LRAA && \
    cd LRAA && \
    git checkout ${LRAA_CO}



ENV PATH="$SRC/LRAA/:$SRC/LRAA/util/:$SRC/LRAA/util/sc/:$SRC/LRAA/util/misc:${PATH}"
