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

RUN pip install argparse pysam networkx lmdb psutil --break-system-packages

WORKDIR $SRC

RUN pip install intervaltree --break-system-packages

## Samtools
ENV SAMTOOLS_VERSION=1.13
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure && make && make install



## R-stuff
RUN apt-get install -y r-base r-base-dev

RUN R -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'

RUN apt-get update && \
    apt-get install -y --fix-missing \
        libharfbuzz-dev \
        libfribidi-dev \
        libtiff5-dev

    
RUN apt-get install -y libharfbuzz-dev libfribidi-dev libtiff5-dev 

RUN R -e 'BiocManager::install("tidyverse")'
RUN R -e 'BiocManager::install("cowplot")'
RUN R -e 'BiocManager::install("argparse")'
RUN R -e 'BiocManager::install("devtools")'
RUN R -e 'BiocManager::install("Matrix")'
RUN R -e 'BiocManager::install("FNN")'
RUN R -e 'BiocManager::install("gtable")'
RUN R -e 'BiocManager::install("pak")'
RUN R -e 'pak::pak("igraph/rigraph")'
RUN R -e 'BiocManager::install("leiden")'

RUN R -e 'remotes::install_github("satijalab/seurat", "seurat5", quiet = FALSE)'

RUN R -e 'BiocManager::install("clustermole")'

RUN pip install numpy --break-system-packages

RUN pip install pytest --break-system-packages

RUN python -m pip install \
    "numpy>=1.26" \
    "pandas>=2.2" \
    "pyarrow>=15" \
    "scipy>=1.11" \
    "matplotlib>=3.8" --break-system-packages

        
WORKDIR $SRC

ENV LRAA_VERSION="v0.9.0-postdev2"
ENV LRAA_CO=b82ef247b0fa99c71bf431cc2576d02f9c9aa3a6


# update VERSION and LRAA_CO files with version info and checkout

RUN git clone --recursive https://github.com/MethodsDev/LongReadAlignmentAssembler.git LRAA && \
    cd LRAA && \
    git checkout ${LRAA_CO}



ENV PATH="$SRC/LRAA/:$SRC/LRAA/util/:$SRC/LRAA/util/sc/:$SRC/LRAA/util/misc:${PATH}"


