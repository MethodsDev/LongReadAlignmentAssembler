# Single-cell, differential-usage and plotting work: the R stack (Seurat,
# DropletUtils, tidyverse, edgeR) and the scientific Python stack that the
# util/sc and util/misc helpers need.
#
# Builds from Dockerfile.base rather than from lraa-core, so that bumping the
# pinned LRAA commit does not invalidate the R layers below and force an hour of
# recompiling Seurat.  The LRAA checkout is the final layer instead.
#
# Tasks that need this image rather than the core one:
#   run_seurat_from_gene_sparseM        (WDL/subwdls/LRAA-gene_sparseM_to_seurat_clusters.wdl)
#   run_filter_good_cells               (WDL/subwdls/LRAA-filter_good_cells.wdl)
#   sc_build_sparse_matrices[_from_tracking]
#   LRAA_sqanti_like_reads_eval_task and the two multi-sample summary tasks
#   RunSaturation                       (WDL/FSM_and_isoform_identifiability_saturation.wdl)

ARG LRAA_BASE_IMAGE=lraa-base:latest
FROM ${LRAA_BASE_IMAGE}

ENV DEBIAN_FRONTEND=noninteractive

# R itself, plus the system libraries the CRAN/Bioconductor sources compile
# against.  No GDAL/GEOS/PROJ: no spatial R package is installed here.
RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends \
    ca-certificates \
    dirmngr \
    gnupg \
    software-properties-common \
    wget && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        libbz2-dev \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libfribidi-dev \
        libglpk-dev \
        libharfbuzz-dev \
        libicu-dev \
        libjpeg-dev \
        liblzma-dev \
        libpng-dev \
        libssl-dev \
        libtiff5-dev \
        libudunits2-dev \
        libuv1-dev \
        libxml2-dev \
        pandoc \
        patch \
        pkg-config \
        r-base \
        r-base-dev \
        zlib1g-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)' >> /usr/lib/R/etc/Rprofile.site

# One layer, caches dropped before it is committed.  Seurat with
# dependencies = NA: its Suggests (BPCells, Rfast, Rfast2, Rnanoflann) are ~400MB
# and nothing in this repository loads them.
#
# clustermole is absent, as it is in the current production image: it needs
# GSVA -> SpatialExperiment -> magick, magick needs libmagick++-dev, and without
# those headers BiocManager::install only warns and moves on.  Only
# util/sc/notebook_templates/Seurat_get_cell_clusters.Rmd loads it and no
# workflow renders that notebook; add libmagick++-dev above to bring it back.
#
# The final check turns any such silent skip into a build failure.
RUN Rscript --vanilla -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")' && \
    Rscript --vanilla -e 'install.packages(c("argparse", "tidyverse", "cowplot", "data.table", "Matrix", "gridExtra", "pheatmap", "patchwork", "ggdendro", "ggrepel", "knitr", "future", "remotes"), repos="https://cloud.r-project.org")' && \
    Rscript --vanilla -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); BiocManager::install(c("DropletUtils", "edgeR", "limma"), ask = FALSE, update = FALSE)' && \
    Rscript --vanilla -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); remotes::install_github("satijalab/seurat", ref = "seurat5", upgrade = "never", dependencies = NA)' && \
    Rscript --vanilla -e 'required <- c("argparse", "tidyverse", "cowplot", "data.table", "Matrix", "gridExtra", "pheatmap", "patchwork", "ggdendro", "ggrepel", "knitr", "future", "remotes", "DropletUtils", "edgeR", "limma", "Seurat"); missing <- setdiff(required, rownames(installed.packages())); if (length(missing)) stop("R packages failed to install: ", paste(missing, collapse = ", ")); invisible(lapply(required, function(p) suppressPackageStartupMessages(library(p, character.only = TRUE)))); cat("Seurat", as.character(packageVersion("Seurat")), "\n")' && \
    rm -rf /root/.cache/R /tmp/Rtmp* /tmp/downloaded_packages

# Scientific Python for util/sc and util/misc: sparse matrix construction
# (pandas/scipy), plotting (matplotlib/seaborn), differential usage
# (statsmodels), and pytest, which pylib/SQANTI_like_annotator.py imports at
# module scope.
RUN pip install --no-cache-dir --break-system-packages \
    pandas \
    scipy \
    matplotlib \
    seaborn \
    statsmodels \
    pytest

ARG LRAA_VERSION
ARG LRAA_CO
ENV LRAA_VERSION=${LRAA_VERSION}
ENV LRAA_CO=${LRAA_CO}

# Last layer, so a version bump reuses everything above.
RUN if [ -z "${LRAA_CO}" ]; then echo "build arg LRAA_CO is required; see Docker/LRAA_CO.txt" >&2; exit 1; fi; \
    cd ${SRC} && \
    curl -sSL https://github.com/MethodsDev/LongReadAlignmentAssembler/archive/${LRAA_CO}.tar.gz -o lraa.tar.gz && \
    tar xzf lraa.tar.gz && \
    mv LongReadAlignmentAssembler-${LRAA_CO} LRAA && \
    rm -rf lraa.tar.gz LRAA/testing
