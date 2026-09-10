# Single-cell, differential-usage and plotting work: the LRAA checkout on top of
# the R and scientific-Python tier.
#
# Everything expensive is in Dockerfile.sc-base, published as lraa-sc-base and
# consumed here by tag.  This file adds one layer, so cutting a release costs
# seconds rather than an hour of recompiling Seurat.  Add packages to
# Dockerfile.sc-base, not here.
#
# Tasks that need this image rather than the core one:
#   run_seurat_from_gene_sparseM        (WDL/subwdls/LRAA-gene_sparseM_to_seurat_clusters.wdl)
#   run_filter_good_cells               (WDL/subwdls/LRAA-filter_good_cells.wdl)
#   sc_build_sparse_matrices[_from_tracking]
#   LRAA_sqanti_like_reads_eval_task and the two multi-sample summary tasks
#   RunSaturation                       (WDL/FSM_and_isoform_identifiability_saturation.wdl)

ARG LRAA_SC_BASE_IMAGE=lraa-sc-base:latest
FROM ${LRAA_SC_BASE_IMAGE}

ARG LRAA_VERSION
ARG LRAA_CO
ENV LRAA_VERSION=${LRAA_VERSION}
ENV LRAA_CO=${LRAA_CO}

# Real provenance, and the only provenance: the SHA the checkout below was
# fetched with.  Readable without running the image, unlike an ENV.
LABEL org.opencontainers.image.revision=${LRAA_CO}

# Last layer, so a version bump reuses everything above.
COPY lraa_checkout.tar.gz lraa_checkout.sha /tmp/
RUN if [ -z "${LRAA_CO}" ]; then echo "build arg LRAA_CO is required; the build scripts pass git rev-parse HEAD" >&2; exit 1; fi; \
    if [ "`cat /tmp/lraa_checkout.sha`" != "${LRAA_CO}" ]; then \
        echo "the staged checkout is `cat /tmp/lraa_checkout.sha` but LRAA_CO is ${LRAA_CO}: the build context holds a stale tarball, so this image would be labelled with a commit it does not contain. Re-run the build script, which regenerates it." >&2; \
        exit 1; \
    fi; \
    cd ${SRC} && \
    tar xzf /tmp/lraa_checkout.tar.gz && \
    rm -f /tmp/lraa_checkout.tar.gz /tmp/lraa_checkout.sha
