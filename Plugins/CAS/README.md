# Cellarium CAS plugin

Runs the Cellarium Cell Annotation Service (CAS) on an LRAA single-cell matrix
directory and writes `<prefix>.h5ad` and `<prefix>.tsv` with the top-k most granular
cell type calls per cell.

## Example usage

    # token in a private file (chmod 600); it is never put on the command line
    ./run_CAS_docker.sh \
        --matrix-dir refQuant^gene-sparseM \
        --api-token-file ~/.cas_token \
        --output-prefix cas.test.out

or directly with docker:

    export CAS_API_TOKEN=$(cat ~/.cas_token)
    docker run --rm -e CAS_API_TOKEN -v `pwd`:/data us-central1-docker.pkg.dev/methods-dev-lab/lraa/cas \
         --matrix-dir /data/refQuant^gene-sparseM \
         --output-prefix /data/cas.test.out

In the WDL, prefer `api_token_file` over `api_token`: a String input is written into
the task's command script.

## Notes

- **Client version.** The image pins `cellarium-cas` exactly (see `Dockerfile`). The CAS
  server refuses clients older than its minimum version; as of 2026-11-18 that is 1.9.1.
  When runs fail with `ClientTooOldError`, bump the pin and rebuild with `build_docker.sh`.
- **Model.** The default is `scanvi-classification-v8.hs-general`, the best of the 1.9 human
  models on our single-nucleus brain data. The CAS server's own default is still the older
  `cellarium-pca-v1.1beta`; request it with `--cas-model-name server-default`. Avoid
  `pca-v2.hs-general` on brain, which assigned about a quarter of cortical nuclei to other
  tissues. `CASClient(api_token=...).print_models()` lists the available models.
- **Empty cells.** CAS labels every cell it receives, including cells with no counts, often
  with high scores. Cells below `--min-counts` (default 1) are dropped before the request and
  are absent from the outputs.
- **Single-nucleus data.** Most nuclear reads are unspliced. A quant-only matrix counts only
  reads compatible with annotated transcripts and undercounts long, intron-rich genes (in
  our cortex data by about 5-fold for genes over 1 Mb). Use the gene matrix from a
  ref-guided LRAA run, where intronic transcripts carry their host gene's symbol; it raised
  agreement with marker-based labels by about 3 points for every model.
- **Feature-overlap warning.** Client 1.9.1 warns when under half of the model's
  61,497-gene schema is present. LRAA writes only genes with reads, so typical matrices hold
  about 45% and trigger it. Missing genes are filled with zeros, and the warning can be
  ignored when the matrix is otherwise normal.
