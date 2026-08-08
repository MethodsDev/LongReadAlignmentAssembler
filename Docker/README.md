# LRAA Docker images

Four images are built from this directory, and a fifth name, plain `lraa`, is an
alias for the smallest of them.

| image | Dockerfile | size | contents |
|---|---|---|---|
| `lraa-core` | `Dockerfile.core` | 422 MB | LRAA and what an LRAA run reaches: python with pysam, networkx, intervaltree, tqdm, lmdb, psutil, numpy, igraph and leidenalg; samtools; htslib; minimap2; gffcompare; perl |
| `lraa-orf` | `Dockerfile.orf` | 615 MB | `FROM lraa-core`, plus TransDecoder, diamond, and the `blastp`/`makeblastdb` pair TransDecoder's `--blast_tool blastp` path runs |
| `lraa-sc` | `Dockerfile.sc` | 2.69 GB | `FROM lraa-core`, plus R with Seurat, DropletUtils, tidyverse, edgeR and limma, and pandas, scipy, matplotlib, seaborn, statsmodels, pytest |
| `lraa-combined` | `Dockerfile` | 3.95 GB | everything in one image, as it was before the split |
| `lraa` | none | 422 MB | the same digest as `lraa-core`, pushed under the plain name |

Registry: `us-central1-docker.pkg.dev/methods-dev-lab/lraa/`

We pay egress on every pull, so the plain `lraa` name has to be the cheapest
thing that can run an assembly. Someone who pulls it to assemble isoforms gets
422 MB instead of the 3.95 GB of R and TransDecoder they will never invoke.
`docker tag` gives the two names one digest, so the alias costs no extra
storage and no extra build.

Before v0.18.1 the plain name held the combined image. Anyone who needs that
content should switch to `lraa-combined`; `lraa-core` covers isoform discovery
and quantification on its own.

The specialized images derive from the core one, so a workflow that mixes them
pulls the shared layers once.

## Which image a task runs on

Sorted by what each task's `command` block actually executes. Of the 27 tasks
in `WDL/`, 18 run on core, 8 on sc, 1 on orf.

| image | tasks |
|---|---|
| core | `LRAA_runner_task`; the five tasks in `LRAA.wdl` (`merge_GTFs`, `mergeQuantResults`, `mergeGenomeTxArbSummaries`, `count_bam`, `filterBamToSecondaryRescue`); `partition_by_chromosome_task`; `normalize_bam_by_strand`; `partition_bam_by_cell_cluster`; `run_gffcompare`; `incorporate_gene_symbols_sc`; `merge_bams`; `validate_pre_normalized_inputs`; the tar, merge and pseudobulk tasks in `LRAA-cell_cluster_guided.wdl` |
| sc | `run_seurat_from_gene_sparseM`; `run_filter_good_cells`; `sc_build_sparse_matrices`; `sc_build_sparse_matrices_from_tracking`; `LRAA_sqanti_like_reads_eval_task`; the two `summarize_mult_samples` tasks; `RunSaturation` |
| orf | `LRAA_ORF_prediction_task` |

The core set is everything that runs LRAA, samtools, minimap2, gffcompare, perl
or plain-stdlib python. The sc set is every task that invokes `Rscript`, plus
the two python helpers that import pandas and scipy
(`util/sc/singlecell_tracking_to_sparse_matrix.py` and
`util/misc/read_FSM_and_identifiability_saturation_fit.py`).

Workflows expose the choice through inputs. `docker` names the core image
everywhere. `LRAA-singlecell.wdl` and `LRAA-cell_cluster_guided.wdl` also
declare `docker_sc`, and the former passes its value into the latter so the
nested single-cell task follows the caller rather than a default two levels
down.

## Tags

| tag | set by | used by |
|---|---|---|
| `latest` | `build_docker.latest.sh` | defaults written inside the `.wdl` files |
| `<version>` from `VERSION.txt` | `build_docker.versioned.sh` | release pins |
| `testing` | `build_docker.versioned.sh` | input JSONs under `testing/` |

## Building

```bash
cd Docker
bash build_docker.versioned.sh   # <version> and testing
bash build_docker.latest.sh      # latest
```

Both scripts build core first and pass its tag to the other two as
`--build-arg LRAA_CORE_IMAGE=`, so the derived images build against the core
image just made rather than whatever copy the registry holds. The versioned
script applies the `testing` tag with `docker tag` instead of a second build; a
rebuild costs about an hour of R compilation for a byte-identical result.

Each script finishes by retagging the core image it just built as plain `lraa`
and pushing that too, so the alias never drifts from `lraa-core`.

A cold run of the versioned script takes about an hour, nearly all of it
compiling Seurat and its dependencies from source for `lraa-sc` and for the
combined image. Core and orf take a few minutes between them. The second script
then reuses those layers and only tags and pushes.

## Releasing a new version

1. Update `VERSION` in the top-level `LRAA` script and `Docker/VERSION.txt` to
   the same value.
2. Commit and push the code.
3. Set `ENV LRAA_CO` in both `Dockerfile` and `Dockerfile.core` to the commit
   SHA you just pushed, and `ENV LRAA_VERSION` to the version. The images fetch
   `https://github.com/MethodsDev/LongReadAlignmentAssembler/archive/${LRAA_CO}.tar.gz`
   rather than cloning, so this pin is the only thing that decides which code
   ships.
4. Commit the pin, then run both build scripts.
5. Confirm what shipped:

```bash
docker run --rm <registry>/lraa-core:<version> \
  bash -c 'LRAA --version; minimap2 --version'
```

## Adding a dependency

Put it in the smallest image that needs it. Anything reachable from an `LRAA`
run belongs in `Dockerfile.core` and is paid for by every task; anything an R
script or a pandas/scipy helper needs belongs in `Dockerfile.sc`.

Two dependencies are load-bearing in ways the imports do not show, both found
by running the images rather than reading the source:

- `igraph` and `leidenalg` are imported under `try/except`, but
  `Transcript.recluster_transcripts_to_genes` raises without them on the default
  discovery path. They belong in core.
- `pytest` is imported at module scope by `pylib/SQANTI_like_annotator.py`, so
  the image running the SQANTI-like task needs it.

The R layer in `Dockerfile` and `Dockerfile.sc` ends by loading every package it
claims to install. Keep new packages in that list. `BiocManager::install` only
warns when a dependency fails to build, so without the check a layer can exit 0
having installed nothing. That is how `clustermole` was missing from every image
built up to v0.18.1 while the layer that installed it reported success.

## What is deliberately absent

- `clustermole`, which needs GSVA, SpatialExperiment and magick, and so needs
  `libmagick++-dev`. Only `util/sc/notebook_templates/Seurat_get_cell_clusters.Rmd`
  loads it and no workflow renders that notebook. Add `libmagick++-dev` to the
  apt layer and the package back to the R list to restore it, at a cost of about
  250 MB.
- The GDAL, GEOS and PROJ stack, worth about 400 MB. No spatial R package is
  installed.
- `default-jre`. Nothing in the tree invokes java.
- Seurat's `Suggests`, worth about 620 MB and led by BPCells at 219 MB. Seurat
  installs with `dependencies = NA`.
- Most of NCBI BLAST+. TransDecoder's blastp path runs `blastp` and
  `makeblastdb`; the other thirty programs cost 619 MB and are never invoked.
- `libghc-zlib-dev`, a Haskell binding that drags in the 667 MB GHC toolchain.
  The C builds want `zlib1g-dev`.
- The repository's `testing/` fixtures, about 240 MB with the git pack that
  carries them. Only the repo's own Makefiles read them, from a host checkout.

## minimap2

Transcriptome read rescue realigns reads against local transcript sequences with
minimap2, and rescue is on by default. Every image installs it. If it is missing,
LRAA now exits at startup naming the setting that requires it rather than
skipping the rescue and reporting the smaller read count as the answer.
