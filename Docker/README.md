# LRAA Docker images

Four images are published from this directory, plus the plain `lraa` name as an
alias for the smallest of them. A fifth, `lraa-base`, is a build input that
carries the shared dependencies and is never pushed.

| image | Dockerfile | size | contents |
|---|---|---|---|
| `lraa-base` | `Dockerfile.base` | 400 MB | not published. Python with pysam, networkx, intervaltree, tqdm, lmdb, psutil, numpy, igraph and leidenalg; samtools; htslib; minimap2; gffcompare; perl |
| `lraa-core` | `Dockerfile.core` | 422 MB | `FROM lraa-base`, plus the LRAA checkout |
| `lraa-orf` | `Dockerfile.orf` | 615 MB | `FROM lraa-base`, plus TransDecoder, diamond, the `blastp`/`makeblastdb` pair TransDecoder's `--blast_tool blastp` path runs, and the LRAA checkout |
| `lraa-sc` | `Dockerfile.sc` | 2.69 GB | `FROM lraa-base`, plus R with Seurat, DropletUtils, tidyverse, edgeR and limma, pandas, scipy, matplotlib, seaborn, statsmodels, pytest, and the LRAA checkout |
| `lraa-combined` | `Dockerfile` | 3.95 GB | everything in one image, as it was before the split |
| `lraa` | none | 422 MB | the same digest as `lraa-core`, pushed under the plain name |

Registry: `us-central1-docker.pkg.dev/methods-dev-lab/lraa/`

We pay egress on every pull, so the plain `lraa` name has to be the cheapest
thing that can run an assembly. Someone who pulls it to assemble isoforms gets
422 MB instead of the 3.95 GB of R and TransDecoder they will never invoke.
`docker tag` gives the two names one digest, so the alias costs no extra
storage and no extra build.

The plain name used to hold the combined image, through the 0.18.0 tags. Anyone
who needs that content should switch to `lraa-combined`; `lraa-core` covers
isoform discovery and quantification on its own.

## Layer order, and why lraa-base exists

The LRAA checkout is the final layer of every image. Docker keys a layer's
cache on its parent image id, so anything sitting above the checkout is
rebuilt whenever the pinned commit moves.

`lraa-sc` and `lraa-orf` used to build `FROM lraa-core`, which put the R and
TransDecoder layers above core's checkout. Bumping a version then recompiled
Seurat from source, and a release that changed nothing but a git SHA took an
hour. Building all three from a dependency-only base instead costs one small
layer per image.

The three published images therefore share `lraa-base` rather than deriving
from each other, so a workflow that mixes them still pulls the shared layers
once. What is no longer shared is the checkout itself, a few MB repeated in
each image, which is the price of the arrangement.

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

Both scripts build `lraa-base` first and pass it to the other builds as
`--build-arg LRAA_BASE_IMAGE=`, so the three published images are built against
the base just made rather than a stale local copy. The base is not pushed. The
version and the pinned commit are passed as `--build-arg LRAA_VERSION=` and
`--build-arg LRAA_CO=`, read from `VERSION.txt` and `LRAA_CO.txt`, so the SHA
lives in one file instead of four Dockerfiles that can drift apart.

The versioned script applies the `testing` tag with `docker tag` instead of a
second build; a rebuild would cost an hour of R compilation for a byte-identical
result. Each script finishes by retagging the core image as plain `lraa` and
pushing that too, so the alias never drifts from `lraa-core`.

A cold run takes about an hour, nearly all of it compiling Seurat from source
for `lraa-sc` and for the combined image. A run that changes only `LRAA_CO.txt`
and `VERSION.txt` takes a few minutes, because the checkout is the last layer of
every image.

## Releasing a new version

1. Update `VERSION` in the top-level `LRAA` script and `Docker/VERSION.txt` to
   the same value.
2. Commit and push the code.
3. Put the commit SHA you just pushed in `Docker/LRAA_CO.txt`. The images fetch
   `https://github.com/MethodsDev/LongReadAlignmentAssembler/archive/${LRAA_CO}.tar.gz`
   rather than cloning, so this one line decides which code ships.
4. Commit the pin, then run both build scripts.
5. Confirm what shipped:

```bash
docker run --rm <registry>/lraa-core:<version> \
  bash -c 'LRAA --version; minimap2 --version'
```

## Adding a dependency

Put it in the smallest place that covers everything needing it. Anything
reachable from an `LRAA` run belongs in `Dockerfile.base`, since all three
published images run LRAA and each pays for the base. `Dockerfile.core` adds
nothing but the checkout, so a dependency put there would be missing from
`lraa-sc` and `lraa-orf`, which build from the base rather than from core.
Anything only an R script or a pandas/scipy helper needs belongs in
`Dockerfile.sc`; anything only TransDecoder needs belongs in `Dockerfile.orf`.
Whatever you add, add it to `Dockerfile` too, which is standalone.

Two dependencies are load-bearing in ways the imports do not show, both found
by running the images rather than reading the source:

- `igraph` and `leidenalg` are imported under `try/except`, but
  `Transcript.recluster_transcripts_to_genes` raises without them on the default
  discovery path. They belong in the base, so that every image running LRAA has
  them.
- `pytest` is imported at module scope by `pylib/SQANTI_like_annotator.py`, so
  the image running the SQANTI-like task needs it.

The R layer in `Dockerfile` and `Dockerfile.sc` ends by loading every package it
claims to install. Keep new packages in that list. `BiocManager::install` only
warns when a dependency fails to build, so without the check a layer can exit 0
having installed nothing. That is how `clustermole` was missing from every image
built before the check existed while the layer that installed it reported
success.

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
