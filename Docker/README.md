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
| core | `LRAA_runner_task`; the five tasks in `LRAA.wdl` (`merge_GTFs`, `mergeQuantResults`, `mergeReadAssignmentSummaries`, `count_bam`, `filterBamToSecondaryRescue`); `partition_by_chromosome_task`; `normalize_bam_by_strand`; `partition_bam_by_cell_cluster`; `run_gffcompare`; `incorporate_gene_symbols_sc`; `merge_bams`; `validate_pre_normalized_inputs`; the tar, merge and pseudobulk tasks in `LRAA-cell_cluster_guided.wdl` |
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

`LRAA-ORF-prediction.wdl` names the orf image, under the same `docker` input.
A workflow with two image inputs needs both set to override it: an unset one
keeps its `:latest` default, and half the run tests the last release.  The test
targets under `testing/` set them from one place, `testing/lraa_test_docker.mk`,
which names the core, sc and orf images and the tag they share.

## Tags

| tag | set by | built from | used by |
|---|---|---|---|
| `latest` | `build_docker.latest.sh` | `git rev-parse HEAD` | defaults written inside the `.wdl` files |
| `<version>` from `VERSION.txt` | `build_docker.versioned.sh` | `git rev-parse HEAD` | release pins |
| `testing` | `build_docker.testing.sh` | `git rev-parse HEAD` | local WDL test targets, through `testing/lraa_test_docker.mk` |
| `<version>-testing` | `build_docker.testing.sh` | `git rev-parse HEAD` | testing that outlives the moving tag, e.g. runs dispatched to VMs |

Neither testing tag is a release artifact. Both come out of the same build, so
they cannot drift apart, and only the release scripts write `latest` or a bare
version.

`testing` is a moving pointer: it names the last commit someone validated and is
overwritten as often as anyone runs the script. That is what the local test
targets want, since they are checking the tree in front of you.

`<version>-testing` exists for testing that has to remain identifiable after
`testing` has moved on -- a run dispatched to a VM, or two candidates compared
against each other. The suffix is the point: a bare version tag promises a
published release, and someone pulling `0.19.0` has no way to tell it was cut
mid-development. Production work uses a release tag, never either of these.

### What the registry actually holds, as of v0.20.0

The table above says what each tag means. What is published diverges from it in
two places, deliberately, and both are easy to misread:

| repository | `:latest` serves | backed by a release? |
|---|---|---|
| `lraa` | **v0.17.7** | yes -- `origin/main` is v0.17.7, and its wdls hardcode `lraa:latest` in twenty places |
| `lraa-core`, `lraa-sc`, `lraa-orf`, `lraa-combined` | 0.18.3-era digests | **no** |

`lraa:latest` is the one public users resolve, and it is correct: it points at
the image built for the last published GitHub release. Nothing on `main`
references the four split repositories.

Their `:latest` tags therefore assert something no release backs. They are left
in place rather than deleted or moved, because only devel-branch wdls name them
and those runs override the tag anyway through `testing/lraa_test_docker.mk`.
Moving them to a development version would re-assert the same false claim about
a newer commit; deleting them would break any internal caller that has quietly
come to depend on them. Neither is worth doing until the release that closes the
0.17.7-to-devel gap, which is when they should be pointed at that release and
the claim becomes true again.

`lraa:testing` is a third case: it holds a 0.18.3-era digest and **no current
script writes it**. `build_docker.testing.sh` deliberately avoids the plain name,
so nothing will ever move it. Treat it as stale. `lraa-sc:0.18.3-sklearn` is a
hand-made tag outside this taxonomy from the same period.

## Building

Three paths. All three build the commit you are sitting on; they differ only in
the tag they write:

```bash
cd Docker
bash build_docker.testing.sh     # testing,   from git HEAD
bash build_docker.versioned.sh   # <version>, from git HEAD
bash build_docker.latest.sh      # latest,    from git HEAD
```

Use `build_docker.testing.sh` to validate the commit you are on before it is
released -- running the WDLs against `latest` tests the last release, so it
cannot fail on anything you have changed since. Use the versioned script for a
release someone can pin, and the latest script for the release the `.wdl`
defaults resolve to; a release normally runs both.

All three build `lraa-base` first and pass it to the other builds as
`--build-arg LRAA_BASE_IMAGE=`, so the three published images are built against
the base just made rather than a stale local copy. The base is not pushed. The
version and the commit are passed as `--build-arg LRAA_VERSION=` and
`--build-arg LRAA_CO=`, so the SHA lives in one place instead of four
Dockerfiles that can drift apart.

All three ask git for the SHA, and all three refuse to build a commit GitHub
cannot serve: the Dockerfiles fetch the checkout by SHA, so an unpushed commit
fails several layers in with an error naming a tarball rather than the mistake.
The testing script additionally asserts the version the built image reports,
which catches a stale cached layer or a `VERSION.txt` that has drifted from the
`VERSION` in the `LRAA` script.

There is no pinned-SHA file. One existed, `Docker/LRAA_CO.txt`, and it drifted:
it named a 0.18.3-era commit while `VERSION.txt` said 0.19.0, so the release
scripts would have published 0.19.0 images built from code predating the
release. It also shipped inside every image as part of the checkout, where it
read like provenance while describing a different commit than the surrounding
code. A build now describes the tree it was run from, and the SHA it used is
recorded in each image twice: `ENV LRAA_CO`, and the OCI label
`org.opencontainers.image.revision` for reading without running the image.

Only the release scripts retag the core image as plain `lraa` and push that too,
so the alias never drifts from `lraa-core`. `testing` gets no such alias: the
plain name is what a pipeline that never named an image ends up pulling, and a
pre-release commit does not belong there.

A cold run takes about an hour, nearly all of it compiling Seurat from source
for `lraa-sc` and for the combined image. A run that only moves the commit takes
a few minutes, because the checkout is the last layer of every image. That is
what makes a per-commit testing build cheap enough to be routine.

## Releasing a new version

1. Update `VERSION` in the top-level `LRAA` script and `Docker/VERSION.txt` to
   the same value.
2. Commit and push the code.
3. Validate the commit you just pushed before naming it a release: run
   `bash build_docker.testing.sh`, then the WDL test targets under `testing/`,
   which run against `:testing`. A green run against `:latest` would only tell
   you the previous release still works.
4. Run both release build scripts. They build the commit you just pushed and
   fetch it as
   `https://github.com/MethodsDev/LongReadAlignmentAssembler/archive/${LRAA_CO}.tar.gz`
   rather than cloning, so what ships is the tree you validated in step 3 --
   provided you have not committed anything since.
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
