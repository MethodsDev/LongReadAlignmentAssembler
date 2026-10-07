# Terra / C3D branch: working notes

Status as of 2026-10-07. Everything described here is **uncommitted** in the working tree of
`devel` (nothing staged by the assistant; the user has staged part of
`WDL/LRAA-cell_cluster_guided.wdl` while reviewing, so use `git diff`, not `--cached`).
It is intended for a Terra-specific branch; the repo owner may pull pieces into `devel` later.

| File | What it holds |
|---|---|
| [changes.md](changes.md) | Every change, by area, with the new/renamed WDL inputs and the image requirement |
| [measurements.md](measurements.md) | All benchmarks: numbers, the data they came from, how to re-run them |
| [open-items.md](open-items.md) | Decisions on hold, known caveats, things not verified, ideas not done |

## Read this first
1. **Images must be rebuilt before any of this runs on Terra.** The WDLs now pass flags the
   released `lraa-core`/`lraa-sc` images do not have (`--bam-compression-level`,
   `--intermediate_compression_level`, `--cpu` for `merge_LRAA_GTFs.py`, and the changed
   sparse-matrix converter). Nothing was pushed or built into the registry by the assistant.
2. **ON HOLD (until the repo owner is consulted): the GTF merge change.** `util/merge_LRAA_GTFs.py --cpu`
   (process pool), its WDL wiring in `lraa_merge_gtf_task`, and the tests for it. It is in the working
   tree but must not be treated as settled. Reason: isoform numbering of exact ties is
   process-state dependent (see open-items.md). How to back it out is listed there.
3. Decisions already made by the user: normalize stays at **8 CPUs** (`c3d-highcpu-8`);
   N4/N4D/C4/C4D are not usable (Hyperdisk only, Cromwell `disks:` takes LOCAL/SSD/HDD);
   machine families are C3D (main) and N2D (tiny tasks).
4. Each Terra task runs on its own VM: the old comments about "draining the box" / local-host
   core budgets are not the constraint there. Use the cores the machine has.

## Tools and paths
- womtool: `java -jar ~/tools/womtool-92.jar validate <wdl>` (matches Cromwell 92; also run `miniwdl check`).
- Real runs: local read-only copies of a de novo, a quant-only and a reference-guided single-cell run (not in the repo)
  (read-only; test in scratch copies).
- Reference: a local copy of the 10x GRCh38 reference (`fasta/genome.fa`, `genes/genes.gtf`).
- Images: `us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:latest`, `lraa-sc:latest` (older commit than the repo).
