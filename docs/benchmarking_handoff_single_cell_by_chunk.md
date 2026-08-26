# Handoff: single-cell `by_chunk` benchmarking on a 28-core / 187 GB host

Written for an agent picking this up on a larger machine. Everything below is
current as of `devel` @ `51d37647`, with `lraa-core`/`lraa-sc`/`lraa-orf` pushed
at that same SHA.

## What you are benchmarking, and why it is new

LRAA's chunked mode (`scattering = by_chunk`) partitions each contig at
low-coverage positions between annotated loci, processes chunks independently and
merges. It is now the default for bulk, and the single-cell workflows can select
it per phase. Three things landed recently that change single-cell behaviour and
have NOT been characterised at scale:

1. **One run-wide chunk plan.** Cut selection happens ONCE per run, before phase
   1, and the same geometry serves the init pass, per-cluster discovery and the
   final quant. Previously each phase selected its own; a ref-guided smoke run
   measured 76 distinct geometries across 270 extractions in the discovery phase
   alone. Cross-cluster comparability depends on this.
2. **A shared splice-graph BAM, sliced per chunk.** Every cluster builds its graph
   from ONE depth-normalized merged BAM, sliced to each chunk with geometry
   identical to that chunk's own reads. Verified: 270 final-quant extractions
   across 27 clusters yielded only 10 distinct `(chrom, lend, rend)` triples.
3. **Cluster-local pass-1 priors** (`--bam_for_priors`). Pass 1 estimates the
   theta that pass 2 uses to apportion ambiguous reads. It used to read the
   SHARED graph BAM, so every cluster's apportionment was coupled to every other
   cluster's expression. It now reads that cluster's own normalized reads.
   MEASURED, holding `--bam` and `--bam_for_sg` fixed on one cluster of a chr19
   fixture: own-priors versus pooled differs on 152 of 654 transcripts in read
   assignment and 215 in TPM. Not cosmetic.

So the three roles per cluster are three different files:

    --bam              this cluster's full reads          (pass 2 assignment)
    --bam_for_sg       the shared merged normalized bam   (graph only)
    --bam_for_priors   this cluster's normalized reads    (pass-1 theta)

## Data

Use `pbmc_pbio`. The Terra sample table on the dev workstation is
`/home/unix/bhaas/projects/LRAA_PAPER_Analyses/lraa_bmark_terra.sample_table.tsv`;
the row is `PBMCs_pbio`. Its GCS paths, which you will need to stage:

    inputBAM    gs://fc-6a6064fc-0245-4f36-9fd9-c674ef382a6b/submissions/2da07230-68d7-484f-9247-621446f271ff/Minimap2_LR/3af08a9f-a2b0-4d59-a50e-e5ff7ad1053b/call-minimap2_ubam/attempt-4/PBMCs_pbio.aligned.sorted.bam
    genome      gs://mdl-data/Benchmarking/IsoIDandQuant/Cell_lines_SIRVs/GRCh38_no_alt.fa
    annotation  gs://mdl-refs/GRCh38/GTEx_Resources_GENCODE47/human.refdata-gex-GRCh38-GENCODE47.gtf
    clusters    gs://fc-6a6064fc-0245-4f36-9fd9-c674ef382a6b/PMBCs_for_LRAA_paper/PBMCs_LRAA_init_denovo.genes-cell_cluster_assignments.v0.4.0.tsv

Sizes on the dev host: BAM 7.8 GB, genome 3.1 GB, annotation ~135 MB (a
GENCODE v39 variant), cluster assignments 188 KB holding 9,936 barcodes across 14
clusters. A different assignments file will give a different cluster count, and
cluster count drives the fan-out, so record which one you used.

Barcodes are bare (no `-1` suffix) and match the BAM's `CB` tags directly; 9,933
of 9,936 appear in a chr1 sample alone. `XM` carries the UMI. The BAM has no `XW`,
which is correct -- it is raw, and normalization adds the weights.

**Restrict to a few chromosomes** for anything exploratory. `chr21 chr22 chrM` is
a good shape: two ordinary chromosomes plus one that is oversimplified by default
(`--oversimplify chrM`), and it is what a successful Terra submission used. On the
dev host that yields 11 chunks; the whole genome yields 475.

## Invocation

The single-cell entrypoint is `WDL/LRAA-singlecell.wdl`. A working restricted
inputs file is at
`testing/single_cells/sc_full_pipe/pbmc-chr21_22_M.by_chunk.QUANTONLY.inputs.json`
-- copy it and repoint the four paths. Run it the way the Makefile targets do:

```bash
cd testing/single_cells/sc_full_pipe
_cc=$(../../miniwdl_call_cache_dir.sh \
    us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:testing \
    us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-sc:testing \
    us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-orf:testing)
MINIWDL__CALL_CACHE__DIR="$_cc" miniwdl run --cfg ../../../miniwdl.test.cfg \
    ../../../WDL/LRAA-singlecell.wdl -i YOUR.inputs.json -d OUTDIR \
    docker=...lraa-core:0.28.0-51d3764 docker_sc=...lraa-sc:0.28.0-51d3764 --verbose
```

Bind the COMMIT-PINNED image tags, not `:testing` -- `:testing` moves on every
rebuild, so a run bound to it is not reproducible. Tags exist as
`0.28.0-<short-sha>` for every build.

## Modes worth exploring

Three annotation modes, each with all three scatter phases set to `by_chunk`.
Committed smoke tests exist for all three on a 2 Mb chr19 fixture (`make
test_sc_cg_by_chunk_{qonly,refguided,denovo}`), so the shapes are known to work;
what is unknown is behaviour at real depth and cluster count.

| mode | `quant_only` | `initial_annot_gtf` |
|---|---|---|
| quant-only | `true` | supplied |
| ref-guided | `false` | supplied |
| de novo | `false` | OMITTED (it is `File?`) |

The scatter knobs are independent, so the interesting axis beyond annotation mode
is which phases chunk:

    scattering_init          default by_chunk
    scattering_per_cluster   default by_chromosome
    scattering_final_quant   default by_chromosome

Set all three to `by_chunk` for the new architecture. Comparing
`by_chunk` against `by_chromosome` per phase is the wall-time and peak-memory
question nobody has answered on a real library.

**Chunk width matters.** `approx_MB_per_cut` defaults to 10 Mb, and a contig
shorter than it yields no cuts. The `_init` variants
(`approx_MB_per_cut_init`, `approx_MB_per_cut_wiggle_window_init`) set the init
pass separately from the rest. On small fixtures 0.2/0.05 is needed to get any
fan-out at all; on a real chromosome the defaults are appropriate.

## Resources: what is measured, what is a guess

Six chunk resource knobs are now reachable from the single-cell entrypoints
(they were inert before -- they stopped at `LRAA.wdl` and nothing a caller passed
had any effect):

    chunkMakeChunksCpu   chunkMakeChunksMemoryGB
    chunkCpu             chunkMemoryGB
    chunkMergeCpu        chunkMergeMemoryGB

Current defaults and their basis, from `subwdls/LRAA_chunk_scatter.wdl`:

| task | default | basis |
|---|---|---|
| `make_chunks` | 16 cpu, 32 GiB | **guess.** Only per-subprocess numbers exist (80 MiB for a cut selection, 90 MiB for an extraction, implying ~1.4 GiB at 16 workers) and that metric is a LOWER BOUND -- the same interval sampler understated a leaf's true peak by 4x. 8 GiB is the documented candidate; do not apply it without a task-level number. |
| `process_chunk` (leaf) | 2 cpu, 16 GiB, preemptible 3 | measured max 3,003 MiB over 150 sampled leaves of a whole-genome run (median 1,037, p95 2,026). NOT lowered, because a width-derived value would under-request the pathological chunk: `approx_MB_per_cut` is a TARGET, not a bound, and an uncut contig can be far wider. A whole chr1 unit measures 5,594.7 MiB. |
| `merge_chunks` | 2 cpu, 8 GiB | **measured.** 371 MiB over 556 real units (9.8M tracking rows). Generalizes because `tracking_merge_peak_resident_rows` landed exactly on the external sort's 500,000 cap, so the peak is bounded by the cap rather than by input size. 8 rather than the implied 4 per a project rule: safe minimum 8 GiB wherever a measurement suggests 4 or less. |

**The leaf number is the open one and your host is well suited to settling it.**
The correct fix is for `make_chunks` to emit a per-chunk memory value from each
chunk's ACTUAL span -- it already computes the spans, and
`pylib/ChunkedRun.py:chunk_unit_peak_mib` (`1024 + 22 MiB per genomic Mb`, per
UNIT; a strandless leaf runs two) already exists -- as a file per chunk globbed
in the same order as `chunkTars` so the scatter can zip them.

### Telemetry to collect

Every chunk task now emits its own cost as a delocalized output, so a run answers
the sizing question with a table rather than log archaeology:

    chunkResources          one row per leaf: container_peak_rss_bytes, cpu, memory requested
    unitResourceSummaries   LRAA's per-unit interval series (which stage spent it)
    makeChunksResources     the prep task's peak, with chunk count
    makeChunksTiming        its timing.json: per-unit RSS and wall series
    mergeResources          merge peak, unit count, and peak_resident_rows

`container_peak_rss_bytes` comes from `/sys/fs/cgroup/memory.peak`, the kernel's
high-water mark, so unlike interval sampling it cannot miss a spike. Prefer it
over `unitResourceSummaries` for sizing; use the series to see WHICH stage spent
the peak.

Note miniwdl treats `runtime.memory` as a scheduling RESERVATION and does not
enforce it (`memory_limit_multiplier` defaults to 0), so on 187 GB you will get
real concurrency rather than the serialization a smaller box suffers -- but a task
cannot be OOM-killed for exceeding its request either, so peaks are honest.
**Grep for `exit_code: 137` and `OOMKilled` explicitly**; a `make` exit of 0 does
not rule them out.

## Validation standard

Both validators, on every WDL you touch and every file that imports it (imports
are relative):

```bash
miniwdl check WDL/LRAA.wdl
java -jar ~/BIN/womtool-91.jar validate WDL/LRAA.wdl
python3 -m pytest pylib/test_wdl_runtime_completeness.py -q
```

Rationale in `WDL/README.validation.md`. Short version: miniwdl is what the local
harness accepts, womtool is what Terra accepts, and NEITHER catches a task with no
`docker` in its runtime block -- an absent runtime attribute is legal WDL, but GCP
Batch refuses the job. That shipped once and cost a Terra submission plus several
wrong hypotheses.

## Traps that have already cost time here

- **Image/worktree revision skew.** `make test_wdls` refuses to run when the image
  revision and git HEAD disagree (`check_test_image_revision`). That guard is
  correct; committing WHILE a suite runs is what trips it. Do not pass
  `LRAA_TEST_SKIP_IMAGE_CHECK=1` to get past it -- a mixed-version run cost two
  hours of misleading results. Rebuild instead: `Docker/build_docker.testing.sh`,
  which also refuses to build from an unpushed commit.
- **`testing/**` is gitignored** with explicit exceptions. A new inputs JSON needs
  `git add -f`, or the tracked Makefile target referencing it is dead in a fresh
  clone.
- **`"${}"` on a non-optional input.** A Terra inputs JSON binding a non-optional
  `Int` to an empty `${}` is not the same as omitting it. Five `memoryGB*` fields
  were bound that way in a real submission.
- **Do not compare `quant.expr` files byte-wise.** They carry a `# LRAA CMD:`
  header echoing the invocation, so any flag change differs unconditionally.
  Parse and compare values.

## Open questions worth your cycles, in priority order

1. **Leaf memory from actual span**, per above. Needs the per-chunk emission plus
   a run wide enough to include a pathological chunk.
2. **`make_chunks` task-level peak** on a whole-genome run, to settle 32 vs 8 GiB.
   `makeChunksResources` now reports it; nobody has collected one.
3. **`by_chunk` versus `by_chromosome` per phase** on a real library: wall time,
   peak memory, and whether the two non-preemptible tasks `by_chunk` adds per
   `LRAA_wf` invocation cost more than the leaf parallelism buys for a small
   per-cluster BAM. That tradeoff is why per-cluster phases still default to
   `by_chromosome`.
4. **Whole-genome `by_chunk` end to end.** Largest verified is 7 SIRV contigs, a
   3.4 Mb minigenome, and a 305-chunk init pass that was stopped before its merge.
   The ~300-chunk merge on a full library is unmeasured.
5. **Cluster-locality at scale.** The fixture proves cluster A's numbers move when
   its priors change. Worth confirming on a real library that per-cluster
   `reads_total` values differ and sum to approximately the library size --
   before the fix all 32 clusters reported an identical 94,908.

## Reproducing the state you are starting from

    devel @ 51d37647
    pytest pylib                     1245 passed
    make test_wdls                   passed at 254b3d91 (before later WDL-only commits)
    miniwdl + womtool                clean on all seven workflows at 51d37647
    images                           lraa-core/sc/orf at 0.28.0-51d3764, pushed
    smoke tests                      test_sc_cg_by_chunk_{qonly,refguided,denovo} all green
