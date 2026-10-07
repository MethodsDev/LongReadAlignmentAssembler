# Measurements

All timings from this machine (16 cores, 58 GB RAM, local NVMe) with cores pinned via `--cpuset-cpus`. They do NOT include
GCP persistent-disk throughput: pd-ssd throughput scales with disk size (~0.5 MB/s per GB), so a 150 GB disk could make the
partition I/O-bound on Terra. Time one real Terra run before relying on the partition numbers.

## Real-run reference (task_resources.tsv, miniwdl on a 28-core host, v0.44.0)
Per-cluster `partition_by_chromosome_task` (cpu 7, one call PER CLUSTER, 12-16 clusters): 234-337 s unqueued; initial
whole-library call 336 s (cpu 12). `normalize_cluster_bam` 24-233 s (largest cluster 1.79 GB); `normalize_merged_bam` 544 s (8 cpu);
`merge_bams` 214 s; `merge_cluster_trackings` 5.0/7.9/8.7 min; `lraa_merge_gtf_task` 44.5/46.0 min (de novo/ref-guided), peak 13-14.4 GiB;
`LRAA_runner_task` dominates (27/69/73 h summed over 325/775/675 runs, peak 3.1-5.3 GiB of 16). `build_sc_sparse_matrices` is not in the logs.

## Partition
Cluster BAM 3 (518 MB), all 25 contigs, -@4, output size / wall: default(6) 518 MB 22 s | level 1 620 MB 11 s | level 4 543 MB 14 s | level 6 521 MB 21 s.
Workers x -@ on the 6.6 GB merged proxy (stand-in for the shared splice-graph BAM), 6 cores, level 4: 1x-@3 164 s | 2x-@1 325 s | 3x-@0 171 s | 4x-@0 128 s
(peak total RSS 129-209 MB; a first cold-cache 1x-@3 run read 234 s).
**Shipped defaults:** shared split (6.6 GB BAM + genome.fa 3.1 GB + genes.gtf 1.4 GB) on 16 cores, 14 workers x -@0: **54 s, 974 MB peak**.
Per-cluster call (518 MB cluster BAM) on 4 cores, 2 workers x -@0: **20 s, 114 MB peak**. Before: FASTA+GTF split alone is 26 s per call and
~4.4 GB written per cluster; the shared-BAM split is most of the old 234-337 s.
Single-extraction knee (older measurement, in the WDL): -@2 3.81 s, -@3 2.58, -@4 1.97, -@5 1.63, -@8 1.67 on a 1.30 GiB slice.

## Normalize (6.6 GB merged proxy = `samtools merge` of the 12 quant-only cluster BAMs, 34.5 M records)
| Run | Time |
|---|---|
| image/repo original code, 8 workers | 547 s (real run: 544 s) |
| repo original scripts (like-for-like), 8 workers | 545 s |
| shipped config (intermediate level 1, threaded collapse), 8 workers | **353 s** incl. index; records, header (no @PG) and file size identical to original |
| same, 4 workers pinned to 4 cores | 645 s |
Phases at 8 workers before -> after: strand split 158 -> 100 s, per-contig normalize 245 -> 170 s, merge 59 s, collapse 73 -> 23 s, index 10 s.
Cost (spot, see below): highcpu-8 353 s ~0.75 c, highcpu-4 645 s ~0.68 c; break-even 706 s. User decision: keep 8.
pysam write of one contig: default 5.1 s / 39.7 MB, level 4 3.1 s / 41.9 MB, level 1 2.2 s / 48.5 MB. Collapse scan on 6.6 GB: 125 s (1 thread) vs 47 s (8).
Cluster 3 normalize (518 MB, 8 cores): 42 s (real run 36 s).

## GTF merge (de novo run's 15 per-cluster GTFs; 1.08 M output lines; `--oversimplify chrM`)
Repo original serial: 1,418 s, ~17 GB. Pool: 6 workers 280 s; 8 workers 200 s, PSS peak 22.7 GB; 4 workers 334 s, 17.6 GB. Serial parse of inputs ~43 s.
(Summed RSS overstates because forked workers share pages; use PSS or cgroup.) Image-version script gave a different line count (1.71 M vs 1.08 M)
than the repo code: image and repo differ in version; always benchmark against `git show HEAD:...`.

## merge_trackings / sparse matrices
- 3 x 30 M-line trackings: Python loop 8 m 39 s -> shell/pigz 28 s (16 cores unpinned). Real quant-only 12 inputs (41 M rows), 4 cores: 36 s, output rows identical to the old merge.
  25 inputs pinned to 8 cores: 7 m 04 s. Peak RSS ~0.1 GiB.
- `sc_build_sparse_matrices` on the 41 M-row merged tracking, 3 cores: 9 m 49 s -> 6 m 55 s (stream phase 8 m 44 s -> 5 m 49 s); profile: pandas python-engine parse 72%,
  `_check_comment` a third of the run. Level 1 gzip is fine (matrix 298 MB raw: level 1 88.6 MB 0.95 s, level 6 74.8 MB 3.2 s).

## Spot prices (us-central1, from the user) and derived
c3d-standard-8 $0.09256/h; c3d-highmem-4 $0.06244/h. Derived assuming linear pricing: $0.00753/vCPU-h, $0.00101/GB-h -> highcpu-4 ~$0.0382/h, highcpu-8 ~$0.0764/h (DERIVED, check).
Merge: standard-8 200 s = 0.514 c vs highmem-4 334 s = 0.579 c (standard-8 cheaper and faster).

## End-to-end smoke (miniwdl, overlay images)
`testing/single_cells/sc_full_pipe_scattered` (~5 min) with default scattering and with `scattering_final_quant=by_chromosome` (exercises the shared
splice-graph split): new vs committed WDL+scripts: 25-26 identical outputs; matrices and tracking identical as sorted lines; differences only
rds/pdf timestamps and a `bam_identity` hash in the chunk plan. Per-cluster partition got only `inputBAM`; shards got fasta/gtf/SG by name from the shared task.

## Recipes
- Proxy BAM: `tar -xzf .../<quant-only run>/partitioned_cluster_bams_tar/*.tar.gz --strip-components=1 -C DIR`, then `samtools merge --no-PG -@N merged_all.bam DIR/*.bam`.
- GTF inputs: `tar -xzf .../<de novo run>/prelim_cluster_gtfs/*.tar.gz`; give the GTFs in cluster order (`sort -V`).
- Overlay image for smoke: `FROM lraa-core:latest` + `COPY LRAA, pylib, util` over `/usr/local/src/LRAA` (the top-level `LRAA` script must be copied too);
  `miniwdl run --cfg miniwdl.test.cfg WDL/LRAA-singlecell.wdl docker=... docker_sc=... -i inputs.json -d OUT` (key=value args right after the WDL path).
- Run containers with `--user $(id -u):$(id -g)`, mount real-run dirs `:ro`, scratch under a local temp dir (`mktemp -d -p`). zsh: use `${var}` before `:`, and word-split with `${=var}`.
  Never `pkill -f` a pattern that appears in your own command line.

## Shared split at 4 cores (decision: keep shared_partition_cpu = 16)
6.6 GB proxy BAM + real genome.fa + genes.gtf, pinned to 4 cores, level 4, -@0: 2 workers 272 s | 3 workers 267 s | 4 workers 270 s (CPU-bound, worker count irrelevant)
vs 54 s on 16 cores. Compute cost about equal (~0.29 c vs ~0.23 c), so the small VM only adds minutes. User kept 16 (larger inputs expected later).

## Partition CPU use on the real 8.5 GB input BAM (a merged long-read BAM, 25 contigs, level 4, -@0), mpstat on the pinned cores
16 cores/14 workers 70 s: busy 83% mean, 17% of seconds below 80%; 16 workers 69 s (no gain); 8 workers on 16 cores 89 s (53% busy);
8 cores/6 workers 118 s; 4 cores/2 workers 335 s at only 52% busy (2 of 4 cores) -- the earlier "CPU-bound at 4 cores" inference was NOT supported.
Wall at 16 cores is set by the longest single contig (chr6 70 s, chr1 68 s, chr19 67 s at ~3.3 M alignments/min): more cores/workers cannot beat that floor.
Input: a merged BAM (+.bai) from a Terra workspace, copied locally.

## Incorporate_gene_symbols
gffcompare v0.12.6 cannot read .gz GTFs (parse error) and has no threading option; its real-data run is ~21 s. incorporate_gene_symbols_in_sc_features.py
takes 9 s on the real run; gunzip of the 1.4 GB reference ~5 s. Change: plain inputs are symlinked (not copied), the two gunzips run concurrently; outputs identical.

## Cromwell 92 (Local backend, Docker) end-to-end on testing/single_cells/sc_full_pipe_scattered fixture, images lraa-core/lraa-sc :cg-terra-testing (commit 5870393 + partition read_map fix)
cluster-guided scattered: Succeeded, 639 s, 59 workflow outputs (empty ones = features this fixture disables). basic + quant_only (same fixture): Succeeded, 163 s, 13 non-empty outputs
(no gtf outputs, as expected for quant-only). Only ERROR lines: Cromwell's cost estimate not knowing C3D machine types. Not covered locally: GCP Batch machine types, preemption, real data size.
