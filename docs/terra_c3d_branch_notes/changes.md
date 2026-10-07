# Changes (uncommitted, working tree of `devel`)

## 1. Machine types (first rounds)
- All tasks reachable from `LRAA-singlecell.wdl` use `predefinedMachineType` (Cromwell GCP Batch rejects it
  together with `cpu`/`memory`, so those lines were removed from the converted `runtime` blocks).
  Files: `LRAA.wdl`, `LRAA-cell_cluster_guided.wdl`, `LRAA_quant_by_cluster.wdl`, and `subwdls/` Incorporate_gene_symbols,
  LRAA-build_sparse_matrices_from_tracking, LRAA-filter_good_cells, LRAA-gene_sparseM_to_seurat_clusters,
  LRAA_chunk_scatter, LRAA_runner, Normalize_bam, Partition_data_by_chromosome, partition_bam_by_cell_cluster.
- Sizing logic is copied from another scattered C3D workflow of the team's: cpu tier and memory tier
  (4/8/16/30/60/90/180/360) round up to the first C3D size at or above the old request; highcpu/standard/highmem from memory.
- Disks that were `HDD` are now `SSD` (the example uses SSD; "C3D rejects HDD" was assumed, not verified).
- Fixed tiny tasks (1 cpu, 1-4 GiB): `n2d-highcpu-2` (<=2 GiB) or `n2d-standard-2` (`mergeQuantResults`, 4 GiB).
  `build_cluster_pseudobulk_matrices`: `n2d-standard-2` (previous size 1 cpu / 8 GiB).
- Threading/`--cpu_budget` options now use `c3d_effective_cpu` (the machine's cores) in: `LRAA_runner_task`
  (`--cpu_budget`, samtools -@), `make_chunks`, `process_chunk`, `emit_shared_chunk_plan`, `merge_bams`,
  `normalize_bam_by_strand`, `partition_bam_by_cell_cluster`. Deliberately not changed: tasks with no threaded command.
- `count_bam`: static `c3d-highcpu-4`, `samtools view -@ 3` (measured pinned to 4 cores, 1.5 GB BAM: -@2 3.79 s, -@3 2.61, -@4 2.82, -@5 2.59).
  The `countBamThreads` workflow input and the task's `samtools_threads` input were removed (callers binding `countBamThreads` must drop it).
- `mergeQuantResults` (LRAA.wdl): the single-stream python gzip tracking merge is replaced by the same parallel
  decompress|filter|pigz-per-shard pipeline as `LRAA_merge_trackings` (schema check kept); fixed `c3d-highcpu-4`
  (was n2d-standard-2). Measured on 25 proxy shards cut from the real quant-only tracking (51.8 M rows), 4 pinned
  cores: 293 s -> 43 s, decompressed output byte-identical; also identical on empty-first shard, plain-text shard,
  single input, and the schema-mismatch error still fires. Real runs: up to 379-514 s for the whole-sample merge.
- `validate_scattering` (cpu 1, 1 GiB, 10 HDD), `gather_shard_cut_plans` (cpu 1, 2 GiB), `derive_contigs` (cpu 1, 2 GiB, 20 HDD) and `mergeReadAssignmentSummaries` (cpu 1, 2 GiB): left on their ORIGINAL runtimes -- near-zero work (user decision).
- `LRAA_tar_outputs`: tars in place (`tar -zcvhf ... --transform 's,^.*/,<dir>/,' --files-from=<list>`), no copy; disk 5x -> 2x.
- `LRAA_merge_trackings`: parallel decompress|filter|pigz pipeline per input, parts concatenated in input order
  (valid multi-member gzip). New input `cpu = 4`; `memoryGB = 4` -> `c3d-highcpu-4`. Verified identical to the
  old Python merge on the real quant-only run (inputs must be in the same order, which was lexicographic `0 1 10 11 2 ...`).
- (SUPERSEDED by the devel rebase, upstream's own rewrite is used as is) `util/sc/singlecell_tracking_to_sparse_matrix.py`: `#` lines filtered by a `zcat -f | grep -av '^#'` subprocess
  instead of `read_csv(comment="#")`, and `drop_duplicates()` before the mapping loop; raises if the subprocess fails.
  Identical outputs; 1.5x faster stream phase. Behaviour note: pandas' `comment` also cut at a `#` mid-line; the new
  filter only drops lines that START with `#`.

## 2. Partition by chromosome (this round)
Files: `util/partition_data_by_chromosome.py`, `WDL/subwdls/Partition_data_by_chromosome.wdl`, `LRAA.wdl`,
`LRAA-cell_cluster_guided.wdl`, `LRAA_quant_by_cluster.wdl`, `LRAA-singlecell.wdl`.
- Script: `--bam-compression-level N` (pysam/samtools `--output-fmt-option level=N`); units carry the level as 7th tuple element.
- Task defaults now: `cpu = 4` (C3D tier), `samtools_threads = 1` (one core per extraction; `-@` = threads-1, clamped so one
  worker + the 2 FASTA/GTF jobs fit the machine), `Int? partition_workers` (unset = `(cores-2)/(extra+1)`, min 1),
  `bam_compression_level = 4`, `memoryGB = 4`, `preemptible_tries = 3`. `--reserved-cpu` is now the machine's cores
  (`c3d_effective_cpu`), never the requested number. The tier takes cpu and memory only (no cycle with workers).
- New outputs: `chromosomeFASTAsByName`, `chromosomeGTFsByName`, `chromosomeBAMsForSGByName` (`Map[String,File]` from
  `read_map` of TSVs written in the command).
- **Shared inputs are split once**, by the parent, and handed to each cluster by contig NAME:
  - `LRAA.wdl` new optional inputs `internal_presplit_fastas`, `internal_presplit_gtfs`, `internal_presplit_sg_bams`
    (`Map[String,File]?`), `Int partition_cpu = 4`, `Int? partition_workers` (was `Int = 1`). When a map is set, the matching
    input is not passed to `splitByChr` at all (conditional blocks leave the `File?` undefined); shards read
    `select_first([map])[contig_name]`.
  - `LRAA-cell_cluster_guided.wdl`: parent call `split_discovery_shared_inputs` (fasta + reference annot gtf) for the
    by_chromosome discovery phase, preceded by `LRAA.derive_contigs` when `main_chromosomes == ""`.
  - `LRAA_quant_by_cluster.wdl`: parent call `split_shared_inputs` (merged normalized splice-graph BAM + fasta + consolidated
    gtf) when `scattering == "by_chromosome"`.
  - New/changed workflow inputs: `initial_partition_cpu = 16`, `Int? initial_partition_workers` (was `Int = 2`),
    `cluster_partition_cpu = 4`, `Int? cluster_partition_workers` (was `Int = 1`), `shared_partition_cpu = 16`
    (in `LRAA-singlecell.wdl`, `LRAA-cell_cluster_guided.wdl`, `LRAA_quant_by_cluster.wdl`). No Terra config in the repo bound the old ones.

## 3. Normalize (this round)
Files: `util/normalize_bam_by_strand.py`, `util/separate_bam_by_strand.py`, `WDL/subwdls/Normalize_bam.wdl`.
- Finding: the 34 `@PG` records are INHERITED from the input BAMs (minimap2 + upstream shards), `--no-PG` is already at every
  write site, so a flag cannot remove them; the collapse (`util/misc/collapse_bam_pg_header.py`) is the right tool.
- The collapse call now gets `--threads` (its full-BAM `PG:Z:` scan was single-threaded: 125 s -> 47 s on 6.6 GB; byte-identical).
- `--intermediate_compression_level N` on both scripts (pysam `format_options=[b"level=N"]`): strand-split BAMs and per-contig
  parts only. Never applied to the final output: disabled for `--input_is_single_strand` and for whole-file units
  (`scope is None`). WDL input `intermediate_compression_level = 1`, `preemptible_tries = 3`, `cpu = 8` (unchanged; user keeps 8).
- A level-4 final merge was tried and dropped (saved ~5 s, +4.9% file size).

## 4. GTF merge (this round) - ON HOLD, see open-items.md
- `util/merge_LRAA_GTFs.py`: loop body moved into `_merge_one_group()`; new `--cpu` (default 1 = serial); forked pool, one
  fresh process per (contig, strand) group (`maxtasksperchild=1`), largest first, results written in original order.
- `lraa_merge_gtf_task` (`LRAA-cell_cluster_guided.wdl`): `cpu = 8`, `preemptible_tries = 3`, C3D shape
  (`c3d-standard-8` at 8 cores / 32 GB), passes `--cpu ~{c3d_effective_cpu}`; workflow inputs `cpuMergeGTFs = 8`
  (`LRAA-cell_cluster_guided.wdl`, `LRAA-singlecell.wdl`). `memoryGBmergeGTFs` unchanged (32).

## 5. Tests
- Updated: `pylib/test_partition_contig_fanout.py` (new invariants: workers fill the machine, wiring of cpu/shared splits,
  default fits machine, compression level reaches samtools).
- New: `pylib/test_merge_gtfs_pool_order.py` (ordering/plumbing with a stub group; hold-related), `pylib/test_intermediate_compression.py`.
- Environment notes: images lack pytest (`pip install --target DIR pytest`, mount, `PYTHONPATH`). Partition tests need >= 7
  visible cores and a writable cwd. Baseline on the committed tree: 227 passed / 3 skipped (8 cores); working tree: 234 / 3
  (before the later test additions: 111 passed on the targeted files). Full `pylib` run: identical 13 failures/errors on both trees
  (missing deps in `lraa-core`); the sc-related tests pass 35/35 in `lraa-sc`.

## 6. Not done on purpose
- No `CHANGELOG.txt` entry and no version bump (release decision; v0.30.0 set the precedent of bumping when WDLs need new script contracts).
- No git commits/branch.

## Incorporate_gene_symbols.wdl (run_gffcompare, incorporate_gene_symbols_sc)
- Machine: 2-vCPU N2D picked from `memoryGB` (<=2 highcpu-2, <=8 standard-2 [gffcompare default 8], <=16 highmem-2 [integrate default 16], then highmem-4/8). Not C3D: neither task threads
  (gffcompare v0.12.6 has no thread option; real run 21-33 s, peak 0.87 GiB). No disk type change.
- GTF prep: `.gz` inputs are decompressed (gffcompare cannot read gz), both at the same time; plain inputs are `ln -s "$1"` (path exactly as given, no realpath, so it follows
  whatever symlink/bind the backend set up -- same reasoning as the `tar -h` comment in LRAA_tar_outputs) instead of `cp`. Outputs verified identical on the real run (plain vs gz).
- Not done: the other 4 `cp` lines (id mappings + 3 sparse tarballs); python-side .gz reading (single reader, saves seconds at most).
- Not verified: an actual Apptainer run.

## LRAA-filter_good_cells.wdl / LRAA_chunk_scatter.wdl merge_chunks
- filter_good_cells: the 3 input tarballs are extracted concurrently; the 3 output tarballs are packed concurrently with `tar -cf - | pigz -1 -p 2` (inner files are already gzipped;
  level 1 vs 6 = +0.02% size; serial gzip 11 s -> ~1 s); barcodes zcat overlaps; summary `cp` -> `mv` after the gene tar is packed (local file, not a localized input, so mv is backend-safe).
  Verified on the real initial matrices (gene/isoform/splice): decompressed contents of every file in all 3 tarballs, barcodes and summary identical to the old command (the .gz bytes differ only by gzip headers).
  Real runs of this task: 164-238 s, peak 2.5-3.3 GiB of 32 GB requested (machine is c3d-highmem-4; sizing left to the user).
- merge_chunks: `cp` into staged/ -> `ln -s "$f"` (path as given, like make_chunks). Verified byte-identical merge outputs on the pylib fixtures with symlinked staging
  (throwaway test, not committed). disk factor (3x inputs) not lowered.
- merge_chunks machine: C3D -> N2D (`n2d-standard-2` for mergeMemoryGB<=8, highmem-2/4/8 above). The merge is single-threaded Python; mergeCpu is not used by the command (only logged). Real in-process stage-6
  merges (same code) took 2.4-62 s per shard; no real timing exists for the WDL merge_chunks task itself (the real runs never used it) -- revisit with a real run.

## Preemptible round (user decision)
- New task input `Int preemptible_tries = 3` + `preemptible: preemptible_tries` on every small task reachable from LRAA-singlecell.wdl that lacked it: validate_scattering, count_bam, derive_contigs,
  gather_shard_cut_plans, mergeQuantResults, mergeReadAssignmentSummaries, merge_GTFs, LRAA_merge_trackings, LRAA_tar_outputs, require_annot_gtf, build_cluster_pseudobulk_matrices, sc_build_sparse_matrices,
  merge_bams, collate_read_assignment_summaries, validate_pre_normalized_inputs, emit_shared_chunk_plan (was explicit 0), run_gffcompare, incorporate_gene_symbols_sc, run_filter_good_cells,
  run_seurat_from_gene_sparseM, sc_build_sparse_matrices_from_tracking, merge_chunks (was explicit 0), partition_bam_by_cell_cluster. Defaults only; no workflow-level inputs added.
- LRAA_runner_task: `preemptible_tries = 2` (85% of logged compute; 95% of calls < 13.6 min; max 77 min). Old "absent attribute, unverified default" comment replaced.
- Deliberately NOT changed: make_chunks (explicit `preemptible: 0`; holds the whole input, runtime unmeasured), FSM/ORF/sqanti/saturation WDLs (not reachable from singlecell).
- Runner runtime vs size (no doc exists): per-chunk wall ~= 47 s + 2.2 ms x records_total in the chunk (R2 0.61, MAE 32 s, 8,235 real chunks); median chunk 42 s, p90 131 s, max 4,165 s -- the tail is not explained by record count.

## Rebase onto devel (51be586, v0.44.x) -- conflict resolutions
- util/sc/singlecell_tracking_to_sparse_matrix.py: took upstream's file whole (their "direct" parser + drop_duplicates replace my grep-subprocess/drop_duplicates change).
- util/separate_bam_by_strand.py: upstream's stamped-strand header writer + my `**_writer_options()` (intermediate compression level).
- util/merge_LRAA_GTFs.py: kept upstream's boundary-support filters, `--ignore_TSS_POLYA`, and the new 3'-end annotation (annotate_polyA_signal / filter_internally_primed_transcripts), the latter moved
  INTO `_merge_one_group` (after recluster) so the `--cpu` pool path gets it. Checked on the real prelim cluster GTFs, --contig chr21: --cpu 1 and --cpu 4 byte-identical.
- LRAA.wdl: kept upstream's build_sc_sparse_shards / docker_sc / memoryGBscShardSparse inputs, dropped `countBamThreads`; shard_contig_length_bp uses shard_fasta (presplit maps); merge_GTFs is
  `n2d-highmem-2` because upstream deliberately reserves 16 GiB there (an OOM is unrecoverable on Terra).
- LRAA-cell_cluster_guided.wdl: kept both import sets and both call-input sets.
- LRAA-build_sparse_matrices_from_tracking.wdl: first task stays on my C3D block (c3d_cpu 3 -> 2 to follow upstream's cpu: 2; still c3d-highmem-8 at the 64 GB default). Upstream's NEW tasks
  sc_build_shard_sparse (1,075 real calls, median 32 s) and merge_sc_shard_sparse (134-247 s) converted to 2-vCPU N2D picked from memoryGB (8 -> standard-2, 16 -> highmem-2), SSD, preemptible_tries 3;
  the unused `cpu` input of sc_build_shard_sparse was removed.
- Re-run after resolving: womtool on every WDL ok; pylib test_partition_contig_fanout / test_intermediate_compression / test_merge_gtfs_pool_order: 43 passed, 2 skipped.
