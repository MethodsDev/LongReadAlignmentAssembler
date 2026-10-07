# Open items, holds and caveats

## ON HOLD until the repo owner is consulted: GTF merge pool
Files in play: `util/merge_LRAA_GTFs.py` (`_merge_one_group`, `_iter_group_results`, `--cpu`), `pylib/test_merge_gtfs_pool_order.py`,
the `lraa_merge_gtf_task` block (cpu, preemptible, C3D shape, `--cpu`) and `cpuMergeGTFs` in `LRAA-cell_cluster_guided.wdl` / `LRAA-singlecell.wdl`.
To back out: restore `util/merge_LRAA_GTFs.py` from git (`git checkout -- util/merge_LRAA_GTFs.py`), delete the test file, and revert those WDL hunks
(task back to `cpu 1`, `n2d` machine chosen from `memoryGB`; the previous task block used `memoryGB` only and no `--cpu`).

### The ID issue (why it is on hold)
- Output of the pooled run vs the original serial run (15 GTFs, 48 contig-strand groups): 8 of 1.08 M GTF lines differ. Two isoforms of gene
  `g:chr1:-:comp-2073` swap numbers (`iso-19` <-> `iso-20`); they differ only in which of two alternative 52 bp exons they use. Gene ids, the set of
  ids and the set of structures are identical; the tracking rows follow the new numbering. Tracking row ORDER already differs between any two runs of the
  original code (unordered), GTF is byte-identical between two serial runs of the same group.
- Pool output is identical across worker counts (6 vs 8 workers, byte for byte) and equals a fresh single-group run (`--contig chr1-`), i.e. what a contig-restricted
  or scattered run gives. The serial whole-genome run is the odd one out.
- Mechanism (proven): `Transcript.structural_sort_key` = `(lend, rend, cdna_len, simple_path_str, exons_string)`; `recluster_transcripts_to_genes` numbers isoforms by it.
  The two isoforms tie on the first three fields, so `simple_path_str` decides: node ids joined as TEXT, compared lexicographically, and intron ids come from a
  process-wide counter (`GenomeFeature.Intron.intron_id_counter`, never reset). Pre-setting the counter in a fresh `chr1-` run flips the order: offsets 0 and 99 give one order,
  100000 and 1000000 the other. Hash seed (PYTHONHASHSEED 0/1/2) does not matter. Running `chr1+` first in the same process also flips it.
- Scope: the same function is called by the main `LRAA` discovery per contig-strand in one process, so tied-isoform numbering can depend on the sharding layout
  (off = one process for the genome; by_chromosome = one process per chromosome, + then -; by_chunk = per chunk with per-unit namespacing). The merge itself works
  from structure, not ids, and on this data grouping/collapse decisions were identical. NOT checked: discovery runs across scattering modes; whether the same tie order
  (it also orders vertices for community clustering) can ever change a gene's membership.
- Candidate fix (NOT made; changes numbering of ties in existing outputs): put `exons_string` before `simple_path_str` in `structural_sort_key` (pure structure), or compare the
  simple path numerically. This is in `pylib/Transcript.py`, the repo owner's area.
- Options for the merge task: ship the pool as is (documented tie difference), serial path through the same fork-per-group code so `--cpu 1` equals `--cpu N`, or wait for the sort-key fix.

## Pending decisions / not done
- Image rebuild and release (and the version bump / changelog entry) are the user's call.
- Preemptible: set (`preemptible_tries = 3`) on partition, normalize, merge GTF. Not set on other tasks. The whole-library partition of a very large BAM is the one that could run long;
  lower `preemptible_tries` for it if measured past ~10 min. No retry/`maxRetries` settings were added.
- Normalize stays at 8 CPUs (user). The per-cluster normalizes could use 4 (cheaper by <0.1 c, ~2x slower, off the critical path); the merged one is on the critical path.
- Parallel parse for the GTF merge (43 s serial of 200 s) and parallel parse in the other tasks (`build_cluster_pseudobulk_matrices`, `lraa_merge_gtf_task`) were measured/considered but not done.
- `build_cluster_pseudobulk_matrices` (n2d-standard-2, user said leave as is) and `sc_build_sparse_matrices` (parse of 160 M rows is single-threaded; levels could run concurrently) are further candidates.
- `LRAA_runner_task` was left alone (main process; hardest to optimize).

## Unverified / risks
- Not run on GCP Batch or Cromwell. womtool 92 and miniwdl accept every changed WDL; `Map[String,File]` outputs via `read_map` were exercised only under miniwdl (runtime localization in Cromwell untested).
- Timings are local-NVMe; Terra pd-ssd throughput may bind the partition and normalize (disk size sets throughput).
- "C3D rejects HDD" and whether a C4 `-lssd` shape could work were never verified.
- 2 GiB is treated as 2 GB when picking machine RAM (user agreed); the merge runs 32 GiB requests on a 32 GB machine (peak measured 22.7 GB).
- Spot prices for highcpu shapes are derived, not looked up.
- Quant-only/de novo/ref-guided real-run outputs (local read-only copies) were only read, never modified.

## Idea (not done): have the chunk tool estimate per-chunk / per-shard workload
Chunks are cut at ~equal genomic span (median 10.0 Mb, p10-p90 9.6-10.4 Mb; annotation blocks can stretch one, e.g. 13.3 Mb on chr17), NOT equal work. Real runs (3 single-cell runs, 8,235 chunks):
- median chunk 42 s, p90 131 s, max 4,165 s; the 285 chunks > 300 s (3.5%) hold 30% of all chunk time.
- The slow chunks are the same loci in every shard/cluster (27 each): chr6:29.5-40 Mb (MHC) median 354 s, chr11:60-70 Mb 235 s, chr19:39.8-50.2 Mb and chr19:10-20 Mb ~200 s,
  chr17:70-83 Mb, chr12:50-60 Mb, chr1:150-160 Mb. Worst single chunks: 4,165 s (chr11:60-70, 1.2 M records, 2,621 transcripts), 3,362 s (chr6 MHC), 3,307 s (chr19:39.8-50.2).
- They are dense in records AND annotated transcripts (1,900-2,800 transcripts / 10 Mb); chr17:70-83 Mb has only 438 k records but 2,741 transcripts and took 2,690 s.
- Rough fit, chunk wall ~= 47 s + 2.2 ms x records_total (R2 0.61, MAE 32 s); adding annotated transcript count in the chunk -> R2 0.68 (~0.1 s per transcript). The tail is not explained by either alone.
Why it would help: the chunk tool already knows each chunk's region, record count (split_counts.records_total) and transcript count (gtf_transcripts_emitted) at cut time, so it could emit a predicted
wall/work score per chunk and per shard. Uses: (a) choose `preemptible_tries` per shard (0 for shards holding the known-dense loci, 2 otherwise), (b) pick cores/memory per shard,
(c) cut dense regions finer (equal work instead of equal span), (d) order/schedule the long chunks first. Data: perf/chunk_index.tsv + chunkReports in each real run's perf/ folder.
