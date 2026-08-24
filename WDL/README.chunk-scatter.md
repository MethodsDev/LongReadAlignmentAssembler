# LRAA chunk-scatter workflow

`LRAA-chunk-scatter.wdl` (`LRAA_chunk_scatter_wf`) runs LRAA in three phases,
scattering per CHUNK rather than per chromosome (which is what `LRAA.wdl` does;
that workflow is unchanged and remains the supported chromosome-scatter path).

## The three phases

1. **`make_chunks`** — ONE task running genome-wide cut selection plus chunk
   extraction (`ChunkedRun.py --stop_after_make_chunks`). This is not a serial
   phase: it is a single concurrent pool in which a finished cut selection
   submits its own extractions while other selections still run, so one task
   with `--cpu_budget ~{makeChunksCpu}` already overlaps selection and
   extraction across all contigs. It also resolves the TPM denominator with the
   same `-F 0x904` policy as LRAA's own read count, writes
   `work/chunk_plan.json`, and tars each chunk directory into one archive so
   the scatter carries exactly one file per shard.
2. **`process_chunk`** — one task per chunk (`ChunkedRun.py --only_chunk
   <chunk_id>`): stages 3b (orientation split), 4 (normalization), and 5 (the
   per-unit LRAA run) for exactly one chunk. Each leaf sees only its own mini
   FASTA/BAM/GTF; mode (`discovery`) and the TPM denominator come from
   `chunk_plan.json`, so a leaf cannot disagree with the prep run. Leaf outputs
   are renamed to unit-scoped basenames (`<unit_id>.quant.expr`, ...) so the
   merge can gather every shard's files into one directory without collisions.
3. **`merge_chunks`** — ONE stage-6 merge (`merge_chunk_outputs.py`) over every
   unit of every chunk. The manifest reproduces `ChunkedRun.ordered_units`
   order exactly (`+` before `-`, then extraction order), which is load-bearing
   for the merged tables. It also refuses duplicate `unit_id`s, which is what a
   retried preemptible shard gathered twice would otherwise produce — silent
   double-counting.

## Resource defaults and why

| task | cpu | memory | preemptible | rationale |
|---|---|---|---|---|
| `make_chunks` | 16 | 32 GiB | 0 | The genome-wide selection/extraction pool; cpu IS its `--cpu_budget`. Non-preemptible because it holds the whole input and everything downstream depends on it. |
| `process_chunk` | 2 | 16 GiB | 3 | 16 GiB is provisioned headroom, not a measured peak: ChunkedRun's own scheduler estimates a chunk at a fixed base plus 22 MiB per genomic Mb, a deliberate upper envelope. If a leaf is OOM-killed, raise `chunkMemoryGB` first — do NOT lower `approx_MB_per_cut`, because chunk geometry changes the merge's coordinate offsets and invalidates any equivalence comparison in flight. |
| `merge_chunks` | 2 | 32 GiB | 0 | Sized ABOVE the leaves because stage 6 is where a whole-genome run's peak RSS lives: an external merge sort over every unit's tracking rows. Non-preemptible for the same reason as `make_chunks`: it is a single task holding EVERY shard's output, so a preemption discards the whole gather and re-runs the most expensive non-parallel step in the workflow. |

## `--no_reuse_source_bam` is mandatory for this shape

`make_chunks` always passes it. Without it, a strandless chunk spanning its
whole contig skips writing a mini BAM and its manifest names the SOURCE BAM —
which a leaf task on another machine cannot open. The cost is a mini-BAM copy
for every whole-contig chunk (chrM plus the ~90 read-empty references of a
whole-genome FASTA); the benefit is that every chunk directory is
self-contained and shippable.

## Notes

- Strandless-only: the workflow exposes no `chunk_by_strand` knob. Strand-first
  chunking requires a whole-BAM orientation split (stage 1) that cannot be
  scattered, and strandless is LRAA's default and the faster ordering.
- `quant_only = true` requires `annot_gtf`; `ChunkedRun` refuses the
  combination up front, so there is no WDL-side guard.
- `num_total_reads` may be supplied to skip the counting pass; otherwise it is
  computed inside `make_chunks` and recorded both in `chunk_plan.json` (which
  the leaves read) and as the `numTotalReads` workflow output (provenance).
- Terra config: `TerraWorkflowConfigs/LRAA.chunk_scatter.config.json`. There is
  deliberately no `main_chromosomes` input — the partition is derived from the
  genome FASTA filtered to references present in the BAM header.
