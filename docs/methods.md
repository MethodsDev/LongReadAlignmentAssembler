# LongReadAlignmentAssembler (LRAA) — Methods

This section describes the LRAA algorithm and implementation in sufficient detail to reproduce
results and reason about design choices. LRAA reconstructs transcript isoforms and quantifies
their abundances from long-read RNA-seq alignments.

## Overview

Given a coordinate-sorted, indexed BAM of long-read alignments and a reference genome
(optionally a guiding GTF), LRAA:

1) builds a splice graph per contig/strand from read-supported exon/intron evidence;
2) compiles read-supported simple paths into a MultiPath Directed Acyclic Graph (DAG);
3) reconstructs isoforms via a trellis search with iterative best-path selection;
4) quantifies the resulting transcripts using expectation–maximization (EM) with optional
   regularization and weighted read assignments.

Primary outputs include a GTF of assembled isoforms, per-isoform expression estimates, and a
tracking file linking reads to isoforms. Optional BED and debug artifacts are produced under
debug mode.

## Inputs and preprocessing

- Alignments: BAM/BAI, typically from minimap2 or similar long-read mappers. Reads are filtered
  by mapping quality and alignment percent identity.
- Reference genome: FASTA (used for coordinate context and boundary handling).
- Optional annotation: GTF to guide splice-graph construction and transcript boundary priors.

Read ingestion is implemented in `pylib/Bam_alignment_extractor.py`. LRAA filters:
- Minimum mapping quality during discovery/assembly (config `min_mapping_quality`).
- Minimum mapping quality during final quantification (config
  `min_mapping_quality_for_final_quant`).
- Alignment percent identity derived from the `NM`/`nM` tag and the aligned-base count. Reads
  below `config['min_per_id']` are discarded.
- Primary, non-supplementary alignments only, and this is not configurable. Secondary and
  supplementary records are dropped both where the BAM is split by alignment orientation
  (`util/separate_bam_by_strand.py:62-66`) and again during ingestion
  (`pylib/Bam_alignment_extractor.py:139-145`), as are duplicate, QC-failed, and unmapped
  records. Paired-end flags are respected when present.
- Alignments carrying any intron longer than `config['max_intron_length']` are discarded at the
  strand split (`util/separate_bam_by_strand.py:74-75`) and again during ingestion
  (`pylib/Bam_alignment_extractor.py:147-154`); see the configuration section below.

In quantification-only mode, LRAA also supports an optional transcriptome rescue pass
(`--rescue_unassigned_reads_via_transcriptome_alignment`). This pass retries reads that fail
normal graph/isoform assignment, including `DISCARDED-SPACER` reads, reads anchored to a gene
but incompatible with all candidate isoforms, and reads discarded from genomic assignment
because their genomic alignment percent identity falls below the configured threshold. Rescue
is restricted to the local contig/strand transcript set and is performed before quantification
EM.

For memory efficiency, read objects destined for pretty/realignment inspection are “lightened”
unless the read is a soft-clipping realignment candidate.

## Splice-graph construction

LRAA builds a splice graph per contig and strand (`pylib/Splice_graph.py`). Nodes represent
genomic intervals (exons, introns, and boundary anchors), and edges represent splice junction
connectivity observed in the alignments. Key elements:

- Evidence aggregation: Junctions and exon boundaries are accumulated from CIGAR splice ops and
  read coverage.
- Boundary detection: Transcription start sites (TSS) and polyadenylation (PolyA) sites are
  inferred when supported by read boundaries; annotation can also seed boundary nodes.
- Data structures: The graph uses `networkx` for connectivity and `intervaltree` for interval
  queries. Graph parameters are configured via class-level tunables set by
  `Splice_graph.init_sg_params(...)` from the main CLI.
- Componentization: The splice graph is partitioned into weakly connected components to enable
  independent processing and parallelism.

## MultiPath counts and graph

Read alignments are transformed into read-supported simple paths across splice-graph nodes
(`pylib/MultiPath.py`, `pylib/MultiPathCounter.py`). These are summarized into a MultiPath DAG
(`pylib/MultiPathGraph.py`) where each node corresponds to a simple path segment and edges
preserve permissible concatenations. Distinguishing features:

- Boundary flags: MultiPath nodes carry TSS/PolyA indicators, influencing downstream
  compatibility and sorting.
- Boundary snapping for read paths: During read mapping, terminal TSS/PolyA nodes can be
  attached to read paths when termini fall near existing boundary features, terminal soft
  clipping is within configured thresholds, and splice-graph adjacency is valid. Candidate
  features are searched within half of the configured alternative-site window and ranked by
  nearest distance then support.
- Pruning: Extremely large components are pruned using `config['max_path_nodes_per_component']`
  to keep search tractable.
- Caching/read IDs: Read–path associations are tracked (`MpReadIdStore`, `ReadNameStore`) to
  support later compatibility checks and quantification.

## Isoform reconstruction

Isoforms are reconstructed per component in `pylib/LRAA.py` using a trellis search over the
MultiPath DAG and an iterative best-path selection procedure:

1) Trellis build: The DAG is unrolled into a trellis indexed by topological order. Dynamic
   programming aggregates support scores along candidate paths.
2) Path scoring: Path scores combine read support (counts/coverage), boundary consistency
   (TSS/PolyA), and penalties for weakly supported transitions. Ties are broken using
   boundary-aware sort weights (`Transcript.get_left/right_boundary_sort_weight`).
3) Best-path extraction: The highest-scoring path is materialized as a `Transcript` object
   (simple path representation with boundary flags). Reads strongly explained by the selected
   path are down-weighted or removed.
4) Iteration: Steps 1–3 repeat until no path exceeds minimum support or a maximum count is
   reached. The procedure is applied independently to each component, with large components
   optionally processed via multiprocessing.

This approach balances sensitivity (capturing alternate splicing and variable boundaries) with
specificity (avoiding combinatorial explosion in large graphs) by constraining the search to
read-supported MultiPath segments and using iterative extraction.

## Quantification via EM

After reconstruction, LRAA assigns reads to transcripts and estimates abundances using an EM
algorithm (`pylib/Quantify.py`). Assignment logic:

- Compatibility: A read is considered compatible with a transcript if it matches splicing
  structure and overlaps sufficiently ( `config['fraction_read_align_overlap']`). Exact matches
  may be prioritized; otherwise, compatible matches form a candidate set.
- Assignment cascade behavior: In HiFi-style matching, boundary nodes (TSS/PolyA) are respected
  first; if unresolved, matching falls back to boundary-trimmed non-HiFi-style compatibility
  checks.
- Weighted assignments: If `config['weight_reads_by_3prime_agreement']` is enabled,
  compatibility weights prioritize agreement of read 3' ends with transcript 3' ends.

Transcriptome rescue realigns reads that cannot be placed on the splice graph to local
multi-exon transcript sequences using minimap2 in non-splice mode. Rescue is on or off
(`config['rescue_unassigned_reads_via_transcriptome_alignment']`); there is no mode that
assigns every read by transcript alignment, and none that adjudicates a whole-genome
assignment against a whole-transcriptome one. It runs in two places, not only in
quantification-only mode. During reference-guided discovery, only reads with no path
through the splice graph are eligible: a read that has a path but matches no supplied model is
evidence for a novel isoform, and reshaping it onto a reference structure would suppress that
isoform. During final quantification the eligible reads are those the genomic graph could not
place plus those whose path matched no target, and the rescue targets are the isoforms being
quantified. A rescue is confined to targets the read's own genome alignment overlaps: the
genome alignment is trusted for locality but not for optimality, so rescue may refine
placement within the locus the read occupies and may not move it to another one. Rescue
alignments are accepted only when, after the same
small-indel block-merging logic used for genomic pretty alignments, the alignment collapses to
a single contiguous transcript-coordinate interval and at least
`config['rescue_unassigned_min_aligned_read_frac']` of the read's length aligns to the target.
That coverage bar is aligned length over read length, not matched bases, so it is independent
of platform error rate; percent identity is measured over the aligned portion only and cannot
see a read that aligns locally and clips the remainder. An alignment that explains fewer of
the read's bases than its genome alignment already does, or that agrees with its target worse
than the genome alignment agrees with the genome, is declined -- an accepted rescue detaches
the read from the path it already had, so admitting a worse alignment would withdraw support
from the better structure. Accepted rescue hits are projected back
onto splice-graph node paths; ambiguous transcriptome hits are retained only when all top hits
imply the same node path. These rescued paths are merged into the evidence set before the first
EM iteration, so quantification continues under the same compatibility and weighting model as
native genomic assignments.

Model and updates:

- Let transcripts be indexed by i. The per-locus EM operates over read-compatibility classes --
  reads sharing an identical labeled splice-graph path -- rather than individual reads; index
  those classes by r, each carrying a read count. Define a per-class likelihood L_{ri}
  capturing compatibility/weight.
- Let p_i be the isoform proportion parameters, initialized uniformly. Regularization uses a
  Dirichlet prior whose mass is scaled per isoform by the ambiguous support that isoform
  competes for, with α = `config['EM_alpha']` (see the M-step).

E-step (responsibilities):

$w_{ri} = \frac{p_i \cdot L_{ri}}{\sum_j p_j \cdot L_{rj}}$ for transcripts compatible with
read r, else 0.

M-step (update proportions):

$p_i \leftarrow \frac{\sum_r w_{ri} + \alpha_i}{\sum_k \left(\sum_r w_{rk} + \alpha_k\right)}$,
where $\alpha_i = \alpha \cdot A_i$ and $A_i$ is the read count summed over compatibility
classes that include isoform i but are compatible with more than one isoform. The prior mass
therefore scales with how much ambiguous evidence an isoform competes for; it is not a constant
$\alpha - 1$. An isoform with no ambiguous support receives no prior mass.

Counts and TPM are derived from estimated proportions with effective-length or read-depth
normalization as appropriate for long reads. Iteration stops when the Euclidean norm of the
change in the proportion vector between consecutive iterations falls below
`config['EM_convergence_tol']`, or at the iteration cap (`max_EM_iterations_during_asm` during
assembly, `max_EM_iterations_quant_only` for final quantification and quantification-only
analysis). No log-likelihood is computed.

Iteration stops after an M-step, so a final E-step is then run at the proportions being
returned. Without it the reported responsibilities would come from the previous iteration's
proportions while the reported proportions are the updated ones, and any consumer recombining
them — a streaming assignment pass splitting a read it must resolve itself, or a reader
comparing a `frac_assigned` in `quant.tracking` against the abundance in `quant.expr` — would
be mixing two different estimates with nothing in the output to say so. The gap is bounded by
`config['EM_convergence_tol']` only when EM converged; on a gene that reaches the iteration cap
it is unbounded, and those are the genes with the most ambiguous reads. Proportions are not
re-updated by this step, so abundances are unchanged; only the responsibilities and the counts
derived from them are brought onto the same footing.

## Configuration and tunables

Global configuration lives in `pylib/LRAA_Globals.py` as a `config` dictionary. The CLI updates
many values directly, and additional overrides can be provided via `--config_update` JSON.
Notable keys include:

- Read filtering: `min_mapping_quality`, `min_mapping_quality_for_final_quant`, `min_per_id`.
- Intron length: `max_intron_length` (default 200,000). Any alignment containing an intron
  (CIGAR `N` operation) longer than this is discarded whole, both at the strand split that opens
  depth normalization (`util/separate_bam_by_strand.py:74-75`) and at read ingestion
  (`pylib/Bam_alignment_extractor.py:147-154`); `0` or a negative value disables the filter.
  Both call sites share one implementation, `Util_funcs.has_disqualifying_long_intron`, so the
  splice-graph evidence and the reads offered to assignment cannot disagree about which
  alignments are modellable. Ingestion is on the path of every quantified read, so the discard
  holds under `--no_norm` and `--normalize_max_cov_level 0`, which skip the strand split
  entirely (`LRAA:1385-1440`). This differs in kind from the earlier use of the same threshold,
  which pruned only an overlong intron at an alignment terminus
  (`pylib/Pretty_alignment.py:330,341`) and dropped the overlong intron edge from the graph
  (`pylib/Splice_graph.py:1515`) while the alignment itself still spanned the distance.
- Graph scale: `max_path_nodes_per_component`, thresholds controlling junction/exon evidence.
- Assignment/EM: `fraction_read_align_overlap`, `weight_reads_by_3prime_agreement`, `EM_alpha`.
- Monoexonic isoform confidence: `min_monoexonic_TPM`, plus
  `min_monoexonic_read_span_peak_frac` and `min_monoexonic_adjusted_TPM_ratio`, which require
  the supporting reads of a single-exon model to overlap one another rather than tile a span.
  Either may be set to `0` to disable it.
- Monoexonic isoform confidence, single-cell only: `min_monoexonic_supporting_cells` requires a
  novel single-exon model to be seen in that many distinct cells; reference-matching monoexonic
  models are exempt and bulk data is unaffected. The roster of real cells comes from
  `--cell_list` (per-cluster in cluster-guided runs), so empty droplets and ambient barcodes do
  not count. `0` disables.
- Transcriptome rescue: `rescue_unassigned_reads_via_transcriptome_alignment`,
  `rescue_unassigned_min_aligned_read_frac`, `rescue_unassigned_min_per_id`.
- Parallelism: `CPU`, `min_mpgn_component_size_for_spawn`.
- Coverage normalization: `normalize_max_cov_level` sets the read depth that splice-graph
  evidence is thinned toward, per alignment orientation — the BAM is split on `read.is_forward`
  (`util/separate_bam_by_strand.py:239`), not on an inferred transcribed strand; `0` disables it.
  Coverage below the target is retained in full, and retained reads carry the reciprocal of their
  acceptance probability in the `XW` tag so support stays on the scale of the original BAM
  (inverse-probability weighting, i.e. the Horvitz–Thompson estimator; a weighted sum's precision
  is set by how many reads it retained, roughly `1/sqrt(retained)`). The level is a target rather
  than a ceiling: reads carrying a junction supported by fewer reads than the level are kept
  unconditionally and acceptance is a per-read random draw, so realised depth
  can exceed it (`util/normalize_bam_by_strand.py:275,316-319,332-335`). Quantification reports
  counts derived from the unnormalized BAM; under the default two-pass streaming path the
  abundance-estimating first pass reads the normalized BAM and the second pass streams the
  unnormalized one to assign every read (`LRAA:5702-5703`), while with `--no_stream_reads` a
  single pass reads the unnormalized BAM directly. See `docs/coverage_normalization.md`.
- Streaming quantification: `stream_reads` (on by default since v0.25.0) quantifies in two passes
  — abundances from the coverage-normalized BAM, then a streaming pass over the full BAM that
  assigns each read from a precomputed path-to-isoform table without retaining it, so peak memory
  is bounded by the annotation rather than the library. It requires a first-pass BAM thinner than
  the streamed one and is refused with normalization disabled unless a distinct `--bam_for_sg` is
  supplied (`LRAA:1874-1896`). Related: `stream_reads_rescue_unassigned`,
  `stream_reads_rescue_unassigned_to_targets`. See `docs/streaming_quantification.md`.
- Debug: `--debug` enables extensive intermediate artifacts.

Splice-graph parameters are set via `Splice_graph.init_sg_params(...)` inside `LRAA` to keep
all graph-level thresholds centralized.

## Parallelism and performance

LRAA partitions work by contig and by splice-graph components. Components exceeding
`config['min_mpgn_component_size_for_spawn']` are processed in parallel using a lightweight
multiprocessing manager and picklable objects. Memory use is mitigated by reading alignments
incrementally, discarding low-identity or secondary alignments early, and lightening read objects
unless they are candidates for soft-clip realignment. Beyond that, the default two-pass streaming
path removes per-read state from quantification altogether (`docs/streaming_quantification.md`).

## Outputs

- Assembly: `LRAA.gtf` (and optional `.bed`), with one record per transcript isoform;
  boundaries mark TSS/PolyA when applicable.
- Quantification: `LRAA.quant.expr` (per-transcript counts, final-report TPM, and
  `RPM_total_reads`) and `LRAA.quant.tracking.gz` (read-to-transcript compatibility/assignment
  details). The `TPM` column is normalized over final reported transcripts; `RPM_total_reads`
  scales against `num_total_reads`, the number of genome-mapped reads in the input BAM
  (unmapped, secondary, and supplementary records excluded; `LRAA:5137-5163`), counted from the
  input BAM itself so it does not move when a later filter drops an alignment from splice-graph
  evidence.

  Tracking is always gzipped, including the per-(contig, strand) temporary files the parallel
  merge concatenates. It carries one row per (read, compatible isoform) and is the only output
  that scales with library size rather than with the annotation — roughly 140 bytes per row at
  ~2.2 rows per assigned read, so a billion-read library implies about 0.30 TB uncompressed
  against 0.05 TB compressed (measured 5.57x on ONT chr20). There is no option to disable it:
  a run that could emit either form would leave consumers guessing which of two filenames is
  current, and every in-repo reader already selects its decompressor from the suffix.
- Debug (optional): `__*` files including component descriptions and intermediate GTF/BEDs of
  MultiPath graphs and trellis selections.

## Implementation notes

- Language and libraries: Python; key dependencies are `pysam`, `networkx`, and `intervaltree`.
  Auxiliary tools include `samtools` (BAM indexing), and `gffcompare` used in evaluation
  workflows.
- Code structure: Orchestrated by the `LRAA` script (entrypoint, not a package). Core modules
  include `Splice_graph.py`, `MultiPathGraph.py`, `LRAA.py`, `Quantify.py`, and
  `Transcript.py`. Utilities live in `util/` and are called from the main script as needed.

## Reproducibility

- Deterministic tests: Unit tests in `pylib/test_*.py` cover annotation compatibility,
  MultiPath counting, and assignment logic. Example scenarios are available under `testing/`
  (e.g., SIRVs and quant-only cases), with `testing/Makefile` documenting canonical runs of
  `./LRAA`.
- Environments and containers: A Docker image defined in `Docker/Dockerfile` pins dependencies.
  End-to-end WDL workflows are provided under `WDL/` for smoke tests and evaluation.
- Configuration capture: For publication, record the exact LRAA version (git tag/commit),
  Docker image tag, key `config` overrides, and the read aligner/version/settings.

## Limitations and edge cases

- Extremely complex loci may exceed pruning thresholds; results then reflect the
  highest-confidence subgraph.
- Boundary inference depends on coverage of read starts/ends; weak support can blur TSS/PolyA
  distinctions.
- Read assignment prioritizes splicing compatibility; small indels and alignment artifacts near
  junctions may down-weight compatibility in borderline cases.

## Availability

LRAA is open source under this repository. See `README.md` for installation and quickstart
commands, and `Docker/` for containerized execution.

## Troubleshooting parallel runs

- Per-(contig,strand) worker logs: In parallel mode, each worker writes outputs under
  `__<output_prefix>.contigtmp/<contig>/<strand>/`. On failures, an error log is written to
  `<contig>.<strand>.err.log` in that directory. This file now also captures stderr and Python
  faulthandler dumps for fatal errors (e.g., segmentation faults), so it should contain
  diagnostics even if the process exited abruptly.
- Missing error log: If a job fails and no `.err.log` is present, check the main console output
  and any scheduler logs. Also look at the worker failure line emitted by LRAA showing the
  worker exit code. A negative exit code indicates the signal that terminated the worker (e.g.,
  `-9` → SIGKILL, typically OOM).
- Resuming: On reruns, LRAA will skip completed contig/strand jobs detected by the presence of
  `<contig>.<strand>.ok` and non-empty outputs in the same directory unless
  `--no_resume_parallel` is used.
