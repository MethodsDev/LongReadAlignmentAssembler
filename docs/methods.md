# LongReadAlignmentAssembler (LRAA) — Methods

This section describes the LRAA algorithm and implementation in sufficient detail to reproduce results and reason about design choices. LRAA reconstructs transcript isoforms and quantifies their abundances from long-read RNA-seq alignments.

## Overview

Given a coordinate-sorted, indexed BAM of long-read alignments and a reference genome (optionally a guiding GTF), LRAA:

1) builds a splice graph per contig/strand from read-supported exon/intron evidence;
2) compiles read-supported simple paths into a MultiPath Directed Acyclic Graph (DAG);
3) reconstructs isoforms via a trellis search with iterative best-path selection;
4) quantifies the resulting transcripts using expectation–maximization (EM) with optional regularization and weighted read assignments.

Primary outputs include a GTF of assembled isoforms, per-isoform expression estimates, and a tracking file linking reads to isoforms. Optional BED and debug artifacts are produced under debug mode.

## Inputs and preprocessing

- Alignments: BAM/BAI, typically from minimap2 or similar long-read mappers. Reads are filtered by mapping quality and alignment percent identity.
- Reference genome: FASTA (used for coordinate context and boundary handling).
- Optional annotation: GTF to guide splice-graph construction and transcript boundary priors.

Read ingestion is implemented in `pylib/Bam_alignment_extractor.py`. LRAA filters:
- Minimum mapping quality (config `min_mapping_quality`).
- Alignment percent identity derived from the `NM`/`nM` tag and the aligned-base count. Reads below `config['min_per_id']` are discarded.
- Secondary, duplicate, and QC-failed alignments are excluded. Paired-end flags are respected when present.

For memory efficiency, read objects destined for pretty/realignment inspection are “lightened” unless the read is a soft-clipping realignment candidate.

## Splice-graph construction

LRAA builds a splice graph per contig and strand (`pylib/Splice_graph.py`). Nodes represent genomic intervals (exons, introns, and boundary anchors), and edges represent splice junction connectivity observed in the alignments. Key elements:

- Evidence aggregation: Junctions and exon boundaries are accumulated from CIGAR splice ops and read coverage.
- Boundary detection: Transcription start sites (TSS) and polyadenylation (PolyA) sites are inferred when supported by read boundaries; annotation can also seed boundary nodes.
- Data structures: The graph uses `networkx` for connectivity and `intervaltree` for interval queries. Graph parameters are configured via class-level tunables set by `Splice_graph.init_sg_params(...)` from the main CLI.
- Componentization: The splice graph is partitioned into weakly connected components to enable independent processing and parallelism.

## MultiPath counts and graph

Read alignments are transformed into read-supported simple paths across splice-graph nodes (`pylib/MultiPath.py`, `pylib/MultiPathCounter.py`). These are summarized into a MultiPath DAG (`pylib/MultiPathGraph.py`) where each node corresponds to a simple path segment and edges preserve permissible concatenations. Distinguishing features:

- Boundary flags: MultiPath nodes carry TSS/PolyA indicators, influencing downstream compatibility and sorting.
- Pruning: Extremely large components are pruned using `config['max_path_nodes_per_component']` to keep search tractable.
- Caching/read IDs: Read–path associations are tracked (`MpReadIdStore`, `ReadNameStore`) to support later compatibility checks and quantification.

## Isoform reconstruction

Isoforms are reconstructed per component in `pylib/LRAA.py` using a trellis search over the MultiPath DAG and an iterative best-path selection procedure:

1) Trellis build: The DAG is unrolled into a trellis indexed by topological order. Dynamic programming aggregates support scores along candidate paths.
2) Path scoring: Path scores combine read support (counts/coverage), boundary consistency (TSS/PolyA), and penalties for weakly supported transitions. Ties are broken using boundary-aware sort weights (`Transcript.get_left/right_boundary_sort_weight`).
3) Best-path extraction: The highest-scoring path is materialized as a `Transcript` object (simple path representation with boundary flags). Reads strongly explained by the selected path are down-weighted or removed.
4) Iteration: Steps 1–3 repeat until no path exceeds minimum support or a maximum count is reached. The procedure is applied independently to each component, with large components optionally processed via multiprocessing.

This approach balances sensitivity (capturing alternate splicing and variable boundaries) with specificity (avoiding combinatorial explosion in large graphs) by constraining the search to read-supported MultiPath segments and using iterative extraction.

## Quantification via EM

After reconstruction, LRAA assigns reads to transcripts and estimates abundances using an EM algorithm (`pylib/Quantify.py`). Assignment logic:

- Compatibility: A read is considered compatible with a transcript if it matches splicing structure and overlaps sufficiently (
  `config['fraction_read_align_overlap']`). Exact matches may be prioritized; otherwise, compatible matches form a candidate set.
- Weighted assignments: If `config['use_weighted_read_assignments']` is enabled, compatibility weights reflect alignment quality or overlap extent.

Model and updates:

- Let transcripts be indexed by i, reads by r. Define a per-read likelihood L_{ri} capturing compatibility/weight.
- Let p_i be the isoform proportion parameters (initialized uniformly or proportional to initial counts). With optional Dirichlet prior (regularization) α = `config['EM_alpha']`.

E-step (responsibilities):

$w_{ri} = \frac{p_i \cdot L_{ri}}{\sum_j p_j \cdot L_{rj}}$ for transcripts compatible with read r, else 0.

M-step (update proportions):

$p_i \leftarrow \frac{\sum_r w_{ri} + (\alpha - 1)}{\sum_k \left(\sum_r w_{rk} + (\alpha - 1)\right)}$.

Counts and TPM are derived from estimated proportions with effective-length or read-depth normalization as appropriate for long reads. The algorithm iterates to convergence or a maximum iteration threshold; convergence is detected by small changes in the log-likelihood or parameter deltas.

## Configuration and tunables

Global configuration lives in `pylib/LRAA_Globals.py` as a `config` dictionary. The CLI updates many values directly, and additional overrides can be provided via `--config_update` JSON. Notable keys include:

- Read filtering: `min_mapping_quality`, `min_per_id`.
- Graph scale: `max_path_nodes_per_component`, thresholds controlling junction/exon evidence.
- Assignment/EM: `fraction_read_align_overlap`, `use_weighted_read_assignments`, `EM_alpha`.
- Parallelism: `CPU`, `min_mpgn_component_size_for_spawn`.
- Normalization/debug: flags to normalize BAM by strand; `--debug` enables extensive intermediate artifacts.

Splice-graph parameters are set via `Splice_graph.init_sg_params(...)` inside `LRAA` to keep all graph-level thresholds centralized.

## Parallelism and performance

LRAA partitions work by contig and by splice-graph components. Components exceeding `config['min_mpgn_component_size_for_spawn']` are processed in parallel using a lightweight multiprocessing manager and picklable objects. Memory use is mitigated by streaming reads, discarding low-identity or secondary alignments early, and lightening read objects unless they are candidates for soft-clip realignment.

## Outputs

- Assembly: `LRAA.gtf` (and optional `.bed`), with one record per transcript isoform; boundaries mark TSS/PolyA when applicable.
- Quantification: `LRAA.quant.expr` (per-transcript counts/TPM) and `LRAA.quant.tracking` (read-to-transcript compatibility/assignment details).
- Debug (optional): `__*` files including component descriptions and intermediate GTF/BEDs of MultiPath graphs and trellis selections.

## Implementation notes

- Language and libraries: Python; key dependencies are `pysam`, `networkx`, and `intervaltree`. Auxiliary tools include `samtools` (BAM indexing), and `gffcompare` used in evaluation workflows.
- Code structure: Orchestrated by the `LRAA` script (entrypoint, not a package). Core modules include `Splice_graph.py`, `MultiPathGraph.py`, `LRAA.py`, `Quantify.py`, and `Transcript.py`. Utilities live in `util/` and are called from the main script as needed.

## Reproducibility

- Deterministic tests: Unit tests in `pylib/test_*.py` cover annotation compatibility, MultiPath counting, and assignment logic. Example scenarios are available under `testing/` (e.g., SIRVs and quant-only cases), with `testing/Makefile` documenting canonical runs of `./LRAA`.
- Environments and containers: A Docker image defined in `Docker/Dockerfile` pins dependencies. End-to-end WDL workflows are provided under `WDL/` for smoke tests and evaluation.
- Configuration capture: For publication, record the exact LRAA version (git tag/commit), Docker image tag, key `config` overrides, and the read aligner/version/settings.

## Limitations and edge cases

- Extremely complex loci may exceed pruning thresholds; results then reflect the highest-confidence subgraph.
- Boundary inference depends on coverage of read starts/ends; weak support can blur TSS/PolyA distinctions.
- Read assignment prioritizes splicing compatibility; small indels and alignment artifacts near junctions may down-weight compatibility in borderline cases.

## Availability

LRAA is open source under this repository. See `README.md` for installation and quickstart commands, and `Docker/` for containerized execution.
