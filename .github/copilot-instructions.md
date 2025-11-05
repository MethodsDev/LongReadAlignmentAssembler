# LongReadAlignmentAssembler – AI contributor guide

This repo implements LRAA: isoform discovery and/or quantification from long-read RNA-seq. Use these notes to navigate architecture, workflows, and project-specific patterns quickly.

## Big picture and data flow
- Top-level entrypoint is `LRAA` (Python script, not a package). It orchestrates the run per contig and strand, and writes outputs (`LRAA.gtf`, `LRAA.quant.expr`, `LRAA.quant.tracking`, optional `LRAA.bed`, and prelim debug files prefixed by `__`).
- Core pipeline:
  1) Build splice graph from BAM alignments and genome: `pylib/Splice_graph.py` (see `Splice_graph.build_splice_graph_for_contig`). TSS/PolyA sites may be inferred; exon/intron nodes are constructed using `networkx` + `intervaltree`.
  2) Build MultiPath counts from read alignments, then a MultiPath graph: `pylib/LRAA.py` (class) → `pylib/MultiPathGraph.py`. Nodes represent read-supported simple paths over splice-graph nodes; boundaries use TSS/PolyA where present.
  3) Reconstruct isoforms via trellis search and iterative best-path selection (`pylib/LRAA.py:LRAA.reconstruct_isoforms`, `_reconstruct_isoforms_single_component`) into `Transcript` objects.
  4) Quantify isoforms with EM: `pylib/Quantify.py` assigns reads to transcripts and runs EM (configurable) to produce per-isoform counts/TPM and optional BAM tagging.

## Modes, CLI, and outputs
- CLI lives in `LRAA` (main). Key modes:
  - Discovery + quantification (default): requires `--bam` and `--genome`; optional `--gtf` to guide assembly.
  - Quantification-only: `--quant_only --gtf <targets.gtf>`; writes `*.quant.expr` and `*.quant.tracking`; optional `--tag_bam` annotates reads in the BAM.
  - Restrict by contig/region: `--contig chr2` and/or `--region chr2:12345-56789` (append `+`/`-` to fix strand).
  - Filtering switches: `--ME_only`, `--SE_only`, `--no_filter_isoforms`, `--no_norm` or `--normalize_max_cov_level N`.
  - Performance: `--CPU N` controls multiprocessing; components spawn when size ≥ `config['min_mpgn_component_size_for_spawn']`.
- Outputs per run: `LRAA.gtf` (+ `.bed`), `LRAA.quant.expr`, `LRAA.quant.tracking`; many debug artifacts when `--debug` or `LRAA_Globals.DEBUG` (files/dirs prefixed with `__`).

## Configuration conventions
- Global config is centralized in `pylib/LRAA_Globals.py` (`config` dict). The CLI updates many entries directly; see the section “update config” in `LRAA`.
- For ad-hoc overrides, `--config_update my.json` applies key:value pairs post-CLI parsing. The main script gently casts types to match existing defaults. When adding new config keys, define a default in `LRAA_Globals.config` and (optionally) surface a CLI flag.
- Splice-graph parameters are set via `Splice_graph.init_sg_params(...)` in `LRAA`; update there if adding new SG-level thresholds.

## Core modules and how they interact
- `pylib/Splice_graph.py`: Builds and holds the splice digraph; provides interval queries (exon/intron), node accessors, and component partitioning. Class variables hold tunables (set by `init_sg_params`).
- `pylib/LRAA.py` (class): Bridges splice graph → multipath graph, reconstructs isoforms, manages multiprocessing and scoring. Writes debug GTFs and component dumps when debugging.
- `pylib/MultiPathGraph.py`: Converts read-supported multipaths into a DAG of `MultiPathGraphNode` with boundary flags (TSS/PolyA). Prunes large components using `config['max_path_nodes_per_component']`.
- `pylib/Quantify.py`: Assigns reads to transcripts (exact/compatible tests) and runs EM. See `_assign_path_to_transcript` and `_estimate_isoform_read_support` for assignment logic.
- `pylib/Transcript.py`: Transcript feature container with simple-path representation, TSS/PolyA helpers, and quant fields.

## Project-specific patterns and gotchas
- Boundary nodes: TSS/PolyA affect sorting and compatibility checks; see `Transcript.get_left/right_boundary_sort_weight` and checks in `MultiPathGraph` and `Quantify`.
- Multiprocessing: Only components above a threshold spawn; jobs communicate via `MultiProcessManager.Queue`. Ensure objects passed between processes are picklable.
- Debug mode: Set via `--debug` or `LRAA_Globals.DEBUG`; expect numerous `__*` files (e.g., `__mpgns.*.gtf`, `__MPGN_components_described.*.bed`). Avoid committing these.
- BAM index: `LRAA` auto-builds a `.bai` with `samtools index` if missing. The normalization step uses `util/normalize_bam_by_strand.py` when enabled.
- Read assignment: Default requires substantial overlap (`config['fraction_read_align_overlap']`) and may weight assignments (`config['weight_reads_by_3prime_agreement']`, based on 3' end agreement). EM regularization is `config['EM_alpha']`.

## Dependencies and external tools
- Python libs: `pysam`, `networkx`, `intervaltree`, plus scientific stack for utilities/tests. See `Docker/Dockerfile` for an authoritative list.
- External tools used in workflows/tests: `samtools` (indexing), `gffcompare` (evaluation in `testing/sirvs`), and `miniwdl` for WDL-based smoke tests.
- Utilities live in `util/` (Perl and Python helpers for BAM/GTF transformations) and are invoked by the main script.

## Testing and examples
- End-to-end examples: `testing/Makefile` runs scenarios in subfolders (SIRVs, quant-only, etc.). The `testing/sirvs/Makefile` shows canonical invocations of `../../LRAA` and expected outputs.
- Unit tests: Pytest tests live under `pylib/` (e.g., `test_sqanti_like_annotator.py`, `test_diff_iso_tests.py`). Keep new tests colocated in `pylib/` and prefer small, deterministic fixtures.

## When adding features
- New CLI options: add to `LRAA`, update `LRAA_Globals.config` default, and ensure the option propagates into the module that consumes it (e.g., splice graph params via `init_sg_params`).
- Isoform scoring/selection: changes usually live in `pylib/LRAA.py` (trellis build/score thresholds) and must remain consistent with `Quantify` assignment logic.
- Splice-graph heuristics: prefer class-level tunables in `Splice_graph` and make them controllable via config.
- Maintain logging style (`logging` module) and avoid broad refactors; many types are shared across modules.

If anything above is unclear or you need deeper examples (e.g., trellis details, EM tuning, or testing specific scenarios), ask and we’ll refine this guide.
