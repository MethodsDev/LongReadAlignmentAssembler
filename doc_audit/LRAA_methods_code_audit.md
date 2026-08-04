# LRAA Methods-to-Code Audit

## Audit basis

- Date: 2026-08-03
- Manuscript: `/home/unix/bhaas/projects/LRAA_Isoform_ID_n_Quant_Bmark_Platform/LRAA-draft.Aug032026.pdf`
- Manuscript SHA-256: `69dc771b9adf24d9710af91dabff151b26191197af2ead54a87a04c8f4a50ed3`
- Manuscript modification time: 2026-08-03 19:01:50 -0400
- Repository: `LongReadAlignmentAssembler`, branch `devel`, commit `f20bbc209e2fcf607f4bc9897f992b947e4ff6ac`
- Worktree: the audit used the current dirty worktree. Its tracked diff before this report had SHA-256 `9b0f43e6be5e9297bc153fb7a6e01acc1811baf2729b21c443d261b5bb1e4513`. Relevant modified files included `WDL/LRAA-singlecell.wdl`, `WDL/subwdls/LRAA-gene_sparseM_to_seurat_clusters.wdl`, `util/sc/gene_sparseM_to_seurat_clusters_and_umap.R`, and the Seurat notebook template.
- Methods scope: PDF pages 25–44, extraction lines 288–720. The code audit covered the LRAA algorithm, single-cell workflow, differential isoform usage, ORF/proteoform utilities, and SQANTI-like classification. Wet-lab procedures and benchmark execution were classified separately when this repository could not verify them.
- Exclusions: none carried forward from the prior audit.

## Overall assessment

The revision fixes most of the high-impact algorithm descriptions. The active reconstruction score, one-time trellis construction, complete read exclusion during iterative rescoring, terminal-boundary fallback, gene grouping, nine-step assignment cascade, 3′ distance, EM equations, and EM convergence now match the code.

The Methods are not yet fully code-confirmable. The main remaining discrepancies concern configuration claims, the normalization description, TSS enrichment, several omitted filter conditions, RCC processing order, transcript-filter denominators, contained-model pruning, single-cell resource and merge behavior, and differential isoform usage. One manuscript statement exposes a confirmed implementation defect in the internal-priming exemption.

## Confirmed claims

The following revised claims match active behavior, subject to the mode qualifications stated in the manuscript or below:

- The three main execution modes and their GTF requirements.
- HiFi defaults of 97% identity, TSS/PolyA inference, 1% alternative-intron support, 13-bp terminal-spur limit, and disabled adjacent-boundary aggregation. The default non-HiFi profile uses 80%, no boundary inference, 3%, 14 bp, and enabled 5-nt adjacent-boundary removal.
- Graph construction by contig and strand, canonical strand-consistent junction support, exon mean coverage, and distinct represented-read scoring during reconstruction.
- Removal of the obsolete `<50 bp` exon-segment merge description. Active code merges coordinate-contiguous unbranched exon chains.
- Candidate-site denominators for TSS and PolyA filtering: same-feature support within each initial graph component. The active PolyA default is 1%.
- Terminal-spur topology, reference-exon exemption, and current 14-bp non-HiFi and 13-bp HiFi limits.
- Transcript-boundary graph fracturing in both platform profiles.
- RCC containment stored separately from graph edges, with compatible non-containment edges suppressed by default.
- Distinct-read reconstruction scoring without TSS/PolyA bonuses, feature-support terms, or edge penalties.
- One-time trellis construction, cumulative exclusion of represented read IDs, rescoring, reordering, and default minimum score of one.
- Default monoexonic TPM filtering, terminal percentile adjustment, and the fewer-than-seven-read behavior that retains the existing boundary.
- Preliminary shared-path gene components, 50%-of-shorter Leiden graph, resolution 0.2, and the 1,500-transcript union-find fallback.
- The asymmetric Phase-1 boundary rules, Step-5 boundary behavior, Phase-2 trimming, and the complete nine-step assignment order.
- The 3′ exonic splice-graph distance and softened weight equation.
- The local EM likelihood, E-step, M-step form, Euclidean convergence threshold of `1e-6`, and iteration limits of 1,000 during draft assembly and 250 for final or quantification-only estimation.
- Supplied-model and annotation-free oversimplify branches.
- Assigned-pair tracking output, barcode-tag propagation, mitochondrial regex matching for symbol-prefixed and LRAA contig-derived features, and ranked antisense-aware gene-symbol transfer in the current worktree. Opaque reference gene IDs remain a mitochondrial-QC defect described below.
- emptyDrops, Seurat clustering operations, cluster partitioning by the configured tag, per-cluster discovery, common-graph final quantification, and three sparse-matrix levels.
- The standalone SQANTI-like BAM/GTF classifier.
- The reported ALS/FTD alternative-termini subset is supported by the final artifacts. The splice hash encodes the ordered intron chain. All 1,323 significant rows with matching dominant and alternate hashes paired models that differed at one or both outer transcript boundaries: 365 at TSS only, 919 at TES only, and 39 at both. The notebook does not perform this coordinate check, and the complete GTF contains 15 same-hash groups with duplicate outer coordinates, so future analyses should test terminal coordinates explicitly. Evidence: `pylib/Quantify.py:1346-1352`; `util/sc/notebook_templates/sc_diff_SplicePattern_usage_analysis.Rmd:118-166,250-270`; ALS/FTD `LP00714.withGeneSymbols.gtf` and `LP00714.LRAA.sc^isoform.top_only_w_recip_delta_pi.minCF0.05.diff_iso.tsv`.

## Required manuscript corrections

### 1. “All thresholds are user-configurable” remains false

Claim: PDF page 25 states that all thresholds are user-configurable.

Verdict: Contradicted.

Implemented behavior: many values are configurable, but the secondary-rescue cutoffs are literals. The normalization bin width and seed are not forwarded by LRAA. The late `--config_update` propagation defects found during this audit have been corrected in the current worktree, but these remaining fixed values still contradict the universal claim.

Evidence: `util/filter_bam_to_secondary_rescue.py:83-120`; `util/normalize_bam_by_strand.py:24-39`; `LRAA:66-168,831-847,1041-1201`.

Suggested correction: replace “all thresholds are user-configurable” with “many thresholds are user-configurable,” then identify the supported CLI or JSON controls.

### 2. The literal input BAM need not already be indexed, and eligible primary records face additional filters

Claim: pages 25 and 33 require coordinate-sorted, indexed input BAMs; page 25 states that non-supplementary primary records are retained.

Verdict: Partial/overstated.

Implemented behavior: default secondary processing creates a coordinate-sorted, indexed working BAM. Other paths attempt to create a missing index. Primary records can still be rejected for MAPQ, improper pairing, duplicate status, QC failure, or percent identity. Identity filtering is skipped when neither `NM` nor `nM` exists.

Evidence: `LRAA:1224-1254,1784-1844`; `util/filter_bam_to_secondary_rescue.py:268-289,321-373`; `pylib/Bam_alignment_extractor.py:122-179`.

Suggested correction: say that LRAA consumes a coordinate-sorted, indexed working BAM and can create a missing index. State that eligible primary and rescued secondary alignments remain subject to the listed ingestion filters.

### 3. Normalization does not impose an approximate per-base maximum

Claim: page 26 describes an approximate maximum per-base coverage of 1,000× and then gives the alignment-start-bin procedure.

Verdict: Partial/overstated.

Implemented behavior: the utility samples at most 1,000 records per strand-specific 100-bp alignment-start bin. It never computes or caps per-base depth. Reads from adjacent start bins can overlap the same base, so depth can exceed 1,000×.

Evidence: `util/normalize_bam_by_strand.py:20-40,45-86,93-139`; `LRAA:1709-1759`.

Suggested correction: delete “approximate maximum per-base coverage of 1000×.” Retain the exact strand-specific alignment-start-bin description.

### 4. Soft-clip realignment is default-on and narrower than “optionally examines”

Claim: page 26 presents soft-clip realignment as optional.

Verdict: Partial/overstated.

Implemented behavior: `try_correct_alignments` is true by default on ordinary contigs. Only 5–20-base clips are candidates, and extension across an adjacent intron requires an exact sequence match. The operation changes in-memory alignment segments, not the source BAM. Oversimplified contigs disable it.

Evidence: `pylib/LRAA_Globals.py:21-23`; `pylib/Pretty_alignment_manager.py:179-227`; `pylib/Pretty_alignment.py:263-274,319-508`; `LRAA:2974-2988,3215-3229`.

Suggested correction: begin “By default during read-to-graph path mapping” and state the 5–20-base exact-match rule and oversimplify exception.

### 5. The TSS local-enrichment denominator does not match its definition

Claim: page 26 defines

\[
(C_p+1)/(M_{local}+1) \ge 5
\]

with `M_local` equal to the median of positive raw counts at other positions within 50 bp on either side.

Verdict: Contradicted.

Implemented behavior: code takes the median of `[0, 0]` plus the positive neighboring counts. Its scan covers `[p-50,p+49]`, not a symmetric inclusive interval. A focused check with a peak count of 5 and one neighboring count of 1 retained the peak because the code ratio was 6, while the manuscript definition gives 3 and would reject it.

Evidence: `pylib/Splice_graph.py:2939-2994`.

Suggested correction: either define the two zero baselines and exact half-open window, or change the implementation to the stated positive-count symmetric denominator.

### 6. The splice-graph pruning paragraph omits active filters and an eligibility gate

Claim: pages 26–27 describe the active intron and intron-overlapping exon filters.

Verdict: Partial/overstated.

Implemented behavior: non-reference introns longer than 200 kb are removed. A later exon-island pass removes low-support introns even when they do not share a donor or acceptor. The overlapping-exon ratio is evaluated only when intron support is greater than `1/f`, which is at least 101 reads at the 1% default. That pass runs only in discovery and occurs before TSS/PolyA integration.

Evidence: `pylib/Splice_graph.py:1263-1291,1580-1707,1709-1740,1793-1895`; `LRAA:1304-1316`.

Suggested correction: state that the paragraph lists selected rules, or add the 200-kb and exon-island filters. Add the `>1/f` support gate and discovery-only scope to the overlapping-exon rule.

### 7. RCC capping occurs before RCC relationships and later reconstruction components

Claim: page 28 states that the RCC graph is partitioned into weak components and components over 1,000 RCCs are then reduced.

Verdict: Partial/overstated.

Implemented behavior: the 1,000-path cap is applied during `MultiPathGraph` construction to paths grouped by their splice-graph component. This precedes containment and edge construction. Reconstruction later derives components by transitive sharing of splice-graph path nodes. Sequence length breaks equal-support ties.

Evidence: `pylib/MultiPathGraph.py:116-122,170-193,237-263,455-516`; `pylib/LRAA.py:233-251`; `pylib/LRAA_Globals.py:39-50`.

Suggested correction: place the cap before RCC relationship construction, then describe the later shared-path-node reconstruction components.

### 8. There is not always one trellis vertex per RCC

Claim: page 28 says LRAA creates one trellis vertex for each RCC node.

Verdict: Partial/overstated.

Implemented behavior: an RCC contained by more than ten other RCCs is omitted from the trellis vertex set, although it remains represented evidence through its containers.

Evidence: `pylib/LRAA.py:47-50,469-501,1838-1873`; `pylib/Vertex.py:12-24,43-124`.

Suggested correction: say “one vertex for each eligible RCC; RCCs contained by more than ten other RCCs are omitted as vertices but retained as represented evidence.”

### 9. The unique-read fraction denominator and aggressive-filter logic need correction

Claim: page 29 says isoforms below 1% of gene-level reads or below 1% of uniquely assigned reads are removed, and multiple low-abundance isoforms may be removed when more than ten remain.

Verdict: Partial/overstated.

Implemented behavior: LRAA uses two gene-normalized measures. The isoform fraction is the isoform’s fractionally assigned read count divided by the total assigned read count for the gene. The unique-support fraction is the number of reads in multipaths assigned to the isoform with a fraction of at least 0.9999, divided by that same gene-total read count; it is not divided by the number of uniquely assigned gene reads. While more than ten candidates remain, an aggressive round may remove several non-exempt isoforms, but only when both measures are below their thresholds. If nothing has yet been removed in the round, the standard step removes at most one candidate, the lowest-ranked non-exempt isoform failing either threshold. Assignments and fractions are recalculated after each round. Expressed reference-containing models are exempt under the default retention mode. Novel isoforms also require at least two assigned reads before iterative filtering and two uniquely assigned reads within the iterative filter by default.

Evidence: `pylib/TranscriptFiltering.py:72-250,740-782`; `pylib/Quantify.py:1199-1233,1459-1470`; `pylib/LRAA_Globals.py:96-106`; `LRAA:4529-4567,4601-4616`.

Suggested correction: “Within each gene, LRAA iteratively filters non-exempt isoforms using two gene-normalized measures: the isoform fraction, defined as the fractionally assigned read count for an isoform divided by the total assigned read count for the gene, and the unique-support fraction, defined as the number of reads in RCCs assigned to that isoform with an assignment fraction of at least 0.9999, divided by the same gene-total read count. The default minimum for each measure is 1%. While more than ten isoforms remain, multiple isoforms may be removed in one round, but only if they fall below both thresholds. If no isoform has yet been removed during that round, the lowest-ranked eligible isoform falling below either threshold is removed. Read assignments and isoform fractions are recalculated after each round. Under the default reference-retention mode, models containing an expressed reference transcript are exempt. Novel isoforms additionally require at least two assigned reads and two uniquely assigned reads by default.”

### 10. Contained-transcript pruning is broader and has a reciprocal boundary rule

Claim: page 29 limits eligibility to fully contained models with identical splicing and shorter terminal exons, describes the missing-TSS/same-PolyA case, and states that distinct-boundary models above the 20% threshold are retained.

Verdict: Partial/overstated.

Implemented behavior: within each gene, a model is eligible for pairwise pruning when its terminal-node-stripped path contains no more nodes than the comparison path and either is contained as a contiguous node sequence, including exact equality, or has a compatible gap-free overlap with that path. Full transcript containment and shorter terminal exons are not required. With the default boundary-aware setting, an eligible model is pruned if it lacks both TSS and PolyA annotations, lacks a TSS while sharing the comparison model’s PolyA site, lacks a PolyA site while sharing the comparison model’s TSS, or otherwise contributes less than 20% of the pair’s combined expression. Passing one comparison does not guarantee final retention because the model remains eligible for comparison with other isoforms. The comparison model absorbs the pruned model’s multipath evidence, and LRAA reruns quantification after pruning.

Evidence: `pylib/TranscriptFiltering.py:301-321,336-481`; `pylib/Simple_path_utils.py:46-104,859-885`; `pylib/Transcript.py:421-437`; `pylib/LRAA_Globals.py:87-88,111-112`; `LRAA:4563-4599`.

Suggested correction: “Potential degradation products are evaluated pairwise among isoforms assigned to the same gene. After terminal TSS and PolyA nodes are removed, an isoform is eligible for pruning when its path contains no more splice-graph nodes than the comparison path and is either contained within it or overlaps it compatibly without gaps. Under the default boundary-aware behavior, an eligible isoform is pruned if it lacks both TSS and PolyA annotations, lacks a TSS while sharing the comparison isoform’s PolyA site, or lacks a PolyA site while sharing the comparison isoform’s TSS. Otherwise, it is pruned when it accounts for less than 20% of the pair’s combined expression. An isoform that passes one pairwise comparison may still be pruned in comparison with another isoform. When an isoform is pruned, its supporting read-path evidence is transferred to the comparison isoform, and expression is re-estimated after pruning.”

### 11. The EM prior uses compatible ambiguous support, and draft EM has no regularization

Claim: page 32 says prior mass is scaled by ambiguous RCC support assigned to each isoform and presents `alpha=0.01` without a stage qualification.

Verdict: Partial/overstated.

Implemented behavior: `A_i` sums the full read count of every ambiguous RCC compatible with isoform `i`; it is not posterior-assigned support. Final and quantification-only EM use `alpha=0.01`. Draft assembly/filtering EM explicitly uses `alpha=0`.

Evidence: `pylib/EM.py:72-153,241-326`; `pylib/Quantify.py:1097-1197`; `LRAA:4234-4273,4469-4474,4691-4719`.

Suggested correction: replace “assigned to” with “compatible with,” and state that regularization is disabled during draft assembly EM.

### 12. Single-cell matrices sum fractional assignments, and monoexonic splice features are transcript-specific

Claim: page 34 describes gene, isoform, and splice-pattern read counts and says splice patterns group reads independently of transcript identity.

Verdict: Confirmed with qualification.

Implemented behavior: matrices sum `frac_assigned`; they do not count raw assignment rows, apply `read_weight`, or deduplicate UMIs. Multi-exon splice features use an intron-chain hash. Monoexonic features use `transcript_id`, so their splice-pattern grouping is transcript-specific.

Evidence: `util/sc/singlecell_tracking_to_sparse_matrix.py:85-120,186-266,330-353`; `pylib/Quantify.py:1347-1352`; `pylib/Util_funcs.py:111-117`.

Suggested correction: call the matrix values summed fractional read assignments and restrict transcript-independent splice hashes to multi-exon structures.

### 13. Two Seurat edge descriptions are inaccurate

Claim: page 35 says cells exceeding 20% mitochondrial content are removed and PCA retains 12 components.

Verdict: Partial/overstated.

Implemented behavior: cells are retained only at `percent.mt < 20`, so exactly 20% is removed. `RunPCA` computes its default number of components; components 1–12 are then used for neighbors and UMAP. `ScaleData` also regresses `nFeature_RNA`, `nCount_RNA`, and `percent.mt`.

Evidence: `util/sc/gene_sparseM_to_seurat_clusters_and_umap.R:54-99`; `WDL/subwdls/LRAA-gene_sparseM_to_seurat_clusters.wdl:10-36,75-100`.

Suggested correction: say “20% or greater” and “PCs 1–12 are used for neighbors and UMAP.”

### 14. Cluster-GTF merging does not always retain TSS/PolyA and has no separate merge-stage expression threshold

Claim: page 36 says specialized LRAA merging retains annotated TSS/PolyA sites and that single-cluster models are included if they pass merge expression and support thresholds.

Verdict: Partial/overstated.

Implemented behavior: merge reassembly recognizes TSS/PolyA only when `HiFi=true`; the single-cell default is false. Each input model receives synthetic support during merge. The expression/support filtering occurred during per-cluster discovery; there is no separate merge expression matrix or merge-stage expression threshold.

Evidence: `WDL/LRAA-cell_cluster_guided.wdl:131-140,365-428`; `util/merge_LRAA_GTFs.py:190-207,270-359`; `pylib/LRAA.py:2161-2244`.

Suggested correction: state that boundary annotations are respected in HiFi merge mode and that input models have already passed upstream filters before being represented as synthetic merge evidence.

### 15. Nuclear-genome and resource configurability claims are too broad

Claim: page 36 says discovery is restricted to the nuclear genome and separate memory and CPU parameters exist for discovery, quantification, clustering, and matrix construction.

Verdict: Partial/overstated.

Implemented behavior: normal discovery is bypassed only for contigs listed by `oversimplify`. The single-cell default is `chrM`; alternate names such as `MT` or `M` require configuration. Annotated mitochondrial reads use best overlap, while annotation-free runs create one pseudo-transcript per strand. Memory is configurable for several stages, but Seurat uses a fixed CPU count of 4 and sparse-matrix construction uses 3. Some nested resource inputs are not exposed or propagated.

Evidence: `WDL/LRAA-singlecell.wdl:127-173`; `LRAA:2763-2833,3179-3284`; `WDL/subwdls/LRAA-gene_sparseM_to_seurat_clusters.wdl:102-106`; `WDL/subwdls/LRAA-build_sparse_matrices_from_tracking.wdl:82-86`.

Suggested correction: describe default `chrM` oversimplification rather than a nuclear-genome restriction. Limit the resource claim to the knobs that reach their consumers.

### 16. Differential isoform usage rules do not match the stated defaults and filters

Claim: page 37 says unspliced models were excluded from all analyses; reciprocal testing requires both directional sets to meet `|Delta pi| >= 0.1`; each enriched set must meet 5% cell detection; gene and set depths are at least 25 reads per cluster and per set; and significant results require FDR `< 0.001`.

Verdict: Supported in part for the reported ALS/FTD run, but the threshold and set definitions remain imprecise.

Implemented behavior: the ALS/FTD Makefile explicitly enabled the nondefault `--ignore_unspliced` and `--reciprocal_delta_pi` switches and supplied `--min_cell_fraction 0.05`, `--top_isoforms_each 1`, and `--group_by_feature gene_symbol`. The recorded analysis predates the explicit exon-count field: its `:iso-` test identified 37,320 monoexonic features, and none of the 22,233 reported test rows retained that marker. Unspliced exclusion and reciprocal testing therefore did run for this analysis, even though neither behavior is a CLI default. The current worktree no longer infers splicing status from identifier syntax.

The numerical rules differ from the prose. Gene totals must be at least 25 separately in each cluster, whereas each directional set needs at least 25 reads summed across both clusters. Both directional `Delta pi` sums must be strictly greater than 0.1 before rounding; equality fails. Cell-fraction filtering requires at least one member of each directional set to have detection fraction at least 5% in its enriched cluster; it does not compute the fraction of cells expressing the set union or require every member to pass. Although the command requested one top isoform per cluster, missing-direction expansion produced 2-column contingency tables for 22,105 reported rows and 3-column tables for 128 rows. The final significance flag uses adjusted P value `<= 0.001`, not `< 0.001`.

Artifact checks agree with these rules. No reported row had a gene total below 25 in either cluster or a pooled directional-set total below 25, but thousands had directional-set support below 25 in one cluster. Thirty-six rows contained a multi-isoform directional set with a member below 5% cell detection and passed because another member met the threshold. A focused synthetic check also confirmed that a directional set can pass with zero reads in one cluster when its pooled two-cluster support is 25, and that an exactly representable `|Delta pi|=0.125` fails when the configured threshold is 0.125.

Evidence: ALS/FTD analysis `diff_iso_usage_for_alt_termini_analysis/Makefile:4-5` and `LP00714.LRAA.sc^isoform.top_only_w_recip_delta_pi.minCF0.05.diff_iso.tsv`; `util/sc/diff_iso_usage/sc_pseudobulk_test_isoform_DiffUsage.py:197-237,408-490,540-579`; `pylib/DiffIsoformStatTest.py:150-223,248-325,517-546`.

Suggested correction: identify the nondefault switches used for the ALS/FTD run and define the filters as per-cluster gene depth, pooled two-cluster directional-set depth, at-least-one-member cell detection, strict pre-rounding effect-size gates, selected-feature contingency tables, and adjusted P value `<= 0.001`.


## Suspected or confirmed code defects

### Minus-strand projected rescue clips were mapped to opposite genomic boundaries

Verdict: Confirmed and resolved in the current worktree.

Transcriptome rescue aligns against 5′-to-3′ transcript sequences and obtains soft clips from the low- and high-coordinate ends of that transcript reference. Projection reverses transcript coordinates for minus-strand models and then sorts the genomic segments from low to high. The former implementation passed the transcript-reference clips through unchanged, so the transcript 5′ clip controlled the genomic-left PolyA end and the transcript 3′ clip controlled the genomic-right TSS end. The corrected projection swaps the two clip values for minus-strand models before calling the genomic path mapper; plus-strand behavior is unchanged.

A production-object regression fixture parses asymmetric `95M5S` and `5S95M` transcriptome SAM alignments through `_parse_rescue_alignments`, minus-strand coordinate projection, `LRAA._map_read_to_graph`, and boundary snapping. Before the correction, both rescued paths incorrectly acquired both terminal features. After the correction, a transcript-3′ clip blocks only the genomic-left PolyA snap, and a transcript-5′ clip blocks only the genomic-right TSS snap.

Evidence: `pylib/IsoformReadRescue.py:573-615,789-857`; `pylib/LRAA.py:1394-1582`; `pylib/test_isoform_read_rescue_projection.py:51-139`; active callers in `pylib/LRAA.py:177-228` and `LRAA:4111-4126,4144-4161,4211-4221`.

### Internal-priming annotation exemption is ineffective for monoexonic models

The Methods state that a flagged monoexonic model is retained when its 3′ end matches a known annotation. Current code removes a flagged monoexonic model before consulting the known-end interval tree. With the default `restrict_internal_priming_filter_to_monoexonic=True`, multiexonic models bypass removal before the lookup, while monoexonic models enter the unconditional-removal branch; the interval tree is therefore built but never queried. It only affects behavior when that setting is overridden to `False`, in which case it protects multiexonic models near known 3′ ends from the broader internal-priming filter. A real-object check placed a monoexonic candidate exactly at a supplied known 3′ end and confirmed that it was removed. History identifies this as an omission rather than a regression: commit `5079dbf` introduced the known-end lookup only in the multiexonic branch, and commit `d214266` changed the source of known ends while explicitly retaining the multiexonic-only behavior. No later commit moved the exemption into the monoexonic branch.

Evidence: `pylib/TranscriptFiltering.py:493-597`; `pylib/test_transcript_filtering_internal_priming.py:14-55`; commits `5079dbf`, `166f6bd`, and `d214266`.

Decision needed: move the known-end check before monoexonic removal, or delete the exemption from the Methods.

### Late configuration defects are resolved in the current worktree

Verdict: Resolved.

Named CLI values now seed the authoritative configuration before mode-specific handling, and `--config_update` remains the final precedence layer. Normalization, splice-graph initialization, iterative filtering, and discovery final quantification consume post-update values. Omitted splice-frequency options are tracked with suppressed argparse defaults, so HiFi supplies 0.01 only when `--min_alt_splice_freq` was not given; an explicit `--HiFi --min_alt_splice_freq 0.03` is preserved. Both active alternative-splice consumers use the threshold installed by `Splice_graph.init_sg_params`. `Quantify._assign_reads_to_transcripts` now resolves `fraction_read_align_overlap` from config at call time unless an explicit function argument is supplied.

The ordinary defaults are unchanged: normalization remains 1,000 alignments per strand-specific start bin; alternative-splice thresholds remain 3% in default ONT mode and 1% in HiFi mode; the alternative-unspliced threshold remains 1%; iterative filtering uses the 1% isoform-fraction threshold with EM; relaxed read assignment uses 75%; and final quantification uses `EM_alpha=0.01`.

Evidence: `LRAA:66-168,831-847,1041-1201,1246-1315,1342-1345,1448-1465`; `pylib/Splice_graph.py:1365-1378,1835-1867`; `pylib/Quantify.py:149-162`; `pylib/test_cli_config_resolution.py:58-137`; `pylib/test_splice_graph_alt_splice_threshold.py:37-75`; `pylib/test_quant_read_assignment_arbitration.py:63-116`.

### Top-level single-cell alignment-policy inputs intentionally stop before cluster-guided phases

Decision: retain the current behavior. `LRAA-singlecell.wdl` forwards `allow_secondary_alignments` and `rescue_unassigned_reads_via_transcriptome_alignment` to the initial LRAA call but not to `ClusterGuided`. Per-cluster discovery and final quantification therefore use the nested workflow defaults. Applying these policies to every cluster-guided call would increase the computational cost of an already intensive stage.

Manuscript implication: describe these two inputs as controls for the initial pooled LRAA run, not as global single-cell workflow settings.

Evidence: `WDL/LRAA-singlecell.wdl:190-214,274-301`; `WDL/LRAA-cell_cluster_guided.wdl:20-23,82-105,144-168`.

### Unspliced DIU filtering now uses explicit exon counts

Resolved in the current worktree. LRAA records `num_exons` on every quantification tracking row, including the oversimplify paths. Active Python sparse-matrix construction carries this value into the transcript/splice-hash mapping, gene-symbol incorporation preserves it, and the legacy standalone R converter now does likewise. Both WDL tracking-merge paths reject inputs with differing headers instead of combining old seven-column and new eight-column rows under one schema. When `--ignore_unspliced` is enabled, DIU resolves `num_exons` through original or renamed transcript and splice-hash identifiers and removes features with `num_exons == 1`. Missing, invalid, conflicting, or unmapped metadata is an error; there is no identifier-pattern fallback. The option remains disabled by default.

The recorded ALS/FTD analysis used the former `:iso-` rule and successfully removed all 37,320 marked monoexonic features. Regenerating the mapping with the current workflow makes the filter reliable for quantification-only GENCODE transcript IDs as well.

Evidence: `LRAA:2407-2416,2758-2772,3133-3158,3302-3315`; `pylib/Quantify.py:1354-1404`; `WDL/LRAA.wdl:340-362`; `WDL/LRAA-cell_cluster_guided.wdl:325-349`; `util/sc/singlecell_tracking_to_sparse_matrix.py:199-264`; `util/sc/singlecell_tracking_to_sparse_matrix.R:18-53`; `util/sc/incorporate_gene_symbols_in_sc_features.py:150-190`; `util/sc/diff_iso_usage/sc_pseudobulk_test_isoform_DiffUsage.py:133-219,565-580`; `util/sc/test_splicing_status_propagation.py`; `pylib/test_quantify_tpm_reporting.py:63-113`.

### Mitochondrial QC misses opaque reference gene identifiers

The revised default regex matches conventional `MT-`/`mt-` symbols and LRAA contig-derived identifiers, but quantification-only tracking preserves reference GTF gene IDs such as `ENSG...`. Gene-symbol incorporation occurs after Seurat filtering and clustering. These opaque IDs do not match the regex, so `percent.mt` can remain zero and high-mitochondrial cells can pass QC in the manuscript's GENCODE quantification-only path.

Evidence: `pylib/Transcript.py:879-880`; `pylib/Quantify.py:1377-1383`; `util/sc/gene_sparseM_to_seurat_clusters_and_umap.R:71-72`; `WDL/LRAA-singlecell.wdl:253-267,312-323`.

Decision needed: incorporate symbols before mitochondrial QC or derive mitochondrial features from contig/annotation metadata rather than feature-name regex alone.

### Prior single-cell implementation issues

- Barcode tag: resolved. The configured tag reaches the Python partitioner.
- Default mitochondrial regex: partly resolved. It recognizes `MT-`, `mt-`, and LRAA identifiers on `chrM`, `MT`, or `M`, but not opaque reference GTF gene IDs such as `ENSG...` before symbol incorporation.
- Gene-symbol ranking: resolved in the current dirty worktree by the shared ranked, antisense-aware parser.
- Unmapped partition records: still retained when they carry a recognized barcode. The current Methods no longer claim otherwise.
- Cluster BAM indexes: the partition task still emits no BAIs. Downstream LRAA creates missing indexes, so this is an artifact/efficiency issue rather than a current Methods contradiction.
- Tracking scope: resolved by the revised assigned-pair wording.

Evidence: `util/sc/partition_bam_by_cell_cluster.py:122-182`; `WDL/subwdls/partition_bam_by_cell_cluster.wdl:26-56`; `LRAA:1236-1240`; `pylib/gene_symbol_utils.py:7-126`; `pylib/Transcript.py:879-880`; `pylib/Quantify.py:1377-1383`.

## Unconfirmed implementation risks

The following source-level concerns require dedicated end-to-end fixtures before being reported as defects:

- [INFERENCE] oversimplify can count multiple retained alignment records for one read name while tracking only one final choice (`LRAA:3023-3057,3102-3127,3215-3240,3271-3284`).
- [INFERENCE] default non-HiFi graphs can retain stale component sets after terminal-spur deletion when no TSS/PolyA objects trigger rediscovery (`pylib/Splice_graph.py:313-350,2425-2499`).

## Claims not verifiable from this repository

- Exact study invocations using GENCODE v47, sample-specific mitochondrial thresholds, Seurat marker outputs, GProfiler, and CAS cannot be established without archived commands and run outputs.
- The novel-proteoform CDS-event classifier, protein-sequence deduplication, representative-isoform selection, ColabFold, InterProScan, and ChimeraX execution are absent from this repository. The current ORF wrapper uses TransDecoder 6.0.0, defaults to 100 amino acids, and does not require complete ORFs by default; this does not verify the stated historical TransDecoder 5.7.1 runs at 50, 100, and 200 amino acids.
- Simulation generation, SIRV/MORF processing, truth-set construction, ARD scoring, proxy-truth rules, and multi-tool benchmarking belong to `MethodsDev/LRAA_Isoform_ID_n_Quant_Bmark_Platform`, which the manuscript cites at lines 714–716. They require a separate reproducibility audit against that repository and the deposited data.
- Wet-lab procedures and instrument run statistics require laboratory records or deposited primary data.

## Status of the prior 24 findings

| Prior finding | Current status |
|---|---|
| 1. Normalization mechanism | Partly resolved; exact bin sampling added, per-base characterization remains wrong. |
| 2. Primary-alignment scope | Not resolved; additional ingestion filters remain omitted. |
| 3. Secondary heuristic | Qualitative rule resolved; literal hard-coded thresholds conflict with configurability claim. |
| 4. Quantification BAM wording | Resolved with qualification; it is the filtered unnormalized working BAM. |
| 5. Soft-clip correction | Not resolved; default-on behavior and exact restriction remain omitted. |
| 6. Exon-segment merge | Resolved. |
| 7. Boundary denominators | Resolved; new TSS local-median discrepancy remains. |
| 8. Node support versus reconstruction score | Resolved. |
| 9. Spur and intron-overlap rules | Spur resolved; intron-overlap eligibility and scope remain incomplete. |
| 10. Transcript-boundary fracturing | Resolved. |
| 11. RCC containment versus edges | Resolved. |
| 12. Expression collapse stage | Resolved; later transcript pruning is now separated. |
| 13. Reconstruction scoring | Resolved. |
| 14. Trellis rebuilding and 90% down-weighting | Resolved. |
| 15. Iterative filtering | Partly resolved; denominator and standard/aggressive logic remain wrong. |
| 16. Contained-transcript pruning | Partly resolved; compatible overlap and reciprocal 3′ rule remain omitted. |
| 17. Terminal-boundary fallback | Resolved. |
| 18. Initial gene grouping | Resolved. |
| 19. Transcriptome rescue availability | Mostly resolved; early de novo discovery lacks targets, while final de novo quantification can use reconstructed targets. |
| 20. Phase-1 boundary matching | Resolved. |
| 21. 3′ distance | Resolved. |
| 22. EM convergence | Resolved. |
| 23. Unique-read output threshold | The wrong 0.9999 reporting claim was removed; local exact-1 and cross-gene 0.9995 definitions remain undocumented. |
| 24. Oversimplify modes | Resolved in prose; repeated-alignment counting remains an unconfirmed code risk. |

## Verification performed

- Rendered and inspected PDF pages 26, 32, 33, 37, and 43 where text extraction dropped equations or reordered content.
- Traced CLI defaults through HiFi updates, `--config_update`, class initialization, WDL nesting, and active consumers.
- Ran 12 focused pytest files covering spur thresholds/topology, PolyA defaults, internal priming, boundary-node support-count fallback, assignment arbitration, degradation pruning, gene reclustering, secondary rescue, cross-gene correction, cell-cluster partitioning, mitochondrial regex filtering for symbol/LRAA identifiers, and gene-symbol ranking: `43 passed`. The fewer-than-seven terminal-coordinate rule was source-inspected, not exercised by this suite.
- Ran `miniwdl check WDL/LRAA-singlecell.wdl` and `miniwdl check WDL/LRAA-cell_cluster_guided.wdl`. Both parsed successfully. The checks also reported the unused cluster-guided `diskSizeGB` input and existing shell/static warnings.
- Executed a TSS enrichment discriminator. Code retained a 5-count peak with one 1-count neighbor; the code ratio was 6 and the manuscript-defined ratio was 3.
- Executed the known-3′-end internal-priming case with real `Transcript` objects. The flagged monoexonic candidate was removed despite matching the supplied known end.
- Executed DIU threshold cases. Directional sets with zero reads in one cluster passed when they had 25 reads across both clusters, and an exactly representable effect equal to the threshold failed the strict `>` gate.
- Ran nine focused configuration tests covering named values, JSON precedence, ONT and HiFi defaults, explicit HiFi `0.03`, worker-facing snapshots, both alternative-splice consumers, call-time read-overlap resolution, and explicit overlap precedence: `9 passed`. The changed Python files also passed bytecode compilation, and full CLI help retained the documented splice defaults.
