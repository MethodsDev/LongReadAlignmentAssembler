# LongReadAlignmentAssembler (LRAA) Methods

## Overview

We developed LongReadAlignmentAssembler (LRAA) as a method for transcript isoform discovery and quantification from long-read RNA sequencing data. LRAA reconstructs full-length transcript isoforms by constructing a splice graph from read-supported exons and introns, enumerating read-consistent paths through this graph, and selecting high-confidence isoforms via iterative best-path extraction. Isoform abundances are estimated using an expectation-maximization algorithm that resolves ambiguous read assignments.

LRAA operates in three modes: (1) annotation-free isoform discovery and quantification, requiring only aligned reads and a reference genome; (2) reference annotation-guided isoform discovery and quantification, where known gene models inform splice-graph construction and filtering; and (3) quantification-only mode, where isoform discovery is bypassed and reads are directly assigned to provided reference isoforms for abundance estimation via EM.

LRAA was primarily developed for PacBio HiFi data, which provides high base accuracy (>99%) and enables precise inference of transcription start sites (TSS) and polyadenylation (PolyA) sites from read termini. The method operates in HiFi mode (`--HiFi` flag) for PacBio datasets, activating stringent filtering parameters optimized for high-accuracy reads. LRAA also supports Oxford Nanopore (ONT) data in default mode, which relaxes filtering thresholds to accommodate higher sequencing error rates. Unless otherwise noted, parameter values and algorithmic decisions described below reflect HiFi mode configuration.


## Input and Preprocessing

LRAA requires coordinate-sorted, indexed BAM files of long RNA-seq reads aligned to a reference genome (typically via minimap2), along with the reference genome in FASTA format. Reference gene annotations in GTF format are required for quantification-only or reference annotation-guided isoform discovery modes. Annotation-free isoform identification and quantification only requires the input alignment bam and refernece genome files.

Filtering of alignments: Stringent quality filters are applied during read ingestion. Reads below minimum mapping quality (1) are discarded. In HiFi mode, reads below 98% identity are filtered to leverage high base accuracy and reduce alignment artifacts; for ONT data, this threshold is relaxed to 80%. Secondary, duplicate, and QC-failed alignments are excluded.

Normalization of read coverage for isoform discovery: To reduce computational requirements in ultra-high-coverage regions, alignments are down-sampled to a maximum per-base coverage of 1000× using a strand-aware normalization procedure applied independently per strand before splice-graph construction. The original unnormalized BAM is used for isoform abundance estimation, ensuring quantification accuracy is preserved.

Correction of alignments at soft-clipped termini: Long reads often contain soft-clipped bases at alignment boundaries when short segments fail to align across splice junctions to adjacent exons. Optionally, LRAA examines soft-clipped termini in the splice graph context and realigns them to connected exon segments to incorporate legitimate exonic extensions. Soft-clipped regions at 3' ends containing ≥7 adenines (or thymines for reverse strand) comprising ≥80% of the sequence are identified as post-transcriptional polyA tails and removed, revealing the genomic polyadenylation site.

## Splice-Graph Construction

For each genomic contig and strand, LRAA constructs a directed acyclic graph representing the exon-intron architecture inferred from read alignments. Graph nodes correspond to contiguous genomic intervals supported by read coverage (exon segments) or splice junctions (introns). Adjacent exon segments separated by less than 50 bp with no intervening splicing are merged. Introns supported by spliced read alignments become nodes that connect the corresponding exon nodes in the graph.

Read support is quantified differently for exon and intron nodes. For introns, each read spanning the junction increments the intron's support count. For exons, genomic base coverage is accumulated across all aligned reads, then the mean coverage across the exon segment's coordinates is computed and assigned as the node's support value. These support metrics are used throughout downstream filtering and scoring steps.

In HiFi mode, TSS and PolyA sites are explicitly inferred from read boundary clustering. Reads terminating with minimal soft-clipping (≤0 bases) within a 50 bp window are aggregated to define boundary candidates. A candidate site is retained if supported by at least 5 reads representing ≥5% of gene-level read coverage. TSS and PolyA nodes are integrated as first-class graph features that influence downstream isoform selection. For ONT data, TSS and PolyA inference is disabled; terminal boundaries are instead refined empirically after reconstruction.

To remove spurious features arising from alignment errors, the graph is filtered using several criteria. Introns supported by fewer than 2 reads are removed unless they match known junctions from a provided reference annotation. For competing splice junctions at the same donor or acceptor site, minor junctions representing less than 1% (HiFi mode) or 3% (ONT mode) of total junction support are pruned. Short terminal exon segments (≤13 bp in HiFi mode, ≤20 bp for ONT) not anchored by TSS or PolyA nodes are removed. Unspliced intervals bridging annotated introns are retained only if supported by ≥1% of junction coverage, distinguishing genuine intron retention from unprocessed pre-mRNA. When a reference GTF is provided, known exons and introns are integrated to further guide graph construction.

After filtering, the refined splice graph is partitioned into weakly connected components, each ideally representing a candidate gene locus. In HiFi mode, TSS and PolyA sites undergo additional filtering within each component. Sites representing less than 5% of total boundary read support within their component are removed as likely noise. TSS sites are further examined for potential RNA degradation artifacts: the graph is traversed along linear exon paths from each high-support TSS, and alternative TSS sites along this path with ≤20% of the dominant TSS support are removed as probable 5' degradation products. After this boundary refinement, components are rediscovered to account for any resulting graph topology changes. 

## Labeled Read Path Graph Construction

Each read alignment is represented as a labeled splice graph read path—an ordered sequence of splice-graph nodes encoding the exon-intron structure traversed by the read. Reads with identical splice structure are aggregated, and their labeled read paths are compiled into a directed acyclic graph (Labeled Read Path Graph). Each node in the Labeled Read Path Graph represents a unique consecutive sequence of splice-graph nodes, annotated with genomic coordinates, read support count, and boundary flags indicating whether the path begins or ends at an inferred TSS or PolyA site. Containment relationships are established between nodes where one path is fully encompassed within another. Importantly, even when splice junctions are compatible and one path is contained within another, distinct TSS or PolyA annotations prevent automatic collapsing—such isoforms are retained as separate transcripts representing genuine alternative transcription start sites or polyadenylation sites. Contained paths lacking boundary annotations or sharing identical boundaries with their containing path are evaluated for collapsing based on expression thresholds, as they likely represent incomplete or degraded transcripts rather than biologically distinct isoforms.

The Labeled Read Path Graph is partitioned into weakly connected components using shared splice-graph nodes. Each component represents a set of transcripts sharing at least one exon or intron and is processed independently. To prevent combinatorial explosion in highly complex loci, components exceeding 1000 nodes (default) are reduced to the 1000 most highly-supported nodes based on read count, ensuring computational tractability while retaining the most confidently supported transcript structures.

## Isoform Reconstruction

LRAA reconstructs transcript isoforms from the Labeled Read Path Graph using trellis-based dynamic programming with iterative best-path extraction. The Labeled Read Path Graph is unrolled into a trellis indexed by topological order, where each layer corresponds to a labeled read path node. Dynamic programming aggregates path scores across all possible routes from source to sink nodes.

Path scores integrate multiple factors: read support (primary weight), boundary consistency (bonus scores for paths with explicit TSS and PolyA boundaries), and transition penalties for weak edges or incompatible boundaries. This scoring function prioritizes high-coverage, boundary-anchored paths while permitting discovery of lower-abundance isoforms in subsequent iterations.

The highest-scoring path is identified via backtracking through the trellis and materialized as a transcript. Reads strongly explained by this transcript (≥90% compatibility) are down-weighted or removed from labeled read path counts. The trellis is rebuilt with updated counts, and the process repeats until no remaining path exceeds minimum support (default: 1 read) or all reads are assigned.

## Isoform Filtering

Reconstructed isoforms are filtered to remove very lowly expressed and low-confidence transcripts. Isoforms representing less than 1% of gene-level reads are removed unless they match annotated transcripts (reference-guided mode) with non-zero expression. This filter is applied iteratively within each gene: read assignments and isoform fractions are recalculated via EM after each filtering round, and isoforms are sorted by expression. When a gene has more than 10 isoforms (default threshold for aggressive filtering), multiple low-abundance isoforms may be removed in a single round to accelerate convergence. Otherwise, only the lowest-expressing isoform failing the threshold is removed per round. Iteration continues until all remaining isoforms satisfy the minimum expression criteria or no further isoforms can be removed.

Single-exon transcripts are subject to a minimum expression threshold of 1 TPM to mitigate false positives from transcriptional noise. In HiFi mode, single-exon transcripts additionally require TSS or PolyA annotation, as boundary evidence evident from HiFi reads is useful for distinguishing genuine monoexonic genes from transcriptional noise or alignment artifacts.

When one isoform is fully contained within another (identical splicing but shorter terminal exons), collapsing decisions depend on boundary annotations and relative expression. Contained isoforms lacking both TSS and PolyA annotations are automatically removed as likely degradation products. Contained isoforms lacking TSS but sharing the same PolyA site with their containing isoform are removed as redundant 5' truncations. For contained isoforms with distinct boundary annotations, those representing less than 20% of the combined expression of both isoforms are removed as probable incomplete transcripts rather than genuine alternative termini. Importantly, contained isoforms with distinct TSS or PolyA boundaries that exceed this expression threshold are retained as separate isoforms, enabling the discovery of biologically meaningful alternative transcription start sites and polyadenylation sites that produce overlapping but functionally distinct transcript variants.

Isoforms with 3' termini within 10 bp of genomic A-rich sequences (≥7 consecutive adenines) are flagged as potentially internally primed—artifacts from oligo-dT priming within transcripts rather than at genuine polyA tails. Flagged monoexonic isoforms are removed unless they match known 3' ends from a provided reference annotation.

Novel isoforms (those not matching reference annotation) must be supported by at least 2 uniquely assigned reads.

## Isoform Terminal Exon Boundary Refinement

For isoforms lacking explicit TSS or PolyA annotation, corresponding terminal exon boundaries are refined using the empirical distribution of read end positions within terminal exons of supported isoforms. Terminal exon coordinates are adjusted using the following percentiles: the 10th percentile of read positions for the left boundary and the 90th percentile for the right boundary. A minimum of 7 reads must terminate in a terminal exon for percentile-based adjustment to proceed; if insufficient reads are available, the method uses the extreme (minimum/maximum) observed positions of isoform-supported reads.

## Gene Assignment

Reconstructed isoforms are assigned to gene loci using a two-stage clustering approach. Initially, isoforms are grouped by genomic overlap: transcripts sharing exonic overlap (≥1 bp) are merged into preliminary gene components. Within each initial overlap component, Leiden community detection is applied to further refine gene assignments. A binary graph is constructed where transcripts are nodes and edges connect transcript pairs satisfying the overlap thresholds: overlap length must be ≥50% of the shorter transcript and ≥20% of the longer transcript. Leiden clustering is performed with a resolution parameter of 0.2. For extremely large components (>1500 transcripts), Leiden is bypassed and a computationally efficient union-find algorithm is applied instead to define connected components based on the minimum isoform pairwise overlap criteria. This two-stage approach partitions the overlap graph into tightly connected subgraphs (communities) ideally representing distinct genes where isoforms share overlapping exons. Transcripts within the same partition or community are assigned shared gene identifiers.

## Read-to-Isoform Assignment

To quantify isoforms, labeled read paths must be assigned to compatible isoforms. Recall that reads sharing identical splice structure (the same ordered sequence of splice-graph nodes) are aggregated into a single labeled read path with an associated read count—this aggregation, performed during Labeled Read Path Graph construction, enables efficient assignment by evaluating each unique splice pattern once rather than processing individual reads redundantly. Isoforms (reconstructed or provided as reference annotations) are similarly assigned labeled paths through the splice graph based on their structure representation within the splice graph. Each labeled read path is compared against each isoform's labeled read path to establish compatibility. Assignment proceeds through a cascading hierarchy of stringency levels, attempting the most stringent criteria first and progressively relaxing constraints until a compatible isoform is identified.

The assignment cascade operates in two phases reflecting HiFi and non-HiFi data characteristics. HiFi-style assignment (Phase 1) leverages high-accuracy boundary information by requiring TSS and PolyA nodes to match when present. If no assignment is made under HiFi rules, the algorithm falls back to non-HiFi-style assignment (Phase 2), trimming boundary nodes entirely to accommodate reads lacking reliable terminal information. Within each phase, multiple compatibility tests are applied in order of decreasing stringency:

**Phase 1: HiFi-style assignment (TSS/PolyA boundaries respected)**
1. *Exact match with boundary anchoring*: Read path exactly matches isoform path end-to-end. If read or isoform begins/ends with TSS or PolyA nodes, boundaries must match precisely.
2. *Compatible and contained with boundary anchoring*: Read path is fully contained within isoform path as a consecutive segment. Boundaries must match when present.
3. *Introns contained with boundary anchoring*: All read introns appear in isoform in the same order within the read's genomic span, with no extra isoform introns intervening, and ≥75% read overlap with isoform structure. Boundaries must match when present.
4. *Compatible overlap with boundary anchoring*: Read and isoform are compatible (no conflicting introns) across their overlapping region with no gaps in the overlap, and ≥75% read coverage. Boundaries must match when present.
5. *Compatible overlap without boundary anchoring*: Same as above, but TSS/PolyA boundaries need not match even if present.

**Phase 2: Non-HiFi-style assignment (TSS/PolyA boundaries ignored)**
If no assignment is made in Phase 1, TSS and PolyA nodes are trimmed from both read and isoform paths, and compatibility tests are applied without boundary constraints:
6. *Exact match*: Trimmed paths match exactly.
7. *Compatible and contained*: Trimmed read path fully contained in isoform as a consecutive segment.
8. *Introns contained*: All read introns appear in isoform in the same order within the read's genomic span, with no extra isoform introns intervening, and ≥75% read overlap with isoform structure.
9. *Compatible overlap*: Read and isoform are compatible (no conflicting introns) across their overlapping region with no gaps in the overlap, and ≥75% read coverage.

For labeled read paths compatible with multiple isoforms, assignment ambiguity is resolved probabilistically using relative expression levels estimated via expectation-maximization (described below). Labeled read paths are additionally weighted by their 3' end proximity to isoform 3' termini to account for alternative polyadenylation, with weights decreasing linearly as genomic distance increases.

## Isoform Quantification

Isoform abundances are estimated using an expectation-maximization (EM) algorithm that resolves ambiguous read assignments. Read-isoform compatibility follows the cascade described in `Read-to-Isoform Assignment`: the read's labeled path must agree with the isoform's exon-intron structure, with exact matches attempted first and progressively relaxed overlap tests (≥75% coverage) applied when needed. Reads are additionally weighted by their 3' end proximity to isoform 3' termini to account for alternative polyadenylation. For a read compatible with multiple isoforms, the weight assigned to each isoform index $$i$$ is proportional to $$1 - d_i / D$$, where $$d_i$$ is the genomic distance between the read's 3' end and isoform $$i$$'s 3' terminus, and $$D$$ is the sum of distances across all compatible isoforms; explicitly, $$w_{ri} = 1 - d_i / D$$. This linear weighting scheme prioritizes isoforms whose 3' ends align more closely with the read's observed termination point. When the `--oversimplify` option is active for selected contigs, this EM stage is bypassed and reads are instead assigned to the single transcript with the greatest exonic base-pair overlap.

Let $$\mathbf{p} = (p_1, \ldots, p_K)$$ be the isoform proportions for gene $$g$$ (where $$\sum_{i=1}^K p_i = 1$$), and let $$R_g$$ be the set of read compatibility classes (RCCs) assigned to gene $$g$$. Each RCC $$r \in R_g$$ corresponds to a unique labeled read path and aggregates $$n_r$$ individual reads sharing that structure. The likelihood of observing RCC $$r$$ is:
$$P(r \mid \mathbf{p}) = \sum_{i: C_{ri}=1} p_i \cdot w_{ri}$$
where $$C_{ri} \in \{0,1\}$$ indicates compatibility between RCC $$r$$ and isoform $$i$$. To regularize low-support isoforms, a Dirichlet prior with base concentration $$\alpha = 0.01$$ is applied, with per-isoform prior mass scaled by the ambiguous RCC support assigned to that isoform (see below).

In the expectation step, posterior probabilities (responsibilities) are computed per RCC:
$$\gamma_{ri} = \frac{p_i \cdot w_{ri}}{\sum_{j: C_{rj}=1} p_j \cdot w_{rj}}$$

In the maximization step, isoform proportions are updated with Dirichlet smoothing that scales the prior by the amount of ambiguous RCC support assigned to each isoform. Let $$A_i$$ be the total read count summed across RCCs that map to more than one isoform but include isoform $$i$$; then the effective prior mass is $$\alpha_i = \alpha \cdot A_i$$. The update becomes
$$p_i \leftarrow \frac{\sum_{r \in R_g} \gamma_{ri} \cdot n_r + \alpha_i}{\sum_{j=1}^K \left(\sum_{r \in R_g} \gamma_{rj} \cdot n_r + \alpha_j\right)}$$
where $$n_r$$ is the total read count represented by RCC $$r$$ and $$\alpha_i = 0$$ when isoform $$i$$ has no ambiguous RCC support.

Iteration continues until the relative change in log-likelihood falls below $$10^{-6}$$ or maximum iterations are reached (250 for quantification-only mode, 1000 during assembly). From the final proportions, unique read counts (RCCs assigned exclusively to a single isoform with $$\gamma_{ri} \geq 0.9999$$), total fractional read counts, isoform fractions, and TPM are computed.

