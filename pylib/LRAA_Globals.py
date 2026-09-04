import os

## global vars / constants

SPACER = "???"

DEBUG = False

LRAA_MODE = "unset"  # options ("ID", "QUANT-ONLY", "MERGE")

# Identifies how splice-graph coverage normalization was performed. It appears in
# the names of the normalized bam and its work directory, so a cache produced by
# a different method is never mistaken for a current one.
#
# BUMP THIS whenever normalize_bam_by_strand.py changes which reads it keeps or
# what it records on them. Nothing downstream can detect a stale cache on its
# own: a bam from the read-start-binning era carries no XW tag, and an absent
# tag legitimately means "weight 1", so its distorted counts would be consumed
# in silence.
#
#   startbin1  read starts binned per 100 bp, each bin capped, no weights
#   cov1       depth-targeted sampling, scarce junctions kept whole, XW weights
#   cov2       two independent meanings, because the branches below were developed in
#              parallel and both bumped cov1. On devel: the strand split first drops
#              secondary, supplementary, duplicate, qcfail and unmapped records, and any
#              alignment carrying an intron longer than max_intron_length. On the
#              normalization branch: depth and junction support are measured only over the
#              records the consumer actually reads.
#   cov3       normalization branch only -- as its cov2, plus the consumer's min_per_id
#              floor.
#   cov4       normalization branch only -- adds the mapping-quality floor and the
#              improper-pair, duplicate and qcfail rejections.
#   cov5       both of the above, and measurement no longer reimplements the retention
#              rules: it calls Util_funcs.quant_discard_reason, the one policy
#              quantification itself consumes, which additionally excludes unmapped,
#              unplaced and long-intron records. Counting records the consumer discards
#              measured a coverage level nothing downstream has -- 8% of reads on an ONT
#              chr20 bam from the identity floor alone, and mapping quality bites hardest
#              when set, since multimapping reads carry MAPQ 0 at exactly the paralogous
#              loci where thinning decisions matter most. Neither cov2 nor cov4 describes
#              this, so reusing either token would be a stale hit rather than a miss.
#   cov6       both of the above, and the XW weight now COMPOUNDS with any weight the
#              input record already carried instead of overwriting it. Thinning an
#              already-thinned bam composes two acceptance rates, so a record kept at p1
#              and then at p2 stands for 1/(p1*p2); cov5 wrote 1/p2 and discarded the
#              first factor, under-weighting every such record in the splice graph, which
#              honours the tag unconditionally. A cov5 artifact of an untagged input is
#              byte-identical under cov6 -- absent tag means weight 1 -- but its name
#              cannot say whether its input was tagged, so the token has to change.
SPLICE_GRAPH_NORMALIZATION_METHOD = "cov6"

# The identity floor --HiFi imposes, as ONE definition. Two processes resolve it:
# LRAA's own preset (_apply_hifi_config_overrides), and pylib/ChunkedRun.py when
# it is invoked directly rather than through the driver -- which is what a
# by_chunk WDL task does. Those two disagreeing is not a smaller answer: prep
# selected cuts, priced severed reads and normalized at 80 while every stage-5
# worker filtered at 97, and a shared-cut-plan consumer refused the plan its
# LRAA-driven emitter had selected. So the value lives here rather than as a
# literal in either caller.
HIFI_MIN_PER_ID = 97.0

config = {
    #########################
    # read alignment criteria
    "HiFi": False,  # set to True when --HiFi is used; enables HiFi-specific filtering
    # The non-HiFi floor. --HiFi raises it to HIFI_MIN_PER_ID above.
    "min_per_id": 80,
    "min_mapping_quality": 0,  # used during isoform discovery; lets multi-mapping reads (mapq=0) inform splice-graph and isoform structure (e.g., paralog-cluster genes)
    "min_mapping_quality_for_final_quant": 0,  # default to retaining MAPQ 0 alignments during final quant; callers can raise this threshold if desired
    # CPU budget, and the two things the old "num_threads_per_worker" conflated. The
    # budget is the total for the invocation; "tool_threads" is what a unit worker may
    # pass to a native tool (samtools -@, minimap2 -t); "component_workers" is how many
    # processes it may fork for large multipath-graph components, 0 meaning none. All
    # three are derived from --cpu_budget by pylib/CpuBudget.py, never set independently.
    "cpu_budget": 1,
    "tool_threads": 1,
    "component_workers": 0,
    "try_correct_alignments": True,
    "max_softclip_realign_test": 20,
    "min_softclip_realign_test": 5,
    "min_frac_alignments_pass_per_id_check": 0.9,
    "min_total_alignments_engage_frac_per_id_check": 1000,
    "read_aln_gap_merge_int": 10,
    "max_intron_length": 200000,
    # rDNA-cassette masking (see pylib/RdnaMask.py). Reads whose alignment overlaps
    # a hit of the mask fasta against --genome are excluded everywhere
    # quant_discard_reason is consulted -- coverage normalization, splice-graph
    # construction, and read-to-transcript assignment alike. On by default because
    # a run that hits one of these loci pays the cost regardless of whether the
    # investigator was looking for it; --no_rdna_mask opts out.
    "rdna_mask_enabled": True,
    # None -> RdnaMask.DEFAULT_RDNA_CASSETTE_FASTA (bundled human rDNA repeat unit,
    # resources/human_rDNA_cassette.fa); set from --rdna_mask_fasta for a different
    # organism's cassette.
    "rdna_mask_fasta": None,
    # Bases of clearance added on each side of a cassette-vs-genome alignment hit
    # before it becomes an excluded region, absorbing alignment-boundary slop
    # (indels, soft-clips) at the edge of a real rDNA-homologous span.
    "rdna_mask_pad": 500,
    # A cassette-vs-genome hit must clear BOTH floors to become part of the mask,
    # or a single short, coincidentally-homologous alignment -- unavoidable at
    # genome scale -- would exclude a region indistinguishable from a real
    # rDNA-repeat-unit copy. Real copies observed on the reference genomes tested
    # span 1-44 kb at effectively full identity, so both floors sit far below any
    # genuine hit and only ever reject noise. See RdnaMask._sam_hit_spans.
    "rdna_mask_min_hit_length": 200,
    "rdna_mask_min_per_id": 80,
    # Minimum overlap, in bases, between a read's alignment and the mask before
    # the read is discarded (see RdnaMask.read_overlaps_mask). Every excluded
    # region already carries rdna_mask_pad bp of padding specifically to absorb
    # boundary slop around a genuine hit, so a read that only grazes that padding
    # is far more likely an ordinary read from adjacent unique sequence than one
    # implicated in the locus's multi-mapping ambiguity; a read genuinely inside
    # a masked repeat copy overlaps by its whole aligned length and clears this
    # trivially, so the floor only ever spares boundary-adjacent reads.
    "rdna_mask_min_overlap_bp": 50,
    # NOT a CLI setting. Populated once per LRAA invocation, immediately after
    # --genome is resolved, from RdnaMask.build_rdna_mask_bed +
    # RdnaMask.load_mask_bed: a {contig: IntervalTree} mask, or None when masking
    # is disabled or found nothing for this genome. quant_discard_reason reads it
    # here so every consumer converges on one built-once mask without threading it
    # through every call site by hand.
    "rdna_mask_intervals": None,
    #
    ####################################
    # splice graph construction criteria
    "min_SE_read_ME_exon_overlap_pct": 50,  # min % of SE read length that must overlap with ME exon to filter out SE read
    # default tuned for not HiFi (e.g., ONT); HiFi mode overrides to 0.01 via --HiFi
    "min_alt_splice_freq": 0.03,
    "min_alt_unspliced_freq": 0.01,
    "min_feature_frac_overlap": 0.50,
    "max_exon_spur_length": 14,  # maximum terminal exon spur length; HiFi sets 13
    "aggregate_adjacent_splice_boundaries": True,
    "aggregate_splice_boundary_dist": 5,
    "aggregate_splice_boundary_max_rel_support": 0.2,  # inclusive; collapse only when alt/top support <= this. 1.0 = legacy unconditional collapse
    "fracture_splice_graph_at_input_transcript_bounds": True,
    "max_path_nodes_per_component": 1000,  # max number of path graph nodes per connected component
    # transcript reclustering (gene definition) criteria
    # gene reclustering overlap thresholds
    "min_recluster_overlap_shorter_iso_frac": 0.50,  # (overlap_len / shorter_transcript_len) >= this to connect isoforms in second-stage graph
    "min_recluster_overlap_longer_iso_frac": 0.20,  # also require (overlap_len / longer_transcript_len) >= this to avoid linking large multi-exon to long single-exon with tiny shared portion
    # shared splice junctions as gene evidence. Denominator is the SMALLER intron
    # count of the pair, so a fragmentary model is judged on how much of its OWN
    # splice pattern agrees. Both conditions must hold. Isoforms with identical intron
    # chains are one gene unconditionally and are subject to NEITHER threshold.
    "min_recluster_shared_intron_frac": 0.20,  # shared_introns / min(intron count) > this to connect isoforms
    # Floor on the shared COUNT. The fraction alone is met by ONE shared junction
    # whenever the smaller isoform has <= 4 introns, which is enough for a single
    # fragment to bind two neighbouring genes together. Measured over four contigs,
    # raising this from 1 to 2 gave up 16 of 42 recovered gene ids and avoided 8 of 10
    # added fusions -- an exchange rate of 13:1 rather than 4.2:1.
    "min_recluster_shared_introns": 2,
    # community clustering (Leiden) for transcript→gene reassignment
    "use_community_clustering": True,     # enabled by default; use Leiden communities within initial clusters
    "community_resolution": 0.2,          # Leiden resolution parameter (higher → more, smaller communities)
    "community_random_seed": 42,          # seed for deterministic Leiden partitions
    # safety valve for very large overlap components: skip community clustering when too large
    "max_transcripts_for_community_clustering": 1500,  # if an initial cluster exceeds this size, fall back to lightweight overlap-based DSU reclustering
    #
    ############
    # TSS config
    "infer_TSS": False,  # include TSS feature in read path assignments (HiFi enables)
    "max_dist_between_alt_TSS_sites": 50,
    "min_alignments_define_TSS_site": 5,
    "max_soft_clip_at_TSS": 0,
    # Length of a leading untemplated-G run to strip before judging a TSS. Reverse
    # transcriptase adds these opposite the cap during template switching, so they
    # mark a genuine transcript start; with max_soft_clip_at_TSS at 0 they instead
    # disqualify the read. Measured on chr20: 83.8% of primary alignments are clipped
    # at their 5' end, the first clipped base is G in 99.9%, 96.5% are pure G runs of
    # three or fewer, and none of 257,880 clipped bases matches the reference beyond
    # the alignment. Enabled by default on that evidence.
    #
    # Honest caveat for whoever tunes this next. On chr20 de novo the strip cost 13
    # true chains and gained none; 13 of 17 lost chains had a strict subset emitted in
    # their place, so the mechanism of the loss is truncation. Untemplated G's mark
    # where reverse transcription terminated rather than where the cap is, and for a
    # degraded transcript RT stops internally and still adds them, so some admitted
    # ends are internal. The biology justifies stripping; the chromosome-scale chain
    # count did not. 0 disables the stripping and restores the pre-0.18.3 behaviour.
    "max_untemplated_G_at_TSS": 3,
    "min_TSS_iso_fraction": 0.05,  # during initial TSS definition, require for a 'gene' that a TSS has at least this fraction of TSS-candidate gene reads assigned.
    "TSS_window_read_enrich_len": 50,
    "TSS_window_read_enrich_factor": 5,
    #
    ## - alt TSS isoform pruning
    # during splice graph construction: walking exon segments from a more dominant site, removing less supported sites below fraction of dominant
    # during isoform resolution: comparing isoform i that contains j, j >= this frac of i TSS read support
    "max_frac_alt_TSS_from_degradation": 0.20,
    # to retain j TSS when comparing to i TSS, j TSS must have >= read support fraction of all gene reads
    #
    ####################
    ## polyA site config
    "infer_PolyA": False,  # include PolyA site feature in read path assignments (HiFi enables)
    "max_dist_between_alt_polyA_sites": 50,
    "min_alignments_define_polyA_site": 5,
    "min_frac_alignments_define_polyA_site": 0.1,
    "min_PolyA_ident_length": 7,  # examine softclipped ends of reads, if have polyA with at least this number of bases at terminus, strip it and extended match out
    "min_PolyA_iso_fraction": 0.05,  # during initial TSS definition, require for a 'gene' that a TSS has at least this fraction of polyA-candidate gene reads assigned..
    "max_soft_clip_at_PolyA": 0,  # max amount of softclipping allowed at the end of an alignment to mark it as a candidate boundary
    "min_soft_clip_PolyA_base_frac_for_conversion": 0.8,  # if soft-clipped is at least this frac polyA evidence, then removing soft clipping and marking as candidate polyA read.
    #
    ####################
    ## Terminal boundary definition
    "terminal_boundary_method": "percentile",  # choices: "extreme" (min/max), "mean", "median", "quartile" (Q1/Q3), "percentile" - method for defining terminal coords when TSS/PolyA not annotated
    "terminal_boundary_percentile": 90,  # when terminal_boundary_method is "percentile", use this percentile (e.g., 90 means 10th percentile for left, 90th for right)
    # Minimum read count required for mean/median/quartile/percentile adjustment; below this threshold, the existing boundary is retained.
    # Consider whether a future implementation should instead fall back to the observed extreme (minimum/maximum) read position.
    "min_reads_for_terminal_adjustment": 7,
    #
    # compatible and contained isoform filtering
    "max_rel_frac_expr_alt_compat_contained": 0.2,  # if iso-j contained by iso-i has < this frac of their combined expression, iso-j gets pruned
    #
    ## read assignment to transcript criteria
    "fraction_read_align_overlap": 0.75,  # min fraction of read length that must overlap the compatible transcript isoform structure
    #
    # misc settings
    "min_path_score": 1,  # min number of reads required for reporting isoform
    #
    # transcript filtering criteria
    "min_transcript_length": 200,
    "min_isoform_fraction": 0.01,
    "min_frac_gene_unique_reads": 0.01,  # minimum fraction of all uniquely assigned reads per gene
    "min_monoexonic_TPM": 1.0,
    # A multi-exonic model is kept only if quantification gave it expression above this
    # value. The default of 0 means "any expression at all", which is the judgement
    # ref_trans_filter_mode=retain_expressed asks for: supplied models are selectable
    # from the trellis on their synthetic template read, then this decides whether they
    # were actually expressed. Raise it to demand more than a trace; set it negative to
    # disable the check and report every selected multi-exonic model.
    "min_multiexonic_TPM": 0.0,
    # A monoexonic model has no intron chain to corroborate it, so its only structural
    # evidence is that its reads describe one contiguous thing. Reads that tile a long
    # span without overlapping each other describe a covered region, not a transcript:
    # 500bp reads can never establish a 20kb unspliced isoform. This is the fraction of
    # a model's supporting reads that must mutually overlap at some single base
    # (i.e. peak read depth / supporting read count). Self-calibrating: long reads and
    # genuinely stacked support pass regardless of model length. 0 disables the check.
    "min_monoexonic_read_span_peak_frac": 0.5,
    # Minimum ratio of a monoexonic model's coverage-depth ("adjusted") TPM to its
    # read-count TPM. Equivalently, the mean fraction of the model an individual
    # supporting read covers. Reads that each span the model give a ratio near 1;
    # reads that tile it drive the ratio toward 1/read_count. Scale-free, so it
    # measures agreement rather than abundance. 0 disables the check.
    "min_monoexonic_adjusted_TPM_ratio": 0.20,
    # Single-cell only. Minimum number of distinct cells that must contribute a read
    # to a NOVEL monoexonic model; monoexonic models containing a reference model are
    # exempt under ref_trans_filter_mode retain_expressed, and bulk input carries no
    # barcode in its read names so the check self-disables there. An absolute count
    # rather than a fraction of the cluster: measured across 14 PBMC clusters (122 to
    # 1,506 cells), a fraction's stringency swung with roster size, while recovery of
    # reference-matching monoexons against an absolute bar was stable -- 98% at 3
    # cells, 92% at 5, 73% at 10. 0 disables.
    "min_monoexonic_supporting_cells": 5,
    # Internal priming is rejected during PolyA site identification
    "filter_internal_priming": True,
    "restrict_internal_priming_filter_to_monoexonic": True,
    # When True, a monoexonic transcript that looks internally primed is retained if its
    # 3' end agrees with a reference annotation 3' end -- proximity to a known_transcripts
    # terminus, not to any measured cleavage atlas. Off by default: a monoexonic model has
    # no intron chain corroborating it, so agreement alone is weaker evidence there.
    "spare_monoexonic_internal_priming_with_known_3prime": False,
    # Internal-priming veto at PolyA site identification: when a READ-DERIVED candidate
    # sits at a 3' end the supplied reference annotation also calls, the reference is
    # independent evidence that cleavage happens there, so the A-rich context veto is
    # waived.  Inert without a reference -- ref-free runs have no known 3' ends -- so
    # this only ever loosens a ref-guided run.
    #
    # "Also calls" means within max_dist_between_alt_polyA_sites / 2 (25 nt), not an
    # exact coordinate match: the candidate coordinate is the most-supported read end in
    # a 50 nt aggregation window, so it is not base-precise the way the annotation is,
    # and this is the window the transcript-level reprieve
    # (spare_monoexonic_internal_priming_with_known_3prime, and the multi-exonic rule
    # above it in TranscriptFiltering) already uses for the same question.  Measured
    # exposure on chr20: of 3,138,689 positions both strands where the veto would fire,
    # an exact match waives 193 and +/-25 waives 7,644.
    #
    # Affects graph construction, hence registered in _SPLICE_GRAPH_CONFIG_KEYS.
    "spare_polyA_veto_at_known_3prime": True,
    "ref_trans_filter_mode": "retain_expressed",  # choices ["retain_expressed", "retain_filtered"]
    "min_reads_novel_isoform": 2,
    "min_unique_reads_novel_isoform": 2,
    "min_isoform_count_aggressive_filtering_iso_fraction": 10,  # allow for filtering mult isoforms in a single round if more than this number of isoform candidates.
    # Reads matching an isoform's intron chain exactly -- a full splice match --
    # are direct evidence for that whole structure. Where a splice pattern has at
    # least this many, the last isoform carrying it is kept rather than filtered
    # on isoform fraction or unique-read fraction.
    #
    # Both of those are relative to the gene, so a minor isoform of a deeply
    # sequenced gene must clear a bar that rises with the gene's expression: at
    # 3,324 gene reads the 1% unique-read floor asks for 33, and an annotated
    # SEC11A isoform with 21 reads carrying its exact chain was dropped for having
    # 0.63%. This is an absolute count so that direct evidence of a structure
    # cannot be outvoted by the depth of its neighbours.
    #
    # Only the last carrier is spared, so terminal variants of one chain do not
    # all survive on the strength of the reads they share. 0 disables.
    "min_FSM_reads_retain_isoform": 0,
    # Substitute an absolute count of full-splice-match reads for the relative
    # unique-read fraction when deciding a model is too weakly supported. The
    # default 0 keeps the fraction. On chr20 ref-guided a gate of 2 moved
    # precision 0.329 -> 0.363, dropping 268 false chains for 10 true ones, where
    # tightening the fraction to 0.02 reached the same precision at a cost of 58
    # more true chains -- the gain is in the quantity, not the cut.
    "min_FSM_reads_gate": 0,
    #
    ##########
    # assembly
    "normalize_max_cov_level": 1000,
    "restrict_asm_to_collapse": True,  # if True, no chaining of overlapping/extended paths
    #
    ###############################################
    # chunked parallelism: cutting a contig-strand
    #
    # ON by default: the orientation split runs inside each chunk, concurrently with
    # every other chunk, rather than as a serial pass over the whole bam first.
    # MEASURED 151.2 s against strand-first's 255.2 s on the same input, because the
    # whole-genome split is the single largest serial phase a chunked run has and
    # strandless does not have it at all. The two orientations of an interval then
    # share one extraction -- one mini FASTA, one mini GTF, one pass over the region.
    #
    # Stages 4 and 5 are identical either way: each still receives one
    # orientation-pure bam for one chunk. Strandless changes WHERE the split happens,
    # not whether reads are processed per strand.
    #
    # Opt out with --chunk_by_strand, which restores the strand-first ordering.
    # There is no correctness reason to: chunked-vs-unchunked parity is measured in
    # both modes. It exists so a regression can be bisected against the older path.
    "strandless_chunks": True,
    #
    # A contig-strand is split into chunks that are normalized and processed
    # independently, then merged. Both values below are in MEGABASES.
    #
    # Target cut positions sit at multiples of this across the contig, so a
    # contig gets roughly length / approx_MB_per_cut chunks with no cap: chr20
    # (64.4 Mb) gets 6, chr1 (248.9 Mb) gets 25. Sizing is by SPAN, not by
    # alignment or gene count, so a chunk's coordinates are predictable from the
    # contig length alone.
    "approx_MB_per_cut": 10,
    # TOTAL width of the search window centred on each target, in megabases --
    # the whole window, not the half-width. A target at T is searched over
    # [T - wiggle/2, T + wiggle/2], so the default 1 means T +/- 0.5 Mb. This is
    # the MAXIMUM radius, not the radius searched: the selector starts small and
    # widens progressively, stopping as soon as a compliant position severs
    # nothing.
    #
    # ABSOLUTE, and deliberately not derived from approx_MB_per_cut. THE
    # MISREADING TO PREVENT is making this proportional to the spacing, which
    # looks tidy because 10% of the shipped 10 Mb spacing is exactly this 1.
    # Measured on HG002 PacBio Kinnex at 2 Mb spacing, de novo, counting the
    # retained primary alignments the cuts sever:
    #
    #     contig | 20 kb window | 200 kb window | 1 Mb window
    #     chr21  |          940 |           743 |          0
    #     chr1   |         2598 |            87 |          0
    #
    # A proportional rule would give 200 kb at 2 Mb spacing, which still severs
    # 743 alignments on chr21 where this absolute 1 Mb severs none. And chr1 and
    # chr21 disagree 8.5-fold at identical parameters -- 87 against 743 -- because
    # chr1 offers 6,258 read-free gap runs to chr21's 1,295. The distance a search
    # must travel is a property of the sequence and the library in BASES: the
    # closest zero-cost grid position sits up to 382.8 kb from a chr21 target and
    # 348.3 kb from a chr1 one. It tracks the GENOME, not the geometry, so it
    # cannot be expressed as a fraction of the chunk.
    #
    # Callers testing a finer spacing may of course pass a smaller window, and a
    # run that then severs reads is not a malfunction: the selector still places
    # the best position it can reach and reports what that cost.
    "approx_MB_per_cut_wiggle_window": 1,
    # What a severed MULTI-EXON alignment costs cut selection, against 1 for a
    # monoexonic one. Severing is a cost to minimise and never a veto -- at depth
    # every base is covered, so a rule forbidding it would decline every cut --
    # but the two are not worth the same: a spliced alignment carries junction
    # evidence, which is what the splice graph's edges are built from, while a
    # monoexonic one carries none. 10 makes one severed junction-bearing read
    # outweigh nine severed monoexonic ones.
    "chunk_severed_multiexon_weight": 10,
    #
    # The remaining chunking constants. Each of these previously existed as two
    # to four independent copies across LRAA, ChunkedRun, select_contig_cut_points
    # and normalize_bam_by_strand, all of which happened to agree. A divergence
    # would not have surfaced as a mismatch: the values are baked into the
    # stage-2 cache token, so one copy moving turns a stale cache entry into a
    # HIT that reuses the old geometry while asserting the new parameters.
    #
    # Resolution in bases at which read depth is measured when scoring candidate
    # cut positions, and the grid the normalizer's depth windows sit on. The two
    # must be the same number or normalization thins differently on either side
    # of a boundary.
    "chunk_depth_window": 100,
    # Bases of clearance a cut must leave on both sides of every annotated locus.
    # 4x the largest boundary-snapping distance in this file (50 bp:
    # max_dist_between_alt_TSS_sites, max_dist_between_alt_polyA_sites,
    # TSS_window_read_enrich_len), so no snapping can reach across a cut.
    "chunk_margin": 200,
    # Absolute reference coordinate the depth-window grid is anchored to, so the
    # same locus lands in the same window whether it is normalized whole or as
    # part of a chunk.
    "chunk_grid_origin": 0,
    # Seed for the normalizer's reproducible down-sampling.
    "chunk_random_seed": 42,
    #
    #######
    # quant
    "num_total_reads": None,  # for TPM and filtering - set by CLI or within LRAA by counting bam records
    "run_EM": True,
    "max_EM_iterations_quant_only": 250,  # don't set too high, as even at 1000 small biases get greatly amplified.
    "max_EM_iterations_during_asm": 1000,  # for asm, want higher iterations to amplify small diffs and weed out poorly supported isoforms.
    "aggressively_assign_reads": False,
    "rescue_unassigned_reads_via_transcriptome_alignment": True,
    "rescue_unassigned_minimap2_preset": "auto",
    "rescue_unassigned_minimap2_filter_fraction": 0,
    # Fraction of a read's length that must align to the target transcript before the
    # alignment can be accepted as rescue evidence. Measured as aligned length over
    # read length (clipping excluded from the numerator), not as matched bases, so
    # platform error rates do not make it unsatisfiable -- mismatches are bounded
    # separately by rescue_unassigned_min_per_id. A partial alignment describes a read
    # that only locally resembles the target and must not count as support for the
    # whole isoform. 0 disables the check.
    "rescue_unassigned_min_aligned_read_frac": 0.95,
    "rescue_unassigned_min_per_id": None,
    # Longest insertion or deletion tolerated in a transcriptome rescue alignment.
    # A read that skips or inserts this many bases relative to the target disagrees
    # with it structurally rather than by sequencing error, so the alignment is
    # declined and the read keeps its genome alignment.
    #
    # Calibrated on chr20 reads whose genome intron chain exactly matches an
    # annotated transcript, realigned to that transcript's cDNA: for correctly
    # placed reads the largest deletion observed was 32 (PacBio HiFi, p99.9 = 15)
    # and 45 (ONT cDNA, p99.9 = 37); insertions reach p99.9 of 12 and 13. The
    # defaults sit below those maxima deliberately -- rejecting ~0.3% of correctly
    # placed reads on either platform buys a bar low enough to catch exon-sized
    # disagreement, which error-calibrated caps would not. The HiFi block lowers
    # this to 10. 0 disables the check.
    "rescue_unassigned_max_indel_length": 30,
    # When True, weight ambiguous read assignments by agreement of read 3' ends with transcript 3' ends
    # (previously "use_weighted_read_assignments" which weighted by both 5' and 3' ends)
    "weight_reads_by_3prime_agreement": True,
    # XW coverage-normalization weights are honoured unconditionally and have no setting.
    # A weight is present exactly where thinning happened, and an untagged read weighs 1
    # (Pretty_alignment.get_normalization_weight), so honouring the tag is a no-op on a bam
    # nobody thinned -- which makes weighting a property of the DATA rather than a mode.
    # The input roles are what guarantee that: --bam must be the full library and
    # --bam_for_sg must already be normalized, both checked in LRAA's setup.
    #
    # A single pass may still opt out, via _populate_read_multi_paths(weight_reads=False).
    # Discovery's pre-filter quantification does, because its isoform gates mix EM-derived
    # quantities that follow a weight with integer tallies that cannot.
    #
    # One acceptance probability per read is a precondition, and it holds by
    # construction rather than by a check: alignment intake discards secondary and
    # supplementary records unconditionally, so a read reaches weighting as at most
    # one record.
    # Diagnostic dumps for evaluating a streaming assignment pass, all off unless set to
    # an output prefix via --config_update. They must exist here with defaults or
    # --config_update rejects them as unknown keys.
    #
    # dump_read_path_map      read name -> the canonical path chosen to represent it
    # dump_mp_fraction_table  canonical path -> the fractional split over transcripts
    # dump_rescue_candidates  read name -> why it is a rescue candidate, written by the
    #                         batch path and by the streaming path under distinct names so
    #                         the two populations can be diffed read for read. Counts alone
    #                         would let two different sets of the same size look equal.
    #
    # The first two are keyed on canonical paths (feature type plus genomic coordinates)
    # rather than node or multipath ids, since those are process-global counters and drift
    # between runs.
    "dump_read_path_map": None,
    "dump_mp_fraction_table": None,
    "dump_rescue_candidates": None,
    # Two-pass alternative to the default final quantification, ON by default since
    # v0.25.0 (see --no_stream_reads). The first pass quantifies normally against
    # the coverage-normalized bam; the second streams the full bam, looks each
    # read's path up in the table the first pass produced, writes its tracking row
    # and forgets it. Nothing per-read is retained, which is what makes a
    # billion-read library tractable -- the non-streaming path holds each shard's
    # alignments and every read name it will report.
    #
    # This reports one expectation step at the first pass's abundances, where the
    # non-streaming path re-estimates them, so its counts are close to but not
    # identical with the non-streaming path's. --no_stream_reads reverts to the
    # pre-v0.25.0 single-pass, in-memory behaviour.
    "stream_reads": True,
    # No key here bounds how much of a streamed unit the first pass's table has to answer.
    # How much it DID answer -- the served fraction -- is reported per contig-strand by
    # StreamingQuant and gated on by nothing. A max-unseen-path-read-fraction tripwire used
    # to sit here at 0.25 and was enforced after the streaming loop returned, by which point
    # every read had been mapped, looked up, written to the tracking file and dropped: it
    # refused a complete and correct output for having been slow, and no threshold value
    # makes that right. It was also blind to the failure that matters -- a miss resolves
    # with the same theta the first pass would have used, so if pass 1 was starved theta is
    # unreliable for the cache HITS too, and a high hit rate over unusable abundances passed
    # in silence. A correctness gate belongs on pass 1, asking whether its abundance
    # estimates are usable, BEFORE pass 2 spends the work; not implemented here.
    # Rescue candidates against the local transcriptome from inside the streaming pass,
    # using a resident mappy index instead of the batch path's minimap2 subprocess. This
    # bare key is only the FALLBACK baseline for callers that read it directly (e.g.
    # ChunkedRun.py's own standalone parser); LRAA's own CLI does not read it as a static
    # default any more. Since --stream_reads is on by default and transcriptome rescue is
    # on by default, forcing this flat False would make every default invocation refuse
    # itself (--stream_reads requires transcriptome rescue turned off, unless
    # --stream_reads_rescue_unassigned is given -- see LRAA's guard). So the CLI resolves
    # --stream_reads_rescue_unassigned's default dynamically, to whatever transcriptome
    # rescue itself resolves to, unless the caller states either flag explicitly. See
    # --no_stream_reads_rescue_unassigned to opt out even though rescue stays on
    # elsewhere.
    #
    # The candidate population is exactly the one the batch path collects at its three
    # gated sites -- reads the extractor discarded for low_perID, reads whose graph path
    # contains a spacer, and reads with no graph path -- so the two paths target the same
    # reads. Measured identical on ONT chr20, read for read: 14,455 of 120,370 records.
    # The batch path's fourth category is deliberately NOT included here; see the key
    # below. Outcomes may still differ: mappy exposes no equivalent of minimap2's -f, and
    # no alignment score, so best-hit ranking falls back to matched-minus-NM. See
    # pylib/StreamingRescue.py.
    "stream_reads_rescue_unassigned": False,
    # Extend streaming rescue to the fourth candidate category: reads that DID map to a
    # graph path, but whose path matched no target. Off by default, and behind its own
    # flag rather than the one above, because this category is the one place the two
    # paths cannot target the same reads. The batch path derives it from its own first
    # pass; under --stream_reads that first pass reads the coverage-normalized bam while
    # the stream reads the full one, so the streaming population is a strict superset.
    # Measured on ONT chr20: batch 3,442, streaming 11,196, batch-only 0. Enabling it
    # therefore rescues against a larger candidate set than the batch path would, which
    # is a deliberate extension rather than a reproduction of it.
    "stream_reads_rescue_unassigned_to_targets": False,
    "EM_alpha": 0.01,  # regularization
    "EM_convergence_tol": 1e-6,  # L2 change in normalized abundances; shared by both EM passes
    # assignment fraction at or above which a read counts as uniquely assigned.
    # Reporting requires a whole read; isoform filtering tolerates EM rounding.
    "unique_read_report_min_frac": 1.0,
    "unique_read_filter_min_frac": 0.9995,
    # low-memory tuning knobs (now implicit defaults: always avoid in-memory read-name storage; always track spans)
    #
    ######
    # single cell
    "cell_barcode_tag": "CB",
    "read_umi_tag": "XM",
    ######
    # parallelization
    #
    "min_mpgn_component_size_for_spawn": 150,
    "no_cleanup": False,
    ######
    # resource monitoring
    "resource_monitor_enabled": True,
    # How often a row is WRITTEN. RSS is sampled far more often than this and
    # each row carries the high-water mark over its interval, so the peak does
    # not depend on this value; it only decides how much shape the time series
    # has. 60 s used to be the sampling rate too, against a median work unit of
    # ~53 s, which left most units with a single spot reading.
    "resource_monitor_interval": 15.0,  # seconds
    "resource_monitor_include_children": True,
    ######
    # progress monitoring
    # read mapping to graph stage
    "show_progress_mapping": True,       # emit progress while mapping read alignments to the splice graph
    "mapping_update_every_n": 10000,     # fallback: update every N reads processed
    "mapping_update_interval_sec": 2.0,  # fallback: or at least this often in seconds
    # logging cadence for mapping stage (separate from stderr progress); set None or 0 to disable
    "mapping_log_progress_interval_sec": 30.0,
    # splice-graph population logging cadence (coverage + intron scan); set None or 0 to disable
    "splice_graph_log_progress_interval_sec": 30.0,
    # component/timing instrumentation (post itree validation)
    "log_splice_graph_component_timing": True,  # emit timing/memory stats around connected component discovery
    "log_splice_graph_merge_progress_interval_sec": 120.0,  # optional interval (sec) for progress during exon segment merging (0 disables)
    "log_splice_graph_debug_counts": True,  # log node/edge counts at key refinement checkpoints
    # coverage reset progress (recompute base coverage from pretty alignments)
    "show_progress_cov_reset": True,          # show progress while recomputing base coverage
    "cov_reset_update_every_n": 5000,         # fallback: update every N alignments processed
    "cov_reset_update_interval_sec": 2.0,     # fallback: or at least this often in seconds
    # input transcript integration progress
    "show_progress_integrate_transcripts": True,
    # quant: assign reads to transcripts stage
    "show_progress_quant_assign": True,  # emit periodic progress updates during read->transcript assignment
    "use_tqdm_progress": True,           # if tqdm is available, prefer tqdm-based progress bar
    "progress_update_every_n": 1000,     # update every N multipath-count pairs processed (set None to disable count-based updates)
    "progress_update_interval_sec": 5.0, # or at least this often in seconds (set None to disable time-based updates)
    # isoform reconstruction progress (selection of best transcript paths within large components)
    # emit periodic INFO logs while iterating scored paths; set interval <=0 to disable
    "iso_recon_progress_interval_sec": 120.0,
    # emit progress every N path iterations (in addition to time-based interval); set <=0 to disable
    "iso_recon_progress_every_n": 250,
    # pruning phases (splice graph refinement) optional progress intervals (sec); set <=0 to disable
    "prune_introns_progress_interval_sec": 60.0,
    "prune_unspliced_exons_progress_interval_sec": 60.0,
    # finalize splice graph (interval tree + node indexing) progress logging
    # emit periodic logs while populating interval trees for very large graphs
    # time-based interval (sec); set <=0 to disable
    "finalize_splice_graph_progress_interval_sec": 60.0,
    # count-based logging: every N nodes processed (set <=0 to disable)
    "finalize_splice_graph_progress_every_n": 10000,
    # connected component discovery progress (second pass included)
    # time-based interval (sec); set <=0 to disable
    "cc_discovery_progress_interval_sec": 60.0,
    # count-based logging: every N nodes considered (set <=0 to disable)
    "cc_discovery_progress_every_n": 25000,
    # TSS pruning progress controls
    "tss_prune_progress_interval_sec": 60.0,
    "tss_prune_progress_every_n": 500,
    # PolyA pruning progress controls
    "polya_prune_progress_interval_sec": 60.0,
    "polya_prune_progress_every_n": 500,
    # node→component assignment progress controls
    "component_assign_progress_interval_sec": 60.0,
    "component_assign_progress_every_n": 20000,
    # multipath graph build progress
    "mp_graph_build_progress_interval_sec": 60.0,
    "mp_graph_build_progress_every_n": 10000,
    # multipath component discovery progress
    "mp_component_discovery_progress_interval_sec": 60.0,
    "mp_component_discovery_progress_every_n": 5000,
    # multipath pruning progress (large component removal, node pruning)
    "mp_prune_progress_interval_sec": 60.0,
    "mp_prune_progress_every_n": 10000,
    ######
    # disk-backed storage backend for read tracking stores
    # choices: 'auto' (prefer lmdb if available, else sqlite), 'lmdb', 'sqlite', 'memory'
    # default changed to 'memory' for faster runs when persistence is unnecessary
    "store_backend": "memory",
    ######
    # oversimplify (best-overlap) mode
    # When enabled via CLI --oversimplify <contig[,contig2,...]>, specified contigs in quant-only runs
    # will bypass graph/EM and assign each read to the single best-overlapping reference transcript.
    "oversimplify_enabled": False,
    "oversimplify_contigs": [],  # list of contig names (e.g., ["chrM", "MT"]) to treat with simplified assignment
    # Polyadenylation signal annotation.  Defaults are human: the two canonical hexamers
    # and the transcript-sense window LRAA's own PAS analyses used.  Both are settable
    # for other organisms -- plant and many invertebrate signals are more degenerate and
    # sit at different spacings -- via --polyA_signal_motifs / --polyA_signal_window.
    # All motifs must share one length, since the containment bound is derived from it.
    # Purely annotation: these affect the PAS and PAS_offset GTF attributes and nothing
    # else, so they are not part of the splice-graph cache key.
    "polyA_signal_motifs": ["AATAAA", "ATTAAA"],
    "polyA_signal_window": [-40, -10],
}


def resolve_min_polya_iso_fraction(
    min_isoform_fraction, min_polya_iso_fraction_override=None
):
    if min_polya_iso_fraction_override is None:
        return min_isoform_fraction

    return min_polya_iso_fraction_override

# Default read-store backend: favor in-memory unless caller overrides later (CLI/env).
if "LRAA_READSTORE_BACKEND" not in os.environ:
    try:
        os.environ["LRAA_READSTORE_BACKEND"] = str(config.get("store_backend", "memory"))
    except Exception:
        pass

# Global, per-run external stores for read tracking (set at runtime by entry script)
# When set, MultiPath.get_read_names() can stream read names via these stores even when
# in-memory read name retention is disabled.
READ_NAME_STORE = None  # type: ignore
MP_READ_ID_STORE = None  # type: ignore

# Per-run map from compact read ID to the coverage-normalization weight the read
# carries -- the reciprocal of its acceptance probability, taken from the bam's XW
# tag. Summing these rather than counting reads is what recovers the support an
# unnormalized bam would have shown, so quantification must consult it wherever a
# multipath's read tally stands in for abundance.
#
# Keyed on the ID rather than the name because MultiPath._coerce_read_identifier is
# a pure function of the read name: every path that rebuilds a multipath -- the
# genome pass, transcriptome rescue after minimap2 has discarded the tag, and any
# split or clone -- arrives at the same key without having to carry the weight
# along. A read absent from the registry weighs 1, so an unnormalized bam and a bam
# predating the tag both quantify exactly as they did before.
READ_WEIGHT_REGISTRY = {}  # type: ignore


def reset_read_weight_registry():
    """Drop every recorded weight. Call when a quant pass starts reading a bam.

    Weights belong to one bam. A run can pass over several -- the splice graph reads
    the normalized bam while quantification reads the original, and discovery quants
    twice -- and a weight surviving from one into another would be applied to reads
    that never carried it, silently scaling support in the arm that is supposed to be
    the control.
    """
    READ_WEIGHT_REGISTRY.clear()


def register_read_weight(read_id, weight):
    """Record one read's normalization weight; the last write wins.

    Acceptance probability is a property of an alignment record, not of a read, so a
    read with several records has several candidate weights. The caller resolves that
    by writing exactly one authoritative value: the weight of the record whose path
    was actually chosen to represent the read. Overwriting rather than combining
    keeps this honest -- a maximum or a sum over competing records would describe a
    read that was never observed that way.

    Under DEBUG a repeat write must agree with what is already recorded. The only
    legitimate repeat is the provisional write for a read whose chosen record turns
    out to be the one already used, so a disagreement means a read is being described
    two different ways.
    """
    if read_id is None:
        return
    try:
        w = float(weight)
    except (TypeError, ValueError):
        return
    if not (w > 0.0):
        return
    key = int(read_id)
    if DEBUG:
        prior = READ_WEIGHT_REGISTRY.get(key)
        if prior is not None and abs(prior - w) > 1e-9:
            raise RuntimeError(
                "conflicting normalization weights for read id {}: {} then {}".format(
                    key, prior, w
                )
            )
    READ_WEIGHT_REGISTRY[key] = w


def read_weight_for_id(read_id):
    """This read's weight, or 1 when it was never thinned."""
    try:
        return READ_WEIGHT_REGISTRY.get(int(read_id), 1.0)
    except (TypeError, ValueError):
        return 1.0


# Read IDs of the synthetic multipaths injected for input transcripts. They give the
# reference structure a template in the graph, but they are not observations, so path
# scoring excludes them: a path supported by nothing but these has no read evidence and
# must not be selected as a candidate. Populated per splice graph while incorporating
# input transcripts, and cleared when a new multipath graph is built.
SYNTHETIC_READ_IDS = set()

# Barcodes accepted as real cells, from --cell_list. Empty means no list was given,
# in which case the supporting-cell filter trusts every barcode in the BAM, which
# is only safe when the BAM was already restricted to called cells (as the
# cluster-guided partitioner does).
CELL_ROSTER = set()
