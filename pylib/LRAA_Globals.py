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
SPLICE_GRAPH_NORMALIZATION_METHOD = "cov5"

config = {
    #########################
    # read alignment criteria
    "HiFi": False,  # set to True when --HiFi is used; enables HiFi-specific filtering
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
    "fracture_splice_graph_at_input_transcript_bounds": True,
    "max_path_nodes_per_component": 1000,  # max number of path graph nodes per connected component
    # transcript reclustering (gene definition) criteria
    # gene reclustering overlap thresholds
    "min_recluster_overlap_shorter_iso_frac": 0.50,  # (overlap_len / shorter_transcript_len) >= this to connect isoforms in second-stage graph
    "min_recluster_overlap_longer_iso_frac": 0.20,  # also require (overlap_len / longer_transcript_len) >= this to avoid linking large multi-exon to long single-exon with tiny shared portion
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
    # [T - wiggle/2, T + wiggle/2], so the default 1 means T +/- 0.5 Mb. Within
    # that window the cut is placed where it severs the fewest retained primary
    # alignments; the window is never widened to find a zero-crossing position.
    "approx_MB_per_cut_wiggle_window": 1,
    #
    #######
    # quant
    "num_total_reads": None,  # for TPM and filtering - set by CLI or within LRAA by counting bam records
    "run_EM": True,
    "max_EM_iterations_quant_only": 250,  # don't set too high, as even at 1000 small biases get greatly amplified.
    "max_EM_iterations_during_asm": 1000,  # for asm, want higher iterations to amplify small diffs and weed out poorly supported isoforms.
    "aggressively_assign_reads": False,
    "quant_read_assignment_mode": "rescue_unassigned",
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
    # Honor the XW coverage-normalization weight during quantification, so a
    # multipath's support is the weight its reads stand for rather than a count of
    # the reads that survived thinning. Off by default: quantification normally
    # reads the unnormalized bam, where every weight is 1 and this changes nothing,
    # and the paths that reconstruct multipaths from more than one alignment record
    # per read have not been validated under weighting. Requires primary-only
    # alignments; multi-record groups raise rather than produce a weighted number
    # from a configuration that was never checked.
    "use_XW_read_weights_for_quant": False,
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
    # Two-pass alternative to the default final quantification, off by default. The
    # first pass quantifies normally against the coverage-normalized bam; the second
    # streams the full bam, looks each read's path up in the table the first pass
    # produced, writes its tracking row and forgets it. Nothing per-read is retained,
    # which is what makes a billion-read library tractable -- the default path holds
    # each shard's alignments and every read name it will report.
    #
    # This reports one expectation step at the first pass's abundances, where the
    # default path re-estimates them, so its counts are close to but not identical
    # with the default's. Left off so the standard behaviour is untouched.
    "stream_reads": False,
    # Refuse the run when more than this share of streamed reads land on paths the table
    # never saw. Not a correctness bound: those reads are resolved in-stream by the same
    # cascade and theta the first pass would have used, so they are assigned, not dropped.
    # It is a "the table is not doing its job" tripwire -- a graph mismatch, or a first pass
    # thinned so far it precomputed nothing useful, shows up here as a rate near 100%.
    #
    # Measured on ONT chr20 with the default target of 1000: 7.39% of reads on chr20+ and
    # 7.55% on chr20- landed on unseen paths, across 985 and 998 distinct paths against
    # tables of 4724 and 4621. Resolution is cached per path, so that is ~4.6 reads per
    # resolve and roughly 21% extra cascade work rather than 7.5%. The earlier 2% default
    # was calibrated against a different quantity and fired on healthy data.
    "stream_reads_max_unseen_path_read_frac": 0.25,
    # Rescue candidates against the local transcriptome from inside the streaming pass,
    # using a resident mappy index instead of the batch path's minimap2 subprocess. Off
    # by default, and refused unless --stream_reads_rescue_unassigned is given: without
    # it --stream_reads still requires --quant_read_assignment_mode genome, because a
    # streaming pass that silently skipped rescue would report the first pass's rescue
    # summary as if it covered the whole bam.
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
    "resource_monitor_interval": 60.0,  # seconds
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
