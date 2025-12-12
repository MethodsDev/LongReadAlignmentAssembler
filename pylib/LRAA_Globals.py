import os

## global vars / constants

SPACER = "???"

DEBUG = False

LRAA_MODE = "unset"  # options ("ID", "QUANT-ONLY", "MERGE")

config = {
    #########################
    # read alignment criteria
    "HiFi": False,  # set to True when --HiFi is used; enables HiFi-specific filtering
    "min_per_id": 80,
    "min_mapping_quality": 1,
    "try_correct_alignments": True,
    "max_softclip_realign_test": 20,
    "min_softclip_realign_test": 5,
    "min_frac_alignments_pass_per_id_check": 0.9,
    "min_total_alignments_engage_frac_per_id_check": 1000,
    "min_terminal_splice_exon_anchor_length": 15,
    "read_aln_gap_merge_int": 10,
    "max_intron_length": 100000,
    #
    ####################################
    # splice graph construction criteria
    # default tuned for not HiFi (e.g., ONT); HiFi mode overrides to 0.01 via --HiFi
    "min_alt_splice_freq": 0.03,
    "min_alt_unspliced_freq": 0.01,
    "min_feature_frac_overlap": 0.50,
    "max_exon_spur_length": 20,  # exon spurs not tied to TSS or PolyA and at most this length get pruned (HiFi sets 13)
    "aggregate_adjacent_splice_boundaries": True,
    "aggregate_splice_boundary_dist": 5,
    "fracture_splice_graph_at_input_transcript_bounds": False,  # enabled under HiFi mode
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
    "terminal_boundary_method": "extreme",  # choices: "extreme" (min/max), "mean", "median", "quartile" (Q1/Q3) - method for defining terminal coords when TSS/PolyA not annotated
    "min_reads_for_terminal_adjustment": 7,  # minimum number of reads terminating in terminal exon required for mean/median/quartile adjustment (falls back to extreme if below threshold)
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
    "filter_internal_priming": True,
    "restrict_internal_priming_filter_to_monoexonic": True,
    "ref_trans_filter_mode": "retain_expressed",  # choices ["retain_expressed", "retain_filtered"]
    "min_reads_novel_isoform": 2,
    "min_unique_reads_novel_isoform": 2,
    "min_isoform_count_aggressive_filtering_iso_fraction": 10,  # allow for filtering mult isoforms in a single round if more than this number of isoform candidates.
    #
    ##########
    # assembly
    "normalize_max_cov_level": 1000,
    "restrict_asm_to_collapse": True,  # if True, no chaining of overlapping/extended paths
    "collapse_alt_TSS_and_PolyA": False,  # if True, collapses paths that are overlapping and contained but differ in TSS or PolyA
    #
    #######
    # quant
    "num_total_reads": None,  # for TPM and filtering - set by CLI or within LRAA by counting bam records
    "run_EM": True,
    "max_EM_iterations_quant_only": 250,  # don't set too high, as even at 1000 small biases get greatly amplified.
    "max_EM_iterations_during_asm": 1000,  # for asm, want higher iterations to amplify small diffs and weed out poorly supported isoforms.
    "aggressively_assign_reads": False,
    # When True, weight ambiguous read assignments by agreement of read 3' ends with transcript 3' ends
    # (previously "use_weighted_read_assignments" which weighted by both 5' and 3' ends)
    "weight_reads_by_3prime_agreement": True,
    "EM_alpha": 0.01,  # regularization
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
}

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

