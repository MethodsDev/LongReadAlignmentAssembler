## global vars / constants

SPACER = "???"

DEBUG = False

LRAA_MODE = "unset"  # options ("ID", "QUANT-ONLY", "MERGE")

config = {
    #########################
    # read alignment criteria
    "min_per_id": 98,
    "min_mapping_quality": 0,
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
    "min_alt_splice_freq": 0.01,
    "min_alt_unspliced_freq": 0.01,
    "min_feature_frac_overlap": 0.50,
    "max_exon_spur_length": 5,  # exon spurs not tied to TSS or PolyA and at most this length get pruned
    "aggregate_adjacent_splice_boundaries": False,
    "aggregate_splice_boundary_dist": 5,
    "fracture_splice_graph_at_input_transcript_bounds": True,  # disabled under LowFi mode
    "max_path_nodes_per_component": 1000,  # max number of path graph nodes per connected component
    #
    ############
    # TSS config
    "infer_TSS": True,  # include TSS feature in read path assignments
    "max_dist_between_alt_TSS_sites": 0,
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
    "min_frac_gene_alignments_define_TSS_site": 0.1,
    #
    ####################
    ## polyA site config
    "infer_PolyA": True,  # include PolyA site feature in read path assignments
    "max_dist_between_alt_polyA_sites": 50,
    "min_alignments_define_polyA_site": 5,
    "min_frac_alignments_define_polyA_site": 0.1,
    "min_PolyA_ident_length": 7,  # examine softclipped ends of reads, if have polyA with at least this number of bases at terminus, strip it and extended match out
    "min_PolyA_iso_fraction": 0.05,  # during initial TSS definition, require for a 'gene' that a TSS has at least this fraction of polyA-candidate gene reads assigned..
    "max_soft_clip_at_PolyA": 3,
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
    "ref_trans_filter_mode": "retain_expressed",  # choices ["retain_expressed", "retain_filtered"]
    "min_reads_novel_isoform": 2,
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
    "use_weighted_read_assignments": True,
    "EM_implementation_use": "CGPT",  # choices: "BJH" or "CGPT"
    "EM_alpha": 0.01,  # regularization
    #
    ######
    # parallelization
    #
    "min_mpgn_component_size_for_spawn": 150,
}
