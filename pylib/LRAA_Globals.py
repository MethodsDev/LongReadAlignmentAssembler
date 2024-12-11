## global vars / constants

SPACER = "???"

DEBUG = False

config = {
    #########################
    # read alignment criteria
    "min_per_id": 98,
    "min_mapping_quality": 20,
    "try_correct_alignments": True,
    "max_softclip_realign_test": 20,
    "min_frac_alignments_pass_per_id_check": 0.9,
    "min_total_alignments_engage_frac_per_id_check": 1000,
    #
    ####################################
    # splice graph construction criteria
    "min_alt_splice_freq": 0.01,
    "min_alt_unspliced_freq": 0.01,
    "min_feature_frac_overlap": 0.50,
    "max_exon_spur_length": 5,  # exon spurs not tied to TSS or PolyA and at most this length get pruned
    "max_aggregate_splice_boundary_per_id": 90,  # if min_per_id <= this val, alt splice w/ shared bundary have the other boundary aggregated if within aggregate_splice_boundaary_dist
    "aggregate_splice_boundary_dist": 5,
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
    "min_frac_gene_alignments_define_TSS_site": 0.05,
    #
    ####################
    ## polyA site config
    "infer_PolyA": True,  # include PolyA site feature in read path assignments
    "max_dist_between_alt_polyA_sites": 50,
    "min_alignments_define_polyA_site": 3,
    "min_frac_alignments_define_polyA_site": 0.1,
    "min_PolyA_ident_length": 7,  # examine softclipped ends of reads, if have polyA with at least this number of bases at terminus, strip it and extended match out
    "min_PolyA_iso_fraction": 0.05,  # during initial TSS definition, require for a 'gene' that a TSS has at least this fraction of TSS-candidate gene reads assigned..
    "max_soft_clip_at_PolyA": 3,
    #
    ## read assignment to transcript criteria
    "fraction_read_align_overlap": 0.75,  # min fraction of read length that must overlap the compatible transcript isoform structure
    #
    # misc settings
    "min_path_score": 2,  # min number of reads required for reporting isoform
    #
    # transcript filtering criteria
    "min_transcript_length": 200,
    "min_isoform_fraction": 0.01,
    "min_frac_gene_unique_reads": 0.01,  # minimum fraction of all uniquely assigned reads per gene
    "min_monoexonic_TPM": 1.0,
    "filter_internal_priming": True,
    #
    ##########
    # assembly
    "normalize_max_cov_level": 10000,
    "restrict_asm_to_collapse": True,  # if True, no chaining of overlapping/extended paths
    "collapse_alt_TSS_and_PolyA": False,  # if True, collapses paths that are overlapping and contained but differ in TSS or PolyA
    #
    #######
    # quant
    "num_total_reads": None,  # for TPM and filtering
}
