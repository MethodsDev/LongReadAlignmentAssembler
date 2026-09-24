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
    # ON for every platform as of v0.37.0, paired with
    # strip_polyA_on_proximal_window below. The two are coupled: a read carrying
    # a polyA tail cannot define a polyA site while the tail is still attached,
    # because max_soft_clip_at_PolyA is 0, so inference without the strip fix is
    # inert on untrimmed data and the strip fix without inference has no consumer.
    #
    # MEASURED on SGNex MCF7 (real untrimmed human ONT, chr1/2/12), which is the
    # authoritative substrate for ONT decisions here. Counting polyA sites that
    # land within 10 bp of a GENCODE annotated 3' end, against a background of
    # the same sites shifted 20-60 kb:
    #
    #   inference off (previous default)     0 sites   -- none are built at all
    #   inference on, strip off            178 accurate of 206   (86.4%, bg 0.5%)
    #   inference on, strip on           2,339 accurate of 3,222 (72.6%, bg 0.4%)
    #
    # 13x more correctly placed sites. The per-site rate falls only because 15x
    # more sites are admitted; the absolute count of accurate ones rises, and at
    # 72.6% against a 0.4% background these are overwhelmingly real cleavage
    # sites rather than noise. Transcript level on the same runs: exact GENCODE
    # matches 2,009 -> 2,206 (+9.8%), sensitivity flat at 7.0%, precision
    # 32.2% -> 30.4%.
    #
    # The five-dataset Sequins/SIRV benchmark scores this -0.026 median F1strict
    # and is NOT the basis for the default. Those FASTQs are pychopper-processed,
    # so an alignment's 3' end is set by upstream trimming rather than by the
    # molecule, the strip fix is inert on them by construction (whole-clip and
    # proximal-window agree to within 20 reads in 20,000), and their read
    # boundaries are dirtier than the spike-in design suggests. A whole-structure
    # metric that requires exact termini therefore penalises any terminus change
    # there whether or not it is an improvement.
    #
    # infer_TSS was previously along for the ride: with only the whole-clip rule
    # it moved 3 of 2,206 exact matches on MCF7, because a raw ONT 5' clip is
    # adapter plus primer and only 0.03% of reads could define a TSS at all.
    # min_proximal_untemplated_G_at_TSS below is the 5' counterpart to the polyA
    # strip and closes that: eligibility on raw MCF7 goes 0.02% -> 29.91% and the
    # run calls 371 TSS sites where it previously called zero.
    "infer_TSS": True,  # include TSS feature in read path assignments
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
    #
    # Companion to the key above, for clips the key above cannot reach. That one
    # asks whether the WHOLE clip is a short pure G run; this one asks whether the
    # run of G's TOUCHING THE ALIGNMENT is at least this long, however much
    # adapter and primer sits beyond it. On ONT cDNA the clip is adapter + primer
    # + G run at a median of 80 bp, so the whole-clip form essentially never
    # fires: 139 of 470,290 MCF7 reads, which call ZERO TSS sites.
    #
    # 3 because THREE is the template-switch signature, not a free constant:
    # reverse transcriptase adds three non-templated C's on reaching the cap and
    # the SQK-DCS109/PCS109 strand-switching primer's GGG anneals to them. A
    # fourth G is the explainable variant rather than stronger evidence -- it
    # arises when only two of the primer's G's anneal to the three C's, leaving
    # one C to template an extra G.
    #
    # Together with the key above the rule reads: admit a read when the proximal
    # run is >= this, whatever follows it, or when the clip is nothing but G's.
    # One principle -- the whole clip must be ACCOUNTED FOR. With the signature
    # present, the rest is adapter; without it, there must be no unexplained
    # sequence at all.
    #
    # The threshold trades recall against terminal-vertex precision, and both
    # ends are defensible. MEASURED on MCF7 chr20 against FANTOM5 CAGE caps
    # (score>=100), which unlike GENCODE can credit a novel start:
    #   >= 4  ->  90 sites, 63 cap-backed (70.0%), background 0.00%
    #   >= 3  -> 371 sites, 109 cap-backed (29.4%), background 0.54%
    # So 3 recovers 46 more genuine starts and admits 235 more sites without cap
    # support. Some of those will be MCF7-specific starts FANTOM5 never sampled;
    # the rising background says not all of them are. A wrong TSS becomes a graph
    # terminal and truncates every model through the locus, so raise this to 4 if
    # precision at the terminus matters more than completeness.
    #
    # Kit-coupled: it means "the primer's own G count", not the literal 3. A
    # different strand-switching primer needs a different value. 0 disables it,
    # leaving only the whole-clip rule above.
    "min_proximal_untemplated_G_at_TSS": 3,
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
    # ON for every platform as of v0.37.0. This is the switch the MCF7 evidence
    # in the TSS block above actually argues for -- see there for the numbers;
    # it must move together with strip_polyA_on_proximal_window below.
    "infer_PolyA": True,  # include PolyA site feature in read path assignments
    "max_dist_between_alt_polyA_sites": 50,
    "min_alignments_define_polyA_site": 5,
    "min_frac_alignments_define_polyA_site": 0.1,
    "min_PolyA_ident_length": 7,  # examine softclipped ends of reads, if have polyA with at least this number of bases at terminus, strip it and extended match out
    "min_PolyA_iso_fraction": 0.05,  # during initial TSS definition, require for a 'gene' that a TSS has at least this fraction of polyA-candidate gene reads assigned..
    "max_soft_clip_at_PolyA": 0,  # max amount of softclipping allowed at the end of an alignment to mark it as a candidate boundary
    "min_soft_clip_PolyA_base_frac_for_conversion": 0.8,  # if soft-clipped is at least this frac polyA evidence, then removing soft clipping and marking as candidate polyA read.
    #
    # The two keys above decide whether to STRIP a clip and move a boundary; the
    # two below decide whether a clip counts as EVIDENCE of a real 3' end, and
    # they are deliberately separate because the questions differ. Stripping is
    # scored over the whole clip, which is safe only when the clip is nothing but
    # tail. An ONT cDNA clip is tail + adapter + barcode, so that test fires on
    # 1.1% of reads at GENCODE 3' ends on SGNex MCF7 where scoring the 20 bases
    # nearest the alignment fires on 65.2%. Evidence is therefore read from the
    # proximal window.
    "polyA_tail_proximal_window": 20,  # bases of the clip nearest the alignment scored for tail evidence
    "min_proximal_tail_base_frac": 0.8,  # min A (fwd) / T (rev) fraction in that window to call an external tail
    #
    # Whether proximal tail evidence may also STRIP the clip, and not merely be
    # recorded. A read may only define a polyA site when its residual clip there
    # is <= max_soft_clip_at_PolyA (0 below), so a tail left in place silences
    # the read that carries the evidence.
    #
    # This is SUPPORT FOR RAW cDNA, not a defect in the behaviour above. The
    # whole-clip test is correct for the input LRAA has assumed: Pychopper and
    # Kinnex remove adapters and primers but leave the tail, so the clip IS the
    # tail and scoring all of it works. A raw ONT cDNA clip is tail + adapter +
    # barcode, and the same test then fires on 1.1% of reads. MEASURED on raw
    # MCF7: of 1,729 reads carrying a genuine proximal tail, 27 (1.6%) were
    # allowed to contribute to a polyA site. Enabling this takes eligibility
    # there from 1.1% to 50.7%.
    #
    # Unlike the evidence slot, this MOVES PolyA vertices, so it changes the
    # assembly and invalidates both the alignment and splice-graph caches.
    #
    # On by default, and inert on input that was already trimmed -- which is
    # what makes it safe to default: it only acts on clips the whole-clip test
    # misses, and only matters at all when polyA sites are being inferred.
    # Verified inert on
    # adapter-trimmed input -- on the pychopper-processed Sequins/SIRV benchmark
    # the two tests agree to within 20 reads in 20,000 and all five datasets
    # produce byte-identical GTFs in both modes. Also inert on Kinnex/Iso-Seq,
    # where 98.4% of reads carry a 3' clip shorter than min_PolyA_ident_length.
    # Where it does act -- untrimmed ONT with inference on -- it recovers 187
    # additional exactly-matching transcripts against GENCODE on SGNex MCF7
    # (+9.3%), for 1.8 points of precision.
    "strip_polyA_on_proximal_window": True,
    #
    # An untrimmed polyA tail can be ALIGNED rather than soft-clipped: minimap2
    # in splice mode will happily place it on a genomic A-run kilobases away,
    # behind a spurious terminal intron.  The soft-clip handling above never sees
    # those bases, because by then they are an exon.  Such an alignment reports a
    # 3' end in the wrong place and invents a junction, so the whole alignment is
    # discarded rather than repaired -- the transcript's remaining reads, which
    # terminate correctly, are what should define its boundary.
    "no_exclude_polyA_terminal_segment": False,
    "max_polyA_terminal_segment_length": 24,  # only a block this short can be a landed tail rather than an exon
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
    # A contained model's terminal feature stops protecting it from absorption
    # once that feature holds this fraction or less of the strongest feature of
    # the same type in its gene. Consulted only by the containment decision in
    # TranscriptFiltering.prune_likely_degradation_products; the site itself is
    # never deleted, it just no longer confers immunity.
    #
    # 0.20 mirrors max_frac_alt_TSS_from_degradation, which asks the same question
    # of the same kind of evidence. MEASURED on five Arena ONT datasets: polyA
    # sites that appear only once infer_PolyA is on hold a median 0.06-0.14 of the
    # gene's strongest site, against exactly 1.000 at annotated 3' ends. The
    # expression test this sits in front of separates those two classes not at all
    # (median 0.350 vs 0.349 on LSK109). Set 0 to require exact feature identity,
    # which is the prior behaviour.
    "max_frac_alt_terminal_feature_absorbable": 0.20,
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
    #########################################################################
    # Whole-genome alignment-mismapping filter (v0.40.0). A POST-MERGE stage --
    # it runs once on the merged genome-wide gtf + quant.expr, never per chunk --
    # that removes isoforms which are alignment/strand-mismapping artifacts of a
    # much-higher-expressed transcript. Two independent detectors, unioned:
    #   (1) MIRROR (coordinate): a multi-exon model on the opposite strand to a
    #       higher-expressed model, with high exonic base overlap and every splice
    #       site within mismap_junction_tolerance bp of that model's exon
    #       boundaries. Catches wrong-strand ("s") near-mirrors and needs no genome.
    #   (2) SEQUENCE (minimap2 cDNA all-vs-all): a model whose spliced cDNA is
    #       >= mismap_min_seq_identity%% identical over >= mismap_min_seq_coverage of
    #       its length to a DIFFERENT-gene, higher-expressed model. Catches run-ons
    #       ("x"), chimeras, and same-strand mismappings.
    # Both gate on mismap_max_expr_fraction: only removed when the model carries
    # < that fraction of the matched model's expression. That low-expression gate
    # is what makes it safe (well-expressed paralogs are retained) and lets the
    # quant be repaired by dropping the rows and renormalizing TPM rather than
    # requantifying. The cross-gene requirement (sequence detector) keeps genuine
    # minor same-gene isoforms. Requires the genome fasta (for cDNA extraction).
    "filter_mismappings": True,  # master switch; --no_filter_mismappings sets False
    "mismap_min_seq_identity": 99.0,  # percent identity of the cDNA-vs-cDNA match
    "mismap_min_seq_coverage": 0.85,  # fraction of the query cDNA the match must cover
    "mismap_max_expr_fraction": 0.01,  # remove only if expr < this fraction of the match's expr
    "mismap_junction_tolerance": 20,  # bp tolerance, mirror splice site vs opp-strand exon boundary
    "mismap_min_base_overlap": 0.5,  # mirror: min exonic base overlap fraction with the opp-strand model
    # No model survives assembly/discovery filtering that quantification could not put a
    # whole read on. A COUNT of assigned reads, applied to every model -- novel and
    # reference-containing, monoexonic and spliced -- before any other filter runs, and
    # again after the isoform-fraction EM moves the counts.
    #
    # Scoped to filtering, not to the reported number. The final quant that follows can
    # still put an all_reads below this, and is deliberately not gated: see the NOTE
    # beside it in LRAA. Quant-only does not apply it either -- a quantification answers
    # about every transcript it was asked about.
    #
    # The thresholds below it are all relative -- TPM against library depth, isoform
    # fraction against the gene, cells against the roster -- so each of them can be
    # cleared by a model holding a hundredth of a read in a quiet neighbourhood. This is
    # the one absolute floor, and it is in the unit the question is actually asked in:
    # "did a read support this". Fractional because EM assignments are fractional, and a
    # total rather than a unique count because an isoform indistinguishable from its
    # neighbours over every individual read still earns its row once a read's worth of
    # mass lands on it.
    #
    # 0 disables it and restores the previous behaviour, where any nonzero EM mass was
    # enough for a multi-exonic model.
    "min_reads_retain_isoform": 1.0,
    "min_monoexonic_TPM": 1.0,
    # Require a single-exon model to show some evidence of a real 3' end. A
    # multi-exonic model is vouched for by its intron chain; a monoexonic one has
    # nothing structural, so we ask whether anything marks where it ends, and
    # accept ANY of three independent channels: an inferred PolyA site, a genomic
    # PAS hexamer upstream of the terminus, or an assigned read whose own clip
    # there looks like a tail.
    #
    # This replaces require_terminal_feature_for_monoexonic, which demanded one
    # SPECIFIC channel (inferred TSS or PolyA) and so measured whether inference
    # succeeded in that run rather than whether a model is real: it deleted
    # 99.98% of monoexonic models on A549 ONT and 55% of all output on MCF7,
    # while 88.7% of HiFi models clear it. The three channels fail on opposite
    # inputs -- Kinnex strips tails, ONT rarely yields callable PolyA sites -- so
    # the disjunction is what makes one rule portable across platforms.
    #
    # MUST be declared here, not merely read with config.get(): --config_update
    # drops keys absent from this dict ("ignoring unknown config key"), so an
    # undeclared key is silently a no-op in every chunked run.
    "require_terminal_evidence_for_monoexonic": True,
    # How close a read's polyA tail must sit to a model's 3' terminus to count as
    # evidence FOR THAT MODEL. A tail further along is evidence about a different
    # overlapping transcript.
    "max_dist_tail_evidence_to_terminus": 100,
    # A multi-exonic model is kept only if quantification gave it expression above this
    # value. The default of 0 means "any expression at all": supplied models are
    # selectable from the trellis on their synthetic template read, and this decides
    # whether they were actually expressed. Raise it to demand more than a trace; set it
    # negative to disable the check and report every selected multi-exonic model.
    #
    # This is NOT the reference reprieve and no longer matches it. retain_expressed asks
    # for min_reads_retain_reference (a whole assigned read); this asks only for nonzero
    # mass, and it applies to novel and reference-containing models alike.
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
    # to a NOVEL monoexonic model. Monoexonic models containing a reference model are
    # exempt via reference_model_reprieved(), so the exemption is conditional on
    # min_reads_retain_reference and is not categorical; bulk input carries no barcode
    # in its read names so the check self-disables there. An absolute count
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
    # Where the internal-priming veto DELETES a read-derived PolyA candidate during
    # site identification (Splice_graph._incorporate_PolyA_objects), rather than letting
    # it through to be judged later at transcript filtering.
    #
    #   "always"       -- delete in every graph. Behaviour up to v0.34.0.
    #   "spliced_only" -- DEFAULT. Delete in the ME (spliced) graph, KEEP in the SE
    #                     (monoexonic) graph built separately by build_SE_transcripts.
    #   "never"        -- keep everywhere. NOT recommended: it rewrites the ME graph
    #                     too, spawning competing spliced 3'-variants. Measured on
    #                     chr22 it triples FSM reference loss (46 vs 15) and inflates
    #                     ISM 12%, for roughly twice the monoexonic reduction.
    #
    # Deleting the site removes a terminus the graph needs: no path can END there, so
    # every model over the locus runs on to the next available 3' vertex. Measured at
    # DGCR2 (chr22:19.117-19.120 Mb, two +-strand monoexonic coverage peaks antisense to
    # the gene): all 21 +-strand candidates were rejected, and the models built then
    # spanned BOTH peaks -- 19,117,546-19,119,735 -- across a valley where coverage falls
    # from ~280 to ~22. The left peak was not absent from the output, it was absorbed
    # into an over-long neighbour.
    #
    # The same placement error costs annotated 3' ends: spare_polyA_veto_at_known_3prime
    # below can only spare a candidate that EXISTS, so deleting it first makes that
    # reprieve unreachable.
    #
    # "spliced_only" exists because keeping the site everywhere has a measured cost on
    # SPLICED models, which is where this rule was never needed: the terminus only has
    # to exist for monoexonic reconstruction. On chr22, "never" moved se_antisense
    # 100 -> 24 but also lost 46 FSM reference isoforms and raised ISM 12% -- the
    # absorbing-vertex truncation the rejection site warns about. The SE graph is built
    # separately (build_SE_transcripts -> Splice_graph(restrict_splice_type="SE")), so
    # the deferral can be confined to it and the spliced graph left untouched.
    #
    # Judgement is deferred to TranscriptFiltering.filter_internally_primed_transcripts,
    # which annotates every emitted 3' terminus and deletes primed monoexonic models.
    #
    # Affects graph construction, hence registered in _SPLICE_GRAPH_CONFIG_KEYS.
    # Measured "spliced_only" vs "always": se_antisense 124 -> 59 (chr20) and
    # 100 -> 51 (chr22); ALL spliced output unchanged to the model on chr22
    # (4,252 both ways), ISM +1%, 8-15 FSM reference isoforms lost of ~1,500.
    "reject_internally_primed_polyA_sites": "spliced_only",
    # Whether the ME (spliced) graph may emit SINGLE-EXON models.
    #
    # The ME graph is fed only reads with an intron in their CIGAR
    # (Pretty_alignment_manager partitions on has_introns()), but that is a GLOBAL
    # property of the alignment. A read whose junction lies far outside a given window,
    # or whose junction the graph rejected for lack of support, contributes only exonic
    # blocks locally. Where such blocks accumulate with no validated junction at their
    # boundaries, the ME graph grows an ISOLATED exon segment and path enumeration can
    # only emit it as a single-exon model.
    #
    # Measured at chr20:35,674,748-35,676,911 (antisense to NFS1): the ME graph held one
    # exon segment E:4928[+] with no incident intron -- nearest + introns 43 kb upstream
    # and 15 kb downstream -- built from 47 spliced reads of which 42 had every intron
    # outside the window. It emitted a 2,163 bp single-exon model that overran a 147-read
    # internally primed 3' stack by 153 nt and swallowed a second 63-read primed site.
    # 19 of 181 ME models (10%) in that region were single-exon.
    #
    # Those models are also load-bearing downstream: ME_transcripts becomes the
    # SE_read_encapsulation_mask (LRAA), so an ME single-exon model masks the very
    # monoexonic reads the SE graph would otherwise have used to place the terminus
    # correctly. That makes the defect self-reinforcing.
    #
    # DEFAULT False: single-exon reconstruction is routed to the SE graph alone, which
    # is the graph that has the monoexonic reads, the monoexonic filters, and -- under
    # reject_internally_primed_polyA_sites="spliced_only" -- the retained PolyA termini.
    # Reference single-exon transcripts are unaffected: the ME builder is handed only the
    # intron-bearing reference subset, so it never carries them.
    #
    # Measured on chr20 on top of "spliced_only": se_antisense 59 -> 40, with FSM
    # reference loss essentially unchanged (8 -> 9 of ~1,500). At the NFS1 locus it is
    # what lets the SE graph see the reads at all and place termini on the two primed
    # sites instead of one model spanning both peaks.
    "ME_graph_emits_monoexonic_models": False,
    # Subtract the intronic coverage floor before the SE graph is segmented.
    #
    # Monoexonic coverage inside introns has a substantial near-uniform floor: measured
    # on PBMC chr20/chr22, 491k/465k monoexonic reads lie wholly outside annotated
    # exons, ~2 reads deep at the median intronic base. Segmenting that directly fuses
    # neighbouring real features into one smeared exon segment, which is where loose
    # single-exon boundaries come from -- measured 20-24% of a surviving monoexonic
    # model's span sits at or below its own host intron's background.
    #
    # SE graph only; the ME graph is untouched. Only bases intronic in the ME transcript
    # set and exonic in none of them are eligible, so nothing the spliced graph called an
    # exon is decremented. Measured with that protection: 33-35% of covered intronic
    # bases fall below the floor, 0 protected exonic bases touched, and coverage islands
    # RISE (3,450 -> 4,436 on chr20, ~1,100-1,400 splits per contig) as smears resolve
    # into separate features. Model counts are not the target and barely move: the
    # antisense monoexons sit 25-54x above their local background.
    #
    # ON BY DEFAULT. Measured cost, chr20, applied on top of spliced_only +
    # ME-monoexonic suppression: novel monoexonic 170 -> 148, se_genic 50 -> 30,
    # se_intergenic 1 -> 0, and 15 FSM reference isoforms lost against 22 gained. The
    # FSM movement is NOT caused in the spliced graph -- __ME_isoforms.gtf is
    # byte-identical with and without this setting (verified md5 6c0505c02b1e at
    # chr20:19.95-20.10 Mb). It arises downstream, where ME and SE models are combined
    # (LRAA) and quantified together: nine fewer SE models shift gene denominators and
    # EM mass, and the threshold filters that run next -- the absolute 1.0-read floor,
    # the weakest-first isoform-fraction filter, the degradation pruner -- are
    # order-sensitive and grant no reference reprieve, so marginal multi-exon models
    # flip. The models lost sit at median isoform fraction 0.031 against 0.061 overall.
    #
    # It removes whole low-coverage models more than it tightens boundaries: of 187
    # matched monoexonic models only 12 shortened (median 236 bp) and 174 were
    # unchanged, because a model peaking 25-54x above its local background does not
    # move its edges when 2 reads are subtracted.
    #
    # Affects graph construction, hence registered in _SPLICE_GRAPH_CONFIG_KEYS.
    "SE_subtract_intronic_background": True,
    # Intron-length percentile above which the subtraction is skipped. Long introns are
    # heterogeneous enough that a single median describes them poorly, and they hold most
    # intronic sequence but few of the fusions this addresses: p80 is 9,720 bp on chr20
    # and 7,558 bp on chr22, covering 80% of introns but only 22% of intronic bp.
    # Computed per contig-strand from the ME transcripts actually present.
    "SE_intronic_background_intron_length_pctile": 80,
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
    # What retain_expressed demands of a reference-containing model before it is
    # exempted from the discovery filters. A COUNT of assigned reads, not a rate:
    # the reprieve previously asked get_TPM() > 0, which is
    # read_counts_assigned / num_total_reads and so fires on any nonzero EM mass at
    # all. On a 52M-read PBMC library that admitted 14,193 multi-exon reference
    # chains whose entire assigned mass was below 0.05 of one read -- 99% of them
    # carried some fractional assignment, none carried a read. write_expr prints
    # all_reads at one decimal, so they surfaced as "0.0" and read as unsupported
    # models being reported on the strength of their annotation.
    #
    # 1.0 means "the reads assigned to this structure sum to at least one read".
    # Deliberately not a UNIQUE read: an isoform indistinguishable from its
    # neighbours over every read still earns its output row once a full read's worth
    # of mass lands on it. Fractional because EM assignments are fractional; set 0
    # to restore the any-nonzero-mass behaviour.
    "min_reads_retain_reference": 1.0,
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
    # "Uniquely assigned" means exclusively assignable: exactly one compatible
    # isoform. One definition feeds every consumer -- the novel-isoform floor and
    # frac_gene_unique_reads. unique_read_filter_min_frac no longer decides
    # uniqueness; it is retained only as a decision-log diagnostic so the divergence
    # from the old threshold stays observable.
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
    # Assignment fraction at or above which the SUPERSEDED criterion counted a read
    # as uniquely assigned. It no longer decides uniqueness anywhere -- exclusive
    # assignability does -- and is read only to record how far the two diverge, as the
    # decision log's near_unique and effectively_unique columns.
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
