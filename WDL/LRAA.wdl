version 1.0

import "subwdls/Partition_data_by_chromosome.wdl" as PartByChr
import "subwdls/LRAA_runner.wdl" as LRAA_runner
import "subwdls/LRAA_chunk_scatter.wdl" as ChunkScatter


workflow LRAA_wf {
     input {
        String sample_id
                 
        File referenceGenome 
        File inputBAM
        # INTERNAL PLUMBING, not a user knob. Splice-graph evidence supplied by the
        # CALLING workflow: the single-cell final quant normalizes every cluster bam,
        # merges them, normalizes the merge once more, and hands that ONE bam to all
        # per-cluster quant jobs (LRAA_quant_by_cluster.wdl), so every cluster is
        # quantified against the same splice graph while counting only its own reads.
        # Supplying it IS the statement "already normalized, do not normalize again":
        # no_norm is derived from it below rather than offered separately, since
        # no_norm alone means "build the splice graph from uncapped depth", which is
        # not a thing to offer. Prefixed `internal_` so a workspace still binding the
        # old top-level `bam_for_sg` / `no_norm` fails on an unknown input rather than
        # quietly changing what the splice graph is built from.
        File? internal_bam_for_sg
        File? internal_bam_for_sg_index
        # INTERNAL PLUMBING, the same shape and the same reason as
        # internal_bam_for_sg above, and a DIFFERENT ROLE. Three bams, three roles,
        # per run:
        #
        #   inputBAM                 the full library. Pass 2 assigns THESE reads.
        #   internal_bam_for_sg      the splice graph. Shared across sibling runs and
        #                            never reconstructed here.
        #   internal_bam_for_priors  pass-1 theta. THIS run's own normalized reads.
        #
        # Unset, pass 1 falls back to the splice-graph bam whenever stream_reads is on
        # (LRAA:_first_pass_assignment_bam), which is the default -- and in the
        # cluster-guided shape that bam is ONE file shared by every cluster, so theta
        # comes from POOLED evidence and each cluster's ambiguous reads are apportioned
        # by every OTHER cluster's expression. It looks like it worked, because each
        # cluster still reports its own totals: observed 32 clusters all reporting
        # reads_total 94,908 while assigning 24,083 / 27,616 / 17,414 / ...
        #
        # Prefixed `internal_` for the reason internal_bam_for_sg is: which reads
        # estimate theta is a role the CALLING workflow fills from files it produced,
        # not a knob a workspace binds. Wanted in ALL THREE scattering modes -- `off`
        # is a single whole-genome invocation, and it still runs a pass 1.
        File? internal_bam_for_priors
        File? internal_bam_for_priors_index
        # INTERNAL PLUMBING, same shape and same reason as internal_bam_for_sg above:
        # ONE chunk plan, produced once by the calling workflow from the WHOLE
        # pre-partition bam and handed to every sibling run, so all of them cut at
        # identical positions. Without it each per-cluster run selects cuts on its own
        # bam, so the shared internal_bam_for_sg gets sliced at a different geometry per
        # cluster and the clusters stop being comparable -- exactly what supplying one
        # splice-graph bam was for. GEOMETRY ONLY: num_total_reads and discovery come
        # from each consuming run, never from the plan.
        #
        # Consumed by by_chunk AND by by_chromosome. by_chromosome is NOT an unchunked
        # path: subwdls/LRAA_runner.wdl passes --chunk unconditionally, so each
        # chromosome shard chunks its contig INSIDE itself and, without a plan, picks
        # those cut positions from that shard's own bam. The plan is genome-wide and
        # per-contig, so a shard restricted with --contig just uses its contig's entry.
        # Do not "simplify" this back to by_chunk-only: doing so silently restores
        # per-cluster chunk boundaries in the mode the single-cell pipeline defaults to.
        #
        # off is the one mode that does not take a plan -- it is a single whole-genome
        # invocation, so there is no sibling run for its geometry to have to match --
        # and supplying one there is refused in validate_scattering rather than dropped.
        File? internal_chunk_plan
        File? annot_gtf
        Boolean HiFi = false
         
        String main_chromosomes = "" # ex. "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"

        # HOW the work is divided across tasks. Chunking itself is unconditional in
        # all three; what differs is where the chunks run.
        #
        #   off             one LRAA task for the whole (optionally filtered)
        #                   genome, chunked internally under its own cpu_budget.
        #   by_chromosome   one task per contig, each chunked internally. The
        #                   makespan floor is the longest contig.
        #   by_chunk        one make-chunks task, a scatter of one task per CHUNK,
        #                   then one merge. The floor is one chunk, and the leaves
        #                   are small enough to run preemptible.
        #
        # main_chromosomes is a CONTIG FILTER in every mode, not the mode switch it
        # used to be: by_chromosome partitions those contigs, by_chunk restricts the
        # partition to them, and off passes a single name through as --contig.
        String scattering = "by_chunk"
        # Preemptibility of the per-chunk leaf tasks in by_chunk mode. The only
        # preemptible knob surfaced here: a chunk leaf is short, independent and
        # cheap to redo, whereas a chromosome shard or a merge is long-running and
        # loses substantial work to a preemption.
        Int chunkPreemptible = 3

        # by_chunk resource knobs, forwarded to subwdls/LRAA_chunk_scatter.wdl.
        # Optional on purpose: unset takes that subworkflow's defaults, which are
        # measured from a whole-genome PBMC run and documented there. Named with a
        # chunk prefix where they would otherwise collide with this workflow's own
        # cpu/memoryGB, which size a chromosome shard rather than a 10 Mb chunk.
        Int? chunkMakeChunksCpu
        Int? chunkMakeChunksMemoryGB
        Int? chunkCpu
        Int? chunkMemoryGB
        Int? chunkMergeCpu
        Int? chunkMergeMemoryGB
        String? region # example: "chr1:100000-200000"; when set, workflow will not split by chromosome and will pass --region to LRAA
        String? oversimplify # comma-separated contig names to run in oversimplify mode
        
        Float? min_per_id
        Boolean no_EM = false
        # Apportion an ambiguous read by how well its 3' end agrees with each
        # candidate transcript, or split it flat. True makes pylib/EM.py use a flat
        # 1.0 weight instead of the computed one, so it CHANGES REPORTED NUMBERS.
        # Reaches all three scattering modes: by_chromosome and off as LRAA's own
        # --no_weight_reads_by_3prime_agreement, by_chunk as a config override on
        # each leaf. Until now it reached none of them -- no workflow declared it,
        # and on a chunked run LRAA parsed the flag and then dropped it, because
        # its config key is the NEGATION of its dest and _NEGATED_CONFIG_FLAGS did
        # not map the pair.
        Boolean no_weight_reads_by_3prime_agreement = false
        Boolean quant_only = false
        # Return the depth-normalized BAM(s) the splice graph was built from. Single-cell
        # workflows that call this one set this false; they never surface the file, and
        # delocalizing it would cost them storage for nothing.
        Boolean retain_normalized_splice_graph_bam = true
        Boolean rescue_unassigned_reads_via_transcriptome_alignment = true
        Int min_mapping_quality = 0
        Int min_mapping_quality_for_final_quant = 0

        # If set, the count_bam step is skipped and this value is passed to LRAA as --num_total_reads
        Int? num_total_reads
        # Barcodes accepted as real cells. Without it the supporting-cell filter
        # trusts every barcode in the BAM, which is only safe when the BAM was already
        # restricted to called cells.
        File? cell_list

        String cell_barcode_tag = "CB"
        String read_umi_tag = "XM"

        # Chunked mode is UNCONDITIONAL -- there is no `chunk` input to turn it off. The
        # scatter below parallelises per CONTIG, so a whole-genome run's makespan is its
        # longest chromosome and more cores cannot shorten it. Chunking is what breaks
        # that floor: it splits each contig-strand at low-coverage positions between
        # annotated loci INSIDE the per-contig task, runs the pieces concurrently under
        # that task's one --cpu_budget, and merges. On chr1 at budget 8: 2.8x end to end,
        # and peak RSS 3.58 GB against 10.03 GB.
        #
        # That memory property is why the toggle is gone rather than merely defaulted on.
        # An unchunked run's peak scales with the contig it holds, so supporting both
        # modes meant every memory number here had to be a function of input size or of
        # cpu; with only one mode, peak follows chunk concurrency and chunk width and the
        # two fixed numbers below say all there is to say. A workspace still binding
        # `chunk` will fail on an unknown input rather than quietly change modes.
        #
        # Works in all three modes, so nothing is lost by forcing it: quant_only with
        # annot_gtf is quant-only, annot_gtf alone is ref-guided discovery, neither is de
        # novo. Chunked discovery agrees with unchunked EXACTLY on chr21 -- 1460 = 1460
        # models, 0 of 11,811 GTF rows differing, in both rescue configurations. LRAA's
        # one refusal is --quant_only without a gtf, which has nothing to quantify and
        # which an unchunked run refuses too.
        Float? approx_MB_per_cut
        Float? approx_MB_per_cut_wiggle_window
        # Chunk each contig-STRAND separately, splitting the whole bam by orientation
        # first as a serial phase. OFF by default: the strandless ordering splits inside
        # each chunk, concurrently with every other chunk, and removes the largest
        # serial phase a chunked run has -- 151.2 s against 255.2 s on the same input.
        #
        # Named for what setting it DOES. It replaces a `strandless_chunks` input that
        # read as a double negative once strandless became the default, and whose false
        # case emitted no flag at all, so asking for strand-first silently got you
        # strandless. A workspace still binding the old name will fail on an unknown
        # input rather than quietly run the other mode.
        #
        # Neither ordering needs num_total_reads from the caller. LRAA counts the library
        # itself before any chunking, with the same -F 0x904 policy as the count_bam
        # task below, so the scattered and direct paths agree and neither depends on
        # the caller getting `samtools view -c` right.
        Boolean chunk_by_strand = false

        # Two-pass streaming quantification; see subwdls/LRAA_runner.wdl's task input
        # for detail. ON by default since v0.25.0, matching LRAA's own default.
        Boolean stream_reads = true

        #  non-scattered runs
        # Cores for the LRAA task: its cpu request AND the --cpu_budget it divides across
        # work units. numThreadsPerWorker and num_parallel_contigs multiplied instead:
        # 5 x 3 asked LRAA for up to fifteen cores on a five-core task.
        Int cpu = 5
        # Override for the whole-genome box below.
        Int? memoryGB
        # The box a WHOLE-GENOME (non-scattered) run gets. One number, not a function of
        # the input BAM: every run chunks, so peak RSS is how many chunks run at once
        # times what one chunk holds, and neither term scales with library size.
        Int memoryGB_whole_genome = 32

        # scattered runs
        # Per-shard, not uniform: chr1 and chrM used to request the SAME 5 cores despite
        # differing 25x in chunk count, so chr1's chunk phase was core-starved while chrM's
        # sat mostly idle on cores it had no chunks to fill. Computed below from each
        # shard's own contig length, one array entry per contig_index; this stays an
        # OVERRIDE (unset = computed), the same shape as memoryGBPerWorkerScattered below --
        # a caller that still wants one fixed value for every shard can ask for it explicitly.
        Int? cpuScattered
        # Ceiling on the COMPUTED per-shard value only; an explicit cpuScattered above
        # ignores it. 16 matches the fleet-standard n2-standard-16 box this pipeline was
        # benchmarked on.
        Int max_cpu_per_chunked_shard = 16
        Int? memoryGBPerWorkerScattered
        # The box every chromosome shard gets. Uniform across shards on purpose: a shard's
        # peak follows its chunk concurrency, which max_cpu_per_chunked_shard already caps,
        # not the length of the contig it was handed.
        Int memoryGB_per_chromosome_shard = 16
        
        
        Int diskSizeGB = 256
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:latest"
        Int countBamThreads = 16
        
        
    }

    # Memory: two fixed numbers, one per shape of run.
    #
    #   whole genome, non-scattered  -> memoryGB_whole_genome,           32 GiB
    #   one chromosome shard         -> memoryGB_per_chromosome_shard,   16 GiB
    #
    # Fixed, and stated HERE rather than derived downstream, because both of the formulas
    # this replaces got the relative sizing wrong. LRAA_runner_task self-sizes at 2
    # GiB/core, which is proportional to cpu and therefore INVERTED the two: the
    # whole-genome run asks for 5 cores and got 16 GiB while the largest shards ask for
    # max_cpu_per_chunked_shard (16) and got 32, backwards for a run that carries every
    # contig's work plus the whole library's serial phases. And the BAM-derived estimates
    # that used to sit here (1.5x the full BAM, floor 64) sized a chunked run off an input
    # whose size does not predict its peak. Both are gone from this workflow's calls; the
    # task keeps its own formula as the fallback for anyone calling
    # subwdls/LRAA_runner.wdl directly at an arbitrary cpu.
    #
    # memoryGB and memoryGBPerWorkerScattered override the respective number.
    Int direct_memoryGB = select_first([memoryGB, memoryGB_whole_genome])
    Int scattered_memoryGB = select_first([memoryGBPerWorkerScattered, memoryGB_per_chromosome_shard])

    # DERIVED, never asked for. The only caller-supplied splice-graph evidence this
    # workflow accepts is internal_bam_for_sg, which is already depth-normalized by
    # the caller (LRAA_quant_by_cluster.wdl normalizes, merges, normalizes again), so
    # normalizing it a second time would cap already-capped depth. The two were only
    # ever meaningful as a pair.
    Boolean no_norm = defined(internal_bam_for_sg)

    # region forces `off`: a region restriction and a genome-wide partition are two
    # ways of saying the same thing, and only the unscattered path passes --region
    # through to LRAA. Validated in validate_scattering rather than silently
    # overridden here, so a caller who asked for both is told.
    Boolean scatter_by_chromosome = (scattering == "by_chromosome")
    Boolean scatter_by_chunk = (scattering == "by_chunk")
    Boolean run_without_splitting = (scattering == "off")
    String LRAA_output_suffix = if !defined(annot_gtf) && !quant_only then "LRAA.ref-free" else if defined(annot_gtf) && !quant_only then "LRAA.ref-guided" else "LRAA.quant-only"
    String LRAA_output_prefix = sample_id + "." + LRAA_output_suffix

    call validate_scattering {
        input:
            scattering = scattering,
            main_chromosomes = main_chromosomes,
            region_given = defined(region),
            chunk_plan_given = defined(internal_chunk_plan),
            docker = docker
    }

    if (scatter_by_chunk) {
        call ChunkScatter.LRAA_chunk_scatter as chunk_scatter {
            input:
                sample_id = sample_id,
                referenceGenome = referenceGenome,
                inputBAM = inputBAM,
                annot_gtf = annot_gtf,
                bam_for_sg = internal_bam_for_sg,
                bam_for_sg_index = internal_bam_for_sg_index,
                bam_for_priors = internal_bam_for_priors,
                bam_for_priors_index = internal_bam_for_priors_index,
                chunk_plan = internal_chunk_plan,
                main_chromosomes = main_chromosomes,
                oversimplify = oversimplify,
                no_weight_reads_by_3prime_agreement = no_weight_reads_by_3prime_agreement,
                quant_only = quant_only,
                HiFi = HiFi,
                min_mapping_quality = min_mapping_quality,
                min_mapping_quality_for_final_quant = min_mapping_quality_for_final_quant,
                approx_MB_per_cut = approx_MB_per_cut,
                approx_MB_per_cut_wiggle_window = approx_MB_per_cut_wiggle_window,
                cell_list = cell_list,
                cell_barcode_tag = cell_barcode_tag,
                read_umi_tag = read_umi_tag,
                stream_reads = stream_reads,
                rescue_unassigned_reads_via_transcriptome_alignment = rescue_unassigned_reads_via_transcriptome_alignment,
                num_total_reads = num_total_reads,
                chunkPreemptible = chunkPreemptible,
                # Without these six the by_chunk path was UNREACHABLE for sizing:
                # cpu/memoryGB stopped at this workflow, so every caller got the
                # subworkflow's defaults whatever it passed. On a 16-core box that
                # meant make_chunks asking 16 cpu and each of ~300 leaves 16 GiB,
                # and a run sitting in "no suitable node" until something else
                # finished. Deliberately NOT wired to this workflow's own
                # cpu/memoryGB: those size a whole-chromosome shard, which is the
                # wrong shape for a 10 Mb chunk. Unset means the measured defaults
                # in subwdls/LRAA_chunk_scatter.wdl.
                makeChunksCpu = chunkMakeChunksCpu,
                makeChunksMemoryGB = chunkMakeChunksMemoryGB,
                chunkCpu = chunkCpu,
                chunkMemoryGB = chunkMemoryGB,
                mergeCpu = chunkMergeCpu,
                mergeMemoryGB = chunkMergeMemoryGB,
                docker = docker
        }
    }

    if (scatter_by_chromosome) {

        if (!defined(num_total_reads)) {
            call count_bam {
                input:
                    bam = inputBAM,
                    samtools_threads = countBamThreads,
                    docker = docker
            }
        }

        # Unset main_chromosomes means "every contig the fasta and the bam agree
        # on", derived rather than left empty -- an empty list would partition
        # nothing.
        if (main_chromosomes == "") {
            call derive_contigs {
                input:
                    referenceGenome = referenceGenome,
                    inputBAM = inputBAM,
                    docker = docker
            }
        }
        String chromosomes_to_partition = if (main_chromosomes != "")
            then main_chromosomes
            else select_first([derive_contigs.contigs])

        Int scatter_num_total_reads = select_first([num_total_reads, count_bam.count])

        
        ## Split inputs by main chromosomes
        
        call PartByChr.partition_by_chromosome as splitByChr {
            input:
                inputBAM = inputBAM,
                bam_for_sg = internal_bam_for_sg,
                genome_fasta = referenceGenome,
                annot_gtf = annot_gtf,
                chromosomes_want_partitioned = chromosomes_to_partition,
                docker = docker,
            }
     
                  
        Int num_chromosomes = length(splitByChr.chromosomeBAMs)

        scatter (contig_index in range(num_chromosomes)) {
            String contig_name = basename(splitByChr.chromosomeBAMs[contig_index], ".bam")

            # Chunk-count ESTIMATE for THIS shard, from its own contig length -- not a
            # second scatter, not a re-read of the BAM. chromosomeFASTAs is already
            # materialised by splitByChr above, and a FASTA's byte size is sequence length
            # plus a header line and one newline per wrapped row, so this reads a few bytes
            # high rather than low. That is the SAFE direction: the true placed-cut count
            # (util/misc/select_contig_cut_points.py) can only be <= this estimate, because
            # annotation-blocked or too-short targets are declined or tail-merged AWAY, never
            # added, so a request sized off this ceiling never under-provisions a shard that
            # turns out to need more cores than guessed.
            #
            # Only the CPU REQUEST is computed here, not memory: every shard gets the same
            # memoryGB_per_chromosome_shard box whatever its contig length, since what a
            # chunked shard holds follows its chunk concurrency -- capped just below --
            # rather than the size of the chromosome it was handed.
            Float shard_contig_length_bp = size(splitByChr.chromosomeFASTAs[contig_index], "B")
            Float effective_approx_mb_per_cut = select_first([approx_MB_per_cut, 10.0])
            Int shard_chunks_estimate = ceil(shard_contig_length_bp / (effective_approx_mb_per_cut * 1000000.0))
            Int shard_chunks_estimate_floored = if shard_chunks_estimate < 1 then 1 else shard_chunks_estimate
            Int shard_cpu_computed = if shard_chunks_estimate_floored > max_cpu_per_chunked_shard
                then max_cpu_per_chunked_shard
                else shard_chunks_estimate_floored
            # Run LRAA separately per chromosome  
            call LRAA_runner.LRAA_runner as LRAA_scatter {
                input:
                    sample_id = sample_id,
                    shardno = contig_index,
                    inputBAM = splitByChr.chromosomeBAMs[contig_index],
                    bam_for_sg = if defined(internal_bam_for_sg) then select_first([splitByChr.chromosomeBAMsForSG])[contig_index] else internal_bam_for_sg,
                    # NOT partitioned per contig, unlike bam_for_sg above:
                    # subwdls/Partition_data_by_chromosome.wdl slices --bam and the sg
                    # bam, and this shard is a full LRAA invocation restricted with
                    # --contig, so it reads only this contig's reads out of the
                    # whole-genome priors bam and its own chunking slices the rest. The
                    # cost is localizing one per-cluster NORMALIZED bam per shard --
                    # small next to the library, and paying it here keeps the partition
                    # subworkflow's surface unchanged.
                    bam_for_priors = internal_bam_for_priors,
                    # Same shared geometry the by_chunk branch gets. This shard chunks
                    # its contig internally, so without the plan it would select cuts
                    # from THIS cluster's reads and slice bam_for_sg above at its own
                    # boundaries -- per-cluster geometry, which is the defect the plan
                    # removes.
                    chunk_plan = internal_chunk_plan,
                    genome_fasta = splitByChr.chromosomeFASTAs[contig_index],
                    annot_gtf = splitByChr.chromosomeGTFs[contig_index],
                    oversimplify = oversimplify,
                    contig = contig_name,
                    num_total_reads = scatter_num_total_reads,
                    cell_list = cell_list,
                    min_per_id = min_per_id,
                    quant_only = quant_only,
                    HiFi = HiFi,
                    no_norm = no_norm,
                    no_EM = no_EM,
                    no_weight_reads_by_3prime_agreement = no_weight_reads_by_3prime_agreement,
                    retain_normalized_splice_graph_bam = retain_normalized_splice_graph_bam,
                    rescue_unassigned_reads_via_transcriptome_alignment = rescue_unassigned_reads_via_transcriptome_alignment,
                    cell_barcode_tag = cell_barcode_tag,
                    read_umi_tag = read_umi_tag,
                    approx_MB_per_cut = approx_MB_per_cut,
                    approx_MB_per_cut_wiggle_window = approx_MB_per_cut_wiggle_window,
                    chunk_by_strand = chunk_by_strand,
                    stream_reads = stream_reads,
                    cpu = select_first([cpuScattered, shard_cpu_computed]),  # explicit override, else per-shard estimate above
                    min_mapping_quality = min_mapping_quality,
                    min_mapping_quality_for_final_quant = min_mapping_quality_for_final_quant,
                    docker = docker,
                    memoryGB = scattered_memoryGB,  # memoryGB_per_chromosome_shard, 16 GiB
                    diskSizeGB = diskSizeGB
            }
        }

        # Always merge quant outputs regardless of quant_only
        call mergeQuantResults {
            input:
                quantExprFiles = LRAA_scatter.LRAA_quant_expr,
                quantTrackingFiles = LRAA_scatter.LRAA_quant_tracking,
                outputFilePrefix = LRAA_output_prefix,
                docker = docker,
        }

        call mergeReadAssignmentSummaries {
            input:
                summaryFiles = select_all(LRAA_scatter.LRAA_read_assignment_summary),
                outputFilePrefix = LRAA_output_prefix,
                docker = docker
        }

        # Only merge GTFs when not in quant-only mode
        if (!quant_only) {
            call merge_GTFs {
                input:
                    gtfFiles = select_all(LRAA_scatter.LRAA_gtf),
                    outputFilePrefix = LRAA_output_prefix,
                    docker = docker,
            }
        }
    }

    if (run_without_splitting) {
            
        # main_chromosomes is a contig FILTER in every mode (see its declaration),
        # and off's only representation of one is LRAA's own --contig -- which
        # subwdls/LRAA_runner.wdl has always declared, forwarded to its task and
        # emitted. It simply was never connected here, so the input VALIDATED --
        # validate_scattering below deliberately accepts exactly one name with
        # scattering=off rather than refusing it -- and was then dropped, and
        # scattering=off with main_chromosomes="chr21" processed the whole genome
        # while reporting nothing amiss.
        #
        # A nested declaration rather than a ternary because an empty filter must
        # leave off unrestricted and WDL 1.0 has no None to select: a declaration
        # inside a false conditional is the only undefined String? it can produce.
        #
        # SUPPRESSED when region is set, mirroring the precedence LRAA already
        # applies internally: --region names its own contig and narrows within it,
        # so a region on one contig plus a --contig naming another would agree on
        # nothing and quantify nothing. region requires scattering=off, so this is
        # the one mode where the pair is reachable at all.
        if (main_chromosomes != "" && !defined(region)) {
            String direct_contig = main_chromosomes
        }

        call LRAA_runner.LRAA_runner as LRAA_direct {
            input:
                sample_id = sample_id,
                inputBAM = inputBAM,
                bam_for_sg = internal_bam_for_sg,
                bam_for_priors = internal_bam_for_priors,
                genome_fasta = referenceGenome,
                annot_gtf = annot_gtf,
                region = region,
                contig = direct_contig,
                oversimplify = oversimplify,
                min_per_id = min_per_id,
                quant_only = quant_only,
                HiFi = HiFi,
                no_norm = no_norm,
                no_EM = no_EM,
                no_weight_reads_by_3prime_agreement = no_weight_reads_by_3prime_agreement,
                retain_normalized_splice_graph_bam = retain_normalized_splice_graph_bam,
                rescue_unassigned_reads_via_transcriptome_alignment = rescue_unassigned_reads_via_transcriptome_alignment,
                cell_barcode_tag = cell_barcode_tag,
                read_umi_tag = read_umi_tag,
                approx_MB_per_cut = approx_MB_per_cut,
                approx_MB_per_cut_wiggle_window = approx_MB_per_cut_wiggle_window,
                chunk_by_strand = chunk_by_strand,
                stream_reads = stream_reads,
                cpu = cpu,
                num_total_reads = num_total_reads,
                cell_list = cell_list,
                min_mapping_quality = min_mapping_quality,
                min_mapping_quality_for_final_quant = min_mapping_quality_for_final_quant,
                docker = docker,
                memoryGB = direct_memoryGB,  # memoryGB_whole_genome, 32 GiB
                diskSizeGB = diskSizeGB
        }

    }
    
    output {
    # Three producers now, so each output names all of them. by_chunk merges
    # inside its subworkflow, by_chromosome merges here, off produces directly.
    File mergedQuantExpr = select_first([chunk_scatter.mergedQuantExpr, mergeQuantResults.mergedQuantExprFile, LRAA_direct.LRAA_quant_expr])
    File mergedQuantTracking = select_first([chunk_scatter.mergedQuantTracking, mergeQuantResults.mergedQuantTrackingFile, LRAA_direct.LRAA_quant_tracking])
    File? mergedGTF = if (!quant_only) then select_first([chunk_scatter.mergedGTF, merge_GTFs.mergedGtfFile, LRAA_direct.LRAA_gtf]) else LRAA_direct.LRAA_gtf
    # by_chunk has no per-SHARD summaries to surface: its units are chunks, and
    # stage 6 merges them inside the subworkflow. The merged table is the one
    # below, and it is the artifact single-cell consumes.
    Array[File] shardReadAssignmentSummaries = if (run_without_splitting) then select_all([LRAA_direct.LRAA_read_assignment_summary]) else select_all(select_first([LRAA_scatter.LRAA_read_assignment_summary, []]))
    File? mergedReadAssignmentSummary = if (scatter_by_chunk) then chunk_scatter.mergedReadAssignmentSummary else if (run_without_splitting) then LRAA_direct.LRAA_read_assignment_summary else mergeReadAssignmentSummaries.mergedSummaryFile
    # The depth-normalized BAM(s) the splice graph -- and therefore isoform identification --
    # was built from. Quantification does not use these; it runs against the unnormalized quant
    # BAM. Empty when internal_bam_for_sg was supplied, since that is already normalized
    # and no second normalization runs. Scattered runs normalize each chromosome shard separately,
    # so this is one BAM per shard rather than one whole-genome BAM. EMPTY in by_chunk: the
    # normalization there is per UNIT inside a chunk, and surfacing those is a separate change.
    Array[File] normalizedSpliceGraphBams = if (run_without_splitting) then select_all([LRAA_direct.LRAA_normalized_splice_graph_bam]) else select_all(select_first([LRAA_scatter.LRAA_normalized_splice_graph_bam, []]))
    Array[File] normalizedSpliceGraphBais = if (run_without_splitting) then select_all([LRAA_direct.LRAA_normalized_splice_graph_bai]) else select_all(select_first([LRAA_scatter.LRAA_normalized_splice_graph_bai, []]))
    # One per shard: the chunk manifests and per-chunk timings. Populated by the
    # two modes that run a chunked LRAA per task; by_chunk reports its partition
    # through chunkPlan and chunkLogs instead.
    Array[File] chunkReports = if (run_without_splitting) then select_all([LRAA_direct.LRAA_chunk_report]) else select_all(select_first([LRAA_scatter.LRAA_chunk_report, []]))
    # by_chunk only, empty otherwise.
    File? chunkPlan = chunk_scatter.chunkPlan
    Array[File] chunkLogs = select_first([chunk_scatter.chunkLogs, []])
    File? mergedTpmAudit = chunk_scatter.mergedTpmAudit
    File? mergeResult = chunk_scatter.mergeResult
    }
}



 

task merge_GTFs {
    input {
        Array[File] gtfFiles
        String outputFilePrefix
        String docker
    }

    command <<<
        set -eo pipefail

        lraa_version="$(LRAA --version | awk '{print $NF}')"
        version_comment="# LRAA version ${lraa_version}"

        gtf_output="~{outputFilePrefix}.gtf"
        gtf_files_str="~{sep=' ' gtfFiles}"

        lraa_merge_header.py \
            --version_comment "$version_comment" \
            --inputs $gtf_files_str \
            --output "$gtf_output"

        # the merged header above stands in for the shard headers, so drop each
        # shard's whole leading comment block, not just its version line
        awk 'FNR == 1 { in_header = 1 } in_header && /^#/ { next } { in_header = 0; print }' \
            $gtf_files_str >> "$gtf_output"
    >>>

    output {
        File mergedGtfFile = "~{outputFilePrefix}.gtf"
    }
    
    runtime {
        docker: docker
        cpu: 1
        memory: "2 GiB"
        disks: "local-disk " + ceil(size(gtfFiles, "GB") * 2.0 + 5) + " SSD"
    }
}

task mergeQuantResults {
    input {
        Array[File] quantExprFiles
        Array[File] quantTrackingFiles
        String outputFilePrefix
        String docker
    }

    Float quantExprGB = size(quantExprFiles, "GB")
    Float quantTrackingGB = size(quantTrackingFiles, "GB")
    # The merge concatenates the shard expr files and streams the shard trackings into one
    # gzip, so disk needs the localized inputs plus a copy of them: 2.2x total input.
    Float mergeDiskRawGB = (quantExprGB + quantTrackingGB) * 2.2 + 5.0

    command <<<
    set -eo pipefail

    lraa_version="$(LRAA --version | awk '{print $NF}')"
    export LRAA_VERSION_COMMENT="# LRAA version ${lraa_version}"

    lraa_merge_header.py \
        --version_comment "$LRAA_VERSION_COMMENT" \
        --inputs ~{sep=' ' quantExprFiles} \
        --output merge_header.txt

    prepend_lraa_merge_header() {
        local path="$1"
        local tmp="${path}.with_lraa_header"
        if [[ -s "$path" ]] && head -n 1 "$path" | grep -qxF "$LRAA_VERSION_COMMENT"; then
            return
        fi
        cat merge_header.txt "$path" > "$tmp"
        mv "$tmp" "$path"
    }

    merge_LRAA_quant_expr.py \
        --output "~{outputFilePrefix}.quant.expr" \
        --quant_files ~{sep=' ' quantExprFiles}

    prepend_lraa_merge_header "~{outputFilePrefix}.quant.expr"

    python <<CODE
import json
import gzip

tracking_files_json = '["' + '~{sep='","' quantTrackingFiles}' + '"]'
tracking_files_list = json.loads(tracking_files_json)
header_lines = open("merge_header.txt", "rt").read()

with gzip.open("~{outputFilePrefix}.quant.tracking.gz", "wt") as ofh:
    ofh.write(header_lines)
    wrote_header = False
    expected_header = None
    for i, tracking_file in enumerate(tracking_files_list):
        openf = gzip.open if tracking_file.split(".")[-1] == "gz" else open
        with openf(tracking_file, "rt") as fh:
            header = None
            for line in fh:
                if line.startswith("#"):
                    continue
                if header is None:
                    header = line
                    if expected_header is None:
                        expected_header = header
                    elif header != expected_header:
                        raise RuntimeError(
                            "Cannot merge tracking files with different schemas: {} differs from the first input".format(
                                tracking_file
                            )
                        )
                    if not wrote_header:
                        print(header, file=ofh, end='')
                        wrote_header = True
                    continue
                print(line, file=ofh, end='')
CODE

    >>>

    output {
        File mergedQuantExprFile = "~{outputFilePrefix}.quant.expr"
        File mergedQuantTrackingFile = "~{outputFilePrefix}.quant.tracking.gz"
    }
    
    runtime {
        docker: docker
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + ceil(mergeDiskRawGB) + " SSD"
    }
}

task mergeReadAssignmentSummaries {
    input {
        Array[File] summaryFiles
        String outputFilePrefix
        String docker
    }

    command <<<
    set -eo pipefail

    python -c '
import csv
from pathlib import Path

summary_files = [Path(p) for p in """~{sep="\n" summaryFiles}""".splitlines() if p.strip()]
out_path = Path("~{outputFilePrefix}.read_assignment.summary.tsv")

fieldnames = None
worker_rows = []
total_keys = [
    "reads_total",
    "reads_kept_genome",
    "reads_selected_tx_total",
    "reads_selected_tx_missing_genome",
    "reads_selected_tx_failed_genome",
    "reads_rescue_requested",
    "reads_rescue_rescued",
    "reads_rescue_unrescued",
    "reads_rescue_requested_failed_genome",
    "reads_rescue_requested_unassigned_quant",
    "reads_rescue_declined_locality",
    "reads_rescue_displaced_locality",
    "alignments_rescue_rejected_locality",
]
totals = {key: 0 for key in total_keys}

for summary_file in summary_files:
    with summary_file.open("rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if fieldnames is None:
            fieldnames = list(reader.fieldnames or [])
        for row in reader:
            if row.get("row_type") == "TOTAL":
                continue
            worker_rows.append(row)
            for key in total_keys:
                totals[key] += int(row.get(key, 0) or 0)

if fieldnames is None:
    fieldnames = [
        "row_type",
        "chunk_id",
        "contig_acc",
        "contig_strand",
        "reads_total",
        "reads_kept_genome",
        "frac_kept_genome",
        "reads_selected_tx_total",
        "frac_selected_tx_total",
        "reads_selected_tx_missing_genome",
        "frac_selected_tx_missing_genome",
        "reads_selected_tx_failed_genome",
        "frac_selected_tx_failed_genome",
        "reads_rescue_requested",
        "frac_rescue_requested",
        "reads_rescue_rescued",
        "frac_rescue_rescued",
        "frac_rescue_rescued_of_requested",
        "reads_rescue_unrescued",
        "frac_rescue_unrescued",
        "frac_rescue_unrescued_of_requested",
        "reads_rescue_requested_failed_genome",
        "frac_rescue_requested_failed_genome",
        "reads_rescue_requested_unassigned_quant",
        "frac_rescue_requested_unassigned_quant",
        "reads_rescue_declined_locality",
        "reads_rescue_displaced_locality",
        "alignments_rescue_rejected_locality",
        "rescue_alignment_rejections",
    ]

total_reads = totals["reads_total"]

def frac(key):
    if total_reads <= 0:
        return "0.000000"
    return f"{float(totals[key]) / float(total_reads):.6f}"

total_row = {field: "" for field in fieldnames}
total_row["row_type"] = "TOTAL"
total_row["contig_acc"] = "TOTAL"
total_row["contig_strand"] = "."
for key in total_keys:
    total_row[key] = str(totals[key])
total_row["frac_kept_genome"] = frac("reads_kept_genome")
total_row["frac_selected_tx_total"] = frac("reads_selected_tx_total")
total_row["frac_selected_tx_missing_genome"] = frac("reads_selected_tx_missing_genome")
total_row["frac_selected_tx_failed_genome"] = frac("reads_selected_tx_failed_genome")
total_row["frac_rescue_requested"] = frac("reads_rescue_requested")
total_row["frac_rescue_rescued"] = frac("reads_rescue_rescued")
total_row["frac_rescue_unrescued"] = frac("reads_rescue_unrescued")
requested = totals["reads_rescue_requested"]
rescued = totals["reads_rescue_rescued"]
unrescued = totals["reads_rescue_unrescued"]
if requested > 0:
    total_row["frac_rescue_rescued_of_requested"] = f"{float(rescued) / float(requested):.6f}"
    total_row["frac_rescue_unrescued_of_requested"] = f"{float(unrescued) / float(requested):.6f}"
else:
    total_row["frac_rescue_rescued_of_requested"] = "0.000000"
    total_row["frac_rescue_unrescued_of_requested"] = "0.000000"
total_row["frac_rescue_requested_failed_genome"] = frac("reads_rescue_requested_failed_genome")
total_row["frac_rescue_requested_unassigned_quant"] = frac("reads_rescue_requested_unassigned_quant")

with out_path.open("wt", newline="") as ofh:
    writer = csv.DictWriter(ofh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for row in worker_rows:
        writer.writerow(row)
    writer.writerow(total_row)
'
    >>>

    output {
        File mergedSummaryFile = "~{outputFilePrefix}.read_assignment.summary.tsv"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "2 GiB"
        disks: "local-disk " + ceil(size(summaryFiles, "GB") * 2.0 + 5) + " SSD"
    }
}

task validate_scattering {
    # WDL 1.0 has neither an enum nor a fail(), so the refusals live in a task
    # whose command exits non-zero. One place, named values in the message, and it
    # runs before anything expensive.
    input {
        String scattering
        String main_chromosomes
        Boolean region_given
        Boolean chunk_plan_given
        String docker
    }

    command <<<
    set -euo pipefail

    case "~{scattering}" in
        off|by_chromosome|by_chunk) ;;
        *)
            echo "Error: scattering must be off, by_chromosome or by_chunk; got '~{scattering}'" >&2
            exit 1
            ;;
    esac

    n=$(echo "~{main_chromosomes}" | wc -w)

    # region restricts to one interval, which only the unscattered path forwards
    # to LRAA. Scattering a partition across a region is two ways of saying the
    # same thing, and the combination would silently ignore one of them.
    if [[ "~{region_given}" == "true" && "~{scattering}" != "off" ]]; then
        echo "Error: region requires scattering=off; got scattering=~{scattering}" >&2
        exit 1
    fi

    # off runs ONE LRAA invocation, whose own restriction is a single --contig or
    # a --region. A list of contigs has no representation there, so it is refused
    # rather than partly honoured.
    if [[ "~{scattering}" == "off" && "${n}" -gt 1 ]]; then
        echo "Error: scattering=off accepts at most one contig in main_chromosomes, got ${n}: ~{main_chromosomes}" >&2
        exit 1
    fi

    # A chunk plan is geometry shared between SIBLING runs, and off has no siblings:
    # it is one whole-genome invocation, so there is nothing for its cut positions to
    # have to agree with. by_chunk and by_chromosome both consume it -- by_chromosome
    # because LRAA_runner passes --chunk unconditionally, so every chromosome shard
    # chunks internally and would otherwise choose its own cuts.
    #
    # Refused rather than dropped: a caller passing a plan is asserting that every
    # sibling run cuts at the SAME positions, and silently ignoring it would give each
    # run its own geometry again -- the defect the plan exists to remove, and invisible
    # in the outputs.
    if [[ "~{chunk_plan_given}" == "true" && "~{scattering}" == "off" ]]; then
        echo "Error: internal_chunk_plan requires scattering=by_chunk or by_chromosome; got scattering=off" >&2
        exit 1
    fi

    echo "scattering=~{scattering} contigs=${n} region=~{region_given} chunk_plan=~{chunk_plan_given}"
    >>>

    output {
        String checked = read_string(stdout())
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "1 GiB"
        disks: "local-disk 10 HDD"
    }
}


task count_bam {
  input {
    File bam
        Int samtools_threads = 16
        String docker
  }

    Float bam_size_gb = size(bam, "GB")
    Float estimated_disk = ceil(bam_size_gb * 2.2 + 20.0)
    Float disk_gb = if estimated_disk > 150.0 then estimated_disk else 150.0
    Int disk_gb_int = ceil(disk_gb)

  command <<<
    set -ex
        # -F 0x904 drops unmapped/secondary/supplementary so this matches
        # count_reads_from_bam() in the LRAA driver: one count per genome-mapped
        # read. Both paths must agree or TPM depends on which one ran.
        samtools view -@ ~{samtools_threads} -c -F 0x904 ~{bam}

  >>>

  runtime {
        docker: docker
        disks: "local-disk " + disk_gb_int + " SSD"
        cpu: samtools_threads
        memory: "8G"
  }
  output {
    Int count = read_int(stdout())
  }
}



task derive_contigs {
    # The contig list by_chromosome scatters over when the caller named none:
    # every reference the genome fasta and the bam header AGREE on, in fasta
    # order. Same rule ChunkedRun.enumerate_prep_contigs applies for chunking,
    # and for the same reason -- a reference absent from the bam header cannot
    # hold a record, so scattering a task for it produces an empty shard, while
    # a whole-genome fasta carries enough unplaced scaffolds to turn that into
    # hundreds of them.
    input {
        File referenceGenome
        File inputBAM
        String docker
    }

    command <<<
    set -euo pipefail

    # faidx writes beside its argument, and the input mount may be read-only
    mkdir -p work
    ln -s ~{referenceGenome} work/genome.fa
    if [[ -f "~{referenceGenome}.fai" ]]; then
        ln -s ~{referenceGenome}.fai work/genome.fa.fai
    else
        samtools faidx work/genome.fa
    fi

    cut -f1 work/genome.fa.fai > fasta_contigs.txt
    samtools view -H ~{inputBAM} | awk '$1=="@SQ"{for(i=2;i<=NF;i++) if ($i ~ /^SN:/) {sub(/^SN:/,"",$i); print $i}}' > bam_contigs.txt

    # fasta ORDER preserved, so the shard order matches what an unfiltered run
    # would have produced
    awk 'NR==FNR{seen[$0]=1; next} ($0 in seen)' bam_contigs.txt fasta_contigs.txt > kept.txt

    if [[ ! -s kept.txt ]]; then
        echo "Error: no reference in ~{referenceGenome} appears in the bam header, so there is nothing to scatter over" >&2
        exit 1
    fi
    paste -sd' ' kept.txt
    >>>

    output {
        String contigs = read_string(stdout())
        Int num_contigs = length(read_lines("kept.txt"))
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "2 GiB"
        disks: "local-disk 20 HDD"
    }
}
