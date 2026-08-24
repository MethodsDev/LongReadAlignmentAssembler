version 1.0

task LRAA_runner_task {
    input {
        String sample_id
        File genome_fasta
        File inputBAM
        File? bam_for_sg

        File? annot_gtf
        Boolean quant_only
        Boolean HiFi = false
        String? region
        String? oversimplify
        # Optional: disable contig-level parallelization inside LRAA.
        # Keep default = false for non-scattered runs; set to true in scatter contexts to avoid oversubscription.
        Boolean no_parallelize_contigs = false
        String? contig
        
        Int? num_total_reads
        File? cell_list
        Float? min_per_id
        Boolean no_EM 
        Boolean no_norm 
        # Export the depth-normalized splice-graph BAM as a task output. Bulk runs keep it;
        # single-cell callers turn it off so they do not pay to delocalize a BAM they never surface.
        Boolean retain_normalized_splice_graph_bam = true
        Boolean rescue_unassigned_reads_via_transcriptome_alignment = true
        Int min_mapping_quality = 0
        Int min_mapping_quality_for_final_quant = 0
        Float? min_isoform_fraction
        Float? min_monoexonic_TPM
        Boolean? no_filter_internal_priming
        Float? min_alt_splice_freq
        Float? min_alt_unspliced_freq

        String cell_barcode_tag = "CB"
        String read_umi_tag = "XM"

        # Chunked mode is UNCONDITIONAL -- there is no `chunk` input, and the command
        # below always passes --chunk. A workspace still binding it will fail on an
        # unknown input rather than quietly change modes.
        #
        # Chunking splits each contig-strand at low-coverage positions between annotated
        # loci, quantifies the chunks concurrently under one --cpu_budget, and merges. It
        # is what lets a task go faster than its longest contig, and it bounds peak memory
        # by the chunk size instead of the contig size -- measured on chr1 at budget 8,
        # 3.58 GB peak against 10.03 GB unchunked, 2.8x wall. That second property is why
        # the toggle is gone: with one mode, memory follows chunk concurrency and chunk
        # width, and nothing here has to be a function of input size.
        #
        # Nothing is lost by forcing it. All three modes chunk: quant_only with annot_gtf
        # is quant-only, annot_gtf without quant_only is ref-guided discovery, neither is
        # de novo (see output_suffix below, which has always distinguished the three).
        # Chunked discovery agrees with unchunked EXACTLY on chr21 -- 1460 = 1460 models,
        # 0 of 11,811 GTF rows differing, in both rescue configurations. Single-cell runs
        # chunk too: cell_list and the barcode/UMI tags are forwarded into the chunk
        # workers. LRAA's one refusal is --quant_only WITHOUT --gtf, which has nothing to
        # quantify and which an unchunked run refuses as well.
        #
        # The per-chunk intermediates go to <output_prefix>.chunked_work in the task's own
        # working directory, which is writable and is not delocalized. They can be several
        # times the input BAM, so diskSizeGB wants the same headroom it always did.

        # Approximate MEGABASES between cuts. Larger chunks, fewer of them, higher peak
        # memory per chunk. LRAA's default is 10; a contig shorter than this yields no
        # cuts at all and chunking degenerates to a single-chunk run of the whole thing.
        Float? approx_MB_per_cut
        # TOTAL width in MEGABASES of the search window centred on each target cut, not
        # the half-width. LRAA's default is 1, so it must come down alongside a small
        # approx_MB_per_cut or neighbouring windows overlap.
        Float? approx_MB_per_cut_wiggle_window
        # Chunk each contig-STRAND separately, splitting the whole BAM by orientation
        # first as a serial phase. OFF by default, because the strandless ordering --
        # split inside each chunk, concurrently with every other chunk -- removes the
        # largest serial phase a chunked run has: 151.2 s against 255.2 s on the same
        # input.
        #
        # Named for what setting it DOES, not for what it turns off. The input this
        # replaces was `strandless_chunks`, which had to be read as a double negative
        # once strandless became the default, and whose false case emitted no flag at
        # all -- so asking for strand-first silently got you strandless.
        #
        # Neither ordering needs num_total_reads from the caller: LRAA counts the
        # library itself, before any chunking, with the -F 0x904 policy that matches
        # the unchunked path.
        Boolean chunk_by_strand = false

        # Two-pass streaming quantification. ON by default, matching LRAA's own
        # v0.25.0 default. Both true and false are emitted explicitly below
        # (--stream_reads / --no_stream_reads), same reasoning as chunk above.
        # Transcriptome-rescue-inside-streaming is left to LRAA's own dynamic
        # default (tracks whether rescue_unassigned_reads_via_transcriptome_alignment
        # is itself on) rather than surfaced as a separate WDL input.
        Boolean stream_reads = true

        Int? shardno
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:latest"
        # TOTAL cores for this task: both the runtime cpu request and the --cpu_budget
        # that LRAA divides across its work units. One declaration, so the cores the
        # task is given and the cores LRAA believes it has cannot disagree. This is what
        # numThreadsPerWorker and num_parallel_contigs used to MULTIPLY into.
        Int cpu
        Int? memoryGB
        Int diskSizeGB = 128
        Int progress_report_interval_seconds = 300
        Int progress_tail_lines = 20
        Int progress_tail_chars = 5000
    }

    # Memory when the caller does not say. Peak is the number of chunks running at once
    # times what one chunk holds, and BOTH sides of that are bounded by cpu: chunk
    # concurrency is derived from --cpu_budget, and a chunk holds only its own extracted
    # mini-contig, whose length is approx_MB_per_cut and not the chromosome's. So this is
    # sized off cpu, not off the input BAM -- BAM size does not predict a chunked run's
    # peak, and sizing off it asked for 40+ GiB to do work measured at 3.58 GB.
    #
    # 2 GiB/core against 0.45 GiB/core measured (3.58 GB at cpu_budget 8, chr1, default
    # 10 Mb cuts) -- 4x, off one corpus, with a 16 GiB floor for the serial phases and
    # the final merge. A caller who raises approx_MB_per_cut far past the default makes
    # each chunk proportionally bigger and should pass memoryGB, which still overrides.
    #
    # THIS IS A FALLBACK for callers of this subworkflow at an arbitrary cpu. LRAA.wdl
    # does not use it: it caps cpu per shard and passes memoryGB explicitly -- 32 GiB
    # whole-genome, 16 GiB per chromosome shard -- because being proportional to cpu made
    # this formula give the whole-genome run (5 cores) LESS than the biggest shards (16
    # cores), which is backwards.
    Float mem_raw = 2.0 * cpu
    Int computed_memoryGB = if mem_raw > 16.0 then ceil(mem_raw) else 16
    Int effective_memoryGB = select_first([memoryGB, computed_memoryGB])

    String no_norm_flag = if (no_norm) then "--no_norm" else ""
    String no_EM_flag = if (no_EM) then "--no_EM" else ""


    String output_prefix_use = if defined(shardno) then "${sample_id}.shardno-${shardno}" else sample_id
    
    String output_suffix = if !defined(annot_gtf) && !quant_only then "LRAA.ref-free" else if defined(annot_gtf) && !quant_only then "LRAA.ref-guided" else "LRAA.quant-only"
    
    command <<<

        set -e

        : > command_output.log

        emit_log_tail() {
            local log_file="$1"
            local tail_lines="$2"

            if [[ ! -s "$log_file" ]]; then
                return 0
            fi

            tail -c ~{progress_tail_chars} "$log_file" | tr '\r' '\n' | tail -n "$tail_lines" >&2 || true
        }

        progress_reporter() {
            local contig_prefix="~{output_prefix_use}.~{output_suffix}"
            local -a contig_tmp_dirs=(
                "__${contig_prefix}.contigtmp"
                "${contig_prefix}.contigtmp"
            )

            while sleep ~{progress_report_interval_seconds}; do
                if [[ -s command_output.log ]]; then
                    echo "----- LRAA progress $(date) -----" >&2

                    if [[ -r /proc/meminfo ]]; then
                        local mem_stats
                        mem_stats=$(awk '
                            /MemTotal:/ {total=$2}
                            /MemAvailable:/ {avail=$2}
                            /MemFree:/ {free=$2}
                            END {
                                if (!avail && free) avail = free
                                if (total && avail) printf "%d %d\n", total, avail
                            }
                        ' /proc/meminfo) || mem_stats=""

                        if [[ -n "$mem_stats" ]]; then
                            local mem_total_kb mem_available_kb
                            read -r mem_total_kb mem_available_kb <<< "$mem_stats"

                            if [[ -n "$mem_total_kb" && -n "$mem_available_kb" ]]; then
                                local mem_used_kb=$((mem_total_kb - mem_available_kb))
                                local mem_total_mb=$((mem_total_kb / 1024))
                                local mem_used_mb=$((mem_used_kb / 1024))
                                if (( mem_total_kb > 0 )); then
                                    local mem_pct=$((mem_used_kb * 100 / mem_total_kb))
                                    echo "RAM usage: ${mem_used_mb} MiB / ${mem_total_mb} MiB (${mem_pct}%)" >&2
                                else
                                    echo "RAM usage: ${mem_used_mb} MiB" >&2
                                fi
                            fi
                        fi
                    fi

                    local printed_worker_progress=0
                    local contig_tmp_dir
                    for contig_tmp_dir in "${contig_tmp_dirs[@]}"; do
                        if [[ ! -d "$contig_tmp_dir" ]]; then
                            continue
                        fi

                        local nullglob_restore="shopt -u nullglob"
                        nullglob_restore=$(shopt -p nullglob 2>/dev/null) || nullglob_restore="shopt -u nullglob"
                        shopt -s nullglob
                        local -a err_logs=("$contig_tmp_dir"/*/*/*.err.log)
                        eval "$nullglob_restore"

                        if (( ${#err_logs[@]} == 0 )); then
                            continue
                        fi

                        if (( !printed_worker_progress )); then
                            echo "----- contig worker progress -----" >&2
                            printed_worker_progress=1
                        fi

                        local err_log
                        for err_log in "${err_logs[@]}"; do
                            echo "=== ${err_log} ===" >&2
                            emit_log_tail "$err_log" ~{progress_tail_lines}
                        done

                        break
                    done

                    if (( printed_worker_progress )); then
                        echo "----- end contig worker progress -----" >&2
                    fi

                    emit_log_tail command_output.log ~{progress_tail_lines}
                    echo "----- end progress -----" >&2
                fi
            done
        }

        progress_reporter &
        progress_pid=$!

        set +e
        (        
        LRAA --genome ~{genome_fasta} \
                                 --bam ~{inputBAM} \
                                 ~{if defined(bam_for_sg) then "--bam_for_sg " + bam_for_sg else ""} \
                                 --output_prefix ~{output_prefix_use}.~{output_suffix} \
                                 ~{if defined(contig) then "--contig " + contig else ""} \
                                 ~{if defined(region) then "--region " + region else ""} \
                                 ~{if defined(oversimplify) then "--oversimplify " + oversimplify else ""} \
                                 ~{if defined(min_per_id) then "--min_per_id " + min_per_id else ""} \
                                 ~{no_norm_flag} \
                                 ~{no_EM_flag} \
                                 --cpu_budget ~{cpu} \
                                 ~{true='' false='--no_rescue_unassigned_reads_via_transcriptome_alignment' rescue_unassigned_reads_via_transcriptome_alignment} \
                                 ~{"--min_mapping_quality " + min_mapping_quality} \
                                 ~{"--min_mapping_quality_for_final_quant " + min_mapping_quality_for_final_quant} \
                                 ~{if defined(min_isoform_fraction) then "--min_isoform_fraction " + min_isoform_fraction else ""} \
                                 ~{if defined(min_monoexonic_TPM) then "--min_monoexonic_TPM " + min_monoexonic_TPM else ""} \
                                 ~{true="--no_filter_internal_priming" false='' no_filter_internal_priming} \
                                 ~{if defined(min_alt_splice_freq) then "--min_alt_splice_freq " + min_alt_splice_freq else ""} \
                                 ~{if defined(min_alt_unspliced_freq) then "--min_alt_unspliced_freq " + min_alt_unspliced_freq else ""} \
                                 ~{if defined(annot_gtf) then "--gtf " + annot_gtf else ""} \
                                 ~{if defined(num_total_reads) then "--num_total_reads " + num_total_reads else ""} \
                                 ~{if defined(cell_list) then "--cell_list " + cell_list else ""} \
                                 ~{true="--quant_only" false='' quant_only} \
                                 ~{true="--HiFi" false='' HiFi} \
                                 ~{true="--no_parallelize_contigs" false='' no_parallelize_contigs} \
                                 --chunk \
                                 ~{true="--chunk_by_strand" false='' chunk_by_strand} \
                                 ~{if defined(approx_MB_per_cut) then "--approx_MB_per_cut " + approx_MB_per_cut else ""} \
                                 ~{if defined(approx_MB_per_cut_wiggle_window) then "--approx_MB_per_cut_wiggle_window " + approx_MB_per_cut_wiggle_window else ""} \
                                 ~{true="--stream_reads" false="--no_stream_reads" stream_reads} \
                                 ~{if (cell_barcode_tag != "CB") then "--cell_barcode_tag " + cell_barcode_tag else ""} ~{if (read_umi_tag != "XM") then "--read_umi_tag " + read_umi_tag else ""} \
                  > command_output.log 2>&1
        )
        cmd_status=$?
        set -e

        kill $progress_pid 2>/dev/null || true
        wait $progress_pid 2>/dev/null || true

        if [[ $cmd_status -ne 0 ]]; then
            echo "Command failed with exit code $cmd_status" >&2
            echo "Last 100 lines of output:" >&2
            emit_log_tail command_output.log 100
            exit $cmd_status
        fi

        if [[ -f ~{output_prefix_use}.~{output_suffix}.quant.tracking ]]; then
            gzip ~{output_prefix_use}.~{output_suffix}.quant.tracking    
        fi

    
        # only create GTF file when not in quant-only mode
        if [[ "~{quant_only}" != "true" ]]; then
            gtf_out="~{output_prefix_use}.~{output_suffix}.gtf"
            if [[ ! -s "$gtf_out" ]]; then
                lraa_version="$(LRAA --version | awk '{print $NF}')"
                printf '# LRAA version %s\n' "$lraa_version" > "$gtf_out"
            fi
        fi


        # The depth-normalized BAM that the splice graph -- and therefore isoform
        # identification -- is built from. LRAA leaves it in its own cache dir, which is not
        # delocalized. Only the top level of the cache holds the finished file; work_*/ holds
        # per-strand intermediates. Moved rather than copied: the run is over, and these can
        # be large enough for a second copy to matter against the task disk.
        if [[ "~{retain_normalized_splice_graph_bam}" == "true" ]]; then
            normalized_sg_bam_out="~{output_prefix_use}.~{output_suffix}.splice_graph_normalized.bam"
            shopt -s nullglob
            normalized_sg_bams=(__*.norm_cache/*.norm_*.bam)
            shopt -u nullglob
            if (( ${#normalized_sg_bams[@]} > 0 )); then
                if [[ -f "${normalized_sg_bams[0]}.bai" ]]; then
                    mv "${normalized_sg_bams[0]}.bai" "${normalized_sg_bam_out}.bai"
                fi
                mv "${normalized_sg_bams[0]}" "$normalized_sg_bam_out"
                if [[ ! -f "${normalized_sg_bam_out}.bai" ]]; then
                    samtools index -@ ~{cpu} "$normalized_sg_bam_out"
                fi
            fi
        fi

        # The one thing worth keeping out of the chunk work directory: the record of what
        # the run actually partitioned into. It carries one entry per chunk under
        # chunk_manifests, plus per-chunk wall times and the concurrency the budget was
        # split into. A few KB per chunk, against a work directory that holds every
        # chunk's BAM, so the report is lifted out and the directory itself is left behind
        # undelocalized.
        chunk_report_out="~{output_prefix_use}.~{output_suffix}.chunk_report.json"
        chunk_timing="~{output_prefix_use}.~{output_suffix}.chunked_work/timing.json"
        if [[ -s "$chunk_timing" ]]; then
            cp "$chunk_timing" "$chunk_report_out"
        else
            echo "WARNING: chunked run left no timing.json at $chunk_timing" >&2
        fi

        
    >>>

    output {
        File? LRAA_gtf = "~{output_prefix_use}.~{output_suffix}.gtf"
        File LRAA_quant_expr = "~{output_prefix_use}.~{output_suffix}.quant.expr"
        File LRAA_quant_tracking = "~{output_prefix_use}.~{output_suffix}.quant.tracking.gz"
        File? LRAA_read_assignment_summary = "~{output_prefix_use}.~{output_suffix}.read_assignment.summary.tsv"
        File? LRAA_normalized_splice_graph_bam = "~{output_prefix_use}.~{output_suffix}.splice_graph_normalized.bam"
        File? LRAA_normalized_splice_graph_bai = "~{output_prefix_use}.~{output_suffix}.splice_graph_normalized.bam.bai"
        # The chunk manifests and per-chunk timings. Optional only because LRAA can leave
        # no timing.json behind; every run chunks.
        File? LRAA_chunk_report = "~{output_prefix_use}.~{output_suffix}.chunk_report.json"
    }


    runtime {
        docker: docker
        bootDiskSizeGb: 30
        cpu: cpu
        memory: "~{effective_memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }

}


workflow LRAA_runner {
    input {
        String sample_id
        File genome_fasta
        File inputBAM
        File? bam_for_sg

        File? annot_gtf
        Boolean quant_only
        Boolean HiFi = false
        String? region
        String? oversimplify
        # Expose toggle to workflow as well; default on to match current LRAA quant defaults.
        Boolean no_parallelize_contigs = false
        String? contig
        
        Int? num_total_reads
        File? cell_list
        Float? min_per_id
        Boolean no_EM
        Boolean no_norm
        Boolean retain_normalized_splice_graph_bam = true
        Boolean rescue_unassigned_reads_via_transcriptome_alignment = true
        Int min_mapping_quality = 0
        Int min_mapping_quality_for_final_quant = 0
        Float? min_isoform_fraction
        Float? min_monoexonic_TPM
        Boolean? no_filter_internal_priming
        Float? min_alt_splice_freq
        Float? min_alt_unspliced_freq

        String cell_barcode_tag = "CB"
        String read_umi_tag = "XM"

        # Chunked mode is unconditional; see the task's inputs. approx_MB_per_cut and
        # chunk_by_strand tune it -- strandless is the default ordering, as in LRAA.
        Float? approx_MB_per_cut
        Float? approx_MB_per_cut_wiggle_window
        Boolean chunk_by_strand = false
        # Two-pass streaming quantification; see the task's input for detail. ON by
        # default since v0.25.0, matching LRAA's own default.
        Boolean stream_reads = true
                    
        Int? shardno
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:latest"
        # TOTAL cores for the run: the task's cpu request and its --cpu_budget.
        Int cpu
    
        Int? memoryGB
        Int diskSizeGB = 128
        Int progress_report_interval_seconds = 300
        Int progress_tail_lines = 20
        Int progress_tail_chars = 5000
    }

    call LRAA_runner_task {
        input:
            sample_id=sample_id,
            genome_fasta=genome_fasta,
            inputBAM=inputBAM,
            bam_for_sg=bam_for_sg,
            annot_gtf=annot_gtf,
            quant_only=quant_only,
            HiFi = HiFi,
            oversimplify = oversimplify,
            no_parallelize_contigs = no_parallelize_contigs,
            contig = contig,
            num_total_reads=num_total_reads,
            cell_list=cell_list,
            min_per_id=min_per_id,
            no_EM=no_EM,
            no_norm=no_norm,
            retain_normalized_splice_graph_bam=retain_normalized_splice_graph_bam,
            rescue_unassigned_reads_via_transcriptome_alignment=rescue_unassigned_reads_via_transcriptome_alignment,
            min_mapping_quality=min_mapping_quality,
            min_mapping_quality_for_final_quant=min_mapping_quality_for_final_quant,
            min_isoform_fraction=min_isoform_fraction,
            min_monoexonic_TPM=min_monoexonic_TPM,
            no_filter_internal_priming=no_filter_internal_priming,
            min_alt_splice_freq=min_alt_splice_freq,
            min_alt_unspliced_freq=min_alt_unspliced_freq,
            cell_barcode_tag = cell_barcode_tag,
            read_umi_tag = read_umi_tag,
            approx_MB_per_cut = approx_MB_per_cut,
            approx_MB_per_cut_wiggle_window = approx_MB_per_cut_wiggle_window,
            chunk_by_strand = chunk_by_strand,
            stream_reads = stream_reads,
            shardno=shardno,
            docker=docker,
            cpu=cpu,
            memoryGB=memoryGB,
            diskSizeGB=diskSizeGB,
            region=region,
            progress_report_interval_seconds=progress_report_interval_seconds,
            progress_tail_lines=progress_tail_lines,
            progress_tail_chars=progress_tail_chars
     }

     output {
        File? LRAA_gtf = LRAA_runner_task.LRAA_gtf
        File LRAA_quant_expr = LRAA_runner_task.LRAA_quant_expr
        File LRAA_quant_tracking = LRAA_runner_task.LRAA_quant_tracking
        File? LRAA_read_assignment_summary = LRAA_runner_task.LRAA_read_assignment_summary
        File? LRAA_normalized_splice_graph_bam = LRAA_runner_task.LRAA_normalized_splice_graph_bam
        File? LRAA_normalized_splice_graph_bai = LRAA_runner_task.LRAA_normalized_splice_graph_bai
        File? LRAA_chunk_report = LRAA_runner_task.LRAA_chunk_report
    }

}
