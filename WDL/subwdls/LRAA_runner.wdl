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
        Int? num_parallel_contigs
        
        Int? num_total_reads
        Float? min_per_id
        Boolean no_EM 
        Boolean no_norm 
        Boolean allow_secondary_alignments = true
        Int min_mapping_quality = 0
        Int min_mapping_quality_for_final_quant = 0
        Float? min_isoform_fraction
        Float? min_monoexonic_TPM
        Boolean? no_filter_internal_priming
        Float? min_alt_splice_freq
        Float? min_alt_unspliced_freq

        String cell_barcode_tag = "CB"
        String read_umi_tag = "XM"

        Int? shardno
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        # CPU cores per contig worker (passed to --num_threads_per_worker)
        Int numThreadsPerWorker
        Int? memoryGB
        Int diskSizeGB = 128
        Int progress_report_interval_seconds = 300
        Int progress_tail_lines = 20
        Int progress_tail_chars = 5000
    }

    # Dynamic memory for per-chromosome shards: use the larger of two estimates.
    #  - 25x BAM size preserves the prior behavior for larger shards.
    #  - A bounded mid-small-shard bump gives extra headroom to shards near 0.7-1.6 GiB
    #    without raising truly small shards above the 32 GiB floor.
    # For direct (non-scattered) runs the caller computes and passes memoryGB explicitly to override.
    Float bam_size_gib = size(inputBAM, "GiB")
    Float mem_raw_size = 25.0 * bam_size_gib
    Float mem_raw_mid_small_uncapped = if bam_size_gib > 0.5 then 32.0 + (40.0 * (bam_size_gib - 0.5)) else 32.0
    Float mem_raw_mid_small = if mem_raw_mid_small_uncapped > 40.0 then 40.0 else mem_raw_mid_small_uncapped
    Float mem_raw = if mem_raw_size > mem_raw_mid_small then mem_raw_size else mem_raw_mid_small
    Int computed_memoryGB = if mem_raw > 32.0 then ceil(mem_raw) else 32
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
                                 --num_threads_per_worker ~{numThreadsPerWorker} \
                                 ~{true='' false='--no_allow_secondary_alignments' allow_secondary_alignments} \
                                 ~{"--min_mapping_quality " + min_mapping_quality} \
                                 ~{"--min_mapping_quality_for_final_quant " + min_mapping_quality_for_final_quant} \
                                 ~{if defined(min_isoform_fraction) then "--min_isoform_fraction " + min_isoform_fraction else ""} \
                                 ~{if defined(min_monoexonic_TPM) then "--min_monoexonic_TPM " + min_monoexonic_TPM else ""} \
                                 ~{true="--no_filter_internal_priming" false='' no_filter_internal_priming} \
                                 ~{if defined(min_alt_splice_freq) then "--min_alt_splice_freq " + min_alt_splice_freq else ""} \
                                 ~{if defined(min_alt_unspliced_freq) then "--min_alt_unspliced_freq " + min_alt_unspliced_freq else ""} \
                                 ~{if defined(annot_gtf) then "--gtf " + annot_gtf else ""} \
                                 ~{if defined(num_total_reads) then "--num_total_reads " + num_total_reads else ""} \
                                 ~{true="--quant_only" false='' quant_only} \
                                 ~{true="--HiFi" false='' HiFi} \
                                 ~{true="--no_parallelize_contigs" false='' no_parallelize_contigs} \
                                 ~{if defined(num_parallel_contigs) then "--num_parallel_contigs " + num_parallel_contigs else ""} \
                                 ~{"--cell_barcode_tag " + cell_barcode_tag} ~{"--read_umi_tag " + read_umi_tag} \
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
            touch ~{output_prefix_use}.~{output_suffix}.gtf
        fi

        
    >>>

    output {
        File? LRAA_gtf = "~{output_prefix_use}.~{output_suffix}.gtf"
        File LRAA_quant_expr = "~{output_prefix_use}.~{output_suffix}.quant.expr"
        File LRAA_quant_tracking = "~{output_prefix_use}.~{output_suffix}.quant.tracking.gz"
    }


    runtime {
        docker: docker
        bootDiskSizeGb: 30
        cpu: "~{numThreadsPerWorker}"
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
        Int? num_parallel_contigs
        
        Int? num_total_reads
        Float? min_per_id
        Boolean no_EM 
        Boolean no_norm 
        Boolean allow_secondary_alignments = true
        Int min_mapping_quality = 0
        Int min_mapping_quality_for_final_quant = 0
        Float? min_isoform_fraction
        Float? min_monoexonic_TPM
        Boolean? no_filter_internal_priming
        Float? min_alt_splice_freq
        Float? min_alt_unspliced_freq

        String cell_barcode_tag = "CB"
        String read_umi_tag = "XM"
                    
        Int? shardno
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        # CPU cores per contig worker (passed to --num_threads_per_worker)
        Int numThreadsPerWorker
    
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
            num_parallel_contigs = num_parallel_contigs,
            num_total_reads=num_total_reads,
            min_per_id=min_per_id,
            no_EM=no_EM, 
            no_norm=no_norm,
            allow_secondary_alignments=allow_secondary_alignments,
            min_mapping_quality=min_mapping_quality,
            min_mapping_quality_for_final_quant=min_mapping_quality_for_final_quant,
            min_isoform_fraction=min_isoform_fraction,
            min_monoexonic_TPM=min_monoexonic_TPM,
            no_filter_internal_priming=no_filter_internal_priming,
            min_alt_splice_freq=min_alt_splice_freq,
            min_alt_unspliced_freq=min_alt_unspliced_freq,
            cell_barcode_tag = cell_barcode_tag,
            read_umi_tag = read_umi_tag,
            shardno=shardno,
            docker=docker,
            numThreadsPerWorker=numThreadsPerWorker,
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
    }

}
