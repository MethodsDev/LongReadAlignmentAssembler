version 1.0

task LRAA_runner_task {
    input {
        String sample_id
        File genome_fasta
        File inputBAM

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
        Int? min_mapping_quality
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
        Int memoryGB = 32
        Int diskSizeGB = 128
        Int progress_report_interval_seconds = 300
        Int progress_tail_lines = 20
    }

    String no_norm_flag = if (no_norm && !quant_only) then "--no_norm" else ""
    String no_EM_flag = if (no_EM) then "--no_EM" else ""

    String output_prefix_use = if defined(shardno) then "${sample_id}.shardno-${shardno}" else sample_id
    
    String output_suffix = if !defined(annot_gtf) && !quant_only then "LRAA.ref-free" else if defined(annot_gtf) && !quant_only then "LRAA.ref-guided" else "LRAA.quant-only"
    
    command <<<

        set -ex

        : > command_output.log

        progress_reporter() {
            while sleep ~{progress_report_interval_seconds}; do
                if [[ -s command_output.log ]]; then
                    echo "----- LRAA progress $(date) -----" >&2
                    tail -n ~{progress_tail_lines} command_output.log >&2 || true
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
                                 --output_prefix ~{output_prefix_use}.~{output_suffix} \
                                 ~{if defined(contig) then "--contig " + contig else ""} \
                                 ~{if defined(region) then "--region " + region else ""} \
                                 ~{if defined(oversimplify) then "--oversimplify " + oversimplify else ""} \
                                 ~{if defined(min_per_id) then "--min_per_id " + min_per_id else ""} \
                                 ~{no_norm_flag} \
                                 ~{no_EM_flag} \
                                 --num_threads_per_worker ~{numThreadsPerWorker} \
                                 ~{if defined(min_mapping_quality) then "--min_mapping_quality " + min_mapping_quality else ""} \
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
            tail -n 100 command_output.log >&2 || true
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
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }

}


workflow LRAA_runner {
    input {
        String sample_id
        File genome_fasta
        File inputBAM

        File? annot_gtf
        Boolean quant_only
        Boolean HiFi = false
        String? region
        String? oversimplify
        # Expose toggle to workflow as well; default false here and override to true from scatter callers.
        Boolean no_parallelize_contigs = false
        String? contig
        Int? num_parallel_contigs
        
        Int? num_total_reads
        Float? min_per_id
        Boolean no_EM 
        Boolean no_norm 
        Int? min_mapping_quality
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
    
        Int memoryGB = 32
        Int diskSizeGB = 128
        Int progress_report_interval_seconds = 300
        Int progress_tail_lines = 20
    }

    call LRAA_runner_task {
        input:
            sample_id=sample_id,
            genome_fasta=genome_fasta,
            inputBAM=inputBAM,
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
            min_mapping_quality=min_mapping_quality,
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
            progress_tail_lines=progress_tail_lines
     }

     output {
        File? LRAA_gtf = LRAA_runner_task.LRAA_gtf
        File LRAA_quant_expr = LRAA_runner_task.LRAA_quant_expr
        File LRAA_quant_tracking = LRAA_runner_task.LRAA_quant_tracking
    }

}

