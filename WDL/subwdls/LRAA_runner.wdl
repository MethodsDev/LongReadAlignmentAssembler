version 1.0

task LRAA_runner_task {
    input {
        String sample_id
        File genome_fasta
        File inputBAM

        File? annot_gtf
        Boolean quant_only
        Boolean LowFi = false
        String? region
        
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
        Int numThreads
        Int memoryGB = 32
        Int diskSizeGB = 128
    }

    String no_norm_flag = if (no_norm && !quant_only) then "--no_norm" else ""
    String no_EM_flag = if (no_EM) then "--no_EM" else ""

    String output_prefix_use = if defined(shardno) then "${sample_id}.shardno-${shardno}" else sample_id
    
    String output_suffix = if !defined(annot_gtf) && !quant_only then "LRAA.ref-free" else if defined(annot_gtf) && !quant_only then "LRAA.ref-guided" else "LRAA.quant-only"
    
    command <<<

        set -ex

        (        
        LRAA --genome ~{genome_fasta} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{output_prefix_use}.~{output_suffix} \
                                ~{"--region " + region} \
                                 ~{"--min_per_id " + min_per_id} \
                                 ~{no_norm_flag} \
                                 ~{no_EM_flag} \
                                 --CPU ~{numThreads} \
                                 ~{"--min_mapping_quality " + min_mapping_quality } \
                                 ~{"--min_isoform_fraction " + min_isoform_fraction} \
                                 ~{"--min_monoexonic_TPM " + min_monoexonic_TPM} \
                                 ~{true="--no_filter_internal_priming" false='' no_filter_internal_priming} \
                                 ~{"--min_alt_splice_freq " + min_alt_splice_freq} \
                                 ~{"--min_alt_unspliced_freq " + min_alt_unspliced_freq} \
                                 ~{"--gtf " + annot_gtf} \
                                 ~{"--num_total_reads " + num_total_reads} \
                                 ~{true="--quant_only" false='' quant_only} \
                                 ~{true="--LowFi" false='' LowFi} \
                                 ~{"--cell_barcode_tag " + cell_barcode_tag} ~{"--read_umi_tag " + read_umi_tag} \
                  > command_output.log 2>&1
        ) || {
             echo "Command failed with exit code $?" >&2
             echo "Last 100 lines of output:" >&2
             tail -n 100 command_output.log >&2
             exit 1
        }

        if [[ -f ~{output_prefix_use}.~{output_suffix}.quant.tracking ]]; then
            gzip ~{output_prefix_use}.~{output_suffix}.quant.tracking    
        fi
    
        # always ensure an output file exists for the wdl output capture.
        touch ~{output_prefix_use}.~{output_suffix}.gtf

        
    >>>

    output {
        File LRAA_gtf = "~{output_prefix_use}.~{output_suffix}.gtf"
        File LRAA_quant_expr = "~{output_prefix_use}.~{output_suffix}.quant.expr"
        File LRAA_quant_tracking = "~{output_prefix_use}.~{output_suffix}.quant.tracking.gz"
    }


    runtime {
        docker: docker
        bootDiskSizeGb: 30
        cpu: "~{numThreads}"
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
        Boolean LowFi = false
        String? region
        
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
        Int numThreads
    
        Int memoryGB = 32
        Int diskSizeGB = 128
    }

    call LRAA_runner_task {
        input:
            sample_id=sample_id,
            genome_fasta=genome_fasta,
            inputBAM=inputBAM,
            annot_gtf=annot_gtf,
            quant_only=quant_only,
            LowFi = LowFi,
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
            numThreads=numThreads,
            memoryGB=memoryGB,
            diskSizeGB=diskSizeGB,
            region=region
     }

     output {
        File LRAA_gtf = LRAA_runner_task.LRAA_gtf
        File LRAA_quant_expr = LRAA_runner_task.LRAA_quant_expr
        File LRAA_quant_tracking = LRAA_runner_task.LRAA_quant_tracking
    }

}

