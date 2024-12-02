version 1.0

# This task uses LRAA latest version
task lraaTask {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
        String dataType
        String ID_or_Quant_or_Both
        Int? LRAA_min_mapping_quality
        Boolean? LRAA_no_norm
        Int cpu = 4
        Int numThreads = 8
        Int memoryGB = 64
        Int diskSizeGB = 2048
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    String OutDir = "LRAA_out"
    String LRAA_min_mapping_quality_flag = if (defined(LRAA_min_mapping_quality)) then "--min_mapping_quality=" + LRAA_min_mapping_quality else ""

    String no_norm_flag = if (defined(LRAA_no_norm) && LRAA_no_norm) then "--no_norm" else ""


    command <<<
        bash ~{monitoringScript} > monitoring.log &
        
        mkdir -p ~{OutDir}/ID ~{OutDir}/ID_reduced ~{OutDir}/Quant ~{OutDir}/Quant_noEM ~{OutDir}/Quant_minMapQ ~{OutDir}/Quant_noEM_minMapQ

        if [[ "~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/ID/LRAA \
                                 ~{no_norm_flag} --CPU 1

        fi

        if [[ ("~{ID_or_Quant_or_Both}" == "ID" || "~{ID_or_Quant_or_Both}" == "Both") && -n "~{referenceAnnotation_reduced}" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/ID_reduced/LRAA_reduced \
                                 ~{no_norm_flag} \
                                 --gtf ~{referenceAnnotation_reduced} --CPU 1 


        fi

        if [[ ("~{ID_or_Quant_or_Both}" == "Quant" || "~{ID_or_Quant_or_Both}" == "Both") && -n "~{referenceAnnotation_full}" && -z "~{LRAA_min_mapping_quality_flag}" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/Quant/LRAA \
                                 --quant_only \
                                 ~{no_norm_flag} \
                                 --gtf ~{referenceAnnotation_full} \
                                 --EM --CPU 1


            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/Quant_noEM/LRAA.noEM \
                                 --quant_only \
                                 ~{no_norm_flag} \
                                 --gtf ~{referenceAnnotation_full} --CPU 1
        fi

        if [[ ("~{ID_or_Quant_or_Both}" == "Quant" || "~{ID_or_Quant_or_Both}" == "Both") && -n "~{referenceAnnotation_full}" && -n "~{LRAA_min_mapping_quality}" ]]; then
            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/Quant_noEM_minMapQ/LRAA.noEM.minMapQ \
                                 --quant_only \
                                 ~{no_norm_flag} \
                                 --gtf ~{referenceAnnotation_full} \
                                 ~{LRAA_min_mapping_quality_flag} --CPU 1

            /usr/local/src/LRAA/LRAA --genome ~{referenceGenome} \
                                 --bam ~{inputBAM} \
                                 --output_prefix ~{OutDir}/Quant_minMapQ/LRAA.minMapQ \
                                 --quant_only \
                                 ~{no_norm_flag} \
                                 --gtf ~{referenceAnnotation_full} \
                                 ~{LRAA_min_mapping_quality_flag} \
                                 --EM --CPU ~{numThreads} --CPU 1
        fi
    >>>

    output {
        File? lraaGTF = "~{OutDir}/ID/LRAA.gtf"
        File? lraaReducedGTF = "~{OutDir}/ID_reduced/LRAA_reduced.gtf"
        File? lraaCounts = "~{OutDir}/Quant/LRAA.quant.expr"
        File? lraaCounts_noEM = "~{OutDir}/Quant_noEM/LRAA.noEM.quant.expr"
        File? lraa_quant_tracking = "~{OutDir}/Quant/LRAA.quant.tracking"
        File? lraa_quant_tracking_noEM = "~{OutDir}/Quant_noEM/LRAA.noEM.quant.tracking"
        File? lraaCounts_noEM_minMapQ = "~{OutDir}/Quant_noEM_minMapQ/LRAA.noEM.minMapQ.quant.expr"
        File? lraa_quant_tracking_noEM_minMapQ = "~{OutDir}/Quant_noEM_minMapQ/LRAA.noEM.minMapQ.quant.tracking"
        File? lraaCounts_minMapQ = "~{OutDir}/Quant_minMapQ/LRAA.minMapQ.quant.expr"
        File? lraa_quant_tracking_minMapQ = "~{OutDir}/Quant_minMapQ/LRAA.minMapQ.quant.tracking"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: "~{cpu}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        errorStrategy: "Continue"
        bootDiskSizeGb: 30
    }
}

workflow lraaWorkflow {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceGenome
        File referenceGenomeIndex
        File? referenceAnnotation_reduced
        File? referenceAnnotation_full
        String dataType
        String ID_or_Quant_or_Both
        Int? LRAA_min_mapping_quality
        Boolean? LRAA_no_norm
    }

    call lraaTask {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            referenceGenome = referenceGenome,
            referenceGenomeIndex = referenceGenomeIndex,
            referenceAnnotation_reduced = referenceAnnotation_reduced,
            referenceAnnotation_full = referenceAnnotation_full,
            dataType = dataType,
            ID_or_Quant_or_Both = ID_or_Quant_or_Both,
            LRAA_min_mapping_quality = LRAA_min_mapping_quality,
            LRAA_no_norm = LRAA_no_norm
    }

    output {
        File? lraaGTF = lraaTask.lraaGTF
        File? lraaReducedGTF = lraaTask.lraaReducedGTF
        File? lraaCounts = lraaTask.lraaCounts
        File? lraaCounts_noEM = lraaTask.lraaCounts_noEM
        File? lraa_quant_tracking = lraaTask.lraa_quant_tracking
        File? lraa_quant_tracking_noEM = lraaTask.lraa_quant_tracking_noEM
        File? lraaCounts_noEM_minMapQ = lraaTask.lraaCounts_noEM_minMapQ
        File? lraa_quant_tracking_noEM_minMapQ = lraaTask.lraa_quant_tracking_noEM_minMapQ
        File? lraaCounts_minMapQ = lraaTask.lraaCounts_minMapQ
        File? lraa_quant_tracking_minMapQ = lraaTask.lraa_quant_tracking_minMapQ
        File monitoringLog = lraaTask.monitoringLog
    }
}
