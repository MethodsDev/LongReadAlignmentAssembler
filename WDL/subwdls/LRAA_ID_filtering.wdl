version 1.0

task FilterTranscripts {
    input {
        File referenceGenome
        File gtf_path
        File expr_file_path
        Float threshold
        Int memoryGB
        Int diskSizeGB
        String docker
        File monitoringScript = "gs://mdl-ctat-genome-libs/terra_scripts/cromwell_monitoring_script2.sh"
    }

    command <<<
        bash ~{monitoringScript} > monitoring.log &
        python /usr/local/bin/FilterTranscripts ~{gtf_path} ~{expr_file_path} ~{referenceGenome} "LRAA_filtered.gtf" ~{threshold}
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 30
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }

    output {
        File filtered_gtf = "LRAA_filtered.gtf"
    }
}

workflow TranscriptFiltering {
    input {
        File referenceGenome
        File gtf_path
        File expr_file_path
        Float threshold = 1.0
        Int memoryGB = 32
        Int diskSizeGB = 128
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/iso-reconstruct-benchmark/filtertranscripts:latest"
    }

    call FilterTranscripts {
        input:
            referenceGenome = referenceGenome,
            gtf_path = gtf_path,
            expr_file_path = expr_file_path,
            threshold = threshold,
            docker = docker,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }

    output {
        File filtered_gtf = FilterTranscripts.filtered_gtf
    }
}
