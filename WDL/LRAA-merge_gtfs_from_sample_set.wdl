version 1.0

import "LRAA-cell_cluster_guided.wdl" as LRAAClusterGuided

workflow LRAA_merge_gtfs_from_sample_set {

    input {
        String sample_set_id
        Array[File] sample_gtfs
        File referenceGenome

        Boolean HiFi = false
        Int memoryGB = 32
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:latest"
    }

    call LRAAClusterGuided.lraa_merge_gtf_task {
        input:
            sample_id = sample_set_id,
            LRAA_cell_cluster_gtfs = sample_gtfs,
            referenceGenome = referenceGenome,
            HiFi = HiFi,
            docker = docker,
            memoryGB = memoryGB
    }

    output {
        File mergedGTF = lraa_merge_gtf_task.mergedGTF
        File? mergedTracking = lraa_merge_gtf_task.mergedTracking
    }
}
