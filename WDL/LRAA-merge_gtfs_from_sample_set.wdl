version 1.0

import "LRAA-cell_cluster_guided.wdl" as LRAAClusterGuided

workflow LRAA_merge_gtfs_from_sample_set {

    input {
        String sample_set_id
        Array[File] sample_gtfs
        File referenceGenome

        # Default: respect TSS/PolyA annotations on the input gtfs. Set true for the
        # former non-HiFi merge that ignores them.
        Boolean ignore_TSS_POLYA = false
        # LRAA_Globals.config overrides (JSON object) for merge_LRAA_GTFs.py --config_update.
        # Unset by default. '{"allow_unflagged_LRAA_guide": true}' merges sample GTFs
        # written before v0.50.0, which lack ref_splice_pattern / ref_supported.
        String? config_update_json
        Int memoryGB = 64
        Int cpu = 8
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:latest"
    }

    call LRAAClusterGuided.lraa_merge_gtf_task {
        input:
            sample_id = sample_set_id,
            LRAA_cell_cluster_gtfs = sample_gtfs,
            referenceGenome = referenceGenome,
            ignore_TSS_POLYA = ignore_TSS_POLYA,
            config_update_json = config_update_json,
            docker = docker,
            memoryGB = memoryGB,
            cpu = cpu
    }

    output {
        File mergedGTF = lraa_merge_gtf_task.mergedGTF
        File? mergedTracking = lraa_merge_gtf_task.mergedTracking
    }
}
