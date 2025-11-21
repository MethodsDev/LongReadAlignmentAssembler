version 1.0

# partition_bam_by_cell_cluster.wdl
# Reusable task for partitioning a BAM file by cell cluster assignments
# Extracts cluster-specific BAMs based on cell barcode clustering information

task partition_bam_by_cell_cluster {
    input {
        String sample_id
        File cell_clusters_info
        File inputBAM
        String docker
        Int cpu = 8
        Int memoryGB = 16
    }

    Int disksize = ceil(5 * size(inputBAM, "GB"))
    
    command <<<
        set -ex

        mkdir partitioned_bams
        cd partitioned_bams/

        (
            partition_bam_by_cell_cluster.py --bam ~{inputBAM} \
                                             --cell_clusters ~{cell_clusters_info} \
                                             --output_prefix ~{sample_id} \
                                             --threads ~{cpu} > command_output.log 2>&1
        ) || {
            echo "Command failed with exit code $?" >&2
            echo "Last 100 lines of output:" >&2
            tail -n 100 command_output.log >&2
            exit 1
        }
      
        ls -1 *.bam
    >>>
     
    output {
        Array[File] partitioned_bams = glob("partitioned_bams/*.bam") 
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{disksize} HDD"
    }
}
