version 1.0

# partition_bam_by_cell_cluster.wdl
# Reusable task for partitioning a BAM file by cell cluster assignments
# Extracts cluster-specific BAMs based on cell barcode clustering information

task partition_bam_by_cell_cluster {
    input {
        String sample_id
        File cell_clusters_info
        File inputBAM
        String cell_barcode_tag = "CB"
        # Restrict the EMITTED cluster bams to these references, space-separated as
        # every caller's main_chromosomes is. Empty means emit every reference.
        #
        # This is the ONLY place the run's contig restriction has to be applied on
        # the bam-preparation side: normalization, the cluster merge and the merged
        # normalization all consume THESE bams, so they inherit it. The LRAA calls
        # were already restricted by main_chromosomes on their own.
        #
        # It does NOT change any TPM denominator on the full input bam: that bam is
        # untouched here, so the initial pass still counts the whole library. Each
        # per-cluster LRAA call resolves its own denominator from the cluster bam it
        # is handed (LRAA:2420-2421), so restricting these bams does make those
        # denominators subset-relative -- per cluster, as they already were.
        String main_chromosomes = ""
        String docker
        Int cpu = 8
        Int? memoryGB
    }

    Int disksize = ceil(5 * size(inputBAM, "GB"))

    # Dynamic memory: 0.3× BAM size, floor 16 GiB
    Float bam_size_gib = size(inputBAM, "GiB")
    Float mem_raw_cluster = 0.3 * bam_size_gib
    Int computed_memoryGB = if mem_raw_cluster > 16.0 then ceil(mem_raw_cluster) else 16
    Int effective_memoryGB = select_first([memoryGB, computed_memoryGB])
    
    command <<<
        set -ex

        BAM="~{inputBAM}"

        # An index is what lets the partitioner FETCH the requested references
        # instead of decompressing the whole library to use a fraction of it. The
        # localized input carries no index of its own -- no caller declares one --
        # so build it here, and only when a restriction is actually in force.
        #
        # Through a symlink because the localized input is read-only and both
        # `samtools index` and faidx write beside their argument. Reuse a .bai that
        # happens to have been localized alongside rather than rebuilding it.
        if [ -n "~{main_chromosomes}" ]; then
            ln -s "~{inputBAM}" restrict_input.bam
            if [ -f "~{inputBAM}.bai" ]; then
                ln -s "~{inputBAM}.bai" restrict_input.bam.bai
            else
                samtools index -@ ~{cpu} restrict_input.bam
            fi
            BAM="$(pwd)/restrict_input.bam"
        fi

        mkdir partitioned_bams
        cd partitioned_bams/

        (
            partition_bam_by_cell_cluster.py --bam "$BAM" \
                                             --cell_clusters ~{cell_clusters_info} \
                                             --output_prefix ~{sample_id} \
                                             --cell_barcode_tag ~{cell_barcode_tag} \
                                             ~{if main_chromosomes != "" then "--restrict_to_chromosomes '" + main_chromosomes + "'" else ""} \
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
        # One bam per cluster, ALWAYS -- the script materializes an empty bam for
        # every cluster before dispatching any read, so a cluster with nothing on
        # the retained references still appears here, header-only. That is what
        # keeps this array the same length as the cell lists below when
        # main_chromosomes is set.
        Array[File] partitioned_bams = glob("partitioned_bams/*.bam")
        # Same cluster names and the same sort order as the BAMs above, so index i
        # of each array refers to the same cluster.
        Array[File] partitioned_cell_lists = glob("partitioned_bams/*.cells.txt")
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{effective_memoryGB} GiB"
        disks: "local-disk ~{disksize} HDD"
    }
}
