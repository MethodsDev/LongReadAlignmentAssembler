version 1.0

# partition_bam_by_cell_cluster.wdl
# Reusable task for partitioning a BAM file by cell cluster assignments
# Extracts cluster-specific BAMs based on cell barcode clustering information

task partition_bam_by_cell_cluster {
    input {
        Int preemptible_tries = 3
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
    # Machine shape (v0.47.1, generated text, identical in every small or HDD task). cpu: and
    # memory: stay as they are, so miniwdl's reservations do not change; this only names the
    # Terra VM (miniwdl ignores predefinedMachineType). N2D rather than C3D here because the
    # task uses at most 2 cores or an HDD disk: C3D has no 2-core shape, and Compute Engine
    # lists no pd-standard for C3D. Disks stay as they were until a probe settles that.
    Int shape_cpu = cpu
    Int shape_memGB = effective_memoryGB
    Int shape_cpu_tier = if shape_cpu <= 2 then 2
        else if shape_cpu <= 4 then 4
        else if shape_cpu <= 8 then 8
        else if shape_cpu <= 16 then 16
        else if shape_cpu <= 32 then 32
        else if shape_cpu <= 48 then 48
        else if shape_cpu <= 64 then 64
        else if shape_cpu <= 80 then 80
        else 96
    Int shape_mem_tier = if shape_memGB <= 16 then 2
        else if shape_memGB <= 32 then 4
        else if shape_memGB <= 64 then 8
        else if shape_memGB <= 128 then 16
        else if shape_memGB <= 256 then 32
        else if shape_memGB <= 384 then 48
        else if shape_memGB <= 512 then 64
        else if shape_memGB <= 640 then 80
        else 96
    Int shape_cores = if shape_cpu_tier >= shape_mem_tier then shape_cpu_tier else shape_mem_tier
    String machine_type = if shape_memGB <= shape_cores then "n2d-highcpu-~{shape_cores}"
        else if shape_memGB <= shape_cores * 4 then "n2d-standard-~{shape_cores}"
        else "n2d-highmem-~{shape_cores}"

    
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

        # The Seurat step writes its cluster assignments WITH a header row, and
        # partition_bam_by_cell_cluster.py's format is headerless: the row has two
        # tab-separated fields, so the parser accepts it and invents a cluster named
        # "cluster" holding one barcode named "cell_barcode". No read carries that
        # barcode, so the phantom cluster's bam is header-only, and the per-cluster
        # LRAA run it reaches aborts on count_reads_from_bam's "no reads counted"
        # assertion. The parser now drops the row itself; this strips it here as well
        # so a run against an image built before that fix still completes.
        #
        # Matched on the EXACT header rather than dropped positionally: this input is
        # documented as headerless and the fixed-cluster fixtures have no header, so an
        # unconditional `tail -n +2` would silently discard a real cell's assignment.
        CLUSTERS="~{cell_clusters_info}"
        if [ "$(head -n 1 "$CLUSTERS")" = "$(printf 'cell_barcode\tcluster')" ]; then
            tail -n +2 "$CLUSTERS" > clusters.noheader.tsv
            CLUSTERS="$(pwd)/clusters.noheader.tsv"
            echo "stripped cell_barcode/cluster header from cell_clusters_info" >&2
        fi

        mkdir partitioned_bams
        cd partitioned_bams/

        (
            partition_bam_by_cell_cluster.py --bam "$BAM" \
                                             --cell_clusters "$CLUSTERS" \
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

        predefinedMachineType: machine_type

        preemptible: preemptible_tries
        docker: docker
        cpu: cpu
        memory: "~{effective_memoryGB} GiB"
        disks: "local-disk ~{disksize} HDD"
    }
}
