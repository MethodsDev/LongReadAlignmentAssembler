version 1.0

import "LRAA.wdl" as LRAA



workflow LRAA_cell_cluster_guided {

    input {
        String sample_id
        File referenceGenome

        File inputBAM
        File cell_clusters_info
        
        File? annot_gtf
        
        Boolean LowFi = false

        String main_chromosomes = "" # ex. "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
        
        Int numThreadsPerLRAA = 4
        Int memoryGBperLRAA = 32
        Int diskSizeGB = 128
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"

     }

      
     call partition_bam_by_cell_cluster {
         input:
            sample_id = sample_id,
            cell_clusters_info = cell_clusters_info,
            inputBAM = inputBAM,
            docker = docker,
     }
      

     scatter (i in range(length(partition_bam_by_cell_cluster.partitioned_bams))) {

         String cluster_sample_id = sub(basename(partition_bam_by_cell_cluster.partitioned_bams[i]), ".bam$", "")
         
         
         call LRAA.LRAA_wf as LRAA_by_cluster {
             input:
               sample_id = cluster_sample_id,
               referenceGenome = referenceGenome,
               annot_gtf = annot_gtf,
               inputBAM = partition_bam_by_cell_cluster.partitioned_bams[i],
               LowFi = LowFi,
               main_chromosomes = main_chromosomes,
               quant_only = false,
               numThreads = numThreadsPerLRAA,
               memoryGB  = memoryGBperLRAA,
               docker = docker
         }

     }

     # merge gtfs from the per-cluster runs.
     call lraa_merge_task {
         input:
            sample_id = sample_id,
            LRAA_cell_cluster_gtfs = LRAA_by_cluster.mergedGTF,
            referenceGenome = referenceGenome,
            docker=docker,
            memoryGB  = memoryGBperLRAA,
     }
     
     # run final quant
     call LRAA.LRAA_wf as LRAA_quant {
         input:
            sample_id = sample_id, 
               referenceGenome = referenceGenome,
               annot_gtf = lraa_merge_task.mergedGTF,
               inputBAM = inputBAM,
               LowFi = LowFi,
               main_chromosomes = main_chromosomes,
               quant_only = true,
               numThreads = numThreadsPerLRAA,
               memoryGB = memoryGBperLRAA,
               docker = docker
     }


     output {
         File LRAA_gtf = lraa_merge_task.mergedGTF
         File LRAA_pseudobulk_expr = LRAA_quant.mergedQuantExpr
         File LRAA_tracking = LRAA_quant.mergedQuantTracking
     }
     
}

task lraa_merge_task {
    input {
        String sample_id
        Array[File] LRAA_cell_cluster_gtfs
        File referenceGenome
        String docker
        Int memoryGB
    }

    command <<<
        set -ex
        
        merge_LRAA_GTFs.py --genome ~{referenceGenome} \
                           --gtf ~{sep=' ' LRAA_cell_cluster_gtfs } \
                           --output_gtf ~{sample_id}.LRAA.sc_merged.gtf

    >>>


    output {
        File mergedGTF = "~{sample_id}.LRAA.sc_merged.gtf"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk 200 HDD"
    }

}


task partition_bam_by_cell_cluster {

    input {
        String sample_id
        File cell_clusters_info
        File inputBAM
        String docker
    }


    Int disksize = ceil(4 * size(inputBAM, "GB") )
    
    command <<<
         set -ex

         mkdir partitioned_bams
         cd partitioned_bams/

         partition_bam_by_cell_cluster.py --bam ~{inputBAM} \
                                          --cell_clusters ~{cell_clusters_info} \
                                          --output_prefix ~{sample_id}


        ls -1 *.bam
        
         
    >>>
     
    output {
         Array[File] partitioned_bams = glob("partitioned_bams/*.bam") 
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "8 GiB"
        disks: "local-disk ~{disksize} HDD"
    }

}

     


