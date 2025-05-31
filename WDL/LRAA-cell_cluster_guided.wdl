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
     call lraa_merge_gtf_task {
         input:
            sample_id = sample_id,
            LRAA_cell_cluster_gtfs = LRAA_by_cluster.mergedGTF,
            referenceGenome = referenceGenome,
            docker=docker,
            memoryGB  = memoryGBperLRAA,
     }
     
     # run final quants
     scatter (i in range(length(partition_bam_by_cell_cluster.partitioned_bams))) {

         String cluster_sample_id_again = sub(basename(partition_bam_by_cell_cluster.partitioned_bams[i]), ".bam$", "")
              
         call LRAA.LRAA_wf as LRAA_quant_final {
             input:
               sample_id = cluster_sample_id_again, 
               referenceGenome = referenceGenome,
               annot_gtf = lraa_merge_gtf_task.mergedGTF,
               inputBAM = partition_bam_by_cell_cluster.partitioned_bams[i],
               LowFi = LowFi,
               main_chromosomes = main_chromosomes,
               quant_only = true,
               numThreads = numThreadsPerLRAA,
               memoryGB = memoryGBperLRAA,
               docker = docker
         }
     }


     call LRAA_tar_exprs {
         input:
             sample_id = sample_id,
             expr_files = LRAA_quant_final.mergedQuantExpr,
             docker = docker
     }

     call LRAA_merge_trackings {
         input:
             sample_id = sample_id,
             tracking_files = LRAA_quant_final.mergedQuantTracking,
             docker = docker
     }
     
     output {
         File LRAA_gtf = lraa_merge_gtf_task.mergedGTF
         File LRAA_cluster_pseudobulk_exprs = LRAA_tar_exprs.LRAA_cluster_pseudobulk_exprs_tar_gz
         File LRAA_merged_tracking = LRAA_merge_trackings.merged_tracking
     }
     
}


task LRAA_tar_exprs {
    input {
        String sample_id
        Array[File] expr_files
        String docker
    }

    Int memoryGB = 8
    Int disksize = 20 + ceil(10 * length(expr_files))
    
    command <<<

        set -ex

        mkdir ~{sample_id}.cluster_pseudobulk.EXPRs

        for file in "~{sep=' ' expr_files}"; do
           mv $file ~{sample_id}.cluster_pseudobulk.EXPRs/
        done
        
        tar -zcvf ~{sample_id}.cluster_pseudobulk.EXPRs.tar.gz ~{sample_id}.cluster_pseudobulk.EXPRs/

    >>>


    output {
        File LRAA_cluster_pseudobulk_exprs_tar_gz = "~{sample_id}.cluster_pseudobulk.EXPRs.tar.gz"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{disksize} HDD"
    }
}


task LRAA_merge_trackings {
    input {
        String sample_id
        Array[File] tracking_files
        String docker
    }

    Int memoryGB = 8
    Int disksize = 20 + ceil(10 * length(tracking_files))
        
    String outputfile = "~{sample_id}.cluster_merged.quant.tracking"
    
    command <<<
        set -ex

        python <<CODE
        import json


        tracking_files_json = '["' + '~{sep='","' tracking_files}' + '"]'
        tracking_files_list = json.loads(tracking_files_json)    # Parse the JSON string into a Python list

        with open("~{outputfile}", "wt") as ofh:
            for i, tracking_file in enumerate(tracking_files_list):
                with open(tracking_file, "rt") as fh:
                    header = next(fh)
                    if i == 0:
                         print(header, file=ofh, end='')
                    for line in fh:
                         print(line, file=ofh, end='')


        CODE
        
     >>>

     output {
         File merged_tracking = "~{outputfile}"
     }
    
     runtime {
        docker: docker
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{disksize} HDD"
     }
} 




task lraa_merge_gtf_task {
    input {
        String sample_id
        Array[File?] LRAA_cell_cluster_gtfs
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

     


