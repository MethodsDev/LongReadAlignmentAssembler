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
        Int memoryGBperLRAA = 16
        Int memoryGBmergeGTFs = 32
        Int memoryGBquantFinal = 32
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
    
    # package them up
    call LRAA_tar_outputs as tar_cluster_prelim_gtf_files {
          input:
             tar_directory_name = sample_id + ".LRAA.prelim.cluster_gtfs",
             input_files = select_all(LRAA_by_cluster.mergedGTF),
             docker = docker
    }

    call LRAA_tar_outputs as tar_cluster_prelim_tracking_files {
          input:
             tar_directory_name = sample_id + ".LRAA.prelim.cluster_read_trackings",
             input_files = LRAA_by_cluster.mergedQuantTracking,
             docker = docker
    }

        
    call LRAA_tar_outputs as tar_cluster_prelim_pseudobulk_expr_files {
          input:
             tar_directory_name = sample_id + ".LRAA.prelim.cluster_pseudobulk.EXPRs",
             input_files = LRAA_by_cluster.mergedQuantExpr,
             docker = docker
    }
        
    # merge gtfs from the per-cluster runs into a single final gtf.
    call lraa_merge_gtf_task {
          input:
            sample_id = sample_id,
            LRAA_cell_cluster_gtfs = select_all(LRAA_by_cluster.mergedGTF),
            referenceGenome = referenceGenome,
            docker=docker,
            memoryGB = memoryGBmergeGTFs ,
    
     }

     
     call LRAA.LRAA_wf as LRAA_quant_final {
             input:
               sample_id = sample_id, 
               referenceGenome = referenceGenome,
               annot_gtf = lraa_merge_gtf_task.mergedGTF,
               inputBAM = inputBAM,
               LowFi = LowFi,
               main_chromosomes = main_chromosomes,
               quant_only = true,
               numThreads = numThreadsPerLRAA,
               memoryGB = memoryGBquantFinal,
               docker = docker
     }

     output {
         # final outputs
         File LRAA_final_gtf = lraa_merge_gtf_task.mergedGTF
         File LRAA_final_expr = LRAA_quant_final.mergedQuantExpr
         File LRAA_final_tracking = LRAA_quant_final.mergedQuantTracking

         # preliminary intermediate outputs.
         File LRAA_prelim_cluster_gtfs = tar_cluster_prelim_gtf_files.tar_gz
         File LRAA_prelim_cluster_read_trackings = tar_cluster_prelim_tracking_files.tar_gz
         File LRAA_prelim_cluster_pseudobulk_exprs = tar_cluster_prelim_pseudobulk_expr_files.tar_gz
    
     }
     
}


task LRAA_tar_outputs {
    input {
        String tar_directory_name
        Array[File] input_files
        String docker
    }

    Int memoryGB = 8
    Int disksize = 20 + ceil(5 * size(input_files, "GiB"))
    
    command <<<

        set -ex

        mkdir ~{tar_directory_name}

        for file in "~{sep=' ' input_files}"; do
           cp $file ~{tar_directory_name}/
        done
        
        tar -zcvf ~{tar_directory_name}.tar.gz ~{tar_directory_name}/

    >>>


    output {
        File tar_gz = "~{tar_directory_name}.tar.gz"
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
        
    String outputfile = "~{sample_id}.cluster_merged.quant.tracking.gz"
    
    command <<<
        set -ex

        python <<CODE
        import json
        import gzip

        tracking_files_json = '["' + '~{sep='","' tracking_files}' + '"]'
        tracking_files_list = json.loads(tracking_files_json)    # Parse the JSON string into a Python list

        with gzip.open("~{outputfile}", "wt") as ofh:
            for i, tracking_file in enumerate(tracking_files_list):
                openf = gzip.open if tracking_file.split(".")[-1] == "gz" else open
                with openf(tracking_file, "rt") as fh:
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
        Array[File] LRAA_cell_cluster_gtfs
        File referenceGenome
        String docker
        Int memoryGB
    }

    
    command <<<
        set -ex


      (
      
        merge_LRAA_GTFs.py --genome ~{referenceGenome} \
                           --gtf ~{sep=' ' LRAA_cell_cluster_gtfs } \
                           --output_gtf ~{sample_id}.LRAA.sc_merged.gtf  > command_output.log 2>&1
      ) || {
             echo "Command failed with exit code $?" >&2
             echo "Last 100 lines of output:" >&2
             tail -n 100 command_output.log >&2
             exit 1
      }

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

        (
         partition_bam_by_cell_cluster.py --bam ~{inputBAM} \
                                          --cell_clusters ~{cell_clusters_info} \
                                          --output_prefix ~{sample_id} > command_output.log 2>&1
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
        cpu: 1
        memory: "8 GiB"
        disks: "local-disk ~{disksize} HDD"
    }

}

     


