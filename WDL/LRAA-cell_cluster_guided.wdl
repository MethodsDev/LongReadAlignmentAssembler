version 1.0

import "LRAA.wdl" as LRAA



workflow LRAA_cell_cluster_guided {

    input {
        String sample_id
        File referenceGenome

        File inputBAM
        File cell_clusters_info
        
        File? annot_gtf
        
        Boolean HiFi = false
        String? oversimplify # comma-separated contig names to simplify (e.g., "chrM" or "chrM,MT")

        String main_chromosomes = "" # ex. "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
        
        Int numThreadsPerWorker = 2
        Int numThreadsPerWorkerScattered = 9
        Int num_parallel_contigs = 3
        Int numThreadsPerLRAA = 5
        Int memoryGB = 64
        Int memoryGBPerWorkerScattered = 32
        Int memoryGBmergeGTFs = 32
        Int memoryGBquantFinal = 32
        Int memoryGBscSparseMatrices = 16
        Int diskSizeGB = 256
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        Boolean quant_only_cluster_guided = false

     }

      
     call partition_bam_by_cell_cluster {
         input:
            sample_id = sample_id,
            cell_clusters_info = cell_clusters_info,
            inputBAM = inputBAM,
            docker = docker,
     }

    # Package the partitioned cluster BAMs for convenient retrieval
    call LRAA_tar_outputs as tar_partitioned_cluster_bams {
        input:
            tar_directory_name = sample_id + ".LRAA.cluster_partitioned_bams",
            input_files = partition_bam_by_cell_cluster.partitioned_bams,
            docker = docker
    }

    # If running in quant-only mode, require annot_gtf to be provided
    if (quant_only_cluster_guided) {
        call require_annot_gtf {
            input:
                annot_gtf = select_first([annot_gtf])
        }
    }


    if (quant_only_cluster_guided == false) {

        scatter (i in range(length(partition_bam_by_cell_cluster.partitioned_bams))) {
            String cluster_sample_id = sub(basename(partition_bam_by_cell_cluster.partitioned_bams[i]), ".bam$", "")
            call LRAA.LRAA_wf as LRAA_by_cluster {
                input:
                    sample_id = cluster_sample_id,
                    referenceGenome = referenceGenome,
                    annot_gtf = annot_gtf,
                    inputBAM = partition_bam_by_cell_cluster.partitioned_bams[i],
                    HiFi = HiFi,
                    oversimplify = oversimplify,
                    main_chromosomes = main_chromosomes,
                    quant_only = false,
                    numThreadsPerWorker = numThreadsPerWorker,
                    numThreadsPerWorkerScattered = numThreadsPerWorkerScattered,
                    num_parallel_contigs = num_parallel_contigs,
                    memoryGB = memoryGB,
                    memoryGBPerWorkerScattered = memoryGBPerWorkerScattered,
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
                docker = docker,
                memoryGB = memoryGBmergeGTFs ,
        }
    }

     
    # Final quantification (single call): prefer merged GTF if produced, else use provided annot_gtf
    call LRAA_quant_bam_list as LRAA_quant_final_bamlist {
        input:
            sample_id = sample_id,
            referenceGenome = referenceGenome,
            annot_gtf = select_first([lraa_merge_gtf_task.mergedGTF, annot_gtf]),
            bam_files = partition_bam_by_cell_cluster.partitioned_bams,
            HiFi = HiFi,
            oversimplify = oversimplify,
            numThreads = numThreadsPerLRAA,
            memoryGB = memoryGBquantFinal,
            docker = docker
    }

    # Tar the per-cluster final quant outputs for convenient retrieval
    call LRAA_tar_outputs as tar_final_cluster_expr_files {
        input:
            tar_directory_name = sample_id + ".LRAA.final.cluster_quant.EXPRs",
            input_files = LRAA_quant_final_bamlist.quant_exprs,
            docker = docker
    }

    call LRAA_tar_outputs as tar_final_cluster_tracking_files {
        input:
            tar_directory_name = sample_id + ".LRAA.final.cluster_quant.trackings",
            input_files = LRAA_quant_final_bamlist.quant_trackings,
            docker = docker
    }

    # Build cluster-level pseudobulk matrices from the per-cluster quant expr files
    call build_cluster_pseudobulk_matrices as build_cluster_pseudobulk {
        input:
            sample_id = sample_id,
            quant_expr_files = LRAA_quant_final_bamlist.quant_exprs,
            docker = docker
    }

     # Merge all per-cluster tracking files into a single tracking.gz for downstream single-cell matrices
    call LRAA_merge_trackings as merge_cluster_trackings {
        input:
            sample_id = sample_id,
            tracking_files = LRAA_quant_final_bamlist.quant_trackings,
            docker = docker
    }

     # Build sparse matrices (Seurat-compatible) from the merged tracking
     call sc_build_sparse_matrices as build_sc_sparse_matrices {
             input:
                 sample_id = sample_id,
                 tracking_file = merge_cluster_trackings.merged_tracking,
                 docker = docker,
                 memoryGB = memoryGBscSparseMatrices
     }

     output {
         # final outputs
         File? LRAA_final_gtf = lraa_merge_gtf_task.mergedGTF
         File? LRAA_final_gtf_tracking = lraa_merge_gtf_task.mergedTracking
         # partitioned cluster BAMs (always produced)
         File LRAA_partitioned_cluster_bams_tar = tar_partitioned_cluster_bams.tar_gz
         # cluster-level final quant outputs (per-cluster/partition) packaged
         File LRAA_final_cluster_exprs_tar = tar_final_cluster_expr_files.tar_gz
         File LRAA_final_cluster_trackings_tar = tar_final_cluster_tracking_files.tar_gz
         # merged tracking across clusters (used for downstream SC matrices)
         File LRAA_final_tracking = merge_cluster_trackings.merged_tracking

         # single-cell sparse matrices and intermediate counts
         File sc_gene_transcript_splicehash_mapping = build_sc_sparse_matrices.mapping_file
         File sc_gene_counts = build_sc_sparse_matrices.gene_counts
         File sc_isoform_counts = build_sc_sparse_matrices.isoform_counts
         File sc_splice_pattern_counts = build_sc_sparse_matrices.splice_pattern_counts
         # gene/isoform/splice-pattern sparse directories as tar.gz
         File sc_gene_sparse_tar_gz = build_sc_sparse_matrices.gene_sparse_dir_tgz
         File sc_isoform_sparse_tar_gz = build_sc_sparse_matrices.isoform_sparse_dir_tgz
         File sc_splice_pattern_sparse_tar_gz = build_sc_sparse_matrices.splice_pattern_sparse_dir_tgz

         # preliminary intermediate outputs (only in discovery mode)
         File? LRAA_prelim_cluster_gtfs = tar_cluster_prelim_gtf_files.tar_gz
         File? LRAA_prelim_cluster_read_trackings = tar_cluster_prelim_tracking_files.tar_gz
         File? LRAA_prelim_cluster_pseudobulk_exprs = tar_cluster_prelim_pseudobulk_expr_files.tar_gz

         # cluster pseudobulk matrices
         File cluster_gene_counts_matrix = build_cluster_pseudobulk.cluster_gene_counts_matrix
         File cluster_gene_TPM_matrix = build_cluster_pseudobulk.cluster_gene_TPM_matrix
         File cluster_isoform_counts_matrix = build_cluster_pseudobulk.cluster_isoform_counts_matrix
         File cluster_isoform_TPM_matrix = build_cluster_pseudobulk.cluster_isoform_TPM_matrix
         File cluster_isoform_counts_forDiffIsoUsage = build_cluster_pseudobulk.cluster_isoform_counts_forDiffIsoUsage
    
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
        File? mergedTracking = "~{sample_id}.LRAA.sc_merged.gtf.tracking.tsv"
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


task LRAA_quant_bam_list {
    input {
        String sample_id
        File referenceGenome
        File annot_gtf
        Array[File] bam_files
    Boolean HiFi = false
        String? oversimplify
        String cell_barcode_tag = "CB"
        String read_umi_tag = "XM"
        Int numThreads = 4
        Int memoryGB = 32
        String docker
    }

    Int disksize = 50 + ceil(2 * size(bam_files, "GB"))

    String output_prefix = "~{sample_id}.LRAA.quant-only.clusters"

    command <<<
        set -ex

        # build bam list file
        : > bam_list.txt
        for f in ~{sep=' ' bam_files}; do
            echo "$f" >> bam_list.txt
        done

        (
          LRAA --genome ~{referenceGenome} \
               --bam_list bam_list.txt \
               --gtf ~{annot_gtf} \
               --quant_only \
               --CPU ~{numThreads} \
         ~{"--oversimplify " + oversimplify} \
               ~{true="--HiFi" false='' HiFi} \
               --cell_barcode_tag ~{cell_barcode_tag} --read_umi_tag ~{read_umi_tag} \
               --output_prefix ~{output_prefix} > command_output.log 2>&1
        ) || {
             echo "Command failed with exit code $?" >&2
             echo "Last 100 lines of output:" >&2
             tail -n 100 command_output.log >&2
             exit 1
        }

        # gzip any per-cluster tracking files produced (handle zero or many)
        shopt -s nullglob
        for f in ~{output_prefix}*.quant.tracking; do
            gzip "$f"
        done
        shopt -u nullglob
    >>>

    output {
        # Capture either single merged outputs (~{output_prefix}.quant.*) or per-cluster/partition outputs (~{output_prefix}.*.quant.*)
        Array[File] quant_exprs = glob("~{output_prefix}*.quant.expr")
        Array[File] quant_trackings = glob("~{output_prefix}*.quant.tracking*")
    }

    runtime {
        docker: docker
        cpu: "~{numThreads}"
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{disksize} HDD"
    }
}


task build_cluster_pseudobulk_matrices {
    input {
        String sample_id
        Array[File] quant_expr_files
        String docker
        Int memoryGB = 8
    }

    Int disksize = 20 + ceil(10 * length(quant_expr_files))

    String output_prefix = "~{sample_id}.LRAA.cluster_pseudobulk"

    command <<<
        set -ex

        build_LRAA_expr_matrices.py \
            --output_prefix ~{output_prefix} \
            --quant_files ~{sep=' ' quant_expr_files}
    >>>

    output {
        File cluster_gene_counts_matrix = "~{output_prefix}.gene.counts.matrix"
        File cluster_gene_TPM_matrix = "~{output_prefix}.gene.TPM.matrix"
        File cluster_isoform_counts_matrix = "~{output_prefix}.isoform.counts.matrix"
        File cluster_isoform_TPM_matrix = "~{output_prefix}.isoform.TPM.matrix"
        File cluster_isoform_counts_forDiffIsoUsage = "~{output_prefix}.isoform.counts.matrix.forDiffIsoUsage"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{disksize} HDD"
    }
}


task sc_build_sparse_matrices {
    input {
        String sample_id
        File tracking_file
        String docker
        Int memoryGB = 16
    }

    Int disksize = 50 + ceil(2 * size(tracking_file, "GB"))

    String output_prefix = "~{sample_id}.LRAA.sc"

    command <<<
        set -ex

        singlecell_tracking_to_sparse_matrix.py \
            --tracking ~{tracking_file} \
            --output_prefix ~{output_prefix}

        # Tar the generated sparse matrix directories for compact output
        tar -zcvf "~{output_prefix}^gene-sparseM.tar.gz" "~{output_prefix}^gene-sparseM" || true
        tar -zcvf "~{output_prefix}^isoform-sparseM.tar.gz" "~{output_prefix}^isoform-sparseM" || true
        tar -zcvf "~{output_prefix}^splice_pattern-sparseM.tar.gz" "~{output_prefix}^splice_pattern-sparseM" || true
    >>>

    output {
        File mapping_file = "~{output_prefix}.gene_transcript_splicehashcode.tsv"
        File gene_counts = "~{output_prefix}.gene_cell_counts.tsv"
        File isoform_counts = "~{output_prefix}.isoform_cell_counts.tsv"
        File splice_pattern_counts = "~{output_prefix}.splice_pattern_cell_counts.tsv"

        # tar.gz of sparse matrix directories
        File gene_sparse_dir_tgz = "~{output_prefix}^gene-sparseM.tar.gz"
        File isoform_sparse_dir_tgz = "~{output_prefix}^isoform-sparseM.tar.gz"
        File splice_pattern_sparse_dir_tgz = "~{output_prefix}^splice_pattern-sparseM.tar.gz"
    }

    runtime {
        docker: docker
        cpu: 2
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{disksize} HDD"
    }
}


task require_annot_gtf {
    input {
        File annot_gtf
    }

    command <<<
        set -euo pipefail
        echo "Annotation GTF provided: ~{annot_gtf}"
    >>>

    output {
        File ok = annot_gtf
    }

    runtime {
        cpu: 1
        memory: "1 GiB"
        disks: "local-disk 10 HDD"
    }
}

     


