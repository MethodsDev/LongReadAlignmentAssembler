version 1.0

import "subwdls/Partition_data_by_chromosome.wdl" as PartByChr
import "subwdls/LRAA_runner.wdl" as LRAA_runner


workflow LRAA_wf {
     input {
        String sample_id
                 
        File referenceGenome 
        File inputBAM
        File? inputBAMindex
        File? annot_gtf
        Boolean LowFi = false
         
        String main_chromosomes = "" # ex. "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
        
        Float? min_per_id
        Boolean no_EM = false
        Boolean quant_only = false
        Boolean no_norm = false
        Int? min_mapping_quality

        String cell_barcode_tag = "CB"
        String read_umi_tag = "XM"

        Int numThreads = 4
        Int memoryGB = 32
        Int diskSizeGB = 128
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        
        
    }

    Boolean run_without_splitting = (main_chromosomes == "")
    
    if (!run_without_splitting) {

        call count_bam {
            input:
              bam = inputBAM, 
              bamindex = inputBAMindex
        }

        
        ## Split inputs by main chromosomes
        
        call PartByChr.partition_by_chromosome as splitByChr {
            input:
                inputBAM = inputBAM,
                inputBAMindex = inputBAMindex,
                genome_fasta = referenceGenome,
                annot_gtf = annot_gtf,
                chromosomes_want_partitioned = main_chromosomes,
                docker = docker,
            }
     
                  
        scatter (i in range(length(splitByChr.chromosomeBAMs))) {
            call LRAA_runner.LRAA_runner as LRAA_scatter {
                input:
                    sample_id = sample_id,
                    shardno = i,
                    inputBAM = splitByChr.chromosomeBAMs[i],
                    genome_fasta = splitByChr.chromosomeFASTAs[i],
                    annot_gtf = splitByChr.chromosomeGTFs[i],
                    num_total_reads = count_bam.count,
                    min_per_id = min_per_id,
                    quant_only = quant_only,
                    LowFi = LowFi,
                    no_norm = no_norm,
                    no_EM = no_EM,
                    cell_barcode_tag = cell_barcode_tag,
                    read_umi_tag = read_umi_tag,
                    numThreads = numThreads,
                    min_mapping_quality = min_mapping_quality,
                    docker = docker,
                    memoryGB = memoryGB,
                    diskSizeGB = diskSizeGB
            }
        }

        call mergeResults {
            input:
                quantExprFiles = LRAA_scatter.LRAA_quant_expr,
                quantTrackingFiles = LRAA_scatter.LRAA_quant_tracking,
                gtfFiles = LRAA_scatter.LRAA_gtf,
                outputFilePrefix = sample_id + ".LRAA",
                docker = docker,
        }
    }

    if (run_without_splitting) {
            
        call LRAA_runner.LRAA_runner as LRAA_direct {
            input:
                sample_id = sample_id,
                inputBAM = inputBAM,
                genome_fasta = referenceGenome,
                annot_gtf = annot_gtf,
                min_per_id = min_per_id,
                quant_only = quant_only,
                LowFi = LowFi,
                no_norm = no_norm,
                no_EM = no_EM,
                cell_barcode_tag = cell_barcode_tag,
                read_umi_tag = read_umi_tag,
                numThreads = numThreads,
                min_mapping_quality = min_mapping_quality,
                docker = docker,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB
        }
        
    }
    
    output {
        File mergedQuantExpr = select_first([mergeResults.mergedQuantExprFile, LRAA_direct.LRAA_quant_expr, ""]) 
        File mergedQuantTracking = select_first([mergeResults.mergedQuantTrackingFile, LRAA_direct.LRAA_quant_tracking, ""])
        File? mergedGTF = if (quant_only) then annot_gtf else select_first([mergeResults.mergedGtfFile, LRAA_direct.LRAA_gtf, ""])
    }
}



task mergeResults {
    input {
        Array[File] quantExprFiles
        Array[File] quantTrackingFiles
        Array[File] gtfFiles
        String outputFilePrefix
        String docker
    }

    command <<<
        set -eo pipefail

        quant_expr_output="~{outputFilePrefix}.quant.expr"
        for file in ~{sep=' ' quantExprFiles}; do
            if [[ ! -f "$quant_expr_output" ]]; then
                cp "$file" "$quant_expr_output"
            else
                tail -n +2 "$file" >> "$quant_expr_output"
            fi
        done

        #quant_tracking_output="~{outputFilePrefix}.quant.tracking"
        #for file in ~{sep=' ' quantTrackingFiles}; do
        #    if [[ ! -f "$quant_tracking_output" ]]; then
        #        cp "$file" "$quant_tracking_output"
        #    else
        #        tail -n +2 "$file" >> "$quant_tracking_output"
        #    fi
        #done

        gtf_output="~{outputFilePrefix}.gtf"
        touch "$gtf_output"
        gtf_files_str="~{sep=' ' gtfFiles}"
        for file in $gtf_files_str; do
           cat "$file" >> "$gtf_output"
        done


        python <<CODE
        import json
        import gzip

        tracking_files_json = '["' + '~{sep='","' quantTrackingFiles}' + '"]'
        tracking_files_list = json.loads(tracking_files_json)    # Parse the JSON string into a Python list

        with gzip.open("~{outputFilePrefix}.quant.tracking.gz", "wt") as ofh:
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
        File mergedQuantExprFile = "~{outputFilePrefix}.quant.expr"
        File mergedQuantTrackingFile = "~{outputFilePrefix}.quant.tracking.gz"
        File mergedGtfFile = "~{outputFilePrefix}.gtf"
    }
    
    runtime {
        docker: docker
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + ceil((size(quantExprFiles, "GB") + size(quantTrackingFiles, "GB") + size(gtfFiles, "GB")) * 2.2 + 5) + " SSD"
    }
}

task count_bam {
    input {
        File bam
        File? bamindex
    }

    # mapped only
    # '" samtools idxstats ~{bamindex} | awk -F'\t' '{m+=$3; u+=$4} END{print m}'"'
    # mapped + unmapped
    # '" samtools idxstats ~{bamindex} | awk -F'\t' '{m+=$3; u+=$4} END{print m+u}'"'

    String cmd_to_run = if defined(bamindex) then
      "samtools idxstats ~{bam} | awk -F'\\t' '{m+=$3; u+=$4} END{print m+u}'"
    else
      "samtools view -c ~{bam}"



    command <<<
        set -ex
        # samtools view -c ~{bam}
        ~{cmd_to_run}

    >>>

    runtime {
        docker: "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        disks: "local-disk " + ceil(2 * size(bam, "GB") ) + " HDD"
        cpu: 1
        memory: "4G"
    }

    output {
        Int count = read_int(stdout())
    }
}
