version 1.0

import "subwdls/Partition_data_by_chromosome.wdl" as PartByChr
import "subwdls/LRAA_runner.wdl" as LRAA_runner


workflow LRAA_wf {
     input {
        String sample_id
                 
        File referenceGenome 
        File inputBAM
        File? bam_for_sg
        File? annot_gtf
        Boolean HiFi = false
         
        String main_chromosomes = "" # ex. "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
        String? region # example: "chr1:100000-200000"; when set, workflow will not split by chromosome and will pass --region to LRAA
        String? oversimplify # comma-separated contig names to run in oversimplify mode
        
        Float? min_per_id
        Boolean no_EM = false
        Boolean quant_only = false
        Boolean no_norm = false
        Boolean allow_secondary_alignments = true
        Boolean rescue_unassigned_reads_via_transcriptome_alignment = true
        Int min_mapping_quality = 0
        Int min_mapping_quality_for_final_quant = 0

        String cell_barcode_tag = "CB"
        String read_umi_tag = "XM"

        #  non-scattered runs
        Int numThreadsPerWorker = 5
        Int? memoryGB
        Int num_parallel_contigs = 3

        # scattered runs
        Int numThreadsPerWorkerScattered = 5
        Int? memoryGBPerWorkerScattered
        
        
        Int diskSizeGB = 256
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        Int countBamThreads = 16
        
        
    }

    # Dynamic memory defaults based on input BAM size.
    # Direct (non-scattered) run: 1.5× full BAM size, floor 64 GiB.
    # Scattered workers self-size at 25× their shard BAM (see LRAA_runner_task); an optional
    # memoryGBPerWorkerScattered override is passed through to the task when set.
    Float inputBAMsizeGiB = size(inputBAM, "GiB")
    Float mem_raw_direct = 1.5 * inputBAMsizeGiB
    Int computed_memoryGB = if mem_raw_direct > 64.0 then ceil(mem_raw_direct) else 64
    Int effective_memoryGB = select_first([memoryGB, computed_memoryGB])

    Boolean run_without_splitting = (main_chromosomes == "" || defined(region))
    
    if (!run_without_splitting) {

        call count_bam {
            input:
          bam = inputBAM,
          samtools_threads = countBamThreads
        }

        
        ## Split inputs by main chromosomes
        
        call PartByChr.partition_by_chromosome as splitByChr {
            input:
                inputBAM = inputBAM,
                bam_for_sg = bam_for_sg,
                genome_fasta = referenceGenome,
                annot_gtf = annot_gtf,
                chromosomes_want_partitioned = main_chromosomes,
                docker = docker,
            }
     
                  
        Int num_chromosomes = length(splitByChr.chromosomeBAMs)

        scatter (contig_index in range(num_chromosomes)) {
            String contig_name = basename(splitByChr.chromosomeBAMs[contig_index], ".bam")
            # Run LRAA separately per chromosome  
            call LRAA_runner.LRAA_runner as LRAA_scatter {
                input:
                    sample_id = sample_id,
                    shardno = contig_index,
                    inputBAM = splitByChr.chromosomeBAMs[contig_index],
                    bam_for_sg = if defined(bam_for_sg) then select_first([splitByChr.chromosomeBAMsForSG])[contig_index] else bam_for_sg,
                    genome_fasta = splitByChr.chromosomeFASTAs[contig_index],
                    annot_gtf = splitByChr.chromosomeGTFs[contig_index],
                    oversimplify = oversimplify,
                    contig = contig_name,
                    num_parallel_contigs = num_parallel_contigs,
                    num_total_reads = count_bam.count,
                    min_per_id = min_per_id,
                    quant_only = quant_only,
                    HiFi = HiFi,
                    no_norm = no_norm,
                    no_EM = no_EM,
                    run_final_cross_gene_EM = false,
                    allow_secondary_alignments = allow_secondary_alignments,
                    rescue_unassigned_reads_via_transcriptome_alignment = rescue_unassigned_reads_via_transcriptome_alignment,
                    cell_barcode_tag = cell_barcode_tag,
                    read_umi_tag = read_umi_tag,
                    numThreadsPerWorker = numThreadsPerWorkerScattered,
                    min_mapping_quality = min_mapping_quality,
                    min_mapping_quality_for_final_quant = min_mapping_quality_for_final_quant,
                    docker = docker,
                    memoryGB = memoryGBPerWorkerScattered,  # Int? — if unset, task self-sizes from shard BAM using size and mid-small-shard estimates (floor 32 GiB)
                    diskSizeGB = diskSizeGB
            }
        }

        # Always merge quant outputs regardless of quant_only
        call mergeQuantResults {
            input:
                quantExprFiles = LRAA_scatter.LRAA_quant_expr,
                quantTrackingFiles = LRAA_scatter.LRAA_quant_tracking,
                outputFilePrefix = sample_id + ".LRAA",
                runCrossGeneEM = allow_secondary_alignments && !no_EM,
                docker = docker,
        }

        call mergeGenomeTxArbSummaries {
            input:
                summaryFiles = select_all(LRAA_scatter.LRAA_genome_tx_arb_summary),
                outputFilePrefix = sample_id + ".LRAA",
                docker = docker
        }

        # Only merge GTFs when not in quant-only mode
        if (!quant_only) {
            call merge_GTFs {
                input:
                    gtfFiles = select_all(LRAA_scatter.LRAA_gtf),
                    outputFilePrefix = sample_id + ".LRAA",
                    docker = docker,
            }
        }
    }

    if (run_without_splitting) {
            
        call LRAA_runner.LRAA_runner as LRAA_direct {
            input:
                sample_id = sample_id,
                inputBAM = inputBAM,
                bam_for_sg = bam_for_sg,
                genome_fasta = referenceGenome,
                annot_gtf = annot_gtf,
                region = region,
                oversimplify = oversimplify,
                min_per_id = min_per_id,
                quant_only = quant_only,
                HiFi = HiFi,
                no_norm = no_norm,
                no_EM = no_EM,
                run_final_cross_gene_EM = true,
                allow_secondary_alignments = allow_secondary_alignments,
                rescue_unassigned_reads_via_transcriptome_alignment = rescue_unassigned_reads_via_transcriptome_alignment,
                cell_barcode_tag = cell_barcode_tag,
                read_umi_tag = read_umi_tag,
                numThreadsPerWorker = numThreadsPerWorker,
                num_parallel_contigs = num_parallel_contigs,
                min_mapping_quality = min_mapping_quality,
                min_mapping_quality_for_final_quant = min_mapping_quality_for_final_quant,
                docker = docker,
                memoryGB = effective_memoryGB,
                diskSizeGB = diskSizeGB
        }

    }
    
    output {
    File mergedQuantExpr = select_first([mergeQuantResults.mergedQuantExprFile, LRAA_direct.LRAA_quant_expr]) 
    File mergedQuantTracking = select_first([mergeQuantResults.mergedQuantTrackingFile, LRAA_direct.LRAA_quant_tracking])
    File? mergedGTF = if (!quant_only) then select_first([merge_GTFs.mergedGtfFile, LRAA_direct.LRAA_gtf]) else LRAA_direct.LRAA_gtf 
    Array[File] tiedSecondariesBams = if (run_without_splitting) then select_all([LRAA_direct.LRAA_tied_secondaries_bam]) else select_all(select_first([LRAA_scatter.LRAA_tied_secondaries_bam, []]))
    Array[File] tiedSecondariesBais = if (run_without_splitting) then select_all([LRAA_direct.LRAA_tied_secondaries_bai]) else select_all(select_first([LRAA_scatter.LRAA_tied_secondaries_bai, []]))
    Array[File] shardGenomeTxArbSummaries = if (run_without_splitting) then select_all([LRAA_direct.LRAA_genome_tx_arb_summary]) else select_all(select_first([LRAA_scatter.LRAA_genome_tx_arb_summary, []]))
    File? mergedGenomeTxArbSummary = if (run_without_splitting) then LRAA_direct.LRAA_genome_tx_arb_summary else mergeGenomeTxArbSummaries.mergedSummaryFile
    }
}



 

task merge_GTFs {
    input {
        Array[File] gtfFiles
        String outputFilePrefix
        String docker
    }

    command <<<
        set -eo pipefail

        gtf_output="~{outputFilePrefix}.gtf"
        touch "$gtf_output"
        gtf_files_str="~{sep=' ' gtfFiles}"
        for file in $gtf_files_str; do
           cat "$file" >> "$gtf_output"
        done
    >>>

    output {
        File mergedGtfFile = "~{outputFilePrefix}.gtf"
    }
    
    runtime {
        docker: docker
        cpu: 1
        memory: "2 GiB"
        disks: "local-disk " + ceil(size(gtfFiles, "GB") * 2.0 + 5) + " SSD"
    }
}

task mergeQuantResults {
    input {
        Array[File] quantExprFiles
        Array[File] quantTrackingFiles
        String outputFilePrefix
        Boolean runCrossGeneEM = false
        String docker
    }

    Float quantMergeInputGB = size(quantExprFiles, "GB") + size(quantTrackingFiles, "GB")
    Float quantMergeInputGiB = size(quantExprFiles, "GiB") + size(quantTrackingFiles, "GiB")
    Float mergeMemoryRawGiB = if runCrossGeneEM then quantMergeInputGiB * 30.0 + 8.0 else 4.0
    Int mergeMemoryGiB = if mergeMemoryRawGiB > 4.0 then ceil(mergeMemoryRawGiB) else 4
    Float mergeDiskRawGB = if runCrossGeneEM then quantMergeInputGB * 4.0 + 20.0 else quantMergeInputGB * 2.2 + 5.0

    command <<<
    set -eo pipefail

    merge_LRAA_quant_expr.py \
        --output "~{outputFilePrefix}.pre_cross_gene_em.quant.expr" \
        --quant_files ~{sep=' ' quantExprFiles}

    python <<CODE
    import json
    import gzip

    tracking_files_json = '["' + '~{sep='","' quantTrackingFiles}' + '"]'
    tracking_files_list = json.loads(tracking_files_json)

    with gzip.open("~{outputFilePrefix}.pre_cross_gene_em.quant.tracking.gz", "wt") as ofh:
        for i, tracking_file in enumerate(tracking_files_list):
            openf = gzip.open if tracking_file.split(".")[-1] == "gz" else open
            with openf(tracking_file, "rt") as fh:
                header = next(fh)
                if i == 0:
                        print(header, file=ofh, end='')
                for line in fh:
                        print(line, file=ofh, end='')
    CODE

    if [[ "~{runCrossGeneEM}" == "true" ]]; then
        reassign_multigene_tracking_reads.py \
            --quant_expr "~{outputFilePrefix}.pre_cross_gene_em.quant.expr" \
            --tracking "~{outputFilePrefix}.pre_cross_gene_em.quant.tracking.gz" \
            --output_expr "~{outputFilePrefix}.quant.expr" \
            --output_tracking "~{outputFilePrefix}.quant.tracking.gz"
    else
        cp "~{outputFilePrefix}.pre_cross_gene_em.quant.expr" "~{outputFilePrefix}.quant.expr"
        cp "~{outputFilePrefix}.pre_cross_gene_em.quant.tracking.gz" "~{outputFilePrefix}.quant.tracking.gz"
    fi
    >>>

    output {
        File mergedQuantExprFile = "~{outputFilePrefix}.quant.expr"
        File mergedQuantTrackingFile = "~{outputFilePrefix}.quant.tracking.gz"
    }
    
    runtime {
        docker: docker
        cpu: 1
        memory: mergeMemoryGiB + " GiB"
        disks: "local-disk " + ceil(mergeDiskRawGB) + " SSD"
    }
}

task mergeGenomeTxArbSummaries {
    input {
        Array[File] summaryFiles
        String outputFilePrefix
        String docker
    }

    command <<<
    set -eo pipefail

    python <<'CODE'
    import csv
    from pathlib import Path

    summary_files = [Path(p) for p in """~{sep='\n' summaryFiles}""".splitlines() if p.strip()]
    out_path = Path("~{outputFilePrefix}.genome_tx_arb.summary.tsv")

    fieldnames = None
    worker_rows = []
    total_keys = [
        "reads_total",
        "reads_kept_genome",
        "reads_selected_tx_total",
        "reads_selected_tx_missing_genome",
        "reads_selected_tx_failed_genome",
        "reads_selected_tx_higher_per_id",
        "reads_tx_present_but_kept_genome",
        "reads_tx_tied_per_id_kept_genome",
        "reads_tx_lower_per_id_kept_genome",
        "reads_tx_per_id_unavailable_kept_genome",
    ]
    totals = {key: 0 for key in total_keys}

    for summary_file in summary_files:
        with summary_file.open("rt", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if fieldnames is None:
                fieldnames = list(reader.fieldnames or [])
            for row in reader:
                if row.get("row_type") == "TOTAL":
                    continue
                worker_rows.append(row)
                for key in total_keys:
                    totals[key] += int(row.get(key, 0) or 0)

    if fieldnames is None:
        fieldnames = [
            "row_type",
            "contig_acc",
            "contig_strand",
            "reads_total",
            "reads_kept_genome",
            "frac_kept_genome",
            "reads_selected_tx_total",
            "frac_selected_tx_total",
            "reads_selected_tx_missing_genome",
            "frac_selected_tx_missing_genome",
            "reads_selected_tx_failed_genome",
            "frac_selected_tx_failed_genome",
            "reads_selected_tx_higher_per_id",
            "frac_selected_tx_higher_per_id",
            "reads_tx_present_but_kept_genome",
            "frac_tx_present_but_kept_genome",
            "reads_tx_tied_per_id_kept_genome",
            "frac_tx_tied_per_id_kept_genome",
            "reads_tx_lower_per_id_kept_genome",
            "frac_tx_lower_per_id_kept_genome",
            "reads_tx_per_id_unavailable_kept_genome",
            "frac_tx_per_id_unavailable_kept_genome",
        ]

    total_reads = totals["reads_total"]

    def frac(key):
        if total_reads <= 0:
            return "0.000000"
        return f"{float(totals[key]) / float(total_reads):.6f}"

    total_row = {field: "" for field in fieldnames}
    total_row["row_type"] = "TOTAL"
    total_row["contig_acc"] = "TOTAL"
    total_row["contig_strand"] = "."
    for key in total_keys:
        total_row[key] = str(totals[key])
    total_row["frac_kept_genome"] = frac("reads_kept_genome")
    total_row["frac_selected_tx_total"] = frac("reads_selected_tx_total")
    total_row["frac_selected_tx_missing_genome"] = frac("reads_selected_tx_missing_genome")
    total_row["frac_selected_tx_failed_genome"] = frac("reads_selected_tx_failed_genome")
    total_row["frac_selected_tx_higher_per_id"] = frac("reads_selected_tx_higher_per_id")
    total_row["frac_tx_present_but_kept_genome"] = frac("reads_tx_present_but_kept_genome")
    total_row["frac_tx_tied_per_id_kept_genome"] = frac("reads_tx_tied_per_id_kept_genome")
    total_row["frac_tx_lower_per_id_kept_genome"] = frac("reads_tx_lower_per_id_kept_genome")
    total_row["frac_tx_per_id_unavailable_kept_genome"] = frac("reads_tx_per_id_unavailable_kept_genome")

    with out_path.open("wt", newline="") as ofh:
        writer = csv.DictWriter(ofh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in worker_rows:
            writer.writerow(row)
        writer.writerow(total_row)
    CODE
    >>>

    output {
        File mergedSummaryFile = "~{outputFilePrefix}.genome_tx_arb.summary.tsv"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "2 GiB"
        disks: "local-disk " + ceil(size(summaryFiles, "GB") * 2.0 + 5) + " SSD"
    }
}

task count_bam {
  input {
    File bam
        Int samtools_threads = 16
  }

    Float bam_size_gb = size(bam, "GB")
    Float estimated_disk = ceil(bam_size_gb * 2.2 + 20.0)
    Float disk_gb = if estimated_disk > 150.0 then estimated_disk else 150.0
    Int disk_gb_int = ceil(disk_gb)

  command <<<
    set -ex
        samtools view -@ ~{samtools_threads} -c ~{bam}

  >>>

  runtime {
        docker: "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        disks: "local-disk " + disk_gb_int + " SSD"
        cpu: samtools_threads
        memory: "8G"
  }
  output {
    Int count = read_int(stdout())
  }
}
