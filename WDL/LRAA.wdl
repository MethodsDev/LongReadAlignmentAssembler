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
        # Return the depth-normalized BAM(s) the splice graph was built from. Single-cell
        # workflows that call this one set this false; they never surface the file, and
        # delocalizing it would cost them storage for nothing.
        Boolean retain_normalized_splice_graph_bam = true
        Boolean rescue_unassigned_reads_via_transcriptome_alignment = true
        Int min_mapping_quality = 0
        Int min_mapping_quality_for_final_quant = 0

        # If set, the count_bam step is skipped and this value is passed to LRAA as --num_total_reads
        Int? num_total_reads
        # Barcodes accepted as real cells. Without it the supporting-cell filter
        # trusts every barcode in the BAM, which is only safe when the BAM was already
        # restricted to called cells.
        File? cell_list

        String cell_barcode_tag = "CB"
        String read_umi_tag = "XM"

        #  non-scattered runs
        # Cores for the LRAA task: its cpu request AND the --cpu_budget it divides across
        # work units. numThreadsPerWorker and num_parallel_contigs multiplied instead:
        # 5 x 3 asked LRAA for up to fifteen cores on a five-core task.
        Int cpu = 5
        Int? memoryGB

        # scattered runs
        Int cpuScattered = 5
        Int? memoryGBPerWorkerScattered
        
        
        Int diskSizeGB = 256
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-core:latest"
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
    String LRAA_output_suffix = if !defined(annot_gtf) && !quant_only then "LRAA.ref-free" else if defined(annot_gtf) && !quant_only then "LRAA.ref-guided" else "LRAA.quant-only"
    String LRAA_output_prefix = sample_id + "." + LRAA_output_suffix
    
    if (!run_without_splitting) {

        if (!defined(num_total_reads)) {
            call count_bam {
                input:
                    bam = inputBAM,
                    samtools_threads = countBamThreads,
                    docker = docker
            }
        }

        Int scatter_num_total_reads = select_first([num_total_reads, count_bam.count])

        
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
                    num_total_reads = scatter_num_total_reads,
                    cell_list = cell_list,
                    min_per_id = min_per_id,
                    quant_only = quant_only,
                    HiFi = HiFi,
                    no_norm = no_norm,
                    no_EM = no_EM,
                    retain_normalized_splice_graph_bam = retain_normalized_splice_graph_bam,
                    rescue_unassigned_reads_via_transcriptome_alignment = rescue_unassigned_reads_via_transcriptome_alignment,
                    cell_barcode_tag = cell_barcode_tag,
                    read_umi_tag = read_umi_tag,
                    cpu = cpuScattered,
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
                outputFilePrefix = LRAA_output_prefix,
                docker = docker,
        }

        call mergeGenomeTxArbSummaries {
            input:
                summaryFiles = select_all(LRAA_scatter.LRAA_genome_tx_arb_summary),
                outputFilePrefix = LRAA_output_prefix,
                docker = docker
        }

        # Only merge GTFs when not in quant-only mode
        if (!quant_only) {
            call merge_GTFs {
                input:
                    gtfFiles = select_all(LRAA_scatter.LRAA_gtf),
                    outputFilePrefix = LRAA_output_prefix,
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
                retain_normalized_splice_graph_bam = retain_normalized_splice_graph_bam,
                rescue_unassigned_reads_via_transcriptome_alignment = rescue_unassigned_reads_via_transcriptome_alignment,
                cell_barcode_tag = cell_barcode_tag,
                read_umi_tag = read_umi_tag,
                cpu = cpu,
                num_total_reads = num_total_reads,
                cell_list = cell_list,
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
    Array[File] shardGenomeTxArbSummaries = if (run_without_splitting) then select_all([LRAA_direct.LRAA_genome_tx_arb_summary]) else select_all(select_first([LRAA_scatter.LRAA_genome_tx_arb_summary, []]))
    File? mergedGenomeTxArbSummary = if (run_without_splitting) then LRAA_direct.LRAA_genome_tx_arb_summary else mergeGenomeTxArbSummaries.mergedSummaryFile
    # The depth-normalized BAM(s) the splice graph -- and therefore isoform identification --
    # was built from. Quantification does not use these; it runs against the unnormalized quant
    # BAM. Empty when no_norm is set. Scattered runs normalize each chromosome shard separately,
    # so this is one BAM per shard rather than one whole-genome BAM.
    Array[File] normalizedSpliceGraphBams = if (run_without_splitting) then select_all([LRAA_direct.LRAA_normalized_splice_graph_bam]) else select_all(select_first([LRAA_scatter.LRAA_normalized_splice_graph_bam, []]))
    Array[File] normalizedSpliceGraphBais = if (run_without_splitting) then select_all([LRAA_direct.LRAA_normalized_splice_graph_bai]) else select_all(select_first([LRAA_scatter.LRAA_normalized_splice_graph_bai, []]))
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

        lraa_version="$(LRAA --version | awk '{print $NF}')"
        version_comment="# LRAA version ${lraa_version}"

        gtf_output="~{outputFilePrefix}.gtf"
        gtf_files_str="~{sep=' ' gtfFiles}"

        lraa_merge_header.py \
            --version_comment "$version_comment" \
            --inputs $gtf_files_str \
            --output "$gtf_output"

        # the merged header above stands in for the shard headers, so drop each
        # shard's whole leading comment block, not just its version line
        awk 'FNR == 1 { in_header = 1 } in_header && /^#/ { next } { in_header = 0; print }' \
            $gtf_files_str >> "$gtf_output"
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
        String docker
    }

    Float quantExprGB = size(quantExprFiles, "GB")
    Float quantTrackingGB = size(quantTrackingFiles, "GB")
    # The merge concatenates the shard expr files and streams the shard trackings into one
    # gzip, so disk needs the localized inputs plus a copy of them: 2.2x total input.
    Float mergeDiskRawGB = (quantExprGB + quantTrackingGB) * 2.2 + 5.0

    command <<<
    set -eo pipefail

    lraa_version="$(LRAA --version | awk '{print $NF}')"
    export LRAA_VERSION_COMMENT="# LRAA version ${lraa_version}"

    lraa_merge_header.py \
        --version_comment "$LRAA_VERSION_COMMENT" \
        --inputs ~{sep=' ' quantExprFiles} \
        --output merge_header.txt

    prepend_lraa_merge_header() {
        local path="$1"
        local tmp="${path}.with_lraa_header"
        if [[ -s "$path" ]] && head -n 1 "$path" | grep -qxF "$LRAA_VERSION_COMMENT"; then
            return
        fi
        cat merge_header.txt "$path" > "$tmp"
        mv "$tmp" "$path"
    }

    merge_LRAA_quant_expr.py \
        --output "~{outputFilePrefix}.quant.expr" \
        --quant_files ~{sep=' ' quantExprFiles}

    prepend_lraa_merge_header "~{outputFilePrefix}.quant.expr"

    python <<CODE
import json
import gzip

tracking_files_json = '["' + '~{sep='","' quantTrackingFiles}' + '"]'
tracking_files_list = json.loads(tracking_files_json)
header_lines = open("merge_header.txt", "rt").read()

with gzip.open("~{outputFilePrefix}.quant.tracking.gz", "wt") as ofh:
    ofh.write(header_lines)
    wrote_header = False
    expected_header = None
    for i, tracking_file in enumerate(tracking_files_list):
        openf = gzip.open if tracking_file.split(".")[-1] == "gz" else open
        with openf(tracking_file, "rt") as fh:
            header = None
            for line in fh:
                if line.startswith("#"):
                    continue
                if header is None:
                    header = line
                    if expected_header is None:
                        expected_header = header
                    elif header != expected_header:
                        raise RuntimeError(
                            "Cannot merge tracking files with different schemas: {} differs from the first input".format(
                                tracking_file
                            )
                        )
                    if not wrote_header:
                        print(header, file=ofh, end='')
                        wrote_header = True
                    continue
                print(line, file=ofh, end='')
CODE

    >>>

    output {
        File mergedQuantExprFile = "~{outputFilePrefix}.quant.expr"
        File mergedQuantTrackingFile = "~{outputFilePrefix}.quant.tracking.gz"
    }
    
    runtime {
        docker: docker
        cpu: 1
        memory: "4 GiB"
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

    python -c '
import csv
from pathlib import Path

summary_files = [Path(p) for p in """~{sep="\n" summaryFiles}""".splitlines() if p.strip()]
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
'
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
        String docker
  }

    Float bam_size_gb = size(bam, "GB")
    Float estimated_disk = ceil(bam_size_gb * 2.2 + 20.0)
    Float disk_gb = if estimated_disk > 150.0 then estimated_disk else 150.0
    Int disk_gb_int = ceil(disk_gb)

  command <<<
    set -ex
        # -F 0x904 drops unmapped/secondary/supplementary so this matches
        # count_reads_from_bam() in the LRAA driver: one count per genome-mapped
        # read. Both paths must agree or TPM depends on which one ran.
        samtools view -@ ~{samtools_threads} -c -F 0x904 ~{bam}

  >>>

  runtime {
        docker: docker
        disks: "local-disk " + disk_gb_int + " SSD"
        cpu: samtools_threads
        memory: "8G"
  }
  output {
    Int count = read_int(stdout())
  }
}

