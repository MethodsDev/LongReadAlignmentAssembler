version 1.0

workflow LRAA_ORF_prediction_wf {
    input {
        String sample_id
        File input_gtf
        File genome_fasta

        # Optional parameters
        Boolean complete_orfs_only = false
        Boolean single_best_only = false
        Boolean no_refine_starts = true
        Int min_prot_length = 100

        # Optional TransDecoder homology-retention parameters
        File? blast_search_pep
        Float blast_evalue = 0.00001
        Int blast_threads = 4
        String blast_tool = "diamond"  # "diamond" or "blastp"
        File? pfam_search_db

        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    }

    call LRAA_ORF_prediction_task {
        input:
            sample_id = sample_id,
            input_gtf = input_gtf,
            genome_fasta = genome_fasta,
            complete_orfs_only = complete_orfs_only,
            single_best_only = single_best_only,
            no_refine_starts = no_refine_starts,
            min_prot_length = min_prot_length,
            blast_search_pep = blast_search_pep,
            blast_evalue = blast_evalue,
            blast_threads = blast_threads,
            blast_tool = blast_tool,
            pfam_search_db = pfam_search_db,
            docker = docker
    }

    output {
        # Main outputs
        File genome_gff3 = LRAA_ORF_prediction_task.genome_gff3
        File genome_bed_gz = LRAA_ORF_prediction_task.genome_bed_gz
        File genome_bed_tbi = LRAA_ORF_prediction_task.genome_bed_tbi

        # TransDecoder outputs
        File transdecoder_gff3 = LRAA_ORF_prediction_task.transdecoder_gff3
        File transdecoder_pep = LRAA_ORF_prediction_task.transdecoder_pep
        File transdecoder_cds = LRAA_ORF_prediction_task.transdecoder_cds
        File transdecoder_bed = LRAA_ORF_prediction_task.transdecoder_bed

        # Optional homology outputs
        File? blastp_outfmt6 = LRAA_ORF_prediction_task.blastp_outfmt6
        File? pfam_domtblout = LRAA_ORF_prediction_task.pfam_domtblout

        # TransDecoder output directory archive
        File transdecoder_dir_tar_gz = LRAA_ORF_prediction_task.transdecoder_dir_tar_gz
    }
}

task LRAA_ORF_prediction_task {
    input {
        String sample_id
        File input_gtf
        File genome_fasta

        Boolean complete_orfs_only
        Boolean single_best_only
        Boolean no_refine_starts
        Int min_prot_length

        File? blast_search_pep
        Float blast_evalue
        Int blast_threads
        String blast_tool
        File? pfam_search_db

        String docker
    }

    command <<<
        set -ex

        # Use TransDecoder from Docker image
        LRAA_isoforms_to_ORFs.py \
            --genome ~{genome_fasta} \
            --gtf ~{input_gtf} \
            --output_prefix ~{sample_id} \
            --td_dir ${TRANSDECODER_DIR} \
            -m ~{min_prot_length} \
            ~{if complete_orfs_only then "--complete_orfs_only" else ""} \
            ~{if single_best_only then "--single_best_only" else ""} \
            ~{if no_refine_starts then "--no_refine_starts" else ""} \
            ~{if defined(blast_search_pep) then "--blast_search_pep ~{blast_search_pep}" else ""} \
            ~{if defined(blast_search_pep) then "--blast_evalue ~{blast_evalue}" else ""} \
            ~{if defined(blast_search_pep) then "--blast_threads ~{blast_threads}" else ""} \
            ~{if defined(blast_search_pep) then "--blast_tool ~{blast_tool}" else ""} \
            ~{if defined(pfam_search_db) then "--pfam_search_db ~{pfam_search_db}" else ""}

        # Gzip blastp output if it exists
        if [ -f "~{sample_id}.blastp.outfmt6" ]; then
            gzip "~{sample_id}.blastp.outfmt6"
        fi

        # Tar the TransDecoder output directory
        tar czf ~{sample_id}.transdecoder_dir.tar.gz ~{sample_id}.cdna.fasta.transdecoder_dir

    >>>

    output {
        # Main genome-coordinate outputs
        File genome_gff3 = "~{sample_id}.cdna.fasta.transdecoder.genome.gff3"
        File genome_bed_gz = "~{sample_id}.cdna.fasta.transdecoder.genome.bed.gz"
        File genome_bed_tbi = "~{sample_id}.cdna.fasta.transdecoder.genome.bed.gz.tbi"

        # TransDecoder transcript-coordinate outputs
        File transdecoder_gff3 = "~{sample_id}.cdna.fasta.transdecoder.gff3"
        File transdecoder_pep = "~{sample_id}.cdna.fasta.transdecoder.pep"
        File transdecoder_cds = "~{sample_id}.cdna.fasta.transdecoder.cds"
        File transdecoder_bed = "~{sample_id}.cdna.fasta.transdecoder.bed"

        # Optional homology outputs
        File? blastp_outfmt6 = "~{sample_id}.blastp.outfmt6.gz"
        File? pfam_domtblout = "~{sample_id}.pfam.domtblout"

        # TransDecoder output directory archive
        File transdecoder_dir_tar_gz = "~{sample_id}.transdecoder_dir.tar.gz"
    }

    runtime {
        docker: docker
        disks: "local-disk " + ceil(4 * size(genome_fasta, "GB") + 4 * size(input_gtf, "GB") + 4 * size(blast_search_pep, "GB") + 4 * size(pfam_search_db, "GB") + 50) + " HDD"
        cpu: blast_threads + 2
        memory: "32G"
    }
}
