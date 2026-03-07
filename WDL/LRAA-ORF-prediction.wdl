version 1.0

workflow LRAA_ORF_prediction_wf {
    input {
        String sample_id
        File input_gtf
        File genome_fasta

        # Optional parameters
        Boolean complete_orfs_only = false
        Boolean single_best_only = false
        Int min_prot_length = 100

        # Optional blastp/diamond parameters
        File? blastp_db
        Float blastp_evalue = 0.00001
        Int blastp_max_target_seqs = 1
        Int blastp_num_threads = 4
        String search_method = "diamond"  # "diamond" or "blastp"

        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    }

    call LRAA_ORF_prediction_task {
        input:
            sample_id = sample_id,
            input_gtf = input_gtf,
            genome_fasta = genome_fasta,
            complete_orfs_only = complete_orfs_only,
            single_best_only = single_best_only,
            min_prot_length = min_prot_length,
            blastp_db = blastp_db,
            blastp_evalue = blastp_evalue,
            blastp_max_target_seqs = blastp_max_target_seqs,
            blastp_num_threads = blastp_num_threads,
            search_method = search_method,
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

        # Optional blastp output
        File? blastp_outfmt6 = LRAA_ORF_prediction_task.blastp_outfmt6
    }
}

task LRAA_ORF_prediction_task {
    input {
        String sample_id
        File input_gtf
        File genome_fasta

        Boolean complete_orfs_only
        Boolean single_best_only
        Int min_prot_length

        File? blastp_db
        Float blastp_evalue
        Int blastp_max_target_seqs
        Int blastp_num_threads
        String search_method

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
            ~{if defined(blastp_db) then "--blastp_db ~{blastp_db}" else ""} \
            ~{if defined(blastp_db) then "--blastp_evalue ~{blastp_evalue}" else ""} \
            ~{if defined(blastp_db) then "--blastp_max_target_seqs ~{blastp_max_target_seqs}" else ""} \
            ~{if defined(blastp_db) then "--blastp_num_threads ~{blastp_num_threads}" else ""} \
            ~{if defined(blastp_db) then "--search_method ~{search_method}" else ""}

        # Gzip blastp output if it exists
        if [ -f "~{sample_id}.blastp.outfmt6" ]; then
            gzip "~{sample_id}.blastp.outfmt6"
        fi

    >>>

    output {
        # Main genome-coordinate outputs
        File genome_gff3 = "~{sample_id}.transcripts.fasta.transdecoder.genome.gff3"
        File genome_bed_gz = "~{sample_id}.transcripts.fasta.transdecoder.genome.bed.gz"
        File genome_bed_tbi = "~{sample_id}.transcripts.fasta.transdecoder.genome.bed.gz.tbi"

        # TransDecoder transcript-coordinate outputs
        File transdecoder_gff3 = "~{sample_id}.cdna.fasta.transdecoder.gff3"
        File transdecoder_pep = "~{sample_id}.cdna.fasta.transdecoder.pep"
        File transdecoder_cds = "~{sample_id}.cdna.fasta.transdecoder.cds"
        File transdecoder_bed = "~{sample_id}.cdna.fasta.transdecoder.bed"

        # Optional blastp output
        File? blastp_outfmt6 = "~{sample_id}.blastp.outfmt6.gz"
    }

    runtime {
        docker: docker
        disks: "local-disk " + ceil(4 * size(genome_fasta, "GB") + 4 * size(input_gtf, "GB") + 4 * size(blastp_db, "GB") + 50) + " HDD"
        cpu: blastp_num_threads + 2
        memory: "32G"
    }
}
