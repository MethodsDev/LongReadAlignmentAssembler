version 1.0

workflow LRAA_sqanti_like_reads_eval_wf {
    input {
        String sample_id
        File ref_annot_GTF
      
        File? input_BAM
        File? input_BAI
        File? input_GTF
        
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    }


    call LRAA_sqanti_like_reads_eval_task {
        input:
          sample_id = sample_id,
          ref_annot_GTF = ref_annot_GTF,
          input_BAM = input_BAM,
          input_BAI = input_BAI,
          input_GTF = input_GTF,
          docker = docker
    }

    output {

        File iso_cats_tsv = LRAA_sqanti_like_reads_eval_task.iso_cats_tsv
        File iso_cats_summary_counts_tsv = LRAA_sqanti_like_reads_eval_task.iso_cats_summary_counts_tsv
        File iso_cats_summary_counts_pdf = LRAA_sqanti_like_reads_eval_task.iso_cats_summary_counts_pdf
        File iso_cats_bam  = LRAA_sqanti_like_reads_eval_task.iso_cats_bam 
    }

}

task LRAA_sqanti_like_reads_eval_task {
    input {
         String sample_id
      
        File? input_BAM
        File? input_BAI

        File? input_GTF

        File ref_annot_GTF
      
        String docker
    }

    String inputs_string = if defined(input_BAM) then "--input_bam ~{input_BAM}" else if defined(input_GTF) then "--input_GTF ~{input_GTF}"  else ""

    
    command <<<
        set -ex

        SQANTI-like_cats_for_reads_or_isoforms.py --ref_gtf ~{ref_annot_GTF} --output_prefix ~{sample_id} ~{inputs_string}
      
    >>>

    output {

        File iso_cats_tsv = "~{sample_id}.iso_cats.tsv"
        File iso_cats_summary_counts_tsv = "~{sample_id}.iso_cats.summary_counts.tsv"
        File iso_cats_summary_counts_pdf = "~{sample_id}.iso_cats.summary_counts.pdf"
        File iso_cats_bam = "~{sample_id}.iso_cats.bam"
    }


  runtime {
    docker: docker
    disks: "local-disk " + ceil(4 * size(input_BAM, "GB") + 4 * size(input_GTF, "GB") ) + " HDD"
    memory: "32G"
  }


}
