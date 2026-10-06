version 1.0

workflow LRAA_sqanti_like_reads_eval_wf {
    input {
        String sample_id
        File ref_annot_GTF
      
        File? input_BAM
        File? input_BAI
        File? input_GTF
        
        # Named docker_sc, not docker: the R plotting step needs lraa-sc. Renamed 2026-10-06
        # so a saved Terra value for the old `docker` input (lraa:latest, the v0.17.7
        # default) cannot carry over and run this on lraa-core, which has no R.
        String docker_sc = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-sc:latest"
        Int? min_disk_GB

        # Cores requested for the task. The script runs with --CPU auto, so with
        # --input_BAM it classifies contigs on as many cores as the backend actually
        # grants (which can be fewer than requested).
        Int cpu = 8
    }


    call LRAA_sqanti_like_reads_eval_task {
        input:
          sample_id = sample_id,
          ref_annot_GTF = ref_annot_GTF,
          input_BAM = input_BAM,
          input_BAI = input_BAI,
          input_GTF = input_GTF,
          docker = docker_sc,
          min_disk_GB = min_disk_GB,
          cpu = cpu
    }

    output {

        File iso_cats_tsv = LRAA_sqanti_like_reads_eval_task.iso_cats_tsv
        File iso_cats_summary_counts_tsv = LRAA_sqanti_like_reads_eval_task.iso_cats_summary_counts_tsv
        File iso_cats_summary_counts_pdf = LRAA_sqanti_like_reads_eval_task.iso_cats_summary_counts_pdf
        File? iso_cats_bam  = LRAA_sqanti_like_reads_eval_task.iso_cats_bam 
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
        Int? min_disk_GB
        Int cpu
    }

    
    command <<<
        set -ex
 
        if [[ "~{input_BAM}" != "" ]]; then

             # The parallel mode reads the bam by contig, so it needs the index beside
             # it; the bam and bai can be localized to different directories. Link
             # both here, and index only if no bai was given.
             ln -s ~{input_BAM} input.bam
             if [[ "~{input_BAI}" != "" ]]; then
                 ln -s ~{input_BAI} input.bam.bai
             else
                 samtools index -@ ~{cpu} input.bam
             fi

             SQANTI-like_cats_for_reads_or_isoforms.py --ref_gtf ~{ref_annot_GTF} --output_prefix ~{sample_id} --input_bam input.bam --CPU auto

        elif [[ "~{input_GTF}" != "" ]]; then
 
             SQANTI-like_cats_for_reads_or_isoforms.py --ref_gtf ~{ref_annot_GTF} --output_prefix ~{sample_id} --input_gtf ~{input_GTF}

       fi


       gzip *.tsv

    >>>

    output {

        File iso_cats_tsv = "~{sample_id}.iso_cats.tsv.gz"
        File iso_cats_summary_counts_tsv = "~{sample_id}.iso_cats.summary_counts.tsv.gz"
        File iso_cats_summary_counts_pdf = "~{sample_id}.iso_cats.summary_counts.pdf"
        File? iso_cats_bam = "~{sample_id}.iso_cats.bam"
    }


  runtime {
    docker: docker
    disks: "local-disk " + (if defined(min_disk_GB) && ceil( (4 * size(input_BAM, "GB")) + (4 * size(input_GTF, "GB")) + (4 * size(ref_annot_GTF, "GB")) + 50 ) < select_first([min_disk_GB]) then select_first([min_disk_GB]) else ceil( (4 * size(input_BAM, "GB")) + (4 * size(input_GTF, "GB")) + (4 * size(ref_annot_GTF, "GB")) + 50 )) + " HDD"
    memory: "32G"
    cpu: cpu
  }


}
