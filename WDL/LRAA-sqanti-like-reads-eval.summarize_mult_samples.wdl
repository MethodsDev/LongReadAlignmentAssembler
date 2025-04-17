version 1.0

workflow LRAA_sqanti_like_multi_sample_summary_wf {
    input {
        String output_prefix
        Array[File] iso_cats_summary_counts_tsv_files

        # pdf page size
        Int width = 12
        Int height = 12
    
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    
    }

    call LRAA_sqanti_like_multi_sample_summary_task {
        input:
            output_prefix=output_prefix,
            iso_cats_summary_counts_tsv_files = iso_cats_summary_counts_tsv_files,
            width=width,
            height=height,
            docker=docker
    }
        
    output {

        File summary_stats_tsv = LRAA_sqanti_like_multi_sample_summary_task.summary_stats_tsv
        File summary_stats_plot_pdf = LRAA_sqanti_like_multi_sample_summary_task.summary_stats_plot_pdf   
    }

}


task LRAA_sqanti_like_multi_sample_summary_task {

    input {
        String output_prefix
        Array[File] iso_cats_summary_counts_tsv_files
        Int width
        Int height
        String docker
    
    }

    Int total_file_size = ceil(size(iso_cats_summary_counts_tsv_files, "GiB") * 2 + 8)
    
    command <<<
      set -ex

      plot_SQANTI_cats.summarize_mult_samples.Rscript \
          --output_prefix ~{output_prefix} \
          --sample_stats '~{sep="," iso_cats_summary_counts_tsv_files}' \
          --width ~{width} \
          --height ~{height}



    >>>


    output {

         File summary_stats_tsv = "~{output_prefix}.tsv"
         File summary_stats_plot_pdf = "~{output_prefix}.pdf"

    }


    runtime {
        docker: docker
        bootDiskSizeGb: 30
        disks: "local-disk " + total_file_size + " HDD"
        cpu: 1
        memory: "4 GiB"
    }
    
}
    
