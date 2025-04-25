version 1.0

workflow LRAA_sqanti_like_multi_sample_summary_wf {
    input {
        String output_prefix
        Array[File] iso_cats_summary_counts_tsv_files
        Array[File] iso_cats_raw_counts_tsv_files

        Int lengths_sample_n = 100000
    
        # pdf page size
        Int width = 12
        Int height = 12
    
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    
    }

    call LRAA_sqanti_like_multi_sample_CATS_summary_task {
        input:
            output_prefix=output_prefix,
            iso_cats_summary_counts_tsv_files = iso_cats_summary_counts_tsv_files,
            width=width,
            height=height,
            docker=docker
    }


    call LRAA_sqanti_like_multi_sample_LENGTHS_summary_task {
        input:
            output_prefix=output_prefix,
            iso_cats_raw_counts_tsv_files = iso_cats_raw_counts_tsv_files,
            lengths_sample_n = lengths_sample_n,
            width=width,
            height=height,
            docker=docker
    }
        
    output {

        File summary_CATS_stats_tsv = LRAA_sqanti_like_multi_sample_CATS_summary_task.summary_CATS_stats_tsv
        File summary_CATS_stats_plot_pdf = LRAA_sqanti_like_multi_sample_CATS_summary_task.summary_CATS_stats_plot_pdf

        File summary_LENGTHS_stats_tsv = LRAA_sqanti_like_multi_sample_LENGTHS_summary_task.summary_LENGTHS_stats_tsv
        File summary_LENGTHS_stats_plot_pdf = LRAA_sqanti_like_multi_sample_LENGTHS_summary_task.summary_LENGTHS_stats_plot_pdf
    
    }

}


task LRAA_sqanti_like_multi_sample_CATS_summary_task {

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
          --output_prefix ~{output_prefix}.CATS \
          --sample_stats ~{sep=" " iso_cats_summary_counts_tsv_files} \
          --width ~{width} \
          --height ~{height}



    >>>


    output {

         File summary_CATS_stats_tsv = "~{output_prefix}.CATS.tsv"
         File summary_CATS_stats_plot_pdf = "~{output_prefix}.CATS.pdf"

    }


    runtime {
        docker: docker
        bootDiskSizeGb: 30
        disks: "local-disk " + total_file_size + " HDD"
        cpu: 1
        memory: "4 GiB"
    }
    
}
    


task LRAA_sqanti_like_multi_sample_LENGTHS_summary_task {

    input {
        String output_prefix
        Array[File] iso_cats_raw_counts_tsv_files
        Int lengths_sample_n
        Int width
        Int height
        String docker
    
    }

    Int total_file_size = ceil(size(iso_cats_raw_counts_tsv_files, "GiB") * 2 + 8)
    
    command <<<
      set -ex

      plot_feature_lengths.summarize_mult_samples.Rscript \
          --output_prefix ~{output_prefix}.LENGTHS \
          --sample_iso_cats ~{sep=" " iso_cats_raw_counts_tsv_files} \
          --sample_n ~{lengths_sample_n} \
          --width ~{width} \
          --height ~{height}

    
    >>>


    output {

         File summary_LENGTHS_stats_tsv = "~{output_prefix}.LENGTHS.tsv"
         File summary_LENGTHS_stats_plot_pdf = "~{output_prefix}.LENGTHS.pdf"

    }


    runtime {
        docker: docker
        bootDiskSizeGb: 30
        disks: "local-disk " + total_file_size + " HDD"
        cpu: 1
        memory: "32 GiB"
    }
    
}
    
