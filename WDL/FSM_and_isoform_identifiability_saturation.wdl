version 1.0

workflow IsoformSaturationRow {
  input {
    File input_tsv
    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"

    String iso_col = "matching_isoforms"
    String cat_col = "sqanti_cat"
    String fsm_value = "FSM"
    Int seed = 42
    Boolean no_shuffle = false
    Int thin = 10000
    Int limit_rows = 0

    Int cpu = 2
    Int memory_gb = 8
    Int preemptible_tries = 0
    Float disk_multiplier = 3.0
    Int disk_min_gb = 20
  }

  call RunSaturation {
    input:
      input_tsv = input_tsv,
      docker = docker,
      iso_col = iso_col,
      cat_col = cat_col,
      fsm_value = fsm_value,
      seed = seed,
      no_shuffle = no_shuffle,
      thin = thin,
      limit_rows = limit_rows,
      cpu = cpu,
      memory_gb = memory_gb,
      preemptible_tries = preemptible_tries,
      disk_multiplier = disk_multiplier,
      disk_min_gb = disk_min_gb
  }

  output {
    File thin_tsv        = RunSaturation.thin_tsv
    File fit_summary_tsv = RunSaturation.fit_summary_tsv
    File plot_png        = RunSaturation.plot_png
  }
}



    
################################################################################
# Task: Run read_FSM_and_identifiability_saturation_fit.py on one input file
################################################################################
task RunSaturation {
  input {
    File input_tsv
    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    String iso_col = "matching_isoforms"
    String cat_col = "sqanti_cat"
    String fsm_value = "FSM"
    Int seed = 42
    Boolean no_shuffle = false
    Int thin = 10000
    Int limit_rows = 0

    # Resources
    Int cpu = 2
    Int memory_gb = 8
    Int preemptible_tries = 0

    # Disk autosize multiplier (working room for parsing + outputs)
    Float disk_multiplier = 3.0
    Int disk_min_gb = 20
  }

  meta {
    description: "Compute Identifiable vs FSM saturation curves + model fits"
  }

  parameter_meta {
    input_tsv: { description: "Input TSV or TSV.GZ with columns: matching_isoforms, sqanti_cat" }
    docker:    { description: "Docker image (Artifact Registry URL)" }
    iso_col:   { description: "Column containing isoform IDs" }
    cat_col:   { description: "Column containing SQANTI category" }
    fsm_value: { description: "SQANTI value treated as FSM" }
    seed:      { description: "Random seed for shuffling" }
    no_shuffle:{ description: "If true, do not shuffle read order" }
    thin:      { description: "Target #points to keep in curve TSV/plot" }
    limit_rows:{ description: "Process only first N rows (0 = all)" }
    cpu:       { description: "vCPUs" }
    memory_gb: { description: "Memory (GB)" }
    preemptible_tries: { description: "Number of preemptible attempts (GCE only)" }
    disk_multiplier: { description: "Autosize factor for disk (= input_size_gb * factor)" }
    disk_min_gb:     { description: "Minimum disk GB" }
  }

  command <<<
    set -euo pipefail

    # Localize input
    IN="~{input_tsv}"

    # Run
    read_FSM_and_identifiability_saturation_fit.py \
      -i "${IN}" \
      --iso-col "~{iso_col}" \
      --cat-col "~{cat_col}" \
      --fsm-value "~{fsm_value}" \
      --seed ~{seed} \
      ~{if no_shuffle then "--no-shuffle" else ""} \
      --thin ~{thin} \
      --limit-rows ~{limit_rows}

    # The script auto-names outputs based on input basename.
    # Capture them via globs.
    echo "Capturing outputs..."
    ls -1 *.saturation_identifiable_vs_FSM.thin.tsv.gz > thin_list.txt
    ls -1 *.saturation_identifiable_vs_FSM.fit_summary.tsv > fit_list.txt
    ls -1 *.saturation_identifiable_vs_FSM.fit.png > png_list.txt
  >>>

  output {
    File thin_tsv = read_string("thin_list.txt")
    File fit_summary_tsv = read_string("fit_list.txt")
    File plot_png = read_string("png_list.txt")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"

    # Autosize disk from input
    disks: "local-disk " + (ceil(size(input_tsv, "GB") * disk_multiplier) + disk_min_gb) + " HDD"

    preemptible: preemptible_tries
  }
}
