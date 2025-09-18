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

    Int cpu = 2
    Int memory_gb = 8
    Int preemptible_tries = 0

    Float disk_multiplier = 3.0
    Int disk_min_gb = 20
  }

  # ---- derive the script's auto-named outputs from input basename ----
  String base = basename(input_tsv)
  String root = sub(base, "\\.(tsv|txt|csv)(\\.gz)?$", "")
  String thin_path = root + ".saturation_identifiable_vs_FSM.thin.tsv.gz"
  String fit_path  = root + ".saturation_identifiable_vs_FSM.fit_summary.tsv"
  String png_path  = root + ".saturation_identifiable_vs_FSM.fit.png"

  command <<<
    set -euo pipefail

    # Run
    read_FSM_and_identifiability_saturation_fit.py \
      -i "~{input_tsv}" \
      --iso-col "~{iso_col}" \
      --cat-col "~{cat_col}" \
      --fsm-value "~{fsm_value}" \
      --seed ~{seed} \
      ~{if no_shuffle then "--no-shuffle" else ""} \
      --thin ~{thin} \
      --limit-rows ~{limit_rows}

    # sanity-check the expected outputs exist (fail fast if not)
    test -s "~{thin_path}" || { echo "Missing ~{thin_path}" >&2; ls -lh; exit 1; }
    test -s "~{fit_path}"  || { echo "Missing ~{fit_path}"  >&2; ls -lh; exit 1; }
    test -s "~{png_path}"  || { echo "Missing ~{png_path}"  >&2; ls -lh; exit 1; }

    # optional: list for logs
    ls -lh "~{thin_path}" "~{fit_path}" "~{png_path}"
  >>>

  output {
    File thin_tsv        = thin_path
    File fit_summary_tsv = fit_path
    File plot_png        = png_path
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + (ceil(size(input_tsv, "GB") * disk_multiplier) + disk_min_gb) + " HDD"
    preemptible: preemptible_tries
  }
}

