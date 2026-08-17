version 1.0

# Probe workflow for testing/call_cache_soundness/test_cache_soundness.sh.
#
# It exists to answer one question: does a miniwdl task run here observe the
# CURRENT content of its container image, or can it be served an output produced
# under a previous build of the same image tag?  So the only thing it does is
# report a marker baked into the image, and its only input is the image
# reference -- the same input, as a string, that miniwdl folds into its call
# cache key and that a moving tag makes useless there.
#
# Deliberately not one of the real LRAA workflows: the property under test is a
# property of the harness, and a probe that takes seconds can be run three times
# per suite invocation while lraa-core cannot.

workflow call_cache_probe {
  input {
    String docker
  }

  call read_image_marker {
    input:
      docker = docker
  }

  output {
    String marker = read_image_marker.marker
  }
}

task read_image_marker {
  input {
    String docker
  }

  command <<<
    set -eu
    cat /probe_marker
  >>>

  output {
    String marker = read_string(stdout())
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "1 GiB"
  }
}
