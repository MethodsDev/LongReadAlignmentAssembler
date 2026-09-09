#!/bin/bash
# Container-level execution test for the streaming quant-only canonical-path
# collision. Same case as test_streaming_quant_canonical_path_collision.py, with
# no pytest dependency -- for CI that tests a built IMAGE rather than a checkout.
#
#   ./run_test.sh <lraa-core sif|docker uri>
#
# Exit 0 = fixed. Exit 1 = the v0.31.0 collision is back. Runs in ~2 s.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMG="${1:?usage: run_test.sh <lraa-core sif or docker uri>}"
W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT
cp "$HERE"/fixture/* "$W"/
echo '{"HiFi": true, "cpu_budget": 2}' > "$W/cfg.json"

# --containall --no-mount hostfs --fakeroot are not cosmetic. An Apptainer install
# missing etc/apptainer/capability.json -- which a relocated/unprivileged install can
# be -- fails a PLAIN `apptainer exec` with
#     FATAL: while opening capability config file: ... capability.json: no such file
# BEFORE the container starts, so LRAA never runs and the test reports a failure
# indistinguishable from the collision it exists to catch. These are the same flags
# miniwdl's singularity backend passes for every task in the production workflows.
# Override with APPTAINER_EXEC_FLAGS="" where a full install makes them unnecessary.
APPTAINER_EXEC_FLAGS="${APPTAINER_EXEC_FLAGS---containall --no-mount hostfs --fakeroot}"
if [[ "$IMG" == docker://* || "$IMG" == *.sif ]]; then
    RUN=(apptainer exec $APPTAINER_EXEC_FLAGS -B "$W:$W" --pwd "$W" "$IMG" python3 /usr/local/src/LRAA/LRAA)
else
    RUN=(docker run --rm -v "$W:$W" -w "$W" "$IMG" python3 /usr/local/src/LRAA/LRAA)
fi

"${RUN[@]}" \
    --genome locus.fa --bam 'locus.strand.+.bam' --no_chunk --gtf locus.gtf \
    --quant_only --bam_for_sg locus.plus.norm.bam --no_norm \
    --num_total_reads 81523164 --cpu_budget 1 --output_prefix q \
    --min_mapping_quality 0 --min_mapping_quality_for_final_quant 0 \
    --HiFi --stream_reads --stream_reads_rescue_unassigned \
    --config_update cfg.json > "$W/run.out" 2> "$W/run.err"
RC=$?

fail() { echo "FAIL: $1" >&2; sed -n '/two multipaths/,+2p' "$W/run.err" >&2; exit 1; }

grep -q 'two multipaths map to canonical path' "$W/run.out" "$W/run.err" \
    && fail "canonical-path collision has returned"
[ $RC -eq 0 ] || fail "LRAA exited $RC$(tail -5 "$W/run.err")"
[ -s "$W/q.LRAA.quant-only.quant.expr" ] || fail "no quant.expr"
[ -s "$W/q.LRAA.quant-only.quant.tracking.gz" ] || fail "no quant.tracking.gz"
grep -q 'comp-1119' "$W/q.LRAA.quant-only.quant.expr" \
    || fail "the component that collided (comp-1119) is absent from quant.expr"

N=$(grep -vc '^#\|^gene_id' "$W/q.LRAA.quant-only.quant.expr")
echo "PASS: quant-only streaming run clean, $N transcript rows, comp-1119 present"
