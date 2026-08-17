#!/bin/bash
# Snapshot every rescue-observable output of the bundled rescue corpora, so the same
# snapshot taken before and after a refactor can be compared byte for byte.
#
# Not a test harness: it supplies no expectations of its own. The per-case Makefiles own
# the invocation and the expected numbers for the txalign cases; the SIRV invocations are
# copied verbatim from testing/sirvs/Makefile. This only runs them and copies what they
# produced.
#
# usage: rescue_parity_snapshot.sh <snapshot_dir>
set -u
snap_dir="${1:?usage: rescue_parity_snapshot.sh <snapshot_dir>}"
here="$(cd "$(dirname "$0")" && pwd)"
lraa="$here/../LRAA"
mkdir -p "$snap_dir"
snap_dir="$(cd "$snap_dir" && pwd)"
rc=0

snapshot_outputs() {
    # $1 = run dir, $2 = destination
    run_dir="$1"; out="$2"
    mkdir -p "$out"
    # Tracking is gzipped: store the decompressed payload, since a gzip header carries a
    # timestamp and would differ on every run while the content was identical. The
    # "# LRAA " comment lines carry the version and the argv, which differ by invocation
    # path rather than by behaviour.
    for f in "$run_dir"/*.quant.expr "$run_dir"/*.gtf "$run_dir"/*.bed; do
        [ -f "$f" ] && grep -v '^# LRAA ' "$f" > "$out/$(basename "$f")"
    done
    for f in "$run_dir"/*.quant.tracking.gz; do
        [ -f "$f" ] && gzip -dc "$f" > "$out/$(basename "$f" .gz)"
    done
    # Stored under a fixed name so two snapshots compare even if the summary is renamed.
    for f in "$run_dir"/*.read_assignment.summary.tsv; do
        [ -f "$f" ] && cp "$f" "$out/rescue_summary.tsv"
    done
}

snapshot_rescue_log() {
    # $1 = log, $2 = destination. Rescue's own accounting: what it was offered, what came
    # back, and which rule declined the rest.
    grep -E 'transcriptome rescue: offered=|rescue summary: total=|early transcriptome rescue: retrying' "$1" \
        | sed -e 's/^[[:space:]]*//' -e 's/^\[[^]]*\] //' > "$2/rescue_log_lines.txt"
}

# ---- rescue execution corpora: one gene each, mouse microexon and plant repetitive ----
cases=$(cd "$here" && ls -d txalign_path_rescue/microexon/*/ txalign_path_rescue/repetitive/*/ 2>/dev/null)
for case_rel in $cases; do
    case_dir="$here/$case_rel"
    [ -f "$case_dir/Makefile" ] || continue
    tag=$(echo "$case_rel" | tr '/' '_' | sed 's/_$//')
    log="$snap_dir/$tag.log"
    ( cd "$case_dir" && make test ) >"$log" 2>&1
    case_rc=$?
    echo "CASE $tag exit=$case_rc"
    [ "$case_rc" -eq 0 ] || rc=1
    snapshot_outputs "$case_dir" "$snap_dir/$tag"
    snapshot_rescue_log "$log" "$snap_dir/$tag"
done

# ---- SIRVs: 7 contigs, so this is the one that exercises the per-work-unit summary
# ---- write and the parent-side merge, and both --HiFi entry points into rescue.
sirv_data="$here/sirvs/data"
run_sirv() {
    tag="$1"; shift
    run_dir="$snap_dir/$tag.run"
    rm -rf "$run_dir"
    mkdir -p "$run_dir"
    log="$snap_dir/$tag.log"
    ( cd "$run_dir" && "$lraa" --genome "$sirv_data/SIRVs1-7.genome.fa" \
        --bam "$sirv_data/sim.fasta.mm2.bam" --HiFi --cpu_budget 4 "$@" ) >"$log" 2>&1
    case_rc=$?
    echo "CASE $tag exit=$case_rc"
    [ "$case_rc" -eq 0 ] || rc=1
    snapshot_outputs "$run_dir" "$snap_dir/$tag"
    snapshot_rescue_log "$log" "$snap_dir/$tag"
    rm -rf "$run_dir"
}

run_sirv sirv_quant_only --quant_only --gtf "$sirv_data/SIRVs1-7.annot.gtf"
run_sirv sirv_ref_guided --norm 0 --gtf "$sirv_data/SIRVs1-7.annot.half.gtf"

exit $rc
