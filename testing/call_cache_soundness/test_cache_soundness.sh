#!/usr/bin/env bash
#
# Asserts that a miniwdl task in this suite cannot be served an output that was
# produced under a different build of its container image.
#
# This is the regression test for the false green of 2026-08-16.  `make
# test_wdls` exited 0 in about two minutes that day without starting a single
# container: miniwdl keys its call cache on the image as the STRING it was
# given, the test images are named by the moving :testing tag, so every entry
# written against the previous build still matched.  The suite validated the
# build that had just been replaced and said so cheerfully.  Neither guard in
# place could see it -- check_test_image_revision proves the image matches HEAD,
# which says nothing about whether the workflow ran against the image.
#
# What is asserted, in order:
#
#   1. a first run against image content A reports A, and is not a cache hit
#   2. a second run against the same content A IS a cache hit, and still
#      reports A -- because a fix that merely deleted the cache every run would
#      pass every other assertion here while throwing away the whole point of
#      the cache, and this test should not accept that
#   3. after rebuilding the SAME TAG to content B, a run reports B and is not a
#      cache hit
#
# (3) is the defect.  Remove the MINIWDL__CALL_CACHE__DIR namespacing from
# testing/lraa_test_docker.mk and it fails twice over: the marker comes back as
# A, and miniwdl logs the hit that produced it.
#
# The probe image is the test image plus one echoed marker file, so "rebuild the
# image under the same tag" costs a second rather than a full LRAA build and the
# whole test costs well under a minute.  Deriving it from the real image rather
# than from alpine is not only fidelity: miniwdl runs a task's command under
# bash, which alpine does not ship, and pulling one in would put a network
# dependency in a test that has no other reason to need one.
set -u
set -o pipefail

cd "$(dirname "$0")" || exit 1

cfg="${1:-}"
base="${2:-}"
image=lraa-call-cache-probe:testing
outputs=probe.outputs.json

fail() {
    echo "" >&2
    echo "FAIL: call cache soundness -- $*" >&2
    echo "" >&2
    exit 1
}

# Docker-only.  Under Apptainer there is no docker daemon to build a throwaway
# image in and no local image whose content could be read; the namespacing this
# tests falls back to HEAD there, which miniwdl_call_cache_dir.sh says on stderr
# every run.  Skipping loudly beats a failure that means nothing.
if [ -n "$cfg" ] && [ -f "$cfg" ] \
   && grep -Eq '^[[:space:]]*container_backend[[:space:]]*=[[:space:]]*(singularity|podman|udocker)' "$cfg"; then
    echo "SKIP: call cache soundness needs the docker backend; $cfg selects another"
    exit 0
fi
if ! command -v docker >/dev/null 2>&1; then
    echo "SKIP: call cache soundness needs docker, which is not on PATH"
    exit 0
fi
if [ -z "$base" ]; then
    fail "no base image given; expected \$(LRAA_TEST_DOCKER) as the second argument"
fi
docker image inspect "$base" >/dev/null 2>&1 \
    || docker pull "$base" >/dev/null \
    || fail "cannot inspect or pull the probe base image $base"

# A fresh marker per invocation, so a pass can never be an artifact of an image
# or cache entry left behind by a previous invocation of this test.
stamp=$(date +%s)-$$
marker_a="content-A-$stamp"
marker_b="content-B-$stamp"

build_probe() {
    docker build -q -t "$image" - >/dev/null <<EOF || fail "could not build probe image"
FROM $base
RUN echo $1 > /probe_marker
EOF
}

probe_run() {
    # $1 log file.  Through make, so what is tested is the MINIWDL_RUN the real
    # test directories use, not a reimplementation of it.
    rm -f "$outputs"
    make --no-print-directory probe_run LRAA_TEST_DOCKER="$image" >"$1" 2>&1
    local rc=$?
    if [ "$rc" -ne 0 ]; then
        echo "--- last 40 lines of $1 ---" >&2
        tail -40 "$1" >&2
        fail "probe workflow failed (rc=$rc); see $1"
    fi
}

reported_marker() {
    # miniwdl -o writes {"dir": ..., "outputs": {...}}, not the bare outputs.
    python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["outputs"]["call_cache_probe.marker"])' \
        "$outputs" 2>/dev/null
}

cache_hit() {
    # $1 log file.  miniwdl logs "call cache hit" at notice level when it serves
    # a task from the cache.
    grep -q 'call cache hit' "$1"
}

rm -rf ./probe_outdir ./probe.run1.log ./probe.run2.log ./probe.run3.log \
    ./__miniwdl_call_cache "$outputs"

echo "1/3 first run against image content $marker_a"
build_probe "$marker_a"
probe_run probe.run1.log
got=$(reported_marker)
[ "$got" = "$marker_a" ] || fail "first run reported '$got', expected '$marker_a'"
! cache_hit probe.run1.log \
    || fail "first run into an empty cache was served from cache; the cache was not cleared"

echo "2/3 second run against the same content, which must hit"
probe_run probe.run2.log
got=$(reported_marker)
[ "$got" = "$marker_a" ] || fail "second run reported '$got', expected '$marker_a'"
cache_hit probe.run2.log || fail \
    "a repeat run against an unchanged image did NOT hit the call cache.
  Namespacing the cache by image content is supposed to separate image builds,
  not defeat caching -- if the cache is being discarded every run, the fix has
  been replaced by something that throws away what the cache is for."

echo "3/3 rebuilding the SAME TAG $image to content $marker_b"
build_probe "$marker_b"
probe_run probe.run3.log
got=$(reported_marker)
if [ "$got" = "$marker_a" ]; then
    fail "the workflow reported '$got' from an image that now contains '$marker_b'.
  The image was rebuilt under the same tag and the call cache replayed the
  PREVIOUS build's output.  This is exactly the false green of 2026-08-16: the
  suite would report success having tested the code that was just replaced.
  MINIWDL__CALL_CACHE__DIR namespacing in testing/lraa_test_docker.mk is
  missing or not reaching miniwdl."
fi
[ "$got" = "$marker_b" ] || fail "run after rebuild reported '$got', expected '$marker_b'"
! cache_hit probe.run3.log \
    || fail "run after the rebuild logged a call cache hit; the cache key does not name the image content"

echo "ok: call cache separates image builds and still serves repeat runs"
