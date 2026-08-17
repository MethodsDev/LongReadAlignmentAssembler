#!/usr/bin/env bash
#
# Print the absolute miniwdl call-cache directory to use for a wdl test run in
# the current working directory, namespaced by the CONTENT of the container
# images that run will use.
#
#     MINIWDL__CALL_CACHE__DIR="$(testing/miniwdl_call_cache_dir.sh IMAGE...)" \
#         miniwdl run --cfg miniwdl.test.cfg ...
#
# lraa_test_docker.mk wires this into $(MINIWDL_RUN); nothing else should need
# to call it, apart from the soundness test that guards it.
#
#
# Why this exists
# ---------------
# miniwdl keys its call cache on the task name, a digest of the task's WDL
# source, and a digest of the task's inputs -- WDL/runtime/task.py:
#
#     cache_key = f"{task.name}/{task.digest}/{Value.digest_env(inputs)}"
#
# The container image reaches a task as an input, so the key does contain it,
# but only as the STRING the caller passed.  The test images are named by a
# moving tag (:testing, from Docker/build_docker.testing.sh).  Rebuild the tag
# and that string is unchanged, so every entry written against the previous
# build still matches the key computed for the new one.  `make test_wdls` then
# replays outputs produced by the code that was just replaced, exits 0 in about
# two minutes, and never starts a container.  That is not hypothetical: it
# happened on 2026-08-16 and the suite reported green.
#
# check_test_image_revision in lraa_test_docker.mk cannot catch this, and it
# passed on the day.  It proves the image was built from HEAD; it says nothing
# about whether the workflow ran against the image or replayed an entry written
# under an older one.  Two guards, and a stale result walked between them.
#
# The cache key must name every determining input, and the image content is a
# determining input.  miniwdl's key construction is not ours to change, so the
# image content is folded into the cache DIRECTORY, which is configurable.
# Directory plus key is the effective key.  An entry written under one image
# content is then unreachable from another, while entries written under an image
# still in use stay reachable -- which is the whole difference between this and
# deleting the cache on every run.  Repeat runs against one build still replay,
# which is what the cache is for; that is also why a target whose assertion is
# about a task having actually run must use $(MINIWDL_RUN_UNCACHED) instead.
#
#
# Why image content and not a digest-pinned reference
# ---------------------------------------------------
# The more idiomatic fix is to pass the image by digest --
# lraa-core@sha256:... -- so the string miniwdl hashes is itself
# content-addressed.  It does not work here, and it fails in precisely the case
# this is defending:
#
#   - A digest reference can only be formed from .RepoDigests, which docker
#     populates from a push or a pull.  A freshly BUILT image has none.  The
#     scenario is "rebuild the tag, the string does not change" -- and a local
#     rebuild is exactly an image with no repo digest, so the pin would fall
#     back to the tag on the one run that most needs it.
#   - .Id (the image config digest) is present for every local image however it
#     arrived, and changes on every rebuild.
#   - A digest reference also only covers the images we remember to name on the
#     miniwdl command line.  A workflow that grows another image input, left at
#     its default, would keep replaying across rebuilds.  A namespace derived
#     from all the test images covers those too.
#
# Content rather than the org.opencontainers.image.revision label, for the same
# reason: rebuilding at one commit is the normal way to test a Docker/ or
# dependency change, and the revision label does not move when you do.
#
#
# Apptainer
# ---------
# With no docker CLI there is no local image to fingerprint -- miniwdl's
# singularity backend converts docker://<tag> straight to a SIF.  The namespace
# falls back to HEAD, which the revision guard already ties the image to, and
# says so on stderr.  That is weaker: a re-push at one commit is invisible to
# it.  Recorded rather than papered over.
set -u
set -o pipefail

if [ "$#" -eq 0 ]; then
    echo "usage: ${0##*/} IMAGE_REF [IMAGE_REF...]" >&2
    exit 2
fi

base="$PWD/__miniwdl_call_cache"

emit() {
    # $1 namespace, $2 provenance.  The provenance file is written into the
    # namespace so `ls __miniwdl_call_cache/*/_images` answers "which image did
    # this entry come from" without re-deriving anything.  miniwdl would create
    # the directory itself; doing it here is what lets the note land.
    local dir="$base/$1"
    mkdir -p "$dir" || exit 1
    printf '%s\n' "$2" >"$dir/_images"
    printf '%s\n' "$dir"
}

# Escape hatch, for pointing a run at a namespace by hand -- bisecting a cache
# entry, or reusing one deliberately across an image change:
#
#     make test_wdls LRAA_TEST_CALL_CACHE_NAMESPACE=scratch
#
if [ -n "${LRAA_TEST_CALL_CACHE_NAMESPACE:-}" ]; then
    emit "$LRAA_TEST_CALL_CACHE_NAMESPACE" "namespace forced by LRAA_TEST_CALL_CACHE_NAMESPACE"
    exit 0
fi

if ! command -v docker >/dev/null 2>&1; then
    head=$(git rev-parse HEAD 2>/dev/null) || head=""
    echo "" >&2
    echo "docker is not on PATH, so image content cannot be read and the call-cache" >&2
    echo "namespace falls back to the commit.  Expected under Apptainer.  A rebuild or" >&2
    echo "re-push of the image tag at this same commit will NOT invalidate cached task" >&2
    echo "outputs; delete __miniwdl_call_cache by hand after one." >&2
    echo "" >&2
    emit "nodocker-${head:-unknown}" "no docker; namespace from git HEAD ${head:-unknown}"
    exit 0
fi

# An absent image is recorded as absent rather than pulled: pulling images the
# caller's targets may not even use would put a multi-gigabyte surprise in front
# of a one-directory test run.  It costs nothing to be wrong here.  The run
# itself pulls whatever it needs, so the image is present on the next pass, the
# namespace changes, and the only consequence is one extra miss.
manifest=$(
    for ref in "$@"; do
        id=$(docker image inspect --format '{{.Id}}' "$ref" 2>/dev/null) || id=""
        printf '%s=%s\n' "$ref" "${id:-absent}"
    done | sort
)

fingerprint=$(printf '%s\n' "$manifest" | sha256sum | cut -c1-16)
if [ -z "$fingerprint" ]; then
    echo "could not fingerprint test images: $*" >&2
    exit 1
fi

emit "img-$fingerprint" "$manifest"
