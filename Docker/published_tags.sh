#!/bin/bash
#
# What is published, and which of it is the newest CODE.
#
#     bash published_tags.sh                 # lraa-core
#     bash published_tags.sh lraa-sc         # any repository
#     bash published_tags.sh lraa-sc 10      # more rows
#
# The registry cannot answer "which is latest" on its own, and the obvious query
# answers a different question convincingly:
#
#   * `UPDATE_TIME` is the last time the ENTRY changed -- including a tag being
#     added or restored -- not when the image was built.  OBSERVED: restoring the
#     `0.30.0` tag moved its entry above a newer build, and `0.30.0-b105a80`
#     reported the same minute as an 0.31.0 build made hours later.  Sorting by it
#     and taking the top row is how you end up naming older code.
#   * A tag is a POINTER.  The digest beside it is the thing nothing can move;
#     pass `<repo>@sha256:...` to anything that must not change underneath it.
#
# So the commit is resolved for every `<version>-<shortsha>` tag and dated from
# git, which is the only ordering that reflects the code.  A tag whose commit this
# clone does not know is reported as such rather than guessed at.
#
# `gcloud` user credentials are subject to Google Cloud session control and stop
# working unattended.  On a GCE VM the attached service account has no such
# policy, so this falls back to a metadata token and keeps answering.

set -uo pipefail

REPO=${1:-lraa-core}
LIMIT=${2:-5}
REGISTRY=us-central1-docker.pkg.dev
PROJECT=methods-dev-lab
AR_REPO=lraa
GIT_DIR_ARG=(-C "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)")

printf '%s/%s/%s/%s\n' "${REGISTRY}" "${PROJECT}" "${AR_REPO}" "${REPO}"
printf 'rows are RECENTLY UPDATED registry entries, not newest code; read the commit column\n\n'

# Annotate a comma-joined tag list with the commit its <version>-<shortsha> names.
describe_commit() {
    local tags=$1 sha="" t
    for t in ${tags//,/ }; do
        case "${t}" in
            *-*)
                sha=${t##*-}
                # a shortsha is hex; -testing and friends are not
                if [[ "${sha}" =~ ^[0-9a-f]{7,40}$ ]]; then break; fi
                sha=""
                ;;
        esac
    done
    if [ -z "${sha}" ]; then
        echo "no commit in tag"
        return
    fi
    local subject
    if subject=$(git "${GIT_DIR_ARG[@]}" log -1 --format='%cd %h %s' --date=short "${sha}" 2>/dev/null); then
        if git "${GIT_DIR_ARG[@]}" merge-base --is-ancestor "${sha}" origin/devel 2>/dev/null; then
            echo "${subject:0:64}"
        else
            echo "${subject:0:52} [NOT on origin/devel]"
        fi
    else
        echo "commit ${sha} unknown to this clone (git fetch?)"
    fi
}

rows=""
if json=$(gcloud artifacts docker images list \
            "${REGISTRY}/${PROJECT}/${AR_REPO}/${REPO}" \
            --include-tags --sort-by=~UPDATE_TIME --limit="${LIMIT}" \
            --format=json 2>/dev/null) && [ -n "${json}" ] && [ "${json}" != "[]" ]; then
    rows=$(printf '%s' "${json}" | python3 -c '
import json, sys
for r in json.load(sys.stdin):
    print("%s\t%s" % (",".join(r.get("tags") or []) or "(untagged)", r["version"][:19]))
')
else
    echo "  (gcloud unauthenticated; using a token and the registry v2 API -- no digests or order)" >&2
    TOKEN=$(gcloud auth print-access-token 2>/dev/null) || TOKEN=""
    if [ -z "${TOKEN}" ]; then
        TOKEN=$(curl -s -m 5 -H "Metadata-Flavor: Google" \
          "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token" 2>/dev/null \
          | python3 -c 'import json,sys; print(json.load(sys.stdin)["access_token"])' 2>/dev/null) || TOKEN=""
    fi
    if [ -z "${TOKEN}" ]; then
        echo "  no usable credentials: run 'gcloud auth login', or run on a VM with an attached service account" >&2
        exit 1
    fi
    rows=$(curl -s -m 30 -H "Authorization: Bearer ${TOKEN}" \
      "https://${REGISTRY}/v2/${PROJECT}/${AR_REPO}/${REPO}/tags/list" \
      | python3 -c '
import json, sys
tags = json.load(sys.stdin).get("tags", [])
for t in sorted(tags, reverse=True):
    print("%s\t-" % t)
')
fi

[ -n "${rows}" ] || { echo "  (nothing published)"; exit 0; }

while IFS=$'\t' read -r tags digest; do
    [ -n "${tags}" ] || continue
    printf '  %-38s %-19s %s\n' "${tags}" "${digest}" "$(describe_commit "${tags}")"
done <<< "${rows}"

echo
echo "  to see what an image CONTAINS rather than what it is called:"
echo "    docker inspect <image> --format '{{index .Config.Labels \"org.opencontainers.image.revision\"}}'"
