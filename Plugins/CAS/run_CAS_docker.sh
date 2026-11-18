#!/bin/bash

set -euo pipefail

# Docker wrapper script for run_CellAnnotationService.py
# This script runs the Cellarium CAS cell type annotation tool via Docker.
# The API token reaches the container as the CAS_API_TOKEN environment variable,
# so it never appears on a command line or in this script's output.

usage() {
    cat <<EOF
Usage: $0 --matrix-dir <path> --output-prefix <prefix> [token option] [options]

Required arguments:
  --matrix-dir <path>          Path to 10x-style matrix directory containing:
                               matrix.mtx[.gz], barcodes.tsv[.gz], features.tsv[.gz]
  --output-prefix <prefix>     Prefix for output files (will create <prefix>.h5ad and <prefix>.tsv);
                               must be under the matrix directory's parent, which is mounted

Cellarium CAS API token, one of (in order of preference):
  --api-token-file <path>      File holding the token
  CAS_API_TOKEN                Environment variable holding the token
  --api-token <token>          The token itself (avoid: visible in the process list
                               and shell history)

Optional arguments:
  --cas-model-name <name>      CAS model name (default: scanvi-classification-v8.hs-general;
                               'server-default' uses the CAS server's own default)
  --min-counts <int>           Drop cells with fewer total counts before CAS (default: 1)
  --chunk-size <int>           Number of cells per CAS chunk (default: 500)
  --min-acceptable-score <f>   Minimum acceptable evidence score (default: 0.2)
  --top-k <int>                Top-k CAS cell type calls per cell (default: 3)
  --obs-prefix <string>        Prefix for CAS cell type columns (default: cas_cell_type)
  --docker-image <image>       Docker image to use (default: us-central1-docker.pkg.dev/methods-dev-lab/lraa/cas)
  --help                       Show this help message

Example:
  $0 \\
    --matrix-dir /path/to/matrix_folder \\
    --api-token-file ~/.cas_token \\
    --output-prefix /path/to/sample1

This will create:
  /path/to/sample1.h5ad
  /path/to/sample1.tsv

EOF
    exit 1
}

# Default values
DOCKER_IMAGE="us-central1-docker.pkg.dev/methods-dev-lab/lraa/cas"
MATRIX_DIR=""
API_TOKEN_FILE=""
API_TOKEN_ARG=""
OUTPUT_PREFIX=""
CAS_MODEL_NAME=""
MIN_COUNTS=""
CHUNK_SIZE="500"
MIN_ACCEPTABLE_SCORE="0.2"
TOP_K="3"
OBS_PREFIX="cas_cell_type"

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --matrix-dir)           MATRIX_DIR="$2"; shift 2 ;;
        --api-token-file)       API_TOKEN_FILE="$2"; shift 2 ;;
        --api-token)            API_TOKEN_ARG="$2"; shift 2 ;;
        --output-prefix)        OUTPUT_PREFIX="$2"; shift 2 ;;
        --cas-model-name)       CAS_MODEL_NAME="$2"; shift 2 ;;
        --min-counts)           MIN_COUNTS="$2"; shift 2 ;;
        --chunk-size)           CHUNK_SIZE="$2"; shift 2 ;;
        --min-acceptable-score) MIN_ACCEPTABLE_SCORE="$2"; shift 2 ;;
        --top-k)                TOP_K="$2"; shift 2 ;;
        --obs-prefix)           OBS_PREFIX="$2"; shift 2 ;;
        --docker-image)         DOCKER_IMAGE="$2"; shift 2 ;;
        --help)                 usage ;;
        *)
            echo "Error: Unknown option: $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [[ -z "$MATRIX_DIR" ]]; then
    echo "Error: --matrix-dir is required"
    usage
fi

if [[ -z "$OUTPUT_PREFIX" ]]; then
    echo "Error: --output-prefix is required"
    usage
fi

# Resolve the token into CAS_API_TOKEN (file, then environment, then argument)
if [[ -n "$API_TOKEN_FILE" ]]; then
    if [[ ! -r "$API_TOKEN_FILE" ]]; then
        echo "Error: cannot read --api-token-file $API_TOKEN_FILE"
        exit 1
    fi
    CAS_API_TOKEN=$(tr -d '[:space:]' < "$API_TOKEN_FILE")
elif [[ -z "${CAS_API_TOKEN:-}" && -n "$API_TOKEN_ARG" ]]; then
    CAS_API_TOKEN="$API_TOKEN_ARG"
fi
if [[ -z "${CAS_API_TOKEN:-}" ]]; then
    echo "Error: no CAS API token; use --api-token-file, CAS_API_TOKEN or --api-token"
    usage
fi
export CAS_API_TOKEN
unset API_TOKEN_ARG

# Convert paths to absolute paths
MATRIX_DIR=$(cd "$(dirname "$MATRIX_DIR")" && pwd)/$(basename "$MATRIX_DIR")

# Verify matrix directory exists
if [[ ! -d "$MATRIX_DIR" ]]; then
    echo "Error: Matrix directory does not exist: $MATRIX_DIR"
    exit 1
fi

# Create output directory if it doesn't exist
OUTPUT_DIR=$(dirname "$OUTPUT_PREFIX")
mkdir -p "$OUTPUT_DIR"
OUTPUT_PREFIX=$(cd "$OUTPUT_DIR" && pwd)/$(basename "$OUTPUT_PREFIX")

# Determine the base directory to mount (parent of matrix dir)
MATRIX_BASE_DIR=$(dirname "$MATRIX_DIR")
MATRIX_SUBDIR=$(basename "$MATRIX_DIR")

# Output prefix relative to the mounted directory; it must lie inside it
OUTPUT_PREFIX_REL=$(realpath -m --relative-to="$MATRIX_BASE_DIR" "$OUTPUT_PREFIX")
if [[ "$OUTPUT_PREFIX_REL" == ../* ]]; then
    echo "Error: --output-prefix must be under $MATRIX_BASE_DIR, the directory mounted into the container"
    exit 1
fi

# Build the Docker command as an array (no eval). '-e CAS_API_TOKEN' with no value
# passes the variable from this environment, so the token is not on the command line.
# The matrix base directory is mounted read-write so .for_CAS can be created.
DOCKER_CMD=(docker run --rm
    -e CAS_API_TOKEN
    -v "${MATRIX_BASE_DIR}:/data"
    "${DOCKER_IMAGE}"
    --matrix-dir "/data/${MATRIX_SUBDIR}"
    --output-prefix "/data/${OUTPUT_PREFIX_REL}"
    --chunk-size "${CHUNK_SIZE}"
    --min-acceptable-score "${MIN_ACCEPTABLE_SCORE}"
    --top-k "${TOP_K}"
    --obs-prefix "${OBS_PREFIX}")

if [[ -n "$CAS_MODEL_NAME" ]]; then
    DOCKER_CMD+=(--cas-model-name "${CAS_MODEL_NAME}")
fi
if [[ -n "$MIN_COUNTS" ]]; then
    DOCKER_CMD+=(--min-counts "${MIN_COUNTS}")
fi

# Print the command for debugging (it holds no token)
echo "Running Docker command (token passed as CAS_API_TOKEN):"
printf '%q ' "${DOCKER_CMD[@]}"
echo
echo ""

"${DOCKER_CMD[@]}"

echo ""
echo "Done! Output files:"
echo "  ${OUTPUT_PREFIX}.h5ad"
echo "  ${OUTPUT_PREFIX}.tsv"
