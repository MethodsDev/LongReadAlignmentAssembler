#!/bin/bash

set -e

# Docker wrapper script for run_CellAnnotationService.py
# This script runs the Cellarium CAS cell type annotation tool via Docker

usage() {
    cat <<EOF
Usage: $0 --matrix-dir <path> --api-token <token> --output-prefix <prefix> [options]

Required arguments:
  --matrix-dir <path>          Path to 10x-style matrix directory containing:
                               matrix.mtx[.gz], barcodes.tsv[.gz], features.tsv[.gz]
  --api-token <token>          Cellarium CAS API token
  --output-prefix <prefix>     Prefix for output files (will create <prefix>.h5ad and <prefix>.tsv)

Optional arguments:
  --cas-model-name <name>      CAS model name (default: use CAS default)
  --chunk-size <int>           Number of cells per CAS chunk (default: 500)
  --min-acceptable-score <f>   Minimum acceptable evidence score (default: 0.2)
  --top-k <int>                Top-k CAS cell type calls per cell (default: 3)
  --obs-prefix <string>        Prefix for CAS cell type columns (default: cas_cell_type)
  --docker-image <image>       Docker image to use (default: us-central1-docker.pkg.dev/methods-dev-lab/lraa/cas)
  --help                       Show this help message

Example:
  $0 \\
    --matrix-dir /path/to/matrix_folder \\
    --api-token YOUR_CAS_TOKEN \\
    --output-prefix /path/to/output/sample1

This will create:
  /path/to/output/sample1.h5ad
  /path/to/output/sample1.tsv

EOF
    exit 1
}

# Default values
DOCKER_IMAGE="us-central1-docker.pkg.dev/methods-dev-lab/lraa/cas"
MATRIX_DIR=""
API_TOKEN=""
OUTPUT_PREFIX=""
CAS_MODEL_NAME=""
CHUNK_SIZE="500"
MIN_ACCEPTABLE_SCORE="0.2"
TOP_K="3"
OBS_PREFIX="cas_cell_type"

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --matrix-dir)
            MATRIX_DIR="$2"
            shift 2
            ;;
        --api-token)
            API_TOKEN="$2"
            shift 2
            ;;
        --output-prefix)
            OUTPUT_PREFIX="$2"
            shift 2
            ;;
        --cas-model-name)
            CAS_MODEL_NAME="$2"
            shift 2
            ;;
        --chunk-size)
            CHUNK_SIZE="$2"
            shift 2
            ;;
        --min-acceptable-score)
            MIN_ACCEPTABLE_SCORE="$2"
            shift 2
            ;;
        --top-k)
            TOP_K="$2"
            shift 2
            ;;
        --obs-prefix)
            OBS_PREFIX="$2"
            shift 2
            ;;
        --docker-image)
            DOCKER_IMAGE="$2"
            shift 2
            ;;
        --help)
            usage
            ;;
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

if [[ -z "$API_TOKEN" ]]; then
    echo "Error: --api-token is required"
    usage
fi

if [[ -z "$OUTPUT_PREFIX" ]]; then
    echo "Error: --output-prefix is required"
    usage
fi

# Convert paths to absolute paths
MATRIX_DIR=$(cd "$(dirname "$MATRIX_DIR")" && pwd)/$(basename "$MATRIX_DIR")
OUTPUT_DIR=$(dirname "$OUTPUT_PREFIX")
OUTPUT_DIR=$(cd "$OUTPUT_DIR" && pwd)
OUTPUT_BASE=$(basename "$OUTPUT_PREFIX")

# Verify matrix directory exists
if [[ ! -d "$MATRIX_DIR" ]]; then
    echo "Error: Matrix directory does not exist: $MATRIX_DIR"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Build Docker command
DOCKER_CMD="docker run --rm"

# Mount the matrix directory (read-only)
DOCKER_CMD="$DOCKER_CMD -v ${MATRIX_DIR}:/data/matrix:ro"

# Mount the output directory (read-write)
DOCKER_CMD="$DOCKER_CMD -v ${OUTPUT_DIR}:/data/output"

# Add the Docker image
DOCKER_CMD="$DOCKER_CMD ${DOCKER_IMAGE}"

# Add required arguments (ENTRYPOINT handles the script invocation)
DOCKER_CMD="$DOCKER_CMD --matrix-dir /data/matrix"
DOCKER_CMD="$DOCKER_CMD --api-token ${API_TOKEN}"
DOCKER_CMD="$DOCKER_CMD --output-prefix /data/output/${OUTPUT_BASE}"

# Add optional arguments
if [[ -n "$CAS_MODEL_NAME" ]]; then
    DOCKER_CMD="$DOCKER_CMD --cas-model-name ${CAS_MODEL_NAME}"
fi

DOCKER_CMD="$DOCKER_CMD --chunk-size ${CHUNK_SIZE}"
DOCKER_CMD="$DOCKER_CMD --min-acceptable-score ${MIN_ACCEPTABLE_SCORE}"
DOCKER_CMD="$DOCKER_CMD --top-k ${TOP_K}"
DOCKER_CMD="$DOCKER_CMD --obs-prefix ${OBS_PREFIX}"

# Print the command for debugging
echo "Running Docker command:"
echo "$DOCKER_CMD"
echo ""

# Execute the Docker command
eval $DOCKER_CMD

echo ""
echo "Done! Output files:"
echo "  ${OUTPUT_DIR}/${OUTPUT_BASE}.h5ad"
echo "  ${OUTPUT_DIR}/${OUTPUT_BASE}.tsv"
