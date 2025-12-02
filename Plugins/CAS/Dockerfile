# Use a lightweight Python image (cellarium-cas supports 3.7–3.12; 3.10+ recommended)
FROM python:3.10-slim

# Ensures Python output is not buffered and pip doesn't cache packages
ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

# Install minimal system dependencies (hdf5 for AnnData/h5ad IO, build tools for any compiled deps)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Python libraries:
# - cellarium-cas: CAS client library
# - scanpy: for read_10x_mtx and AnnData handling
RUN pip install --no-cache-dir \
    cellarium-cas \
    scanpy

# Create a working directory (where you'll typically mount your data)
WORKDIR /work

# Copy your CAS script into /usr/local/bin and make it executable
# Make sure this Dockerfile is in the same directory as run_CellAnnotationService.py
COPY run_CellAnnotationService.py /usr/local/bin/run_CellAnnotationService.py
RUN chmod +x /usr/local/bin/run_CellAnnotationService.py

# By default, run the script; additional arguments can be passed via `docker run`
ENTRYPOINT ["run_CellAnnotationService.py"]
