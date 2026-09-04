#!/usr/bin/env python3
"""Shared on-disk format for per-shard single-cell sparse matrices.

A shard is a (cluster, contig) partition of the library.  Every (feature, cell)
pair belongs to exactly one shard, because feature -> contig and cell -> cluster
are both functions, so the merge in merge_shard_sparse_matrices.py never has to
sum across shards.  See PLAN.sc_sparse_from_shards.md section 0.

Layout of one shard directory:

    barcodes.tsv            this shard's cell barcodes, sorted
    mapping.tsv             gene_id, transcript_id, splice hash, num_exons
    <level>.features.tsv    this level's feature names, sorted
    <level>.csr             the matrix, CSR, rows in features.tsv order and
                            columns in barcodes.tsv order
    manifest.json           shard name, per-level shape and nnz

The .csr file is deliberately not MatrixMarket: 8 bytes per non-zero against
about 30 for text, which is 15 GB of intermediate across a whole run rather
than 56 GB.
"""
import json, os
import numpy as np

MAGIC = b"LRAASPM1"
LEVELS = ("gene", "isoform", "splice_pattern")


def write_csr(path, indptr, indices, data):
    """indptr int64[n_rows+1], indices int32[nnz], data float32[nnz]."""
    indptr = np.ascontiguousarray(indptr, dtype=np.int64)
    indices = np.ascontiguousarray(indices, dtype=np.int32)
    data = np.ascontiguousarray(data, dtype=np.float32)
    with open(path, "wb") as fh:
        fh.write(MAGIC)
        fh.write(np.array([len(indptr) - 1, 0, len(indices)], dtype=np.int64).tobytes())
        fh.write(indptr.tobytes())
        fh.write(indices.tobytes())
        fh.write(data.tobytes())


def write_csr_with_shape(path, matrix):
    """Write a scipy CSR, recording its column count in the header."""
    with open(path, "wb") as fh:
        fh.write(MAGIC)
        fh.write(np.array([matrix.shape[0], matrix.shape[1], matrix.nnz],
                          dtype=np.int64).tobytes())
        fh.write(np.ascontiguousarray(matrix.indptr, dtype=np.int64).tobytes())
        fh.write(np.ascontiguousarray(matrix.indices, dtype=np.int32).tobytes())
        fh.write(np.ascontiguousarray(matrix.data, dtype=np.float32).tobytes())


_HDR = len(MAGIC) + 3 * 8


def read_csr_header(path):
    with open(path, "rb") as fh:
        magic = fh.read(len(MAGIC))
        if magic != MAGIC:
            raise ValueError(f"{path}: not a shard CSR file (magic {magic!r})")
        n_rows, n_cols, nnz = np.frombuffer(fh.read(24), dtype=np.int64)
    return int(n_rows), int(n_cols), int(nnz)


def open_csr(path):
    """Return (n_rows, n_cols, nnz, indptr, indices_memmap, data_memmap).

    indices and data are memmapped so the merge reads them sequentially without
    holding a shard's matrix in memory.
    """
    n_rows, n_cols, nnz = read_csr_header(path)
    off = _HDR
    indptr = np.memmap(path, dtype=np.int64, mode="r", offset=off, shape=(n_rows + 1,))
    off += (n_rows + 1) * 8
    indices = np.memmap(path, dtype=np.int32, mode="r", offset=off, shape=(nnz,))
    off += nnz * 4
    data = np.memmap(path, dtype=np.float32, mode="r", offset=off, shape=(nnz,))
    return n_rows, n_cols, nnz, np.asarray(indptr), indices, data


def read_lines(path):
    with open(path) as fh:
        return [line.rstrip("\n") for line in fh]


def write_lines(path, lines):
    with open(path, "w") as fh:
        for x in lines:
            fh.write(x)
            fh.write("\n")


def write_manifest(outdir, shard, levels_info, n_barcodes):
    with open(os.path.join(outdir, "manifest.json"), "w") as fh:
        json.dump({"shard": shard, "n_barcodes": n_barcodes,
                   "levels": levels_info, "format": MAGIC.decode()}, fh, indent=2)


def read_manifest(outdir):
    with open(os.path.join(outdir, "manifest.json")) as fh:
        return json.load(fh)
