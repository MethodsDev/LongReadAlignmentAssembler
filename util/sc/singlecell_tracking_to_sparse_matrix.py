#!/usr/bin/env python3
import argparse, os, sys, gzip, logging, faulthandler
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite
import numpy as np
try:
    import psutil
except ImportError:  # pragma: no cover - psutil might be missing in minimal envs
    psutil = None
try:
    import resource
except ImportError:  # pragma: no cover - resource is UNIX-only
    resource = None

logger = logging.getLogger(__name__)


def _rss_bytes():
    if psutil is not None:
        return psutil.Process(os.getpid()).memory_info().rss
    if resource is not None:
        usage = resource.getrusage(resource.RUSAGE_SELF)
        rss = usage.ru_maxrss
        if sys.platform == "darwin":
            return rss
        return rss * 1024
    return 0


def format_rss():
    rss = _rss_bytes()
    if rss <= 0:
        return "unknown"
    return f"{rss / (1024 ** 2):.1f} MiB"

def count_lines(filename):
    """Count lines in file (for progress estimation)."""
    opener = gzip.open if filename.endswith(".gz") else open
    with opener(filename, "rt") as f:
        for i, _ in enumerate(f, 1):
            pass
    return i

def stream_group_counts(filename, feature_col, chunksize=1_000_000, engine="python"):
    """
    Stream the tracking file in chunks and aggregate frac_assigned
    by (feature_id, cell_barcode), with progress monitoring.
    Works for .gz or plain text.
    Returns the aggregated COO matrix plus ordered feature and barcode labels.
    """

    feature_to_index = {}
    feature_labels = []
    barcode_to_index = {}
    barcode_labels = []

    row_chunks = []
    col_chunks = []
    val_chunks = []

    is_gz = filename.endswith(".gz")
    opener = gzip.open if is_gz else open

    total_lines = count_lines(filename)
    header_lines = 1
    total_data_lines = max(1, total_lines - header_lines)
    rows_processed = 0

    with opener(filename, "rt") as f:
        reader = pd.read_csv(
            f,
            sep="\t",
            chunksize=chunksize,
            usecols=["read_name", feature_col, "frac_assigned"],
            engine=engine,
        )
        for i, chunk in enumerate(reader, 1):
            rows_processed += len(chunk)
            pct = (rows_processed / total_data_lines) * 100
            logger.info("[chunk %s] processed ~%.1f%% of rows (RSS %s)",
                        i, pct, format_rss())

            split = chunk["read_name"].str.split("^", n=2, expand=True)
            chunk["cell_barcode"] = split[0]

            feature_values = chunk[feature_col]
            for feat in pd.unique(feature_values):
                if feat not in feature_to_index:
                    feature_to_index[feat] = len(feature_labels)
                    feature_labels.append(feat)
            chunk["feature_idx"] = feature_values.map(feature_to_index)

            barcode_values = chunk["cell_barcode"]
            for bc in pd.unique(barcode_values):
                if bc not in barcode_to_index:
                    barcode_to_index[bc] = len(barcode_labels)
                    barcode_labels.append(bc)
            chunk["barcode_idx"] = barcode_values.map(barcode_to_index)

            grouped = chunk.groupby(["feature_idx", "barcode_idx"], observed=True)["frac_assigned"].sum()
            if not grouped.empty:
                row_chunks.append(grouped.index.get_level_values(0).to_numpy(dtype=np.int64))
                col_chunks.append(grouped.index.get_level_values(1).to_numpy(dtype=np.int64))
                val_chunks.append(grouped.to_numpy(dtype=np.float32))

            del chunk

    logger.info("finished 100.0%% of rows (RSS %s)", format_rss())

    if row_chunks:
        rows = np.concatenate(row_chunks)
        cols = np.concatenate(col_chunks)
        vals = np.concatenate(val_chunks)
    else:
        rows = np.array([], dtype=np.int64)
        cols = np.array([], dtype=np.int64)
        vals = np.array([], dtype=np.float32)

    matrix = sparse.coo_matrix(
        (vals, (rows, cols)),
        shape=(len(feature_labels), len(barcode_labels)),
    )

    if matrix.nnz:
        matrix.sum_duplicates()

    matrix = matrix.tocsr()

    feature_arr = np.array(feature_labels, dtype=object)
    barcode_arr = np.array(barcode_labels, dtype=object)

    if feature_arr.size:
        feature_order = np.argsort(feature_arr)
        matrix = matrix[feature_order]
        feature_arr = feature_arr[feature_order]

    if barcode_arr.size:
        barcode_order = np.argsort(barcode_arr)
        matrix = matrix[:, barcode_order]
        barcode_arr = barcode_arr[barcode_order]

    matrix = matrix.tocoo()

    if matrix.nnz:
        counts_df = pd.DataFrame({
            "feature_id": feature_arr[matrix.row],
            "cell_barcode": barcode_arr[matrix.col],
            "UMI_counts": matrix.data.astype(float),
        })
    else:
        counts_df = pd.DataFrame(columns=["feature_id", "cell_barcode", "UMI_counts"])

    return counts_df, matrix, feature_arr, barcode_arr

def make_sparse_matrix_outputs(matrix, feature_labels, barcode_labels, outdirname):
    logger.info("making sparse matrix outputs for: %s (RSS %s)",
                outdirname, format_rss())
    os.makedirs(outdirname, exist_ok=True)

    sparseM = matrix.tocoo()

    mmwrite(os.path.join(outdirname, "matrix.mtx"), sparseM)

    pd.Series(feature_labels).to_csv(
        os.path.join(outdirname, "features.tsv"), index=False, header=False
    )
    pd.Series(barcode_labels).to_csv(
        os.path.join(outdirname, "barcodes.tsv"), index=False, header=False
    )

    # gzip outputs
    for fn in ["matrix.mtx", "features.tsv", "barcodes.tsv"]:
        path = os.path.join(outdirname, fn)
        with open(path, "rb") as f_in, gzip.open(path + ".gz", "wb") as f_out:
            f_out.writelines(f_in)
        os.remove(path)
    logger.info("done with %s (RSS %s)", outdirname, format_rss())

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tracking", required=True, help="input data file quant.tracking[.gz]")
    parser.add_argument("--output_prefix", required=True, help="output prefix")
    parser.add_argument("--chunksize", type=int, default=1_000_000,
                        help="rows per chunk (default 1e6)")
    parser.add_argument("--csv_engine", choices=["python", "c"], default="python",
                        help="pandas CSV engine to use (default python to avoid C-engine segfaults)")
    args = parser.parse_args()

    faulthandler.enable()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        stream=sys.stderr,
    )

    # mapping file (small)
    opener = gzip.open if args.tracking.endswith(".gz") else open
    ids = pd.read_csv(opener(args.tracking, "rt"), sep="\t",
                      usecols=["gene_id","transcript_id","transcript_splice_hash_code"])
    ids.drop_duplicates().to_csv(f"{args.output_prefix}.gene_transcript_splicehashcode.tsv",
                                 sep="\t", index=False)

    for label, col in [
        ("splice_pattern", "transcript_splice_hash_code"),
        ("gene", "gene_id"),
        ("isoform", "transcript_id")
                       ]:
        logger.info("processing %s level counts (RSS %s)", label, format_rss())
        counts, matrix, feature_labels, barcode_labels = stream_group_counts(
            args.tracking, col, chunksize=args.chunksize, engine=args.csv_engine
        )
        counts.to_csv(f"{args.output_prefix}.{label}_cell_counts.tsv", sep="\t", index=False)
        make_sparse_matrix_outputs(
            matrix,
            feature_labels,
            barcode_labels,
            f"{args.output_prefix}^{label}-sparseM",
        )

    logger.info("all done (RSS %s)", format_rss())

if __name__ == "__main__":
    main()

