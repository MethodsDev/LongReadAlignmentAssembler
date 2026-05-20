#!/usr/bin/env python3
# Note: Explicitly specifying dtypes in pd.read_csv() is critical to avoid
# segmentation faults during type inference on large chunked datasets.
import argparse, os, sys, gzip, logging, faulthandler, shutil, subprocess
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

# Global peak memory tracker
peak_rss_bytes = 0


def _rss_bytes():
    global peak_rss_bytes
    if psutil is not None:
        current_rss = psutil.Process(os.getpid()).memory_info().rss
    elif resource is not None:
        usage = resource.getrusage(resource.RUSAGE_SELF)
        rss = usage.ru_maxrss
        if sys.platform == "darwin":
            current_rss = rss
        else:
            current_rss = rss * 1024
    else:
        current_rss = 0
    
    if current_rss > peak_rss_bytes:
        peak_rss_bytes = current_rss
    
    return current_rss


def format_rss():
    rss = _rss_bytes()
    if rss <= 0:
        return "unknown"
    return f"{rss / (1024 ** 2):.1f} MiB"

def format_peak_rss():
    if peak_rss_bytes <= 0:
        return "unknown"
    return f"{peak_rss_bytes / (1024 ** 2):.1f} MiB"


def compress_output_file(path, gzip_level=6):
    gz_path = path + ".gz"
    pigz_path = shutil.which("pigz")

    if pigz_path is not None:
        logger.info(
            "compressing %s -> %s with pigz -%d (RSS %s)",
            path,
            gz_path,
            gzip_level,
            format_rss(),
        )
        subprocess.run([pigz_path, f"-{gzip_level}", "-f", path], check=True)
        logger.info("finished compressing %s with pigz (RSS %s)", gz_path, format_rss())
        return

    logger.info(
        "compressing %s -> %s with python gzip level %d (RSS %s)",
        path,
        gz_path,
        gzip_level,
        format_rss(),
    )
    with open(path, "rb") as f_in, gzip.open(gz_path, "wb", compresslevel=gzip_level) as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(path)
    logger.info("finished compressing %s with python gzip (RSS %s)", gz_path, format_rss())


class LevelAccumulator:
    """Collect sparse matrix coordinates for one feature level across chunks."""

    def __init__(self, label, feature_col):
        self.label = label
        self.feature_col = feature_col
        self.feature_to_index = {}
        self.feature_labels = []
        self.row_chunks = []
        self.col_chunks = []
        self.val_chunks = []

    def consume_chunk(self, chunk, barcode_idx):
        feature_values = chunk[self.feature_col]
        for feat in pd.unique(feature_values):
            if feat not in self.feature_to_index:
                self.feature_to_index[feat] = len(self.feature_labels)
                self.feature_labels.append(feat)

        feature_idx = feature_values.map(self.feature_to_index)
        grouped = (
            pd.DataFrame(
                {
                    "feature_idx": feature_idx,
                    "barcode_idx": barcode_idx,
                    "frac_assigned": chunk["frac_assigned"],
                }
            )
            .groupby(["feature_idx", "barcode_idx"], observed=True)["frac_assigned"]
            .sum()
        )

        if not grouped.empty:
            self.row_chunks.append(grouped.index.get_level_values(0).to_numpy(dtype=np.int64))
            self.col_chunks.append(grouped.index.get_level_values(1).to_numpy(dtype=np.int64))
            self.val_chunks.append(grouped.to_numpy(dtype=np.float32))

    def finalize(self, barcode_labels):
        if self.row_chunks:
            rows = np.concatenate(self.row_chunks)
            cols = np.concatenate(self.col_chunks)
            vals = np.concatenate(self.val_chunks)
        else:
            rows = np.array([], dtype=np.int64)
            cols = np.array([], dtype=np.int64)
            vals = np.array([], dtype=np.float32)

        matrix = sparse.coo_matrix(
            (vals, (rows, cols)),
            shape=(len(self.feature_labels), len(barcode_labels)),
        )
        logger.info(
            "%s: created sparse matrix with %d features x %d barcodes, %d non-zero entries (RSS %s)",
            self.label,
            len(self.feature_labels),
            len(barcode_labels),
            matrix.nnz,
            format_rss(),
        )

        if matrix.nnz:
            matrix.sum_duplicates()

        logger.info("%s: converting to CSR and sorting (RSS %s)", self.label, format_rss())
        matrix = matrix.tocsr()

        feature_arr = np.array(self.feature_labels, dtype=object)
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

        logger.info("%s: creating counts dataframe (RSS %s)", self.label, format_rss())
        if matrix.nnz:
            counts_df = pd.DataFrame(
                {
                    "feature_id": feature_arr[matrix.row],
                    "cell_barcode": barcode_arr[matrix.col],
                    "UMI_counts": matrix.data.astype(float),
                }
            )
        else:
            counts_df = pd.DataFrame(columns=["feature_id", "cell_barcode", "UMI_counts"])

        logger.info(
            "%s: final counts dataframe has %d rows (RSS %s)",
            self.label,
            len(counts_df),
            format_rss(),
        )
        return counts_df, matrix, feature_arr, barcode_arr

def stream_all_counts(filename, chunksize=1_000_000, engine="python"):
    """
    Stream the tracking file once and aggregate mapping, gene, isoform,
    and splice-pattern outputs together.
    """
    logger.info("starting single-pass aggregation from %s (RSS %s)", filename, format_rss())

    levels = [
        LevelAccumulator("gene", "gene_id"),
        LevelAccumulator("isoform", "transcript_id"),
        LevelAccumulator("splice_pattern", "transcript_splice_hash_code"),
    ]

    mapping_columns = ["gene_id", "transcript_id", "transcript_splice_hash_code"]
    mapping_entries = set()
    barcode_to_index = {}
    barcode_labels = []
    rows_processed = 0

    opener = gzip.open if filename.endswith(".gz") else open
    dtype_spec = {
        "read_name": str,
        "gene_id": str,
        "transcript_id": str,
        "transcript_splice_hash_code": str,
        "frac_assigned": np.float32,
    }

    with opener(filename, "rt") as handle:
        read_csv_kwargs = {
            "sep": "\t",
            "chunksize": chunksize,
            "usecols": ["read_name", "gene_id", "transcript_id", "transcript_splice_hash_code", "frac_assigned"],
            "dtype": dtype_spec,
            "engine": engine,
            "comment": "#",
        }
        if engine == "c":
            read_csv_kwargs["low_memory"] = False

        for chunk_num, chunk in enumerate(pd.read_csv(handle, **read_csv_kwargs), 1):
            rows_processed += len(chunk)
            logger.info("[chunk %d] processed %d rows (RSS %s)", chunk_num, rows_processed, format_rss())

            chunk = chunk.fillna(
                {
                    "gene_id": "",
                    "transcript_id": "",
                    "transcript_splice_hash_code": "",
                }
            )

            for entry in chunk[mapping_columns].itertuples(index=False, name=None):
                mapping_entries.add(entry)

            cell_barcodes = chunk["read_name"].str.split("^", n=2, expand=True)[0]
            for bc in pd.unique(cell_barcodes):
                if bc not in barcode_to_index:
                    barcode_to_index[bc] = len(barcode_labels)
                    barcode_labels.append(bc)
            barcode_idx = cell_barcodes.map(barcode_to_index)

            for level in levels:
                level.consume_chunk(chunk, barcode_idx)

    logger.info(
        "finished single-pass aggregation after %d rows, finalizing outputs (RSS %s)",
        rows_processed,
        format_rss(),
    )

    if mapping_entries:
        mapping_df = pd.DataFrame(sorted(mapping_entries), columns=mapping_columns)
    else:
        mapping_df = pd.DataFrame(columns=mapping_columns)

    results = {}
    for level in levels:
        results[level.label] = level.finalize(barcode_labels)

    return mapping_df, results

def make_sparse_matrix_outputs(matrix, feature_labels, barcode_labels, outdirname, gzip_level=6):
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

    # Compress outputs after writing them uncompressed to disk.
    for fn in ["matrix.mtx", "features.tsv", "barcodes.tsv"]:
        path = os.path.join(outdirname, fn)
        compress_output_file(path, gzip_level=gzip_level)
    logger.info("done with %s (RSS %s)", outdirname, format_rss())

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tracking", required=True, help="input data file quant.tracking[.gz]")
    parser.add_argument("--output_prefix", required=True, help="output prefix")
    parser.add_argument("--chunksize", type=int, default=1_000_000,
                        help="rows per chunk (default 1e6)")
    parser.add_argument("--csv_engine", choices=["python", "c"], default="python",
                        help="pandas CSV engine to use (default python; C engine has segfaulted in production runs)")
    parser.add_argument("--parallel", action="store_true",
                        help="process gene/isoform/splice_pattern levels in parallel (requires more RAM)")
    parser.add_argument("--gzip_level", type=int, choices=range(1, 10), default=1,
                        help="gzip compression level for output files (default 1 for faster compression)")
    args = parser.parse_args()

    faulthandler.enable()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        stream=sys.stderr,
    )

    logger.info("=" * 60)
    logger.info("Starting single-cell tracking to sparse matrix conversion")
    logger.info("Input: %s", args.tracking)
    logger.info("Output prefix: %s", args.output_prefix)
    logger.info("Chunksize: %s rows", args.chunksize)
    logger.info("CSV engine: %s", args.csv_engine)
    logger.info("Parallel mode: %s", args.parallel)
    logger.info("Output gzip level: %s", args.gzip_level)
    logger.info("Initial RSS: %s", format_rss())
    logger.info("=" * 60)

    if args.parallel:
        logger.warning("--parallel is ignored by the single-pass refactor")

    output_parent = os.path.dirname(args.output_prefix)
    if output_parent:
        os.makedirs(output_parent, exist_ok=True)

    mapping_df, level_results = stream_all_counts(
        args.tracking,
        chunksize=args.chunksize,
        engine=args.csv_engine,
    )

    mapping_output = f"{args.output_prefix}.gene_transcript_splicehashcode.tsv"
    mapping_df.to_csv(mapping_output, sep="\t", index=False)
    logger.info(
        "wrote mapping file: %s with %d entries (RSS %s)",
        mapping_output,
        len(mapping_df),
        format_rss(),
    )

    for label, (counts, matrix, feature_labels, barcode_labels) in level_results.items():
        counts.to_csv(f"{args.output_prefix}.{label}_cell_counts.tsv", sep="\t", index=False)
        make_sparse_matrix_outputs(
            matrix,
            feature_labels,
            barcode_labels,
            f"{args.output_prefix}^{label}-sparseM",
            gzip_level=args.gzip_level,
        )

    logger.info("all done (current RSS %s, peak RSS %s)", format_rss(), format_peak_rss())

if __name__ == "__main__":
    main()
