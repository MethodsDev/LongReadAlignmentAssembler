#!/usr/bin/env python3
import argparse, os, sys, gzip, logging
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite
from collections import defaultdict
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

def stream_group_counts(filename, feature_col, chunksize=1_000_000):
    """
    Stream the tracking file in chunks and aggregate frac_assigned
    by (feature_id, cell_barcode), with progress monitoring.
    Works for .gz or plain text.
    """
    agg = defaultdict(float)

    is_gz = filename.endswith(".gz")
    opener = gzip.open if is_gz else open

    # count lines for progress estimation
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
            engine="c",  # enforce the C engine; pyarrow backend has triggered segfaults in Terra runtime
        )
        for i, chunk in enumerate(reader, 1):
            rows_processed += len(chunk)
            pct = (rows_processed / total_data_lines) * 100
            logger.info("[chunk %s] processed ~%.1f%% of rows (RSS %s)",
                        i, pct, format_rss())

            # split read_name into cell_barcode
            split = chunk["read_name"].str.split("^", n=2, expand=True)
            chunk["cell_barcode"] = split[0]

            grouped = chunk.groupby([feature_col, "cell_barcode"], observed=True)["frac_assigned"].sum()
            for (feat, bc), val in grouped.items():
                agg[(feat, bc)] += val

    logger.info("finished 100.0%% of rows (RSS %s)", format_rss())

    # convert dict → dataframe
    feature_ids, barcodes, counts = zip(*((f, c, v) for (f, c), v in agg.items()))
    return pd.DataFrame({"feature_id": feature_ids,
                         "cell_barcode": barcodes,
                         "UMI_counts": counts})

def make_sparse_matrix_outputs(counts_data, outdirname):
    logger.info("making sparse matrix outputs for: %s (RSS %s)",
                outdirname, format_rss())
    os.makedirs(outdirname, exist_ok=True)

    features = counts_data["feature_id"].astype("category")
    barcodes = counts_data["cell_barcode"].astype("category")

    rows = features.cat.codes.to_numpy()
    cols = barcodes.cat.codes.to_numpy()
    vals = counts_data["UMI_counts"].to_numpy()

    sparseM = sparse.coo_matrix((vals, (rows, cols)),
                                shape=(features.cat.categories.size,
                                       barcodes.cat.categories.size))

    mmwrite(os.path.join(outdirname, "matrix.mtx"), sparseM)

    # Write single-column features.tsv and barcodes.tsv (first column used by Read10X)
    features.cat.categories.to_series().to_csv(
        os.path.join(outdirname, "features.tsv"), index=False, header=False
    )
    barcodes.cat.categories.to_series().to_csv(
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
    args = parser.parse_args()

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

    for label, col in [("gene", "gene_id"),
                       ("isoform", "transcript_id"),
                       ("splice_pattern", "transcript_splice_hash_code")]:
        logger.info("processing %s level counts (RSS %s)", label, format_rss())
        counts = stream_group_counts(args.tracking, col, chunksize=args.chunksize)
        counts.to_csv(f"{args.output_prefix}.{label}_cell_counts.tsv", sep="\t", index=False)
        make_sparse_matrix_outputs(counts, f"{args.output_prefix}^{label}-sparseM")

    logger.info("all done (RSS %s)", format_rss())

if __name__ == "__main__":
    main()

