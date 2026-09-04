#!/usr/bin/env python3
# The default reader ("direct") uses no CSV parser at all: the tracking file is
# machine-generated TSV with no quoting, so a plain split() loop is enough, and
# on the VILLAGE 1.5 B-read library it measured 4.1x faster than pandas' python
# engine and 1.7x faster than its C engine -- while having no exposure to the
# segfaults that made the C engine unusable here.
#
# --csv_engine python and c still select the pandas readers.  Explicitly
# specifying dtypes there is critical to avoid segmentation faults during type
# inference on large chunked datasets.
import argparse, os, sys, gzip, gc, logging, faulthandler, shutil, subprocess
from sys import intern
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

_INT32_MAX = 2 ** 31

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
        # factorize once, then map only the distinct values through the dict:
        # O(unique) python lookups instead of one per row
        codes, uniques = pd.factorize(feature_values.to_numpy())
        lut = np.empty(len(uniques), dtype=np.int64)
        for k, feat in enumerate(uniques):
            idx = self.feature_to_index.get(feat)
            if idx is None:
                idx = self.feature_to_index[feat] = len(self.feature_labels)
                self.feature_labels.append(feat)
            lut[k] = idx
        if len(self.feature_labels) >= _INT32_MAX:
            raise SystemExit(f"{self.label}: more than 2^31 features; widen the "
                             "accumulator index dtype")
        feature_idx = pd.Series(lut[codes], index=chunk.index)
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
            # int32 rather than int64: 12 bytes per accumulated entry instead of
            # 20.  On the VILLAGE library that is 37.5 GiB -> 22.5 GiB resident
            # through the whole read, and since throughput falls as the resident
            # set grows it buys time as well as memory.
            self.row_chunks.append(grouped.index.get_level_values(0).to_numpy(dtype=np.int32))
            self.col_chunks.append(grouped.index.get_level_values(1).to_numpy(dtype=np.int32))
            self.val_chunks.append(grouped.to_numpy(dtype=np.float32))

    def finalize(self, barcode_labels, build_counts_df=False):
        # release each list as soon as it is concatenated; holding all three
        # levels' chunk lists to the end of the process was 37.5 GiB of the old
        # peak, and concatenating without releasing doubles the level in flight
        if self.row_chunks:
            rows = np.concatenate(self.row_chunks); self.row_chunks = None
            cols = np.concatenate(self.col_chunks); self.col_chunks = None
            vals = np.concatenate(self.val_chunks); self.val_chunks = None
        else:
            rows = np.array([], dtype=np.int32)
            cols = np.array([], dtype=np.int32)
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

        # No explicit sum_duplicates(): tocsr() sums duplicates itself, and the
        # explicit call costs a full lexsort over every non-zero -- measured at
        # 4.24 GiB peak against 1.91 GiB per 100 M entries, for no change to the
        # output.  On this library it collapsed 0.039% of entries.
        logger.info("%s: converting to CSR and sorting (RSS %s)", self.label, format_rss())
        matrix = matrix.tocsr()
        del rows, cols, vals

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

        # The dense counts dataframe costs 22.4 bytes per non-zero against the
        # sparse accumulator's 11.2 -- twice the memory of the data it
        # duplicates -- and its only consumer is the *_cell_counts.tsv write.
        # Those files are declared task outputs in the WDLs but no workflow
        # output block references them, so they were delocalized and dropped.
        # Off unless --emit_dense_counts is passed.
        counts_df = None
        if build_counts_df:
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
            logger.info("%s: counts dataframe has %d rows (RSS %s)", self.label,
                        len(counts_df), format_rss())

        return counts_df, matrix, feature_arr, barcode_arr


WEIGHTED_TRACKING_MARKER = "use_XW_read_weights_for_quant"
WEIGHTED_TRACKING_MARKER_LINE_PREFIX = "# WARNING:"


def _reject_if_weighted_tracking(filename, opener):
    """Refuse a tracking file whose rows cover only normalization-retained reads.

    Counts here are sums of frac_assigned, one row per read, and XW is never consulted.
    A tracking file produced with XW weighting holds rows for the reads coverage
    normalization kept and none for the reads it discarded, so every cell would come out
    short -- unevenly, and worst at the high-coverage loci that weighting exists to
    correct. That is a wrong matrix with nothing about it to look wrong.

    LRAA marks such files in their leading comment block, and this reads it before pandas
    does, because the parser is configured to treat '#' as a comment and would drop the
    warning silently.

    Matched only against the dedicated "# WARNING:"-prefixed marker line, not any
    leading comment containing the substring: every tracking file also carries a
    "# LRAA CMD: ..." echo of its own command line, which contains this same substring
    whenever --use_XW_read_weights_for_quant was passed at all -- including a
    --stream_reads run, whose merged file is complete and does not carry the marker
    line. A blind substring scan would refuse that file too, on the strength of its own
    command-line echo rather than anything describing incompleteness.
    """
    try:
        with opener(filename, "rt") as handle:
            for line in handle:
                if not line.startswith("#"):
                    break
                if (
                    line.startswith(WEIGHTED_TRACKING_MARKER_LINE_PREFIX)
                    and WEIGHTED_TRACKING_MARKER in line
                ):
                    sys.exit(
                        "Error, {} was produced with --{}.\n"
                        "Its rows cover only the reads retained by coverage "
                        "normalization, and frac_assigned is not weighted by XW, so "
                        "per-cell counts summed from it understate every cell -- most at "
                        "high-coverage loci. Re-run quantification without that flag to "
                        "build single-cell matrices.".format(
                            filename, WEIGHTED_TRACKING_MARKER
                        )
                    )
    except OSError:
        # unreadable input is the caller's problem to report, not this check's
        return


USECOLS = ["read_name", "gene_id", "transcript_id",
           "transcript_splice_hash_code", "num_exons", "frac_assigned"]
MAPPING_COLUMNS = ["gene_id", "transcript_id", "transcript_splice_hash_code",
                   "num_exons"]


def _chunk_frame(read_name_bc, gene, tx, hsh, nexon, frac):
    return pd.DataFrame({
        "cell_barcode": np.asarray(read_name_bc, dtype=object),
        "gene_id": np.asarray(gene, dtype=object),
        "transcript_id": np.asarray(tx, dtype=object),
        "transcript_splice_hash_code": np.asarray(hsh, dtype=object),
        "num_exons": np.asarray(nexon, dtype=np.int64),
        "frac_assigned": np.asarray(frac, dtype=np.float32),
    })


def _iter_chunks_direct(filename, chunksize):
    """No CSV parser: split on tabs, take the six columns by header position.

    The cell barcode is pulled from read_name in the same pass, which removes a
    separate .str.split(expand=True) that measured 2.44 s per million rows -- on
    its own more than the whole sparse aggregation.

    The three id columns and the barcode are interned.  pandas' parsers
    deduplicate repeated strings, so their dict lookups hit the fast
    pointer-equality path; a plain split() allocates a fresh string per field and
    would pay a full hash and compare every time.  read_name is NOT interned: it
    is unique per row.
    """
    opener = gzip.open if filename.endswith(".gz") else open
    with opener(filename, "rt") as handle:
        while True:
            line = handle.readline()
            if not line:
                sys.exit(f"Error, {filename} has no header line")
            if not line.startswith("#"):
                break
        header = line.rstrip("\n").split("\t")
        missing = [c for c in USECOLS if c not in header]
        if missing:
            sys.exit(f"Error, {filename} lacks required column(s): {missing}")
        i_rn, i_g, i_t, i_h, i_n, i_f = (header.index(c) for c in USECOLS)

        while True:
            bc = []; g = []; t = []; h = []; n = []; f = []
            bca, ga, ta, ha, na, fa = (bc.append, g.append, t.append,
                                       h.append, n.append, f.append)
            for _ in range(chunksize):
                raw = handle.readline()
                if not raw:
                    break
                if raw.startswith("#"):
                    continue
                parts = raw.split("\t")
                bca(intern(parts[i_rn].partition("^")[0]))
                ga(intern(parts[i_g])); ta(intern(parts[i_t]))
                ha(intern(parts[i_h])); na(parts[i_n]); fa(parts[i_f])
            if not t:
                return
            yield _chunk_frame(bc, g, t, h, n, f)
            if len(t) < chunksize:
                return


def _iter_chunks_pandas(filename, chunksize, engine):
    """The previous readers, kept so --csv_engine python|c still works."""
    dtype_spec = {
        "read_name": str,
        "gene_id": str,
        "transcript_id": str,
        "transcript_splice_hash_code": str,
        "num_exons": np.int64,
        "frac_assigned": np.float32,
    }
    opener = gzip.open if filename.endswith(".gz") else open
    with opener(filename, "rt") as handle:
        kwargs = {"sep": "\t", "chunksize": chunksize, "usecols": USECOLS,
                  "dtype": dtype_spec, "engine": engine, "comment": "#"}
        if engine == "c":
            kwargs["low_memory"] = False
        for chunk in pd.read_csv(handle, **kwargs):
            chunk = chunk.fillna({"gene_id": "", "transcript_id": "",
                                  "transcript_splice_hash_code": ""})
            out = _chunk_frame(
                chunk["read_name"].str.split("^", n=1, expand=True)[0].to_numpy(),
                chunk["gene_id"].to_numpy(), chunk["transcript_id"].to_numpy(),
                chunk["transcript_splice_hash_code"].to_numpy(),
                chunk["num_exons"].to_numpy(), chunk["frac_assigned"].to_numpy())
            yield out


def stream_all_counts(filename, chunksize=1_000_000, engine="direct"):
    """Stream the tracking file once, accumulating all three feature levels.

    Returns (mapping_df, levels, barcode_labels).  The levels are NOT finalized
    here: main() finalizes and writes them one at a time so that no level's
    matrices are held while another is being built.
    """
    logger.info("starting single-pass aggregation from %s (RSS %s)", filename, format_rss())

    levels = [
        LevelAccumulator("gene", "gene_id"),
        LevelAccumulator("isoform", "transcript_id"),
        LevelAccumulator("splice_pattern", "transcript_splice_hash_code"),
    ]
    mapping_entries = set()
    barcode_to_index = {}
    barcode_labels = []
    rows_processed = 0

    opener = gzip.open if filename.endswith(".gz") else open
    _reject_if_weighted_tracking(filename, opener)

    if engine == "direct":
        chunks = _iter_chunks_direct(filename, chunksize)
    else:
        chunks = _iter_chunks_pandas(filename, chunksize, engine)

    for chunk_num, chunk in enumerate(chunks, 1):
        rows_processed += len(chunk)
        logger.info("[chunk %d] processed %d rows (RSS %s)", chunk_num,
                    rows_processed, format_rss())

        # distinct 4-tuples per chunk, not one set.add per row.  Same result as
        # iterating every row; the mapping table has 327,187 entries on a file
        # with 1.78 billion rows.
        for entry in chunk[MAPPING_COLUMNS].drop_duplicates().itertuples(
                index=False, name=None):
            mapping_entries.add(entry)

        cell_barcodes = chunk["cell_barcode"]
        codes, uniques = pd.factorize(cell_barcodes.to_numpy())
        lut = np.empty(len(uniques), dtype=np.int64)
        for k, bc in enumerate(uniques):
            idx = barcode_to_index.get(bc)
            if idx is None:
                idx = barcode_to_index[bc] = len(barcode_labels)
                barcode_labels.append(bc)
            lut[k] = idx
        if len(barcode_labels) >= _INT32_MAX:
            sys.exit("Error, more than 2^31 cell barcodes; widen the accumulator "
                     "index dtype in LevelAccumulator.consume_chunk")
        barcode_idx = pd.Series(lut[codes], index=chunk.index)

        for level in levels:
            level.consume_chunk(chunk, barcode_idx)

    logger.info(
        "finished single-pass aggregation after %d rows (RSS %s)",
        rows_processed,
        format_rss(),
    )

    if mapping_entries:
        mapping_df = pd.DataFrame(sorted(mapping_entries), columns=MAPPING_COLUMNS)
    else:
        mapping_df = pd.DataFrame(columns=MAPPING_COLUMNS)

    return mapping_df, levels, barcode_labels


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
    parser.add_argument("--csv_engine", choices=["direct", "python", "c"], default="direct",
                        help="reader to use. 'direct' (default) uses no CSV parser at all: "
                             "measured 4.1x faster than the python engine and 1.7x faster than "
                             "the C engine, with no segfault exposure. 'python' and 'c' select "
                             "the pandas readers and are kept for fallback.")
    parser.add_argument("--emit_dense_counts", action="store_true",
                        help="also write the dense *_cell_counts.tsv files. Off by default: "
                             "111 GiB on a 1.5 B-read library, and no workflow output block "
                             "references them.")
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

    mapping_df, levels, barcode_labels = stream_all_counts(
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
    del mapping_df

    # One level at a time: finalize, write, free.  Holding all three levels'
    # finalized results at once was the largest single component of the old
    # peak, and nothing needs them together.
    for level in levels:
        counts, matrix, feature_labels, barcodes = level.finalize(
            barcode_labels, build_counts_df=args.emit_dense_counts)
        if counts is not None:
            counts.to_csv(f"{args.output_prefix}.{level.label}_cell_counts.tsv",
                          sep="\t", index=False)
        del counts
        make_sparse_matrix_outputs(
            matrix,
            feature_labels,
            barcodes,
            f"{args.output_prefix}^{level.label}-sparseM",
            gzip_level=args.gzip_level,
        )
        del matrix, feature_labels, barcodes
        level.feature_to_index = None
        level.feature_labels = None
        gc.collect()
        logger.info("%s: done (RSS %s)", level.label, format_rss())

    logger.info("all done (current RSS %s, peak RSS %s)", format_rss(), format_peak_rss())

if __name__ == "__main__":
    main()
