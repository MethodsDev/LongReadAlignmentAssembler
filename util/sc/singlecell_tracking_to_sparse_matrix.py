#!/usr/bin/env python3
import argparse, os, gzip, sys
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite
from collections import defaultdict

def stream_group_counts(filename, feature_col, chunksize=1_000_000):
    """
    Stream the tracking file in chunks, aggregate frac_assigned
    by (feature_id, cell_barcode), with progress monitoring.
    """
    file_size = os.path.getsize(filename)
    bytes_read = 0
    agg = defaultdict(float)

    with open(filename, "rb") as f:
        reader = pd.read_csv(f, sep="\t", chunksize=chunksize,
                             usecols=["read_name", feature_col, "frac_assigned"])
        for i, chunk in enumerate(reader, 1):
            # track progress by file pointer
            bytes_read = f.tell()
            pct = (bytes_read / file_size) * 100
            sys.stderr.write(f"\r  [chunk {i}] processed ~{pct:.1f}% of file")
            sys.stderr.flush()

            # split read_name
            split = chunk["read_name"].str.split("^", n=2, expand=True)
            chunk["cell_barcode"] = split[0]

            # aggregate within this chunk
            grouped = chunk.groupby([feature_col, "cell_barcode"], observed=True)["frac_assigned"].sum()
            for (feat, bc), val in grouped.items():
                agg[(feat, bc)] += val

    sys.stderr.write("\r  finished 100.0%\n")

    # convert dict → dataframe
    feature_ids, barcodes, counts = zip(*((f, c, v) for (f, c), v in agg.items()))
    df = pd.DataFrame({"feature_id": feature_ids,
                       "cell_barcode": barcodes,
                       "UMI_counts": counts})
    return df

def make_sparse_matrix_outputs(counts_data, outdirname):
    print(f"-making sparse matrix outputs for: {outdirname}")
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
    features.cat.categories.to_series().to_csv(os.path.join(outdirname, "features.tsv"),
                                               index=False, header=False)
    barcodes.cat.categories.to_series().to_csv(os.path.join(outdirname, "barcodes.tsv"),
                                               index=False, header=False)

    # gzip outputs
    for fn in ["matrix.mtx", "features.tsv", "barcodes.tsv"]:
        path = os.path.join(outdirname, fn)
        with open(path, "rb") as f_in, gzip.open(path + ".gz", "wb") as f_out:
            f_out.writelines(f_in)
        os.remove(path)
    print(f"done with {outdirname}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tracking", required=True, help="input data file quant.tracking")
    parser.add_argument("--output_prefix", required=True, help="output prefix")
    parser.add_argument("--chunksize", type=int, default=1_000_000,
                        help="rows per chunk (default 1e6)")
    args = parser.parse_args()

    # write mapping file (small)
    ids = pd.read_csv(args.tracking, sep="\t",
                      usecols=["gene_id","transcript_id","transcript_splice_hash_code"])
    ids.drop_duplicates().to_csv(f"{args.output_prefix}.gene_transcript_splicehashcode.tsv",
                                 sep="\t", index=False)

    for label, col in [("gene", "gene_id"),
                       ("isoform", "transcript_id"),
                       ("splice_pattern", "transcript_splice_hash_code")]:
        print(f"-processing {label} level counts")
        counts = stream_group_counts(args.tracking, col, chunksize=args.chunksize)
        counts.to_csv(f"{args.output_prefix}.{label}_cell_counts.tsv", sep="\t", index=False)
        make_sparse_matrix_outputs(counts, f"{args.output_prefix}^{label}-sparseM")

    print("all done.")

if __name__ == "__main__":
    main()

