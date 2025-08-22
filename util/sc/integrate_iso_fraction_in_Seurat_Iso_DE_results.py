#!/usr/bin/env python3
"""
Augment transcript DE results with:
  - gene_id from the fractions file (matched by transcript_id)
  - isoform fractions for cluster1/cluster2
  - their difference

Guarantees:
- All original columns from --DE are preserved unmodified.
- Output columns are reordered so that gene_id and transcript_id come first.
- New columns are appended at the end:
    isoform_fraction_cluster1,
    isoform_fraction_cluster2,
    isoform_fraction_diff

Assumptions:
- --fractions (TSV) has 'transcript_id', one or more 'Cluster_<int>' columns, and
  (optionally) 'gene_id' (not strictly required—if missing, we add a blank column).
- --DE (TSV) has 'transcript_id', 'cluster1', 'cluster2'.
- cluster1/cluster2 directly map to 'Cluster_<int>'.
"""

import sys
import argparse
import pandas as pd
import numpy as np


def parse_args():
    ap = argparse.ArgumentParser(
        description="Append gene_id and isoform fractions to DE results."
    )
    ap.add_argument(
        "--fractions", required=True, help="TSV: isoform fractions matrix (file 1)."
    )
    ap.add_argument("--DE", required=True, help="TSV: transcript DE results (file 2).")
    ap.add_argument("-o", "--out", help="Output TSV (default: stdout).")
    return ap.parse_args()


def _to_cluster_int(x, colname):
    if pd.isna(x):
        raise ValueError(f"NA cluster index in column '{colname}'")
    s = str(x).strip()
    try:
        return int(round(float(s)))
    except Exception:
        raise ValueError(f"Invalid cluster index {x!r} in column '{colname}'")


def main():
    args = parse_args()

    try:
        f1 = pd.read_csv(args.fractions, sep="\t", dtype=str).replace({"": np.nan})
    except Exception as e:
        sys.exit(f"ERROR: Failed to read fractions file: {e}")
    try:
        f2 = pd.read_csv(args.DE, sep="\t", dtype=str).replace({"": np.nan})
    except Exception as e:
        sys.exit(f"ERROR: Failed to read DE file: {e}")

    if "transcript_id" not in f1.columns:
        sys.exit(
            "ERROR: fractions file (--fractions) missing required column: transcript_id"
        )
    for col in ("transcript_id", "cluster1", "cluster2"):
        if col not in f2.columns:
            sys.exit(f"ERROR: DE file (--DE) missing required column: {col}")

    cluster_cols = [c for c in f1.columns if c.startswith("Cluster_")]
    if not cluster_cols:
        sys.exit("ERROR: fractions file has no 'Cluster_<int>' columns.")
    f1_num = f1.copy()
    for c in cluster_cols:
        f1_num[c] = pd.to_numeric(f1_num[c], errors="coerce")

    f1_lookup = f1_num.set_index("transcript_id", drop=True)

    gene_id_from_f1 = None
    if "gene_id" in f1.columns:
        tmp = f1[["transcript_id", "gene_id"]].dropna(subset=["transcript_id"]).copy()
        dup_counts = tmp.groupby("transcript_id")["gene_id"].nunique()
        ambiguous = dup_counts[dup_counts > 1]
        if len(ambiguous) > 0:
            print(
                f"WARNING: {len(ambiguous)} transcript_id(s) map to multiple gene_id values in --fractions; "
                f"using the first non-null per transcript_id.",
                file=sys.stderr,
            )
        gene_id_from_f1 = (
            tmp.dropna(subset=["gene_id"])
            .drop_duplicates(subset=["transcript_id"], keep="first")
            .set_index("transcript_id")["gene_id"]
        )

    try:
        c1_idx = f2["cluster1"].map(lambda x: _to_cluster_int(x, "cluster1"))
        c2_idx = f2["cluster2"].map(lambda x: _to_cluster_int(x, "cluster2"))
    except ValueError as e:
        sys.exit(f"ERROR: {e}")

    c1_cols = c1_idx.map(lambda i: f"Cluster_{i}")
    c2_cols = c2_idx.map(lambda i: f"Cluster_{i}")

    needed_cols = set(c1_cols).union(set(c2_cols))
    missing_needed = [c for c in sorted(needed_cols) if c not in f1_lookup.columns]
    if missing_needed:
        sys.exit(
            "ERROR: fractions file is missing required columns referenced by DE: "
            + ", ".join(missing_needed)
        )

    def get_fraction(tid, col):
        if pd.isna(tid) or tid not in f1_lookup.index:
            return np.nan
        return f1_lookup.at[tid, col]

    iso1 = [get_fraction(tid, col) for tid, col in zip(f2["transcript_id"], c1_cols)]
    iso2 = [get_fraction(tid, col) for tid, col in zip(f2["transcript_id"], c2_cols)]

    out_df = f2.copy()

    gene_col_name = "gene_id"
    if gene_col_name in out_df.columns:
        gene_col_name = "gene_id_from_fractions"
    if gene_id_from_f1 is not None:
        out_df[gene_col_name] = (
            out_df["transcript_id"].map(gene_id_from_f1).astype("string")
        )
    else:
        out_df[gene_col_name] = pd.Series([pd.NA] * len(out_df), dtype="string")

    out_df["isoform_fraction_cluster1"] = pd.to_numeric(iso1, errors="coerce")
    out_df["isoform_fraction_cluster2"] = pd.to_numeric(iso2, errors="coerce")
    out_df["isoform_fraction_diff"] = (
        out_df["isoform_fraction_cluster1"] - out_df["isoform_fraction_cluster2"]
    )

    # Reorder so gene_id and transcript_id are first
    cols = list(out_df.columns)
    if gene_col_name in cols and "transcript_id" in cols:
        reordered = [gene_col_name, "transcript_id"] + [
            c for c in cols if c not in (gene_col_name, "transcript_id")
        ]
        out_df = out_df.loc[:, reordered]

    if args.out:
        try:
            out_df.to_csv(args.out, sep="\t", index=False, na_rep="")
        except Exception as e:
            sys.exit(f"ERROR: Failed to write output: {e}")
    else:
        out_df.to_csv(sys.stdout, sep="\t", index=False, na_rep="")


if __name__ == "__main__":
    main()
