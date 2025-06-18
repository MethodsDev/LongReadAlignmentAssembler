#!/usr/bin/env python3


import pandas as pd
import glob

# Path to your input files (adjust as needed)
input_files = glob.glob("*expr")

# Read and concatenate all input files
df = pd.concat((pd.read_csv(f, sep="\t") for f in input_files), ignore_index=True)

# Group by (gene_id, transcript_id) to sum relevant fields and keep constant ones
grouped = df.groupby(["gene_id", "transcript_id"], as_index=False).agg(
    {"uniq_reads": "sum", "all_reads": "sum", "exons": "first", "introns": "first"}
)

# Recompute TPM (transcripts per million)
grouped["TPM"] = grouped["all_reads"] / grouped["all_reads"].sum() * 1e6

# Recompute isoform_fraction and unique_gene_read_fraction per gene_id
grouped["isoform_fraction"] = grouped.groupby("gene_id")["all_reads"].transform(
    lambda x: x / x.sum()
)
grouped["unique_gene_read_fraction"] = grouped.groupby("gene_id")[
    "uniq_reads"
].transform(lambda x: x / grouped.loc[x.index, "all_reads"].sum())

# Rearranging columns to match original format
final_cols = [
    "gene_id",
    "transcript_id",
    "uniq_reads",
    "all_reads",
    "isoform_fraction",
    "unique_gene_read_fraction",
    "TPM",
    "exons",
    "introns",
]
result = grouped[final_cols]

# Save to file
result.to_csv("merged_and_recomputed.tsv", sep="\t", index=False)
