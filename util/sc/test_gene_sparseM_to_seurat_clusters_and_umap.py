#!/usr/bin/env python3

import csv
import gzip
import shutil
import subprocess
from pathlib import Path

import pytest


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is unavailable")
def test_default_pattern_filters_lraa_and_symbol_mitochondrial_features(tmp_path):
    mitochondrial_features = [
        "MT-ND1",
        "mt-Nd1",
        "g:chrM:+:comp-1",
        "g:MT:-:comp-2",
        "g:M:+:comp-3",
    ]
    nuclear_features = [f"g:chr1:+:comp-{index}" for index in range(1, 31)]
    feature_names = mitochondrial_features + nuclear_features
    filtered_barcodes = [f"high-mt-{index}" for index in range(1, 6)]
    retained_barcodes = [f"cell-{index}" for index in range(1, 51)]
    barcodes = filtered_barcodes + retained_barcodes

    sparse_dir = tmp_path / "gene-sparseM"
    sparse_dir.mkdir()
    with gzip.open(sparse_dir / "features.tsv.gz", "wt") as features:
        features.write("".join(f"{feature}\n" for feature in feature_names))
    with gzip.open(sparse_dir / "barcodes.tsv.gz", "wt") as barcode_file:
        barcode_file.write("".join(f"{barcode}\n" for barcode in barcodes))

    entries = []
    for cell_index in range(len(barcodes)):
        if cell_index < len(mitochondrial_features):
            entries.append((cell_index + 1, cell_index + 1, 1000))
        else:
            entries.append((1, cell_index + 1, cell_index % 5 + 1))
        for gene_index in range(len(nuclear_features)):
            count = ((gene_index + 1) * (cell_index + 1)) % 7 + 1
            entries.append(
                (len(mitochondrial_features) + gene_index + 1, cell_index + 1, count)
            )

    matrix_lines = [
        "%%MatrixMarket matrix coordinate integer general\n",
        "%\n",
        f"{len(feature_names)} {len(barcodes)} {len(entries)}\n",
    ]
    matrix_lines.extend(f"{row} {column} {value}\n" for row, column, value in entries)
    with gzip.open(sparse_dir / "matrix.mtx.gz", "wt") as matrix_file:
        matrix_file.write("".join(matrix_lines))

    script = Path(__file__).with_name("gene_sparseM_to_seurat_clusters_and_umap.R")
    output_prefix = tmp_path / "result"
    subprocess.run(
        [
            "Rscript",
            str(script),
            "--sparseM_dir",
            str(sparse_dir),
            "--output_prefix",
            str(output_prefix),
            "--min_cells",
            "1",
            "--min_features",
            "1",
            "--percent_mt_max",
            "20",
            "--npcs",
            "2",
            "--resolution",
            "0.2",
            "--n_variable_features",
            "15",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    with open(f"{output_prefix}-cell_cluster_assignments.tsv", newline="") as handle:
        observed_barcodes = {
            row["cell_barcode"] for row in csv.DictReader(handle, delimiter="\t")
        }

    assert observed_barcodes == set(retained_barcodes)
