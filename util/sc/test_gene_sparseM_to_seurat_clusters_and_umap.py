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


def _write_10x(sparse_dir, features, barcodes, counts):
    """counts: {(feature_index, barcode_index): value}, 0-based."""
    sparse_dir.mkdir()
    with gzip.open(sparse_dir / "features.tsv.gz", "wt") as fh:
        fh.write("".join(f"{f}\n" for f in features))
    with gzip.open(sparse_dir / "barcodes.tsv.gz", "wt") as fh:
        fh.write("".join(f"{b}\n" for b in barcodes))
    lines = [
        "%%MatrixMarket matrix coordinate integer general\n",
        "%\n",
        f"{len(features)} {len(barcodes)} {len(counts)}\n",
    ]
    lines.extend(f"{i + 1} {j + 1} {v}\n" for (i, j), v in sorted(counts.items()))
    with gzip.open(sparse_dir / "matrix.mtx.gz", "wt") as fh:
        fh.write("".join(lines))


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is unavailable")
def test_nucseq_ranks_variable_genes_within_span_bins(tmp_path):
    # Short genes, ten of them markers of one cell group, and long genes whose counts
    # are large and overdispersed, as intronic reads make them in gene_full. Ranked
    # within span bins, short genes are chosen that plain vst ranking passes over.
    import numpy as np

    rng = np.random.default_rng(7)
    n_cells = 120
    short = [f"SHORT_{k}^ENSG{k:05d}" for k in range(60)]
    long_ = [f"LONG_{k}^ENSG{100 + k:05d}" for k in range(60)]
    features = short + long_
    barcodes = [f"cell-{j}" for j in range(n_cells)]
    counts = {}
    for j in range(n_cells):
        group_a = j < n_cells // 2
        for i in range(len(short)):
            mu = 12.0 if (i < 10 and group_a) else 3.0
            counts[(i, j)] = int(rng.poisson(mu))
        for i in range(len(long_)):
            counts[(len(short) + i, j)] = int(rng.poisson(rng.gamma(0.5, 100.0)))
    counts = {k: v for k, v in counts.items() if v > 0}
    sparse_dir = tmp_path / "gene_full-sparseM"
    _write_10x(sparse_dir, features, barcodes, counts)

    spans = tmp_path / "spans.tsv"
    with open(spans, "wt") as fh:
        fh.write("gene_id\tcontig\tstrand\tlend\trend\tspan\n")
        for k, g in enumerate(short):
            fh.write(f"{g}\tchr1\t+\t1\t{2000 + 10 * k}\t{2000 + 10 * k}\n")
        for k, g in enumerate(long_):
            fh.write(f"{g}\tchr1\t+\t1\t{900000 + 10 * k}\t{900000 + 10 * k}\n")

    script = Path(__file__).with_name("gene_sparseM_to_seurat_clusters_and_umap.R")
    output_prefix = tmp_path / "nuc"
    common = [
        "Rscript", str(script),
        "--sparseM_dir", str(sparse_dir),
        "--min_cells", "1", "--min_features", "1", "--percent_mt_max", "100",
        "--npcs", "3", "--resolution", "0.3", "--n_variable_features", "20",
    ]
    subprocess.run(
        common + ["--output_prefix", str(output_prefix), "--nucseq",
                  "--gene_spans", str(spans), "--span_bins", "2"],
        check=True, capture_output=True, text=True,
    )
    with open(f"{output_prefix}-variable_features.length_aware.tsv", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    # Seurat names features with "-" for "_"; the span join must have followed it
    assert all(r["span"] != "NA" for r in rows)
    assert {r["span_bin"] for r in rows} == {"1", "2"}
    selected = {r["gene"] for r in rows if r["selected"] == "TRUE"}
    assert len(selected) == 20

    # score_in_bin is the z-score of the standardized variance within the span bin,
    # and the selection is the top 20 by it
    by_bin = {}
    for r in rows:
        by_bin.setdefault(r["span_bin"], []).append(float(r["variance_standardized"]))
    for r in rows:
        v = np.array(by_bin[r["span_bin"]])
        z = (float(r["variance_standardized"]) - v.mean()) / v.std(ddof=1)
        assert abs(z - float(r["score_in_bin"])) < 1e-6
    top = sorted(rows, key=lambda r: -float(r["score_in_bin"]))[:20]
    assert {r["gene"] for r in top} == selected

    # plain ranking by standardized variance takes mostly the long genes
    plain = sorted(rows, key=lambda r: -float(r["variance_standardized"]))[:20]
    n_short_plain = sum(r["gene"].startswith("SHORT") for r in plain)
    n_short_binned = sum(g.startswith("SHORT") for g in selected)
    assert n_short_plain < 5
    assert n_short_binned >= 8
    assert (tmp_path / "nuc-cell_cluster_assignments.tsv").exists()

    # --nucseq without spans is refused
    failed = subprocess.run(
        common + ["--output_prefix", str(tmp_path / "nospans"), "--nucseq"],
        capture_output=True, text=True,
    )
    assert failed.returncode != 0


def _shallow_and_deep_matrix(sparse_dir, n_deep, n_shallow, n_genes=300, seed=3):
    """n_deep barcodes detect ~250 genes, n_shallow detect ~60 (two expression groups)."""
    import numpy as np

    rng = np.random.default_rng(seed)
    features = [f"G{i}^ENSG{i:05d}" for i in range(n_genes)]
    barcodes = [f"deep-{j}" for j in range(n_deep)] + [f"shallow-{j}" for j in range(n_shallow)]
    counts = {}
    for j, b in enumerate(barcodes):
        n_detect = 250 if b.startswith("deep") else 60
        genes = rng.choice(n_genes, size=n_detect, replace=False)
        for i in genes:
            mu = 6.0 if (i < 30 and j % 2 == 0) else 2.0
            counts[(int(i), j)] = int(rng.poisson(mu)) + 1
    _write_10x(sparse_dir, features, barcodes, counts)
    return barcodes


def _run_seurat(sparse_dir, prefix, extra, npcs="5"):
    script = Path(__file__).with_name("gene_sparseM_to_seurat_clusters_and_umap.R")
    return subprocess.run(
        ["Rscript", str(script), "--sparseM_dir", str(sparse_dir),
         "--output_prefix", str(prefix), "--min_cells", "1", "--percent_mt_max", "100",
         "--npcs", npcs, "--resolution", "0.5", "--n_variable_features", "100"] + extra,
        check=True, capture_output=True, text=True,
    )


def _read_settings(path):
    with open(path, newline="") as fh:
        return {r["setting"]: r["value"] for r in csv.DictReader(fh, delimiter="\t")}


def _assigned(prefix):
    with open(f"{prefix}-cell_cluster_assignments.tsv", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is unavailable")
def test_adaptive_min_features_lowers_floor_when_too_few_pass(tmp_path):
    # 10 deep + 90 shallow barcodes: 10 pass min_features 200 (< 200 cells, < 50%), so
    # the floor drops to the 10th percentile of detected genes (about 60).
    _shallow_and_deep_matrix(tmp_path / "m", n_deep=10, n_shallow=90)
    prefix = tmp_path / "adapt"
    _run_seurat(tmp_path / "m", prefix, ["--min_features", "200", "--adaptive_min_features",
                                          "--min_features_lower_bound", "20"])
    rec = _read_settings(f"{prefix}-min_features.tsv")
    assert rec["min_features_requested"] == "200"
    assert rec["adjusted"] == "TRUE"
    assert 20 <= int(rec["min_features_used"]) <= 60
    assert rec["barcodes_passing_requested"] == "10"
    assert int(rec["barcodes_passing_used"]) >= 90
    assert len(_assigned(prefix)) >= 90

    # the lower bound wins over a lower quantile
    prefix_b = tmp_path / "bounded"
    _run_seurat(tmp_path / "m", prefix_b, ["--min_features", "200", "--adaptive_min_features",
                                            "--min_features_lower_bound", "150"])
    rec_b = _read_settings(f"{prefix_b}-min_features.tsv")
    assert rec_b["min_features_used"] == "150"


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is unavailable")
def test_adaptive_min_features_keeps_floor_when_enough_pass(tmp_path):
    # 80 deep + 20 shallow: 80 pass (80%), but fewer than min_passing_cells 200 would
    # trigger, so lower that to 50; then the requested floor stands and the run matches
    # one without the flag.
    _shallow_and_deep_matrix(tmp_path / "m", n_deep=80, n_shallow=20)
    plain = tmp_path / "plain"
    adapt = tmp_path / "adapt"
    _run_seurat(tmp_path / "m", plain, ["--min_features", "200"])
    _run_seurat(tmp_path / "m", adapt, ["--min_features", "200", "--adaptive_min_features",
                                         "--min_passing_cells", "50"])
    rec = _read_settings(f"{adapt}-min_features.tsv")
    assert rec["adjusted"] == "FALSE" and rec["min_features_used"] == "200"
    assert _assigned(plain) == _assigned(adapt)
    assert not (tmp_path / "plain-min_features.tsv").exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is unavailable")
def test_fewer_cells_than_pcs_gives_one_cluster(tmp_path):
    # 4 cells pass with npcs 5: no PCA; all cells in cluster 0, UMAP coordinates NA.
    _shallow_and_deep_matrix(tmp_path / "m", n_deep=4, n_shallow=30)
    prefix = tmp_path / "few"
    _run_seurat(tmp_path / "m", prefix, ["--min_features", "200"], npcs="5")
    rows = _assigned(prefix)
    assert len(rows) == 4 and {r["cluster"] for r in rows} == {"0"}
    with open(f"{prefix}-cell_cluster_assignments.wUMAP.tsv", newline="") as fh:
        umap = list(csv.DictReader(fh, delimiter="\t"))
    assert len(umap) == 4 and all(r["umap_1"] == "NA" for r in umap)
    for suffix in ("-seurat_obj.rds", "-seurat_obj.initial.rds", "-umap.pdf"):
        assert Path(f"{prefix}{suffix}").exists()


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is unavailable")
def test_few_cells_above_npcs_still_cluster(tmp_path):
    # 12 cells (above npcs + 1, below the default neighbour counts 20 / 30) cluster.
    _shallow_and_deep_matrix(tmp_path / "m", n_deep=12, n_shallow=10)
    prefix = tmp_path / "twelve"
    _run_seurat(tmp_path / "m", prefix, ["--min_features", "200"], npcs="5")
    assert len(_assigned(prefix)) == 12
    with open(f"{prefix}-cell_cluster_assignments.wUMAP.tsv", newline="") as fh:
        umap = list(csv.DictReader(fh, delimiter="\t"))
    assert all(r["umap_1"] != "NA" for r in umap)
