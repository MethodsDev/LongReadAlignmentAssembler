#!/usr/bin/env python3
# Checks a sc_nucseq run: the gene_full / spliced / unspliced / ambiguous matrices exist
# for the init pass (and the final quant with --final), share one feature and one barcode
# list, and S + U + A equals gene_full entry for entry; the per-barcode table covers every
# barcode and agrees with the matrices; clustering ran on gene_full with length-aware
# variable genes.

import csv, glob, gzip, io, os, sys, tarfile

from scipy.io import mmread

run_dir = sorted(glob.glob(os.path.join(sys.argv[1], "*_LRAA_singlecell_wf")))[-1]
out = os.path.join(run_dir, "out")
phases = ["init", "final"] if "--final" in sys.argv[2:] else ["init"]


def one(name):
    hits = glob.glob(os.path.join(out, name, "*"))
    assert len(hits) == 1, f"expected one {name} output, found {hits}"
    return hits[0]


def matrix(name):
    with tarfile.open(one(name)) as tar:
        member = {os.path.basename(m.name): tar.extractfile(m).read() for m in tar.getmembers() if m.isfile()}
    m = mmread(io.BytesIO(gzip.decompress(member["matrix.mtx.gz"]))).tocsc()
    lines = lambda k: gzip.decompress(member[k]).decode().split("\n")[:-1]
    return m, lines("features.tsv.gz"), lines("barcodes.tsv.gz")


for phase in phases:
    full, features, barcodes = matrix(f"{phase}_sc_gene_full_sparse_tar_gz")
    parts = {}
    for level in ("spliced", "unspliced", "ambiguous"):
        m, f, b = matrix(f"{phase}_sc_gene_{level}_sparse_tar_gz")
        assert f == features and b == barcodes, f"{phase} {level}: feature or barcode list differs from gene_full"
        parts[level] = m
    resid = abs(parts["spliced"] + parts["unspliced"] + parts["ambiguous"] - full).max()
    assert resid < 1e-6, f"{phase}: S + U + A differs from gene_full by {resid}"

    with gzip.open(one(f"{phase}_gene_body_cell_stats"), "rt") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert [r["cell_barcode"] for r in rows] == barcodes, f"{phase}: cell stats do not cover the barcodes"
    col_sums = full.sum(axis=0).A1
    worst = max(abs(float(r["gene_full"]) - s) for r, s in zip(rows, col_sums))
    assert worst < 1e-2, f"{phase}: cell stats gene_full off by {worst}"
    with open(one(f"{phase}_gene_body_strand_check")) as fh:
        (strand,) = list(csv.DictReader(fh, delimiter="\t"))
    assert int(strand["units"]) > 0 and int(strand["spliced_reads_checked"]) > 0, \
        f"{phase}: strand check counted no spliced reads"
    assert strand["verdict"] in ("stranded", "too_few_spliced_reads"), \
        f"{phase}: the oriented fixture fails the strand check: {strand}"
    print(f"{phase}: strand check {strand['spliced_reads_antisense']} of "
          f"{strand['spliced_reads_checked']} spliced reads antisense, {strand['verdict']}")
    std, _, std_bc = matrix(f"{phase}_sc_gene_sparse_tar_gz")
    print(f"{phase}: gene_full {full.shape[0]} genes x {full.shape[1]} barcodes, "
          f"{full.sum():.0f} reads ({full.sum() / std.sum():.3f}x the gene matrix); "
          f"S {parts['spliced'].sum() / full.sum():.1%}, U {parts['unspliced'].sum() / full.sum():.1%}, "
          f"A {parts['ambiguous'].sum() / full.sum():.1%}")

with open(one("seurat_variable_features_tsv")) as fh:
    hvg = list(csv.DictReader(fh, delimiter="\t"))
assert any(r["selected"] == "TRUE" for r in hvg), "no variable genes selected"
assert sum(r["span"] != "NA" for r in hvg) > 0.9 * len(hvg), "most clustered genes have no span: gene ids do not match the GTF"
if glob.glob(os.path.join(out, "seurat_min_features", "*")):
    with open(one("seurat_min_features")) as fh:
        rec = {r["setting"]: r["value"] for r in csv.DictReader(fh, delimiter="\t")}
    print(f"min_features: requested {rec['min_features_requested']}, used {rec['min_features_used']}")
with open(one("seurat_cluster_assignments")) as fh:
    cells = list(csv.DictReader(fh, delimiter="\t"))
print(f"clustering: {len(cells)} cells in {len(set(c['cluster'] for c in cells))} clusters, "
      f"{sum(r['span'] != 'NA' for r in hvg)} of {len(hvg)} genes with a span")
