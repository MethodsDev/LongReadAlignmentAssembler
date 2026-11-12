#!/usr/bin/env python3
"""tracking_to_cluster_pseudobulk_exprs.py: per-cluster quant.expr from one pooled quant."""

import csv
import gzip
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent / "tracking_to_cluster_pseudobulk_exprs.py"

EXPR_HEADER = ["gene_id", "transcript_id", "uniq_reads", "all_reads", "isoform_fraction",
               "unique_gene_read_fraction", "TPM", "exons", "introns", "splice_hash_code",
               "RPM_total_reads", "uniq_FSM_reads", "has_FSM_read"]
TRACK_HEADER = ["gene_id", "transcript_id", "transcript_splice_hash_code", "num_exons", "mp_id",
                "read_name", "frac_assigned", "read_weight", "is_unique", "is_FSM"]

# (transcript, barcode, frac, is_unique, is_FSM); g1 has t1 and t2, g2 has t3
ROWS = [
    ("t1", "AAA", 1.0, 1, 1),      # cluster 0, unique FSM
    ("t1", "AAA", 0.75, 0, 0),     # cluster 0, ambiguous read r2 ...
    ("t2", "AAA", 0.25, 0, 1),     # ... its other share, FSM to t2
    ("t2", "CCC", 1.0, 1, 0),      # cluster 1, unique non-FSM
    ("t3", "CCC", 1.0, 1, 1),      # cluster 1
    ("t3", "GGG", 1.0, 1, 1),      # cluster 1
    ("t1", "TTT", 1.0, 1, 1),      # barcode in no cluster
]
GENE = {"t1": "g1", "t2": "g1", "t3": "g2", "t4": "g2"}


def _pooled_expr_rows(rows, tamper=None):
    out = []
    for tx in ("t1", "t2", "t3", "t4"):
        mine = [r for r in rows if r[0] == tx]
        uniq = sum(r[3] for r in mine)
        frac = sum(r[2] for r in mine)
        ufsm = sum(1 for r in mine if r[3] and r[4])
        has = int(any(r[4] for r in mine))
        if tamper == tx:
            uniq += 1
        out.append([GENE[tx], tx, str(uniq), f"{frac:.1f}", "0.000", "0.000", "0.000",
                    f"chr1:(+)[[1, 10]]", "", f"hash_{tx}", "0.000", str(ufsm), str(has)])
    return out


def _write_inputs(tmp_path, tamper=None):
    expr = tmp_path / "pooled.quant.expr"
    with open(expr, "w") as fh:
        fh.write("# LRAA version vTEST\n")
        fh.write("\t".join(EXPR_HEADER) + "\n")
        for r in _pooled_expr_rows(ROWS, tamper):
            fh.write("\t".join(r) + "\n")
    tracking = tmp_path / "pooled.quant.tracking.gz"
    with gzip.open(tracking, "wt") as fh:
        fh.write("# LRAA version vTEST\n")
        fh.write("\t".join(TRACK_HEADER) + "\n")
        for i, (tx, bc, frac, uniq, fsm) in enumerate(ROWS):
            fh.write("\t".join([GENE[tx], tx, f"hash_{tx}", "2", f"MP{i}", f"{bc}^UMI{i}^m/{i}",
                                f"{frac:.6f}", "1.000000", str(uniq), str(fsm)]) + "\n")
    clusters = tmp_path / "clusters.tsv"
    clusters.write_text("cell_barcode\tcluster\nAAA\t0\nCCC\t1\nGGG\t1\nEMPTY\t2\n")
    return expr, tracking, clusters


def _run(tmp_path, expr, tracking, clusters):
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--tracking", str(tracking), "--quant_expr", str(expr),
         "--cell_clusters", str(clusters), "--output_prefix", "S", "--output_dir",
         str(tmp_path / "out"), "--partition_summary", str(tmp_path / "partition.tsv")],
        capture_output=True, text=True)


def _read(path):
    with open(path) as fh:
        lines = [l for l in fh if not l.startswith("#")]
    return {r["transcript_id"]: r for r in csv.DictReader(lines, delimiter="\t")}


def test_cluster_files_sum_their_cells(tmp_path):
    result = _run(tmp_path, *_write_inputs(tmp_path))
    assert result.returncode == 0, result.stderr

    names = sorted(p.name for p in (tmp_path / "out").iterdir())
    # every cluster gets a file, the empty one too; the header row is not a cluster
    assert names == [f"S.{c}.LRAA.quant-only.quant.expr" for c in ("0", "1", "2")]

    c0 = _read(tmp_path / "out" / "S.0.LRAA.quant-only.quant.expr")
    assert list(c0) == ["t1", "t2", "t3", "t4"]   # every pooled isoform, sorted
    assert c0["t1"]["uniq_reads"] == "1" and c0["t1"]["all_reads"] == "1.8"
    assert c0["t1"]["uniq_FSM_reads"] == "1" and c0["t1"]["has_FSM_read"] == "1"
    # FSM but not unique: has_FSM_read without uniq_FSM_reads
    assert c0["t2"]["uniq_FSM_reads"] == "0" and c0["t2"]["has_FSM_read"] == "1"
    assert c0["t1"]["isoform_fraction"] == "0.875"          # 1.75 / 2.0 in g1
    assert c0["t1"]["unique_gene_read_fraction"] == "0.500"  # 1 / 2.0
    assert c0["t1"]["TPM"] == "875000.000" == c0["t1"]["RPM_total_reads"]
    assert c0["t3"]["all_reads"] == "0.0" and c0["t3"]["TPM"] == "0.000"
    # non-count columns come from the pooled file
    assert c0["t3"]["splice_hash_code"] == "hash_t3"

    c1 = _read(tmp_path / "out" / "S.1.LRAA.quant-only.quant.expr")
    assert c1["t3"]["uniq_reads"] == "2" and c1["t3"]["uniq_FSM_reads"] == "2"
    assert c1["t2"]["has_FSM_read"] == "0" and c1["t1"]["all_reads"] == "0.0"

    c2 = _read(tmp_path / "out" / "S.2.LRAA.quant-only.quant.expr")
    assert all(r["all_reads"] == "0.0" for r in c2.values())


def test_unclustered_reads_are_in_no_cluster_and_reported(tmp_path):
    result = _run(tmp_path, *_write_inputs(tmp_path))
    assert result.returncode == 0, result.stderr
    t1 = sum(float(_read(tmp_path / "out" / f"S.{c}.LRAA.quant-only.quant.expr")["t1"]["all_reads"])
             for c in ("0", "1", "2"))
    assert t1 == pytest.approx(1.8)   # pooled 2.8 minus the unclustered barcode's 1.0

    lines = [l for l in (tmp_path / "partition.tsv").read_text().splitlines() if not l.startswith("#")]
    part = {r["cluster"]: r for r in csv.DictReader(lines, delimiter="\t")}
    assert part["__unclustered__"]["assigned_reads"] == "1.0"
    assert part["__unclustered__"]["barcodes_with_reads"] == "1"
    assert part["0"]["assigned_reads"] == "2.0" and part["1"]["assigned_reads"] == "3.0"
    assert sum(float(r["assigned_reads"]) for r in part.values()) == pytest.approx(6.0)


def test_a_tracking_from_another_run_is_refused(tmp_path):
    expr, tracking, clusters = _write_inputs(tmp_path, tamper="t3")
    result = _run(tmp_path, expr, tracking, clusters)
    assert result.returncode != 0
    assert "same quantification" in result.stderr
    assert not (tmp_path / "out").exists()
