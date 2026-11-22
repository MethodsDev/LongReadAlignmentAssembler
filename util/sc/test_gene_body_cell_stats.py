#!/usr/bin/env python3
"""Single-nucleus helpers: gene-body matrices -> per-barcode table, GTF gene spans,
and the per-shard sidecar concatenation."""

import csv
import gzip
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
MISC = HERE.parent / "misc"

HEADER = "gene_id\tread_name\tfrac_assigned\tsplice_class\tsource\tn_genes\taligned_bp\texonic_bp\tintronic_bp\n"


def _sidecar(path, rows, comment="# LRAA version test\n", trailer=""):
    with gzip.open(path, "wt") as fh:
        fh.write(comment + HEADER)
        for gene, read, frac, cls in rows:
            fh.write(f"{gene}\t{read}\t{frac:.6f}\t{cls}\tisoform\t1\t100\t100\t0\n")
        fh.write(trailer)


def test_cell_stats_from_gene_body_matrices(tmp_path):
    sidecar = tmp_path / "x.quant.gene_body.tracking.gz"
    _sidecar(
        sidecar,
        [
            ("A^ENSG1", "AAA^u1^r1", 1.0, "S"),
            ("A^ENSG1", "AAA^u2^r2", 1.0, "U"),
            ("B^ENSG2", "AAA^u3^r3", 1.0, "U"),
            ("B^ENSG2", "AAA^u4^r4", 0.5, "A"),
            ("A^ENSG1", "AAA^u4^r4", 0.5, "A"),
            ("B^ENSG2", "CCC^u5^r5", 1.0, "U"),
            ("A^ENSG1", "CCC^u6^r6", 0.0, "U"),  # dropped tie: counted nowhere
            ("A^ENSG1", "untagged_read", 1.0, "S"),  # no cell tags: no column
        ],
    )
    prefix = tmp_path / "P"
    subprocess.run(
        [sys.executable, str(HERE / "gene_body_tracking_to_sparse_matrix.py"),
         "--gene_body_tracking", str(sidecar), "--output_prefix", str(prefix)],
        check=True, capture_output=True,
    )
    out = tmp_path / "stats.tsv"
    subprocess.run(
        [sys.executable, str(HERE / "gene_body_cell_stats.py"),
         "--sparseM_prefix", str(prefix), "--output", str(out)],
        check=True, capture_output=True,
    )
    rows = {r["cell_barcode"]: r for r in csv.DictReader(open(out), delimiter="\t")}
    assert set(rows) == {"AAA", "CCC"}
    a = rows["AAA"]
    assert float(a["gene_full"]) == 4.0
    assert (float(a["spliced"]), float(a["unspliced"]), float(a["ambiguous"])) == (1.0, 2.0, 1.0)
    assert int(a["n_genes_full"]) == 2
    assert float(a["spliced_frac"]) == 0.25 and float(a["unspliced_frac"]) == 0.5
    c = rows["CCC"]
    assert (float(c["gene_full"]), float(c["unspliced_frac"]), int(c["n_genes_full"])) == (1.0, 1.0, 1)


def test_gtf_gene_spans_use_output_gene_ids(tmp_path):
    gtf = tmp_path / "a.gtf"
    gtf.write_text(
        'chr1\tX\texon\t100\t200\t.\t+\t.\tgene_id "ENSG1"; transcript_id "T1"; gene_name "AAA";\n'
        'chr1\tX\texon\t900\t1000\t.\t+\t.\tgene_id "ENSG1"; transcript_id "T2"; gene_name "AAA";\n'
        'chr1\tX\texon\t50\t60\t.\t-\t.\tgene_id "g:chr1:-:comp-3"; transcript_id "t3";\n'
        'chr1\tX\ttranscript\t1\t5000\t.\t+\t.\tgene_id "ENSG1"; transcript_id "T1"; gene_name "AAA";\n'
    )
    out = tmp_path / "spans.tsv"
    subprocess.run(
        [sys.executable, str(HERE / "gtf_gene_spans.py"), "--gtf", str(gtf), "--output", str(out)],
        check=True, capture_output=True,
    )
    rows = {r["gene_id"]: r for r in csv.DictReader(open(out), delimiter="\t")}
    assert set(rows) == {"AAA^ENSG1", "g:chr1:-:comp-3"}
    assert (rows["AAA^ENSG1"]["lend"], rows["AAA^ENSG1"]["rend"], rows["AAA^ENSG1"]["span"]) == ("100", "1000", "901")
    assert rows["g:chr1:-:comp-3"]["span"] == "11"


def test_concat_gene_body_tracking(tmp_path):
    a, b = tmp_path / "a.gz", tmp_path / "b.gz"
    _sidecar(a, [("G1", "x^u^r1", 1.0, "U")], comment="# LRAA version 1\n# LRAA CMD: shard0\n",
             trailer="#strand_check\tchr1+\t40\t1\n")
    _sidecar(b, [("G2", "x^u^r2", 1.0, "S")], comment="# LRAA version 1\n# LRAA CMD: shard1\n",
             trailer="#strand_check\tchr2-\t60\t2\n")
    out = tmp_path / "merged.gz"
    subprocess.run(
        [sys.executable, str(MISC / "concat_gene_body_tracking.py"),
         "--inputs", str(a), str(b), "--output", str(out)],
        check=True, capture_output=True,
    )
    lines = gzip.open(out, "rt").read().splitlines()
    assert lines[:2] == ["# LRAA version 1", "# LRAA CMD: shard0"]
    assert lines[2] + "\n" == HEADER
    assert [l.split("\t")[0] for l in lines[3:] if not l.startswith("#")] == ["G1", "G2"]
    # each shard's strand-check line is kept; its other comment lines are not
    assert [l for l in lines if l.startswith("#strand_check")] == [
        "#strand_check\tchr1+\t40\t1", "#strand_check\tchr2-\t60\t2"]
    assert "# LRAA CMD: shard1" not in lines

    bad = tmp_path / "bad.gz"
    with gzip.open(bad, "wt") as fh:
        fh.write("gene_id\tread_name\n")
    r = subprocess.run(
        [sys.executable, str(MISC / "concat_gene_body_tracking.py"),
         "--inputs", str(a), str(bad), "--output", str(tmp_path / "x.gz")],
        capture_output=True,
    )
    assert r.returncode != 0


def test_matrix_builder_writes_strand_check(tmp_path):
    sidecar = tmp_path / "x.quant.gene_body.tracking.gz"
    _sidecar(sidecar, [("A^ENSG1", "AAA^u1^r1", 1.0, "S")],
             trailer="#strand_check\tchr1+\t300\t6\n#strand_check\tchr1-\t200\t4\n")
    prefix = tmp_path / "P"
    tsv = tmp_path / "P.strand_check.tsv"
    subprocess.run(
        [sys.executable, str(HERE / "gene_body_tracking_to_sparse_matrix.py"),
         "--gene_body_tracking", str(sidecar), "--output_prefix", str(prefix),
         "--strand_check_tsv", str(tsv)],
        check=True, capture_output=True,
    )
    (row,) = list(csv.DictReader(open(tsv), delimiter="\t"))
    assert (row["units"], row["spliced_reads_checked"], row["spliced_reads_antisense"]) == ("2", "500", "10")
    assert row["verdict"] == "stranded" and row["antisense_frac"] == "0.020000"
    # the strand-check lines do not reach the matrices
    full = gzip.open(f"{prefix}^gene_full-sparseM/features.tsv.gz", "rt").read().split()
    assert full == ["A^ENSG1"]
