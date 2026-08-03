#!/usr/bin/env python3

import gzip
import subprocess
import sys
from pathlib import Path


def test_single_cell_symbol_mapping_uses_ranked_strand_aware_gffcompare_parser(
    tmp_path,
):
    script = Path(__file__).with_name("incorporate_gene_symbols_in_sc_features.py")
    ref_gtf = tmp_path / "reference.gtf"
    id_mappings = tmp_path / "id_mappings.tsv"
    tracking = tmp_path / "gffcmp.tracking"
    sparse_matrix_dir = tmp_path / "isoform-sparseM"
    sparse_matrix_dir.mkdir()

    ref_gtf.write_text(
        'chr1\tref\ttranscript\t1\t100\t.\t+\t.\tgene_id "REF_C"; transcript_id "REF_C_TX"; gene_name "GENE_C";\n'
        'chr1\tref\ttranscript\t1\t100\t.\t+\t.\tgene_id "REF_J"; transcript_id "REF_J_TX"; gene_name "GENE_J";\n'
        'chr1\tref\ttranscript\t1\t100\t.\t-\t.\tgene_id "REF_X"; transcript_id "REF_X_TX"; gene_name "GENE_X";\n'
    )
    id_mappings.write_text(
        "gene_id\ttranscript_id\ttranscript_splice_hash_code\n"
        "LRAA_GENE\tLRAA_TX\tLRAA_SPLICE\n"
        "LRAA_ANTI_GENE\tLRAA_ANTI_TX\tLRAA_ANTI_SPLICE\n"
    )
    tracking.write_text(
        "TCONS_1\tXLOC_1\tREF_C|REF_C_TX\tc\tq1:LRAA_GENE|LRAA_TX|1|1.0|1.0|1.0\n"
        "TCONS_2\tXLOC_2\tREF_J|REF_J_TX\tj\tq1:LRAA_GENE|LRAA_TX|1|1.0|1.0|1.0\n"
        "TCONS_3\tXLOC_3\tREF_X|REF_X_TX\tx\tq1:LRAA_ANTI_GENE|LRAA_ANTI_TX|1|1.0|1.0|1.0\n"
    )
    with gzip.open(sparse_matrix_dir / "features.tsv.gz", "wt") as features:
        features.write("LRAA_TX\nLRAA_ANTI_TX\n")

    result = subprocess.run(
        [
            sys.executable,
            str(script),
            "--ref_gtf",
            str(ref_gtf),
            "--id_mappings",
            str(id_mappings),
            "--sparseM_dirs",
            str(sparse_matrix_dir),
            "--gffcompare_tracking",
            str(tracking),
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    with gzip.open(sparse_matrix_dir / "features.tsv.gz", "rt") as features:
        assert features.read().splitlines() == ["GENE_C^LRAA_TX", "LRAA_ANTI_TX"]
