#!/usr/bin/env python3

import importlib.util
import gzip
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]


def _load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SPARSE_CONVERTER = _load_module(
    "singlecell_tracking_to_sparse_matrix_test_module",
    REPO_ROOT / "util" / "sc" / "singlecell_tracking_to_sparse_matrix.py",
)
SYMBOL_ANNOTATOR = _load_module(
    "incorporate_gene_symbols_test_module",
    REPO_ROOT / "util" / "sc" / "incorporate_gene_symbols_in_sc_features.py",
)
DIU = _load_module(
    "sc_pseudobulk_diff_usage_test_module",
    REPO_ROOT
    / "util"
    / "sc"
    / "diff_iso_usage"
    / "sc_pseudobulk_test_isoform_DiffUsage.py",
)
TRACKED_READ_EXTRACTOR = _load_module(
    "extract_tracked_reads_test_module",
    REPO_ROOT / "util" / "extract_tracked_reads_from_bam.py",
)


def _mapping_df():
    return pd.DataFrame(
        {
            "gene_id": ["G1", "G1"],
            "transcript_id": ["ENST000001", "LRAA:iso-spliced"],
            "transcript_splice_hash_code": ["ENST000001", "HASH2"],
            "num_exons": [1, 3],
            "new_gene_id": ["GENE^G1", "GENE^G1"],
            "new_transcript_id": [
                "GENE^ENST000001",
                "GENE^LRAA:iso-spliced",
            ],
            "new_transcript_splice_hash_code": [
                "GENE^ENST000001",
                "GENE^HASH2",
            ],
        }
    )


def test_tracking_mapping_and_symbol_annotation_preserve_num_exons(tmp_path):
    tracking = tmp_path / "sample.quant.tracking"
    tracking.write_text(
        "gene_id\ttranscript_id\ttranscript_splice_hash_code\tnum_exons\t"
        "mp_id\tread_name\tfrac_assigned\tread_weight\n"
        "G1\tENST000001\tENST000001\t1\tmp1\tBC1^U1^r1\t1.000\t1.000\n"
        "G1\tTX2\tHASH2\t3\tmp2\tBC2^U2^r2\t1.000\t1.000\n"
    )

    mapping, _ = SPARSE_CONVERTER.stream_all_counts(str(tracking), chunksize=1)

    assert mapping.columns.tolist() == [
        "gene_id",
        "transcript_id",
        "transcript_splice_hash_code",
        "num_exons",
    ]
    assert mapping.set_index("transcript_id")["num_exons"].to_dict() == {
        "ENST000001": 1,
        "TX2": 3,
    }

    mapping_path = tmp_path / "mapping.tsv"
    mapping.to_csv(mapping_path, sep="\t", index=False)
    SYMBOL_ANNOTATOR.write_annotated_id_mappings_file(str(mapping_path), {}, {})
    annotated = pd.read_csv(str(mapping_path) + ".wAnnotIDs", sep="\t")

    assert annotated["num_exons"].tolist() == [1, 3]


def test_unspliced_filter_uses_explicit_status_for_all_identifier_forms():
    num_exons_by_id = DIU.build_num_exons_mapping(_mapping_df())
    counts = pd.DataFrame(
        {
            "gene_id": ["G1", "G1"],
            "transcript_id": ["ENST000001", "LRAA:iso-spliced"],
            "A": [5, 7],
            "B": [4, 8],
        }
    )
    fractions = pd.DataFrame(
        {
            "gene_id": ["G1", "G1"],
            "transcript_id": ["GENE^ENST000001", "GENE^HASH2"],
            "A": [0.5, 0.7],
            "B": [0.4, 0.8],
        }
    )

    retained_counts = DIU.filter_unspliced_isoforms(
        counts, num_exons_by_id, "counts matrix"
    )
    retained_fractions = DIU.filter_unspliced_isoforms(
        fractions, num_exons_by_id, "fraction matrix"
    )

    assert retained_counts["transcript_id"].tolist() == ["LRAA:iso-spliced"]
    assert retained_fractions["transcript_id"].tolist() == ["GENE^HASH2"]


@pytest.mark.parametrize("invalid_num_exons", [None, 0, -1, 1.5, "two"])
def test_num_exons_mapping_rejects_invalid_values(invalid_num_exons):
    mapping = _mapping_df().iloc[[0]].copy()
    mapping["num_exons"] = mapping["num_exons"].astype(object)
    mapping.loc[mapping.index[0], "num_exons"] = invalid_num_exons

    with pytest.raises(ValueError, match="integer >= 1"):
        DIU.build_num_exons_mapping(mapping)


def test_num_exons_mapping_rejects_missing_column_and_conflicts():
    with pytest.raises(ValueError, match="requires a 'num_exons' column"):
        DIU.build_num_exons_mapping(_mapping_df().drop(columns="num_exons"))

    conflicting = pd.concat(
        [
            _mapping_df().iloc[[0]],
            _mapping_df().iloc[[0]].assign(num_exons=2),
        ],
        ignore_index=True,
    )
    with pytest.raises(ValueError, match="Conflicting num_exons metadata"):
        DIU.build_num_exons_mapping(conflicting)


def test_unspliced_filter_rejects_unmapped_features():
    with pytest.raises(ValueError, match="could not map num_exons metadata"):
        DIU.filter_unspliced_isoforms(
            pd.DataFrame({"transcript_id": ["missing"]}),
            DIU.build_num_exons_mapping(_mapping_df()),
            "counts matrix",
        )


def test_tracking_reader_uses_named_column_in_compressed_commented_file(tmp_path):
    tracking = tmp_path / "sample.quant.tracking.gz"
    with gzip.open(tracking, "wt") as handle:
        handle.write(
            "# LRAA version 0.17.7\n"
            "gene_id\ttranscript_id\ttranscript_splice_hash_code\tnum_exons\t"
            "mp_id\tread_name\tfrac_assigned\tread_weight\n"
            "G1\tTX1\tHASH1\t2\tmp1\tBC1^U1^read-one\t1.000\t1.000\n"
        )

    assert TRACKED_READ_EXTRACTOR.get_tracked_read_names(str(tracking)) == {
        "read-one"
    }


def test_legacy_mapping_remains_usable_when_filter_is_disabled(tmp_path):
    counts = tmp_path / "counts.tsv"
    counts.write_text(
        "gene_id\ttranscript_id\tA\tB\n"
        "G1\tENST000001\t40\t10\n"
        "G1\tTX2\t10\t40\n"
    )
    mapping = tmp_path / "mapping.tsv"
    mapping.write_text(
        "gene_id\ttranscript_id\ttranscript_splice_hash_code\n"
        "G1\tENST000001\tENST000001\n"
        "G1\tTX2\tHASH2\n"
    )
    output_prefix = tmp_path / "legacy"
    script = (
        REPO_ROOT
        / "util"
        / "sc"
        / "diff_iso_usage"
        / "sc_pseudobulk_test_isoform_DiffUsage.py"
    )

    result = subprocess.run(
        [
            sys.executable,
            str(script),
            "--sc_cluster_counts_matrix",
            str(counts),
            "--splice_hashcode_id_mappings",
            str(mapping),
            "--output_prefix",
            str(output_prefix),
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
