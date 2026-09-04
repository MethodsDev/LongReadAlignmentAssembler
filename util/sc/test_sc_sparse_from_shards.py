#!/usr/bin/env python3
"""V1 tests for the per-shard sparse build and streaming merge.

Every case is checked against an independent reference implementation that
tallies the tracking rows directly with pandas, so the tests do not encode the
same assumptions the code under test makes.
"""
import os, subprocess, sys
import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = os.path.join(HERE, "tracking_to_shard_sparse.py")
MERGE = os.path.join(HERE, "merge_shard_sparse_matrices.py")
COLS = ["gene_id", "transcript_id", "transcript_splice_hash_code", "num_exons",
        "mp_id", "read_name", "frac_assigned", "read_weight"]
LEVEL_COL = {"gene": "gene_id", "isoform": "transcript_id",
             "splice_pattern": "transcript_splice_hash_code"}


def write_tracking(path, rows, comments=("# LRAA version test",)):
    with open(path, "w") as fh:
        for c in comments:
            fh.write(c + "\n")
        fh.write("\t".join(COLS) + "\n")
        for gene, tx, hsh, nex, bc, umi, frac in rows:
            fh.write("\t".join([gene, tx, hsh, str(nex), "MP1",
                                f"{bc}^{umi}^read", f"{frac:.6f}", "1.000000"]) + "\n")


def reference(all_rows):
    """Independent tally: the answer the merged matrices must equal."""
    df = pd.DataFrame(all_rows, columns=["gene_id", "transcript_id",
                                         "transcript_splice_hash_code",
                                         "num_exons", "bc", "umi", "frac"])
    out = {}
    for level, col in LEVEL_COL.items():
        g = df.groupby([col, "bc"])["frac"].sum()
        feats = sorted(df[col].unique())
        bcs = sorted(df["bc"].unique())
        fi = {f: i for i, f in enumerate(feats)}
        bi = {b: i for i, b in enumerate(bcs)}
        trip = {(fi[f], bi[b]): float(v) for (f, b), v in g.items()}
        out[level] = (feats, bcs, trip)
    return out


def read_mtx(path):
    with open(path) as fh:
        assert fh.readline().startswith("%%MatrixMarket")
        fh.readline()
        nr, nc, nnz = (int(x) for x in fh.readline().split())
        trip = {}
        n = 0
        for line in fh:
            r, c, v = line.split()
            trip[(int(r) - 1, int(c) - 1)] = float(v)
            n += 1
    assert n == nnz, f"{path}: header says {nnz} entries, body has {n}"
    return nr, nc, trip


def run_case(tmp_path, shards, shared_features=False):
    """shards: {name: [row, ...]}.  Returns (merged results, reference)."""
    dirs = []
    for name, rows in shards.items():
        trk = tmp_path / f"{name}.tracking"
        write_tracking(str(trk), rows)
        outdir = tmp_path / f"shard_{name}"
        subprocess.run([sys.executable, BUILD, "--tracking", str(trk),
                        "--outdir", str(outdir), "--shard", name],
                       check=True, capture_output=True)
        dirs.append(str(outdir))

    prefix = str(tmp_path / "merged")
    cmd = [sys.executable, MERGE, "--shard_dirs", *dirs,
           "--output_prefix", prefix, "--no_gzip"]
    if shared_features:
        cmd.append("--shared-features")
    subprocess.run(cmd, check=True, capture_output=True)

    got = {}
    for level in LEVEL_COL:
        d = f"{prefix}^{level}-sparseM"
        nr, nc, trip = read_mtx(os.path.join(d, "matrix.mtx"))
        feats = [l.rstrip("\n") for l in open(os.path.join(d, "features.tsv"))]
        bcs = [l.rstrip("\n") for l in open(os.path.join(d, "barcodes.tsv"))]
        got[level] = (nr, nc, feats, bcs, trip)

    all_rows = [r for rows in shards.values() for r in rows]
    return got, reference(all_rows)


def check(got, ref):
    for level in LEVEL_COL:
        nr, nc, feats, bcs, trip = got[level]
        r_feats, r_bcs, r_trip = ref[level]
        assert feats == r_feats, f"{level}: feature list differs"
        assert bcs == r_bcs, f"{level}: barcode list differs"
        assert (nr, nc) == (len(r_feats), len(r_bcs)), f"{level}: shape differs"
        assert set(trip) == set(r_trip), (
            f"{level}: support differs; "
            f"missing {sorted(set(r_trip) - set(trip))[:5]} "
            f"extra {sorted(set(trip) - set(r_trip))[:5]}")
        for k, v in r_trip.items():
            assert abs(trip[k] - v) <= 1e-5 * max(1.0, abs(v)), \
                f"{level}: value at {k} is {trip[k]}, expected {v}"


# g, tx, hash, num_exons, barcode, umi, frac
def test_basic_disjoint_features_shared_cells(tmp_path):
    """Basic mode: contigs share cells, features are disjoint."""
    shards = {
        "chr1": [("GA", "TA1", "HA1", 3, "AAAA", "u1", 1.0),
                 ("GA", "TA1", "HA1", 3, "AAAA", "u2", 0.5),   # same pair: must sum
                 ("GA", "TA2", "HA2", 2, "CCCC", "u3", 1.0),
                 ("GB", "TB1", "HB1", 1, "AAAA", "u4", 0.25)],
        "chr2": [("GC", "TC1", "HC1", 4, "AAAA", "u5", 2.0),   # cell seen on both
                 ("GC", "TC1", "HC1", 4, "GGGG", "u6", 1.0)],  # cell unique to chr2
    }
    check(*run_case(tmp_path, shards))


def test_gene_and_splice_are_row_aggregations(tmp_path):
    """Two transcripts of one gene, and two transcripts sharing a splice hash."""
    shards = {
        "chr1": [("GA", "TA1", "SHARED", 3, "AAAA", "u1", 1.0),
                 ("GA", "TA2", "SHARED", 3, "AAAA", "u2", 2.0),
                 ("GA", "TA2", "SHARED", 3, "CCCC", "u3", 4.0),
                 ("GB", "TB1", "HB1", 1, "AAAA", "u4", 8.0)],
    }
    got, ref = run_case(tmp_path, shards)
    check(got, ref)
    # gene GA row for cell AAAA must be 1.0 + 2.0; splice SHARED likewise
    _, _, feats, bcs, trip = got["gene"]
    assert trip[(feats.index("GA"), bcs.index("AAAA"))] == pytest.approx(3.0)
    _, _, feats, bcs, trip = got["splice_pattern"]
    assert trip[(feats.index("SHARED"), bcs.index("AAAA"))] == pytest.approx(3.0)


def test_cluster_guided_shared_features_disjoint_cells(tmp_path):
    """Cluster-guided: same features across clusters, cells partition."""
    shards = {
        "c0_chr1": [("GA", "TA1", "HA1", 3, "AAAA", "u1", 1.0),
                    ("GA", "TA2", "HA2", 2, "AAAA", "u2", 2.0)],
        "c1_chr1": [("GA", "TA1", "HA1", 3, "CCCC", "u3", 4.0)],
        "c0_chr2": [("GC", "TC1", "HC1", 1, "AAAA", "u4", 8.0)],
        "c1_chr2": [("GC", "TC1", "HC1", 1, "CCCC", "u5", 16.0)],
    }
    check(*run_case(tmp_path, shards, shared_features=True))


def test_singleton_and_empty_shard(tmp_path):
    shards = {
        "chr1": [("GA", "TA1", "HA1", 1, "AAAA", "u1", 1.0)],
        "chrY": [("GY", "TY1", "HY1", 1, "TTTT", "u2", 0.125)],
    }
    check(*run_case(tmp_path, shards))


def test_shared_feature_without_flag_is_refused(tmp_path):
    """The invariant that makes summing unnecessary is enforced, not assumed."""
    shards = {"a": [("GA", "TA1", "HA1", 1, "AAAA", "u1", 1.0)],
              "b": [("GA", "TA1", "HA1", 1, "CCCC", "u2", 1.0)]}
    dirs = []
    for name, rows in shards.items():
        trk = tmp_path / f"{name}.tracking"
        write_tracking(str(trk), rows)
        outdir = tmp_path / f"shard_{name}"
        subprocess.run([sys.executable, BUILD, "--tracking", str(trk),
                        "--outdir", str(outdir), "--shard", name],
                       check=True, capture_output=True)
        dirs.append(str(outdir))
    p = subprocess.run([sys.executable, MERGE, "--shard_dirs", *dirs,
                        "--output_prefix", str(tmp_path / "m"), "--no_gzip"],
                       capture_output=True, text=True)
    assert p.returncode != 0
    assert "appears in 2 shards" in p.stderr


def test_column_indices_ascending_within_each_row(tmp_path):
    """The monotonic-remap property the merge relies on, checked on output."""
    shards = {
        "chr1": [("GA", "TA1", "H1", 1, bc, f"u{i}", 1.0)
                 for i, bc in enumerate(["TTTT", "AAAA", "GGGG", "CCCC"])],
        "chr2": [("GB", "TB1", "H2", 1, bc, f"v{i}", 1.0)
                 for i, bc in enumerate(["ACAC", "TGTG"])],
    }
    got, ref = run_case(tmp_path, shards)
    check(got, ref)
    d = tmp_path / "merged^isoform-sparseM" / "matrix.mtx"
    seen = {}
    with open(str(d)) as fh:
        fh.readline(); fh.readline(); fh.readline()
        for line in fh:
            r, c, _ = line.split()
            r, c = int(r), int(c)
            assert c > seen.get(r, 0), f"row {r}: column {c} not ascending"
            seen[r] = c
