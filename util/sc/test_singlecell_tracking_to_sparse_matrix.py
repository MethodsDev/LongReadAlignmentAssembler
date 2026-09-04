#!/usr/bin/env python3
"""Regression tests for singlecell_tracking_to_sparse_matrix.py.

The reader was replaced (no CSV parser by default) and the finalize path was
restructured, so these pin the behaviour that must not change: all three
--csv_engine choices must agree with each other and with an independent tally.
"""
import gzip, os, subprocess, sys
import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "singlecell_tracking_to_sparse_matrix.py")
COLS = ["gene_id", "transcript_id", "transcript_splice_hash_code", "num_exons",
        "mp_id", "read_name", "frac_assigned", "read_weight"]
LEVEL_COL = {"gene": "gene_id", "isoform": "transcript_id",
             "splice_pattern": "transcript_splice_hash_code"}


def make_rows(n_cells=40, n_genes=12, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(600):
        gi = int(rng.integers(n_genes))
        ti = int(rng.integers(3))
        bc = "".join(rng.choice(list("ACGT"), 16))
        rows.append((f"G{gi:03d}", f"G{gi:03d}.T{ti}", f"H{gi:03d}_{ti}",
                     1 + ti, bc, f"u{i}", round(float(rng.random()), 6)))
    # a repeated (feature, cell) pair that must be summed
    rows.append(rows[0])
    return rows


def write_tracking(path, rows):
    with open(path, "w") as fh:
        fh.write("# LRAA version test\n")
        fh.write("\t".join(COLS) + "\n")
        for g, t, h, n, bc, umi, frac in rows:
            fh.write("\t".join([g, t, h, str(n), "MP1", f"{bc}^{umi}^rd",
                                f"{frac:.6f}", "1.000000"]) + "\n")


def run(tmp, tracking, engine, dense=False):
    prefix = str(tmp / f"out_{engine}{'_dense' if dense else ''}")
    cmd = [sys.executable, SCRIPT, "--tracking", str(tracking),
           "--output_prefix", prefix, "--csv_engine", engine, "--chunksize", "97"]
    if dense:
        cmd.append("--emit_dense_counts")
    subprocess.run(cmd, check=True, capture_output=True)
    return prefix


def read_out(prefix, level):
    d = f"{prefix}^{level}-sparseM"
    def lines(name):
        p = os.path.join(d, name)
        op = gzip.open if os.path.exists(p + ".gz") else open
        p = p + ".gz" if os.path.exists(p + ".gz") else p
        with op(p, "rt") as fh:
            return [l.rstrip("\n") for l in fh]
    feats, bcs = lines("features.tsv"), lines("barcodes.tsv")
    mm = lines("matrix.mtx")
    nr, nc, nnz = (int(x) for x in mm[2].split())
    trip = {}
    for line in mm[3:]:
        r, c, v = line.split()
        trip[(int(r) - 1, int(c) - 1)] = float(v)
    assert len(trip) == nnz
    return feats, bcs, (nr, nc), trip


def reference(rows, level):
    df = pd.DataFrame(rows, columns=["gene_id", "transcript_id",
                                     "transcript_splice_hash_code", "num_exons",
                                     "bc", "umi", "frac"])
    col = LEVEL_COL[level]
    feats = sorted(df[col].unique()); bcs = sorted(df["bc"].unique())
    fi = {f: i for i, f in enumerate(feats)}; bi = {b: i for i, b in enumerate(bcs)}
    g = df.groupby([col, "bc"])["frac"].sum()
    return feats, bcs, {(fi[f], bi[b]): float(v) for (f, b), v in g.items()}


@pytest.fixture(scope="module")
def tracking(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("trk")
    rows = make_rows()
    p = tmp / "shard.tracking"
    write_tracking(str(p), rows)
    return p, rows, tmp


@pytest.mark.parametrize("engine", ["direct", "python", "c"])
def test_engine_matches_independent_tally(tracking, engine):
    path, rows, tmp = tracking
    prefix = run(tmp, path, engine)
    for level in LEVEL_COL:
        feats, bcs, shape, trip = read_out(prefix, level)
        r_feats, r_bcs, r_trip = reference(rows, level)
        assert feats == r_feats
        assert bcs == r_bcs
        assert shape == (len(r_feats), len(r_bcs))
        assert set(trip) == set(r_trip)
        for k, v in r_trip.items():
            assert trip[k] == pytest.approx(v, rel=1e-6, abs=1e-9)


def test_engines_agree_with_each_other(tracking):
    path, rows, tmp = tracking
    outs = {e: run(tmp, path, e) for e in ("direct", "python", "c")}
    base = read_out(outs["direct"], "isoform")
    for e in ("python", "c"):
        other = read_out(outs[e], "isoform")
        assert base[0] == other[0] and base[1] == other[1] and base[2] == other[2]
        assert set(base[3]) == set(other[3])
        for k in base[3]:
            assert base[3][k] == pytest.approx(other[3][k], rel=1e-6, abs=1e-9)


def test_dense_counts_off_by_default_and_match_when_on(tracking):
    path, rows, tmp = tracking
    off = run(tmp, path, "direct")
    for level in LEVEL_COL:
        assert not os.path.exists(f"{off}.{level}_cell_counts.tsv")

    on = run(tmp, path, "direct", dense=True)
    for level in LEVEL_COL:
        p = f"{on}.{level}_cell_counts.tsv"
        assert os.path.exists(p)
        dense = pd.read_csv(p, sep="\t")
        feats, bcs, _, trip = read_out(on, level)
        assert len(dense) == len(trip)
        fi = {f: i for i, f in enumerate(feats)}; bi = {b: i for i, b in enumerate(bcs)}
        for f, b, v in dense.itertuples(index=False):
            assert trip[(fi[f], bi[b])] == pytest.approx(v, rel=1e-6, abs=1e-9)


def test_mapping_file(tracking):
    path, rows, tmp = tracking
    prefix = run(tmp, path, "direct")
    got = pd.read_csv(f"{prefix}.gene_transcript_splicehashcode.tsv", sep="\t")
    want = pd.DataFrame(rows, columns=["gene_id", "transcript_id",
                                       "transcript_splice_hash_code", "num_exons",
                                       "bc", "umi", "frac"])[
        ["gene_id", "transcript_id", "transcript_splice_hash_code", "num_exons"]
    ].drop_duplicates().sort_values(["gene_id", "transcript_id",
                                     "transcript_splice_hash_code", "num_exons"])
    assert len(got) == len(want)
    assert set(map(tuple, got.to_numpy())) == set(map(tuple, want.to_numpy()))
