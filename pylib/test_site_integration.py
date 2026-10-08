"""Integrated TSS/PolyA sites: every cluster-guided site, plus basic sites it misses.

A basic site is only a supplement when no cluster-guided site of the same type on the
same contig and strand lies within the window -- max_dist_between_alt_{TSS,polyA}_sites,
inclusive, the distance over which LRAA's site definition absorbs ends into one site. The
edges of that rule are what these tests pin: 50 nt is covered and 51 nt is not at the
default window of 50, and strand, contig and site type each separate otherwise-coincident
sites.
"""
import os
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import SiteIntegration
from SiteIntegration import SiteBedFormatError
from SplicePatternCollapse import TSS_BED_COLUMNS, POLYA_BED_COLUMNS

SCRIPT = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "..",
    "util",
    "integrate_TSS_PolyA_sites.py",
)


def _row(site_type, pos, strand="+", contig="chr1", support=10, tids="t1"):
    fields = [
        contig,
        str(pos - 1),
        str(pos),
        "{}:{}:{}:{}".format(site_type, contig, pos, strand),
        str(support),
        strand,
        str(support),
        str(len(tids.split(","))),
        tids,
    ]
    if site_type == "PolyA":
        fields += ["AATAAA", "-20", "False"]
    return fields


def _write_bed(path, site_type, rows):
    columns = TSS_BED_COLUMNS if site_type == "TSS" else POLYA_BED_COLUMNS
    with open(path, "wt") as ofh:
        ofh.write("# provenance comment\n")
        ofh.write("#" + "\t".join(columns) + "\n")
        for row in rows:
            ofh.write("\t".join(row) + "\n")
    return str(path)


def _data_rows(path):
    return [l.rstrip("\n").split("\t") for l in open(path) if not l.startswith("#")]


def _integrate(tmp_path, site_type, primary, supplement, window=50):
    p = _write_bed(tmp_path / "p.bed", site_type, primary)
    s = _write_bed(tmp_path / "s.bed", site_type, supplement)
    out = str(tmp_path / "out.bed")
    counts = SiteIntegration.integrate_site_beds(p, s, site_type, out, window=window)
    return counts, _data_rows(out)


def test_window_edge_is_inclusive(tmp_path):
    counts, rows = _integrate(
        tmp_path,
        "PolyA",
        [_row("PolyA", 1000)],
        [_row("PolyA", 1050), _row("PolyA", 950), _row("PolyA", 1051), _row("PolyA", 949)],
    )
    kept = sorted((int(r[2]), r[-1]) for r in rows)
    assert kept == [(949, "basic"), (1000, "cluster_guided"), (1051, "basic")]
    assert counts["supplement_dropped"] == 2
    assert counts["supplement_kept"] == 2


def test_exact_duplicate_is_dropped(tmp_path):
    counts, rows = _integrate(tmp_path, "TSS", [_row("TSS", 500)], [_row("TSS", 500)])
    assert [r[-1] for r in rows] == ["cluster_guided"]
    assert counts["supplement_dropped"] == 1


def test_strand_and_contig_separate_sites(tmp_path):
    counts, rows = _integrate(
        tmp_path,
        "TSS",
        [_row("TSS", 500, strand="+", contig="chr1")],
        [_row("TSS", 500, strand="-", contig="chr1"), _row("TSS", 500, strand="+", contig="chr2")],
    )
    assert sorted((r[0], r[5], r[-1]) for r in rows) == [
        ("chr1", "+", "cluster_guided"),
        ("chr1", "-", "basic"),
        ("chr2", "+", "basic"),
    ]
    assert counts["supplement_dropped"] == 0


def test_window_scales(tmp_path):
    # window 20: 20 nt covered, 21 nt kept.
    counts, rows = _integrate(
        tmp_path, "TSS", [_row("TSS", 100)], [_row("TSS", 120), _row("TSS", 121)], window=20
    )
    assert sorted(int(r[2]) for r in rows if r[-1] == "basic") == [121]
    assert counts["window"] == 20


def test_half_window_gap_no_longer_kept(tmp_path):
    # 26-50 nt from a cluster-guided site: closer than LRAA ever calls two sites, so the
    # same site called twice, not a supplement.
    counts, rows = _integrate(
        tmp_path, "TSS", [_row("TSS", 1000)], [_row("TSS", 1026), _row("TSS", 1040), _row("TSS", 960)]
    )
    assert [r[-1] for r in rows] == ["cluster_guided"]
    assert counts["supplement_dropped"] == 3


def test_supplement_near_duplicates_collapse_to_the_strongest(tmp_path):
    # one basic site written once per transcript, a few nt apart (DDAH2-like)
    counts, rows = _integrate(
        tmp_path,
        "TSS",
        [],
        [_row("TSS", 260, tids="a"), _row("TSS", 262, support=30, tids="b"),
         _row("TSS", 268, tids="c"), _row("TSS", 400, tids="d")],
    )
    assert [int(r[2]) for r in rows] == [262, 400]
    assert rows[0][8] == "a,b,c" and rows[0][7] == "3"
    assert rows[0][6] == "30"
    assert counts["supplement_collapsed"] == 2 and counts["supplement_kept"] == 2


def test_supplement_collapse_does_not_chain(tmp_path):
    # 0 and 100 are both within 50 of 50, but not of each other: only the strongest
    # (50) absorbs; 0 and 100 are absorbed by it, nothing chains beyond its window
    counts, rows = _integrate(
        tmp_path, "PolyA", [],
        [_row("PolyA", 100), _row("PolyA", 150, support=20), _row("PolyA", 201)],
    )
    assert [int(r[2]) for r in rows] == [150, 201]
    assert counts["supplement_collapsed"] == 1


def test_columns_preserved_and_sorted(tmp_path):
    primary = [_row("PolyA", 900, tids="a,b")]
    supplement = [_row("PolyA", 100, tids="c"), _row("PolyA", 5000, contig="chr0")]
    _, rows = _integrate(tmp_path, "PolyA", primary, supplement)
    assert [(r[0], int(r[2])) for r in rows] == [("chr0", 5000), ("chr1", 100), ("chr1", 900)]
    # every original column is carried unchanged, source appended
    assert rows[2][:-1] == primary[0]
    assert all(len(r) == len(POLYA_BED_COLUMNS) + 1 for r in rows)
    header = [l for l in open(tmp_path / "out.bed") if l.startswith("#")][-1]
    assert header.rstrip("\n")[1:].split("\t") == POLYA_BED_COLUMNS + ["source"]


def test_swapped_site_type_is_refused(tmp_path):
    tss = _write_bed(tmp_path / "tss.bed", "TSS", [_row("TSS", 1)])
    polya = _write_bed(tmp_path / "polya.bed", "PolyA", [_row("PolyA", 1)])
    with pytest.raises(SiteBedFormatError):
        SiteIntegration.integrate_site_beds(tss, polya, "PolyA", str(tmp_path / "o.bed"))


def test_cli_writes_both_beds_and_summary(tmp_path):
    pt = _write_bed(tmp_path / "pt.bed", "TSS", [_row("TSS", 100)])
    pp = _write_bed(tmp_path / "pp.bed", "PolyA", [_row("PolyA", 900)])
    st = _write_bed(tmp_path / "st.bed", "TSS", [_row("TSS", 100), _row("TSS", 400)])
    sp = _write_bed(tmp_path / "sp.bed", "PolyA", [_row("PolyA", 920)])
    prefix = str(tmp_path / "sample")
    subprocess.check_call(
        [sys.executable, SCRIPT, "--primary_TSS_bed", pt, "--primary_PolyA_bed", pp,
         "--supplement_TSS_bed", st, "--supplement_PolyA_bed", sp, "--output_prefix", prefix]
    )
    assert [r[-1] for r in _data_rows(prefix + ".integrated.TSS.bed")] == ["cluster_guided", "basic"]
    assert [r[-1] for r in _data_rows(prefix + ".integrated.PolyA.bed")] == ["cluster_guided"]
    summary = [l.rstrip("\n").split("\t") for l in open(prefix + ".integrated_sites.summary.tsv")]
    assert summary[0] == ["site_type", "window", "cluster_guided_sites",
                          "basic_sites", "basic_supplement_kept",
                          "basic_within_window_dropped", "basic_near_duplicates_collapsed",
                          "integrated_sites"]
    assert summary[1] == ["TSS", "50", "1", "2", "1", "1", "0", "2"]
    assert summary[2] == ["PolyA", "50", "1", "1", "0", "1", "0", "1"]
