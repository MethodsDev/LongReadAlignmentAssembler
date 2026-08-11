#!/usr/bin/env python3

"""Coverage normalization must thin depth without distorting relative support.

The splice graph decides whether to keep a junction by comparing its support
against the strongest junction in its exon island. That comparison is only
meaningful if every feature's reads survived normalization at the same rate, or
if the rates are recorded so they can be divided back out. Sampling read starts
per bin did neither: a junction at a crowded TSS drew its reads from the one bin
sampled hardest while a junction further into the gene drew from many sparse
bins that survived whole, and nothing recorded the difference.
"""

import sys
from pathlib import Path

import pysam
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
for _sub in ("util", "pylib"):
    _p = str(REPO_ROOT / _sub)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import normalize_bam_by_strand as norm


CONTIG = "chr1"
CONTIG_LEN = 200000
READ_LEN = 300


def _header():
    return pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"},
         "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}]}
    )


def _read(header, name, start, cigar):
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = 60
    a.cigarstring = cigar
    span = sum(l for op, l in a.cigartuples if op in (0, 1))
    a.query_sequence = "A" * span
    a.query_qualities = pysam.qualitystring_to_array("I" * span)
    return a


def _write_bam(path, reads):
    header = _header()
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for r in sorted(reads, key=lambda x: x.reference_start):
            fh.write(r)
    pysam.index(str(path))


def _load(path):
    """read name -> (weight, junctions)"""
    out = {}
    with pysam.AlignmentFile(str(path), "rb") as fh:
        for a in fh.fetch(until_eof=True):
            out[a.query_name] = (
                float(a.get_tag("XW")) if a.has_tag("XW") else None,
                tuple(norm._read_junctions(a)),
            )
    return out


def _corpus(header):
    """A deep locus and a shallow one, with two junctions inside the deep one.

    deep    : 4000 reads stacked over 10,100-10,450 -- far above any cap here
    shallow : 40 reads stacked at 100,000 -- below every cap here
    scarce  : 30 reads sharing one junction, spliced out of the deep stack
    common  : 900 reads sharing another junction out of the same stack

    Both spliced sets start at staggered positions but pad the leading match so
    every read in a set reports the identical junction; support is a property of
    the junction, so a set that spread across ten near-identical junctions would
    be ten scarce ones rather than the single common one intended.
    """
    reads = []
    for i in range(4000):
        reads.append(_read(header, f"deep{i}", 10100 + (i % 50), f"{READ_LEN}M"))
    for i in range(40):
        reads.append(_read(header, f"shallow{i}", 100000 + (i % 20), f"{READ_LEN}M"))
    for i in range(30):
        k = i % 10
        reads.append(_read(header, f"scarce{i}", 10150 - k, f"{50 + k}M2000N100M"))
    for i in range(900):
        k = i % 10
        reads.append(_read(header, f"common{i}", 10150 - k, f"{50 + k}M5000N100M"))
    return reads


@pytest.fixture
def normalized(tmp_path):
    header = _header()
    src = tmp_path / "in.bam"
    _write_bam(src, _corpus(header))
    out = tmp_path / "out.bam"
    norm.sift_bam(str(src), str(out), 100, 100, 42)
    return _load(out)


def test_reads_below_the_target_are_all_kept_and_weigh_one(normalized):
    """Normalization, not downsampling: shallow coverage must pass through whole."""
    shallow = {k: v for k, v in normalized.items() if k.startswith("shallow")}
    assert len(shallow) == 40
    assert all(w == pytest.approx(1.0) for w, _ in shallow.values())


def test_depth_above_the_target_is_thinned(normalized):
    deep = [k for k in normalized if k.startswith("deep")]
    assert len(deep) < 4000 / 4, "4000x coverage was not thinned toward a target of 100"


def test_a_scarce_junction_keeps_every_read(normalized):
    """Support for a junction below the target is exact, not estimated.

    These reads sit inside the deep locus, so a purely depth-driven rule would
    thin them to a handful and leave the frequency test to a coin flip.
    """
    scarce = {k: v for k, v in normalized.items() if k.startswith("scarce")}
    assert len(scarce) == 30
    assert all(w == pytest.approx(1.0) for w, _ in scarce.values())


def test_weights_recover_the_support_of_a_thinned_junction(normalized):
    """The common junction is thinned, but its weights must add back up."""
    kept = [(w, j) for w, j in normalized.values() if j and j[0][1] - j[0][0] == 5000]
    assert len(kept) < 900, "the common junction was not thinned"
    assert sum(w for w, _ in kept) == pytest.approx(900, rel=0.25)


def test_relative_support_survives_normalization(normalized):
    """The quantity the splice graph actually decides on.

    True ratio is 30/900. Counting surviving reads instead of summing weights
    gives roughly 30/(900*p) -- inflated by 1/p, which is the whole defect.
    """
    scarce = sum(w for w, j in normalized.values() if j and j[0][1] - j[0][0] == 2000)
    common = sum(w for w, j in normalized.values() if j and j[0][1] - j[0][0] == 5000)
    assert scarce / common == pytest.approx(30 / 900, rel=0.25)


def test_sampling_is_reproducible_and_independent_of_file_order(tmp_path):
    """A read's fate follows its name, so shuffling the input cannot change it.

    Coordinate order is not a property of the evidence, and a scheme that let
    position decide would reintroduce exactly the positional bias being removed.
    """
    header = _header()
    reads = _corpus(header)

    a_bam = tmp_path / "a.bam"
    _write_bam(a_bam, reads)
    a_out = tmp_path / "a.norm.bam"
    norm.sift_bam(str(a_bam), str(a_out), 100, 100, 42)

    # same reads, written in the opposite order within each start position
    b_bam = tmp_path / "b.bam"
    _write_bam(b_bam, list(reversed(reads)))
    b_out = tmp_path / "b.norm.bam"
    norm.sift_bam(str(b_bam), str(b_out), 100, 100, 42)

    assert set(_load(a_out)) == set(_load(b_out))


def test_a_different_seed_selects_a_different_sample(tmp_path):
    header = _header()
    src = tmp_path / "in.bam"
    _write_bam(src, _corpus(header))
    out1, out2 = tmp_path / "s1.bam", tmp_path / "s2.bam"
    norm.sift_bam(str(src), str(out1), 100, 100, 1)
    norm.sift_bam(str(src), str(out2), 100, 100, 2)
    assert set(_load(out1)) != set(_load(out2))


def test_disabling_the_target_keeps_everything(tmp_path):
    header = _header()
    src = tmp_path / "in.bam"
    _write_bam(src, _corpus(header))
    out = tmp_path / "off.bam"
    norm.sift_bam(str(src), str(out), 0, 100, 42)
    assert len(_load(out)) == 4970
