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


def _flagged(header, name, start, cigar, *, secondary=False, supplementary=False):
    a = _read(header, name, start, cigar)
    a.is_secondary = secondary
    a.is_supplementary = supplementary
    return a


def test_supplementary_records_do_not_inflate_depth(tmp_path):
    """Depth must be measured over the records the consumer reads, and no others.

    Every consumer discards supplementary alignments, so counting them raises the
    measured depth above anything downstream sees. Acceptance probability is
    target/depth, so the reads that remain are thinned harder than the target asks
    and each surviving read's 1/p overstates what it stands for. Here the shallow
    locus sits under the target on real records alone; padding it with supplementary
    records must not push it over.
    """
    header = _header()
    reads = [_read(header, f"real{i}", 50000 + (i % 20), f"{READ_LEN}M") for i in range(40)]
    reads += [
        _flagged(header, f"suppl{i}", 50000 + (i % 20), f"{READ_LEN}M", supplementary=True)
        for i in range(400)
    ]
    src = tmp_path / "sup.bam"
    _write_bam(src, reads)
    out = tmp_path / "sup.norm.bam"
    norm.sift_bam(str(src), str(out), 100, 100, 42, "primary")

    kept = _load(out)
    assert not any(k.startswith("suppl") for k in kept), "supplementary records are not evidence"
    real = {k: v for k, v in kept.items() if k.startswith("real")}
    assert len(real) == 40, "the real locus is under target and must pass through whole"
    assert all(w == pytest.approx(1.0) for w, _ in real.values())


def test_secondary_records_count_only_when_the_consumer_keeps_them(tmp_path):
    """The filter mode must decide, because the consumer's policy decides.

    Under primary-only the secondaries are neither counted nor written, so the locus
    stays under target. With secondaries retained they are part of the evidence and
    the same locus is over target, so it thins.
    """
    header = _header()
    reads = [_read(header, f"real{i}", 50000 + (i % 20), f"{READ_LEN}M") for i in range(40)]
    reads += [
        _flagged(header, f"sec{i}", 50000 + (i % 20), f"{READ_LEN}M", secondary=True)
        for i in range(400)
    ]
    src = tmp_path / "sec.bam"
    _write_bam(src, reads)

    prim = tmp_path / "sec.primary.bam"
    norm.sift_bam(str(src), str(prim), 100, 100, 42, "primary")
    kept_prim = _load(prim)
    assert not any(k.startswith("sec") for k in kept_prim)
    assert len([k for k in kept_prim if k.startswith("real")]) == 40
    assert all(w == pytest.approx(1.0) for w, _ in kept_prim.values())

    both = tmp_path / "sec.with.bam"
    norm.sift_bam(str(src), str(both), 100, 100, 42, "with_secondary")
    kept_both = _load(both)
    assert len(kept_both) < 440, "counted as evidence, the locus is over target and thins"
    assert any(w > 1.0 for w, _ in kept_both.values()), "thinned reads must carry a weight"


def test_records_the_consumer_rejects_do_not_inflate_depth(tmp_path):
    """Depth must be measured over exactly the records quantification consumes.

    Bam_alignment_extractor rejects supplementary, secondary (primary-only), improper-pair,
    low-MAPQ, duplicate, qcfail and sub-identity alignments. Every one of those counted
    into depth here would raise measured coverage above anything downstream sees, lower the
    acceptance probability, and inflate every surviving read's XW weight at that locus.

    The clean reads alone sit under the target; padding with rejects of each kind must not
    push the locus over it, and none of the rejects may appear in the output.
    """
    header = _header()
    reads = []
    for i in range(40):
        a = _read(header, f"clean{i}", 50000 + (i % 20), f"{READ_LEN}M")
        a.set_tag("NM", 1)
        reads.append(a)

    def _pad(name, n, **flags):
        for i in range(n):
            a = _read(header, f"{name}{i}", 50000 + (i % 20), f"{READ_LEN}M")
            a.set_tag("NM", 1)
            for k, v in flags.items():
                setattr(a, k, v)
            reads.append(a)

    _pad("suppl", 100, is_supplementary=True)
    _pad("sec", 100, is_secondary=True)
    _pad("dup", 100, is_duplicate=True)
    _pad("qcf", 100, is_qcfail=True)
    for i in range(100):                     # improper pair
        a = _read(header, f"pair{i}", 50000 + (i % 20), f"{READ_LEN}M")
        a.set_tag("NM", 1)
        a.is_paired = True
        a.is_proper_pair = False
        reads.append(a)
    for i in range(100):                     # below the MAPQ floor
        a = _read(header, f"lowq{i}", 50000 + (i % 20), f"{READ_LEN}M")
        a.set_tag("NM", 1)
        a.mapping_quality = 3
        reads.append(a)
    for i in range(100):                     # below the identity floor
        a = _read(header, f"lowid{i}", 50000 + (i % 20), f"{READ_LEN}M")
        a.set_tag("NM", 30)
        reads.append(a)

    src = tmp_path / "rejects.bam"
    _write_bam(src, reads)
    out = tmp_path / "rejects.norm.bam"
    norm.sift_bam(str(src), str(out), 100, 100, 42, "primary", 97.0, 10)

    kept = _load(out)
    for prefix in ("suppl", "sec", "dup", "qcf", "pair", "lowq", "lowid"):
        assert not any(k.startswith(prefix) for k in kept), (
            f"{prefix} records are not evidence and must not be written"
        )
    clean = {k: v for k, v in kept.items() if k.startswith("clean")}
    assert len(clean) == 40, "700 rejected records must not push a 40x locus over target"
    assert all(w == pytest.approx(1.0) for w, _ in clean.values())


def test_low_identity_reads_do_not_inflate_depth(tmp_path):
    """Depth must be measured at the consumer's per-identity floor, not below it.

    Bam_alignment_extractor discards alignments under config['min_per_id'] -- 97 in HiFi
    mode -- before quantification sees them. Counting them here measures coverage that
    does not exist downstream, so acceptance probability comes out too small and every
    surviving read's 1/p too large. On an ONT chr20 bam at the default 80 this is 8% of
    reads. Here the clean reads alone sit under target; noisy reads must not push the
    locus over it.
    """
    header = _header()
    reads = []
    for i in range(40):
        a = _read(header, f"clean{i}", 50000 + (i % 20), f"{READ_LEN}M")
        a.set_tag("NM", 1)                     # ~99.7% identity
        reads.append(a)
    for i in range(400):
        a = _read(header, f"noisy{i}", 50000 + (i % 20), f"{READ_LEN}M")
        a.set_tag("NM", 30)                    # 90% identity, below a 97 floor
        reads.append(a)
    src = tmp_path / "perid.bam"
    _write_bam(src, reads)

    out = tmp_path / "perid.norm.bam"
    norm.sift_bam(str(src), str(out), 100, 100, 42, "primary", 97.0)
    kept = _load(out)
    assert not any(k.startswith("noisy") for k in kept), "sub-floor reads are not evidence"
    clean = {k: v for k, v in kept.items() if k.startswith("clean")}
    assert len(clean) == 40, "the clean locus is under target and must pass through whole"
    assert all(w == pytest.approx(1.0) for w, _ in clean.values())

    # with the floor disabled the noisy reads count, and the locus thins
    out0 = tmp_path / "perid.nofloor.bam"
    norm.sift_bam(str(src), str(out0), 100, 100, 42, "primary", 0)
    kept0 = _load(out0)
    assert len(kept0) < 440, "counted as evidence, the locus is over target"
    assert any(w > 1.0 for w, _ in kept0.values()), "thinned reads must carry a weight"



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


@pytest.fixture
def source_bam(tmp_path):
    p = tmp_path / "sample.bam"
    p.write_bytes(b"A" * 100)
    return p


def _stem(source, level=1000, filter_mode="with_secondary", min_per_id=0,
          min_mapping_quality=0, depth_window=100, random_seed=42):
    import Util_funcs

    return Util_funcs.splice_graph_norm_cache_stem(
        "s.quant", level, str(source), filter_mode, min_per_id,
        min_mapping_quality, depth_window, random_seed,
    )


def test_cache_name_distinguishes_every_parameter_that_changes_the_bam(source_bam):
    """The driver returns on this stem's checkpoint, so it never sees the utility's token.

    Anything that changes which reads the utility keeps but is absent from this name is
    invisible: the run inherits a bam that does not match its settings, and no consumer
    can tell. The mapping-quality floor, the depth window and the seed are all in that
    category -- the driver holds the last two fixed today, which is exactly why leaving
    them out would go unnoticed until someone changed a default.
    """
    base = _stem(source_bam)
    assert _stem(source_bam, min_mapping_quality=30) != base
    assert _stem(source_bam, depth_window=50) != base
    assert _stem(source_bam, random_seed=7) != base
    assert "mapq30" in _stem(source_bam, min_mapping_quality=30)
    assert "w50" in _stem(source_bam, depth_window=50)
    assert "s7" in _stem(source_bam, random_seed=7)


def test_cache_name_distinguishes_the_per_id_floor(source_bam):
    """HiFi's floor of 97 admits a different read population than the default 80.

    The driver returns as soon as it sees this stem's checkpoint, so it never runs the
    utility and never consults the utility's own token. If the floor is absent here, a
    HiFi run silently inherits the bam a default run normalized, whose depth was
    measured over reads HiFi discards.
    """
    assert _stem(source_bam, min_per_id=97) != _stem(source_bam, min_per_id=80)
    assert "pid97" in _stem(source_bam, min_per_id=97)


def test_cache_name_distinguishes_the_alignment_filter_mode(source_bam):
    """Two filter modes yield different bams from one source, so one key cannot serve both.

    The mode decides which records count toward depth and which get written, so a
    primary-only run and a with-secondary run at the same target disagree on every
    acceptance probability. Sharing a key would let whichever ran first hand its bam
    to the other, and nothing downstream could tell.
    """
    assert _stem(source_bam, filter_mode="primary") != _stem(
        source_bam, filter_mode="with_secondary"
    )
    assert "primary" in _stem(source_bam, filter_mode="primary")


def test_cache_name_distinguishes_the_method(source_bam):
    """A cache from a superseded method must not be reachable by the current key.

    Nothing downstream can detect one. A bam produced by read-start binning
    carries no XW tag, and an untagged read correctly weighs 1, so its distorted
    counts would be consumed in silence.
    """
    import LRAA_Globals

    current = _stem(source_bam)
    previous = LRAA_Globals.SPLICE_GRAPH_NORMALIZATION_METHOD
    try:
        LRAA_Globals.SPLICE_GRAPH_NORMALIZATION_METHOD = "startbin1"
        superseded = _stem(source_bam)
    finally:
        LRAA_Globals.SPLICE_GRAPH_NORMALIZATION_METHOD = previous

    assert current != superseded
    assert previous in current


def test_cache_name_distinguishes_the_target_depth(source_bam):
    """Also scopes the work directory, whose checkpoints are otherwise shared."""
    assert _stem(source_bam, 1000) != _stem(source_bam, 5000)


def test_cache_name_distinguishes_inputs_sharing_a_basename(tmp_path, source_bam):
    """The stem carries only a basename, so identity cannot rest on the name."""
    other_dir = tmp_path / "other"
    other_dir.mkdir()
    other = other_dir / "sample.bam"
    other.write_bytes(b"A" * 100)  # same name, same size, different file

    assert _stem(source_bam) != _stem(other)


def test_cache_name_follows_an_input_replaced_in_place(source_bam):
    """Regenerating a bam at the same path must not hit the old cache."""
    import time

    before = _stem(source_bam)
    time.sleep(0.01)
    source_bam.write_bytes(b"B" * 100)  # identical size, new contents

    assert _stem(source_bam) != before


def test_cache_name_is_stable_for_an_untouched_input(source_bam):
    """The point of the cache: an unchanged input must still hit it.

    Keying on identity is only worth it if it does not also defeat reuse across
    reruns, which is the entire reason the normalized bam is kept.
    """
    assert _stem(source_bam) == _stem(source_bam)


def test_every_retained_read_is_tagged(tmp_path):
    """The weight is what makes a cache self-describing, so none may be missing."""
    header = _header()
    src = tmp_path / "in.bam"
    _write_bam(src, _corpus(header))
    out = tmp_path / "out.bam"
    norm.sift_bam(str(src), str(out), 100, 100, 42)

    with pysam.AlignmentFile(str(out), "rb") as fh:
        untagged = [a.query_name for a in fh.fetch(until_eof=True) if not a.has_tag("XW")]
    assert untagged == []
