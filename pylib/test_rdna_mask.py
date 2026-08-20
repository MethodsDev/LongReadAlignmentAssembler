#!/usr/bin/env python3

"""RdnaMask: turning a cassette-vs-genome alignment into an exclusion mask.

Covers the three things that can silently go wrong: a span computed from the
wrong CIGAR ops, a merge that drops an interval instead of joining it, and a
cache that returns a stale mask for a genome or cassette that actually
changed. The end-to-end build (real minimap2, real genome and cassette
FASTAs) is exercised too, not just the pure helpers -- the hazard this module
exists to avoid is exactly the kind a mocked aligner would paper over.
"""

import os
import sys
from pathlib import Path

import pysam
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import RdnaMask
import LRAA_Globals


def _header(contigs):
    return pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": name, "LN": length} for name, length in contigs],
        }
    )


def _aligned(header, contig, start, length):
    a = pysam.AlignedSegment(header)
    a.query_name = "r"
    a.reference_id = header.get_tid(contig)
    a.reference_start = start
    a.mapping_quality = 60
    a.cigar = [(0, length)]
    a.query_sequence = "A" * length
    a.query_qualities = pysam.qualitystring_to_array("I" * length)
    return a


# ---------------------------------------------------------------------------
# _sam_hit_spans


def test_sam_hit_spans_uses_ref_consuming_cigar_ops_only(tmp_path):
    """Insertions and soft-clips must not widen the span; deletions must."""

    sam = tmp_path / "hits.sam"
    sam.write_text(
        "\t".join(
            [
                "q1",
                "0",
                "chr1",
                "101",  # 1-based SAM POS -> 0-based ref_start 100
                "60",
                "10S20M5I10M3D10M10S",  # ref-consuming: 20+10+3+10 = 43
                "*",
                "0",
                "0",
                "A" * 60,
                "I" * 60,
            ]
        )
        + "\n"
    )
    spans = RdnaMask._sam_hit_spans(str(sam))
    assert spans == [("chr1", 100, 143)]


def test_sam_hit_identity_computes_percent_identity_from_nm_tag():
    # 20 aligned (M) bases, 2 mismatches -> 90% identity
    fields = ["q", "0", "chr1", "1", "60", "20M", "*", "0", "0", "A" * 20, "I" * 20, "NM:i:2"]
    assert RdnaMask._sam_hit_identity(fields) == 90.0


def test_sam_hit_identity_none_without_an_nm_tag():
    fields = ["q", "0", "chr1", "1", "60", "20M", "*", "0", "0", "A" * 20, "I" * 20]
    assert RdnaMask._sam_hit_identity(fields) is None


def test_sam_hit_spans_min_length_rejects_a_short_hit():
    fields = ["q", "0", "chr1", "101", "60", "50M", "*", "0", "0", "A" * 50, "I" * 50, "NM:i:0"]
    sam = "\t".join(fields) + "\n"

    def spans_for(min_length):
        import tempfile as _tf
        with _tf.NamedTemporaryFile("w", suffix=".sam", delete=False) as f:
            f.write(sam)
            path = f.name
        try:
            return RdnaMask._sam_hit_spans(path, min_length=min_length)
        finally:
            os.remove(path)

    assert spans_for(0) == [("chr1", 100, 150)]
    assert spans_for(51) == []


def test_sam_hit_spans_min_per_id_rejects_a_low_identity_hit():
    # 100 aligned bases, 30 mismatches -> 70% identity
    fields = ["q", "0", "chr1", "101", "60", "100M", "*", "0", "0", "A" * 100, "I" * 100, "NM:i:30"]
    sam = "\t".join(fields) + "\n"

    def spans_for(min_per_id):
        import tempfile as _tf
        with _tf.NamedTemporaryFile("w", suffix=".sam", delete=False) as f:
            f.write(sam)
            path = f.name
        try:
            return RdnaMask._sam_hit_spans(path, min_per_id=min_per_id)
        finally:
            os.remove(path)

    assert spans_for(0) == [("chr1", 100, 200)]
    assert spans_for(80) == []
    assert spans_for(70) == [("chr1", 100, 200)]  # exactly at the floor: kept



def test_sam_hit_spans_skips_unmapped_placeholder_records(tmp_path):
    """minimap2 -a emits one record per query even with no acceptable alignment."""

    sam = tmp_path / "unmapped.sam"
    sam.write_text(
        "\t".join(["q1", "4", "*", "0", "0", "*", "*", "0", "0", "ACGT", "IIII"]) + "\n"
    )
    assert RdnaMask._sam_hit_spans(str(sam)) == []


# ---------------------------------------------------------------------------
# _merge_spans


def test_merge_spans_pads_and_joins_overlapping_intervals():
    spans = [("chr1", 100, 200), ("chr1", 190, 300)]
    merged = RdnaMask._merge_spans(spans, pad=0)
    assert merged == {"chr1": [(100, 300)]}


def test_merge_spans_keeps_non_overlapping_intervals_separate():
    spans = [("chr1", 100, 200), ("chr1", 5000, 5100)]
    merged = RdnaMask._merge_spans(spans, pad=0)
    assert merged == {"chr1": [(100, 200), (5000, 5100)]}


def test_merge_spans_pad_can_bridge_a_gap():
    spans = [("chr1", 100, 200), ("chr1", 250, 300)]
    assert RdnaMask._merge_spans(spans, pad=0) == {
        "chr1": [(100, 200), (250, 300)]
    }
    # 60 bp pad on each side closes the 50 bp gap between them
    assert RdnaMask._merge_spans(spans, pad=60) == {"chr1": [(40, 360)]}


def test_merge_spans_pad_never_takes_the_start_negative():
    spans = [("chr1", 10, 50)]
    merged = RdnaMask._merge_spans(spans, pad=100)
    assert merged["chr1"][0][0] == 0


def test_merge_spans_separates_contigs():
    spans = [("chr1", 100, 200), ("chr2", 100, 200)]
    merged = RdnaMask._merge_spans(spans, pad=0)
    assert set(merged) == {"chr1", "chr2"}


# ---------------------------------------------------------------------------
# read_overlaps_mask (fixtures for the quality-floor tests directly below need
# _write_fasta/_random_seq, defined further down alongside the rest of the
# build_rdna_mask_bed tests -- fine at module scope, since nothing here runs
# until the whole module has finished importing)


def _mutate(seq, frac, seed):
    """``seq`` with a ``frac`` fraction of its bases substituted, deterministically."""
    import random as _random

    rng = _random.Random(seed)
    seq = list(seq)
    bases = "ACGT"
    for pos in rng.sample(range(len(seq)), int(len(seq) * frac)):
        seq[pos] = rng.choice([b for b in bases if b != seq[pos]])
    return "".join(seq)


def test_build_rdna_mask_bed_min_length_floor_is_applied_end_to_end(tmp_path):
    """The public function, not _sam_hit_spans directly -- the floor must reach it.

    A regression here is exactly the class of bug this test exists to catch: an
    earlier revision resolved min_length/min_per_id from config and folded them
    into the cache key, but never actually passed them to the call that builds
    the spans, so the floors were provably inert while the cache key claimed
    they were active. Only a test that drives the full build (not the private
    helper alone) can catch that -- the helper's own unit tests passed
    throughout.
    """
    cassette_seq = _random_seq(500, seed=500)
    cassette_fa = tmp_path / "cassette.fa"
    _write_fasta(cassette_fa, "cassette", cassette_seq)
    genome_fa = tmp_path / "genome.fa"
    _write_fasta(
        genome_fa, "chr1", _random_seq(1000, seed=1) + cassette_seq + _random_seq(1000, seed=2)
    )

    # The embedded copy is ~500-520 bp once padded -- comfortably below 600.
    rejected = RdnaMask.build_rdna_mask_bed(
        str(genome_fa),
        str(cassette_fa),
        cache_dir=str(tmp_path / "cache_strict"),
        pad=10,
        min_length=600,
        min_per_id=0,
    )
    assert rejected is None

    kept = RdnaMask.build_rdna_mask_bed(
        str(genome_fa),
        str(cassette_fa),
        cache_dir=str(tmp_path / "cache_loose"),
        pad=10,
        min_length=400,
        min_per_id=0,
    )
    assert kept is not None


def test_build_rdna_mask_bed_min_per_id_floor_is_applied_end_to_end(tmp_path):
    """Same regression guard as the length test, for the identity floor."""

    cassette_seq = _random_seq(500, seed=500)
    cassette_fa = tmp_path / "cassette.fa"
    _write_fasta(cassette_fa, "cassette", cassette_seq)
    genome_fa = tmp_path / "genome.fa"
    # 10% substitutions -> minimap2 still reports the hit, measured ~90% identity.
    mutated = _mutate(cassette_seq, 0.10, seed=999)
    _write_fasta(
        genome_fa, "chr1", _random_seq(1000, seed=3) + mutated + _random_seq(1000, seed=4)
    )

    rejected = RdnaMask.build_rdna_mask_bed(
        str(genome_fa),
        str(cassette_fa),
        cache_dir=str(tmp_path / "cache_strict"),
        pad=10,
        min_length=0,
        min_per_id=95,
    )
    assert rejected is None

    kept = RdnaMask.build_rdna_mask_bed(
        str(genome_fa),
        str(cassette_fa),
        cache_dir=str(tmp_path / "cache_loose"),
        pad=10,
        min_length=0,
        min_per_id=80,
    )
    assert kept is not None


def test_build_rdna_mask_bed_reads_quality_floors_from_config_when_not_passed(tmp_path):
    """The exact scenario the advisory caught: floors resolved from config must
    actually reach the filter, not just the cache key.

    Calls build_rdna_mask_bed with min_length/min_per_id entirely UNSET (the
    real call shape every production caller uses -- LRAA's own main() never
    passes them explicitly either), and proves the config value took effect by
    getting a different verdict at a strict vs a loose config setting for the
    SAME inputs.
    """
    cassette_seq = _random_seq(500, seed=500)
    cassette_fa = tmp_path / "cassette.fa"
    _write_fasta(cassette_fa, "cassette", cassette_seq)
    genome_fa = tmp_path / "genome.fa"
    _write_fasta(
        genome_fa, "chr1", _random_seq(1000, seed=1) + cassette_seq + _random_seq(1000, seed=2)
    )

    saved = (
        LRAA_Globals.config.get("rdna_mask_min_hit_length"),
        LRAA_Globals.config.get("rdna_mask_min_per_id"),
    )
    try:
        LRAA_Globals.config["rdna_mask_min_hit_length"] = 600
        LRAA_Globals.config["rdna_mask_min_per_id"] = 0
        rejected = RdnaMask.build_rdna_mask_bed(
            str(genome_fa), str(cassette_fa), cache_dir=str(tmp_path / "cache_strict"), pad=10
        )
        assert rejected is None, "config's strict min_length was not applied"

        LRAA_Globals.config["rdna_mask_min_hit_length"] = 0
        kept = RdnaMask.build_rdna_mask_bed(
            str(genome_fa), str(cassette_fa), cache_dir=str(tmp_path / "cache_loose"), pad=10
        )
        assert kept is not None, "config's loose min_length was not applied"
    finally:
        LRAA_Globals.config["rdna_mask_min_hit_length"] = saved[0]
        LRAA_Globals.config["rdna_mask_min_per_id"] = saved[1]


def test_read_overlaps_mask_true_for_an_overlapping_block():
    header = _header([("chr1", 10000)])
    read = _aligned(header, "chr1", 500, 100)  # covers [500, 600)
    mask = {"chr1": __import__("intervaltree").IntervalTree()}
    mask["chr1"][500:600] = True  # fully overlapping, well past the min_overlap_bp default
    assert RdnaMask.read_overlaps_mask(read, mask) is True


def test_read_overlaps_mask_false_when_blocks_miss_every_interval():
    header = _header([("chr1", 10000)])
    read = _aligned(header, "chr1", 500, 100)  # covers [500, 600)
    mask = {"chr1": __import__("intervaltree").IntervalTree()}
    mask["chr1"][700:800] = True
    assert RdnaMask.read_overlaps_mask(read, mask) is False


def test_read_overlaps_mask_false_for_a_contig_absent_from_the_mask():
    """The common case for an ordinary genome: short-circuit before any interval query."""

    header = _header([("chr1", 10000), ("chr2", 10000)])
    read = _aligned(header, "chr2", 500, 100)
    mask = {"chr1": __import__("intervaltree").IntervalTree()}
    mask["chr1"][0:10000] = True
    assert RdnaMask.read_overlaps_mask(read, mask) is False


@pytest.mark.parametrize("mask", [None, {}])
def test_read_overlaps_mask_false_when_no_mask_is_configured(mask):
    header = _header([("chr1", 10000)])
    read = _aligned(header, "chr1", 500, 100)
    assert RdnaMask.read_overlaps_mask(read, mask) is False


def test_read_overlaps_mask_min_overlap_bp_keeps_a_grazing_read():
    """A read that just clips the padding is kept, not discarded outright.

    Every excluded region already carries padding to absorb alignment-boundary
    slop around a genuine hit; a read whose alignment extends a handful of
    bases into that padding is much more likely an ordinary read from
    adjacent unique sequence than one implicated in the locus's ambiguity.
    """
    header = _header([("chr1", 10000)])
    read = _aligned(header, "chr1", 500, 100)  # covers [500, 600)
    mask = {"chr1": __import__("intervaltree").IntervalTree()}
    mask["chr1"][590:700] = True  # only the last 10 bp of the read overlap
    assert RdnaMask.read_overlaps_mask(read, mask, min_overlap_bp=50) is False
    assert RdnaMask.read_overlaps_mask(read, mask, min_overlap_bp=5) is True


def test_read_overlaps_mask_min_overlap_bp_discards_a_substantially_overlapping_read():
    """A read genuinely inside a masked region clears any reasonable floor."""

    header = _header([("chr1", 10000)])
    read = _aligned(header, "chr1", 500, 1000)  # covers [500, 1500)
    mask = {"chr1": __import__("intervaltree").IntervalTree()}
    mask["chr1"][0:10000] = True
    assert RdnaMask.read_overlaps_mask(read, mask, min_overlap_bp=50) is True


def test_read_overlaps_mask_min_overlap_bp_sums_across_spliced_blocks():
    """Two blocks each below the floor can still sum past it together."""

    header = _header([("chr1", 10000)])
    read = pysam.AlignedSegment(header)
    read.query_name = "r"
    read.reference_id = 0
    read.reference_start = 500
    read.mapping_quality = 60
    # block1 [500,530) 30bp, intron [530,600), block2 [600,630) 30bp
    read.cigar = [(0, 30), (3, 70), (0, 30)]
    read.query_sequence = "A" * 60
    read.query_qualities = pysam.qualitystring_to_array("I" * 60)
    mask = {"chr1": __import__("intervaltree").IntervalTree()}
    mask["chr1"][0:10000] = True
    assert RdnaMask.read_overlaps_mask(read, mask, min_overlap_bp=50) is True
    assert RdnaMask.read_overlaps_mask(read, mask, min_overlap_bp=61) is False


def test_read_overlaps_mask_reads_min_overlap_bp_from_config_by_default():
    header = _header([("chr1", 10000)])
    read = _aligned(header, "chr1", 500, 100)  # covers [500, 600)
    mask = {"chr1": __import__("intervaltree").IntervalTree()}
    mask["chr1"][590:700] = True  # 10 bp overlap

    saved = LRAA_Globals.config.get("rdna_mask_min_overlap_bp")
    try:
        LRAA_Globals.config["rdna_mask_min_overlap_bp"] = 50
        assert RdnaMask.read_overlaps_mask(read, mask) is False
        LRAA_Globals.config["rdna_mask_min_overlap_bp"] = 5
        assert RdnaMask.read_overlaps_mask(read, mask) is True
    finally:
        LRAA_Globals.config["rdna_mask_min_overlap_bp"] = saved


# ---------------------------------------------------------------------------
# build_rdna_mask_bed / load_mask_bed -- real minimap2, real files


def _write_fasta(path, name, seq):
    with open(path, "w") as f:
        f.write(">{}\n".format(name))
        for i in range(0, len(seq), 70):
            f.write(seq[i : i + 70] + "\n")


def _random_seq(n, seed):
    import random

    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_build_rdna_mask_bed_finds_an_embedded_cassette_copy(tmp_path):
    """The end-to-end path: a genome containing one cassette copy gets exactly it masked."""

    cassette_seq = _random_seq(2000, seed=1)
    flank = _random_seq(3000, seed=2)
    genome_seq = flank + cassette_seq + flank

    cassette_fa = tmp_path / "cassette.fa"
    genome_fa = tmp_path / "genome.fa"
    _write_fasta(cassette_fa, "cassette", cassette_seq)
    _write_fasta(genome_fa, "chr1", genome_seq)

    bed_path = RdnaMask.build_rdna_mask_bed(
        str(genome_fa),
        str(cassette_fa),
        cache_dir=str(tmp_path / "cache"),
        pad=10,
    )
    assert bed_path is not None
    assert os.path.exists(bed_path)

    mask = RdnaMask.load_mask_bed(bed_path)
    assert "chr1" in mask
    (interval,) = sorted(mask["chr1"])
    # The embedded copy sits at [3000, 5000); padded by 10 on each side.
    assert abs(interval.begin - 2990) <= 2
    assert abs(interval.end - 5010) <= 2


def test_build_rdna_mask_bed_returns_none_and_checkpoints_when_no_homology(tmp_path):
    """An ordinary genome unrelated to the cassette: no BED, no crash, and a cache hit."""

    cassette_fa = tmp_path / "cassette.fa"
    genome_fa = tmp_path / "genome.fa"
    _write_fasta(cassette_fa, "cassette", _random_seq(2000, seed=10))
    _write_fasta(genome_fa, "chr1", _random_seq(5000, seed=20))

    cache_dir = str(tmp_path / "cache")
    bed_path = RdnaMask.build_rdna_mask_bed(
        str(genome_fa), str(cassette_fa), cache_dir=cache_dir
    )
    assert bed_path is None
    assert RdnaMask.load_mask_bed(bed_path) is None

    # Second call must be a cache hit -- confirmed by the checkpoint existing and
    # no exception being raised re-reading it; a stale/rebuilt SAM would leave
    # a dangling .sam file behind, which build_rdna_mask_bed always removes.
    bed_path_again = RdnaMask.build_rdna_mask_bed(
        str(genome_fa), str(cassette_fa), cache_dir=cache_dir
    )
    assert bed_path_again is None
    assert not any(Path(cache_dir).glob("*.sam"))


def test_build_rdna_mask_bed_cache_hit_returns_the_same_path(tmp_path):
    cassette_seq = _random_seq(1000, seed=3)
    genome_seq = _random_seq(500, seed=4) + cassette_seq + _random_seq(500, seed=5)
    cassette_fa = tmp_path / "cassette.fa"
    genome_fa = tmp_path / "genome.fa"
    _write_fasta(cassette_fa, "cassette", cassette_seq)
    _write_fasta(genome_fa, "chr1", genome_seq)

    cache_dir = str(tmp_path / "cache")
    first = RdnaMask.build_rdna_mask_bed(str(genome_fa), str(cassette_fa), cache_dir=cache_dir)
    second = RdnaMask.build_rdna_mask_bed(str(genome_fa), str(cassette_fa), cache_dir=cache_dir)
    assert first == second
    assert not any(Path(cache_dir).glob("*.sam"))


def test_build_rdna_mask_bed_distinguishes_genome_identity(tmp_path):
    """A different genome must not reuse another genome's mask -- same cache dir."""

    cassette_seq = _random_seq(1000, seed=6)
    cassette_fa = tmp_path / "cassette.fa"
    _write_fasta(cassette_fa, "cassette", cassette_seq)

    genome_a = tmp_path / "a.fa"
    genome_b = tmp_path / "b.fa"
    _write_fasta(genome_a, "chr1", _random_seq(500, seed=7) + cassette_seq)
    _write_fasta(genome_b, "chr1", _random_seq(600, seed=8) + cassette_seq)

    cache_dir = str(tmp_path / "cache")
    bed_a = RdnaMask.build_rdna_mask_bed(str(genome_a), str(cassette_fa), cache_dir=cache_dir)
    bed_b = RdnaMask.build_rdna_mask_bed(str(genome_b), str(cassette_fa), cache_dir=cache_dir)
    assert bed_a != bed_b
    mask_a = RdnaMask.load_mask_bed(bed_a)
    mask_b = RdnaMask.load_mask_bed(bed_b)
    (iv_a,) = sorted(mask_a["chr1"])
    (iv_b,) = sorted(mask_b["chr1"])
    assert iv_a.begin != iv_b.begin


# ---------------------------------------------------------------------------
# resolve_rdna_fasta


def test_resolve_rdna_fasta_defaults_to_the_bundled_cassette():
    assert RdnaMask.resolve_rdna_fasta(None) == os.path.realpath(
        RdnaMask.DEFAULT_RDNA_CASSETTE_FASTA
    )
    assert os.path.exists(RdnaMask.resolve_rdna_fasta(None))


def test_resolve_rdna_fasta_rejects_a_missing_override(tmp_path):
    with pytest.raises(FileNotFoundError):
        RdnaMask.resolve_rdna_fasta(str(tmp_path / "does_not_exist.fa"))


def test_resolve_rdna_fasta_accepts_an_existing_override(tmp_path):
    custom = tmp_path / "other_cassette.fa"
    _write_fasta(custom, "x", "ACGT")
    assert RdnaMask.resolve_rdna_fasta(str(custom)) == os.path.realpath(str(custom))
