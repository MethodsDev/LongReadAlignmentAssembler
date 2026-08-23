#!/usr/bin/env python3

"""Tests for Util_funcs.alignment_per_id, the single identity definition.

Everything that filters on identity shares this function, so an error here does
not surface as an error -- it surfaces as reads quietly missing from a result.

The case these tests exist for: a CIGAR conventionally uses either M or the =/X
pair, but mixing them in one record is legal and minimap2 emits it. Measured on
XP132160.ucsc.bam over a 1-in-50 sample of 46,123,848 records, 0.32% of
alignments carry a few M ops among many =/X ops. Taking M alone as the aligned
base count made the denominator the M portion only, scoring `627=1M25=` with
NM:i:1 -- 653 aligned bases, 99.85% identity -- as 0.0%, so every consumer
applying a floor discarded it. Values ran to -800%.
"""

import pysam
import pytest

from pylib.Util_funcs import alignment_per_id, alignment_passes_per_id


@pytest.fixture
def header():
    return pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrT", "LN": 100000}]}
    )


def _read(header, cigar, nm, name="r"):
    aln = pysam.AlignedSegment(header)
    aln.query_name = name
    aln.reference_id = 0
    aln.reference_start = 100
    aln.mapping_quality = 60
    aln.cigarstring = cigar
    n_query = sum(
        length for op, length in aln.cigartuples if op in (0, 1, 4, 7, 8)
    )
    aln.query_sequence = "A" * n_query
    aln.query_qualities = pysam.qualitystring_to_array("I" * n_query)
    if nm is not None:
        aln.set_tag("NM", nm)
    return aln


def test_pure_m_cigar(header):
    """The conventional M-only form: 99 of 100 bases match."""
    assert alignment_per_id(_read(header, "100M", 1)) == pytest.approx(99.0)


def test_pure_eqx_cigar(header):
    """The --eqx form, where M never appears: denominator is = plus X."""
    assert alignment_per_id(_read(header, "95=5X", 5)) == pytest.approx(95.0)


def test_mixed_cigar_counts_all_aligned_bases(header):
    """The regression: a few M ops among many =/X must not become the denominator.

    `627=1M25=` is 653 aligned bases. With NM:i:1 that is 99.847% identity.
    Counting only the single M base yields 0.0% and the read is discarded by any
    floor -- which is what happened to 0.32% of a real corpus.
    """
    got = alignment_per_id(_read(header, "627=1M25=", 1))
    assert got == pytest.approx(100 - (1 / 653) * 100, abs=1e-9)
    assert got > 99.8, f"mixed-CIGAR identity collapsed to {got}"


def test_mixed_cigar_high_identity_reads_score_high(header):
    """Identity below zero is arithmetically impossible and was observed at -800%.

    It can only arise from a denominator smaller than the true aligned length.
    Cases drawn from the records actually measured on XP132160.ucsc.bam.
    """
    for cigar, nm in (
        ("7=1X36=1X23=1X73=1D31=4D627=1M604=", 9),
        ("9=1M1=1M5=2D3=1X3=1X5=1X6=1X66=1X742=133N61=", 9),
        ("93=1615N41=1X175=1125N35=1M166=2648N67=11107N163=1893N64=921=", 2),
    ):
        got = alignment_per_id(_read(header, cigar, nm))
        assert got is not None
        assert 0.0 <= got <= 100.0, f"{cigar} gave {got}"
        assert got > 98.0, f"{cigar} gave {got}, expected a high-identity read"


def test_indels_and_skips_are_not_aligned_bases(header):
    """Only M/=/X count. N (skip) and D (deletion) consume reference, not aligned bases.

    Pins the denominator's definition, which is what the regression got wrong: a
    spliced read must not have its intron length counted as aligned sequence.
    """
    spliced = _read(header, "50=1000N50=", 1)
    unspliced = _read(header, "100=", 1)
    assert alignment_per_id(spliced) == pytest.approx(alignment_per_id(unspliced))


def test_missing_mismatch_tag_returns_none(header):
    """No NM and no nM means undeterminable, and callers treat that as passing."""
    assert alignment_per_id(_read(header, "100=", None)) is None
    assert alignment_passes_per_id(_read(header, "100=", None), 80) is True


def test_nm_tag_fallback_to_nm_lowercase(header):
    """STAR-style nM is accepted when NM is absent."""
    aln = _read(header, "100M", None)
    aln.set_tag("nM", 2)
    assert alignment_per_id(aln) == pytest.approx(98.0)


def test_floor_gate_uses_the_same_definition(header):
    """The gate must agree with the metric, since the bug's damage was via the gate."""
    read = _read(header, "627=1M25=", 1)
    assert alignment_passes_per_id(read, 80) is True
    assert alignment_passes_per_id(read, 99.9) is False


@pytest.mark.parametrize(
    "cigar,nm",
    [
        ("100M", 1),
        ("95=5X", 5),
        ("627=1M25=", 1),
        ("50=1000N50=", 1),
        ("10M5I10M", 15),          # insertions dominate: NM charges all 5 I
        ("10M5D10M", 15),          # deletions dominate
        ("1M100I1M", 101),         # pathological: almost all inserted
        ("2=9X", 9),               # every aligned base a mismatch -> exactly 0
        ("5=100D5=", 100),
        ("7=1X36=1X23=1X73=1D31=4D627=1M604=", 9),
    ],
)
def test_identity_is_bounded_for_nm(header, cigar, nm):
    """0 <= identity <= 100 holds by arithmetic for NM, not by luck of the input.

    Provable only for NM, and only because the denominator is the exact set of
    columns NM charges edits against: mismatches <= M+=+X, insertions = I,
    deletions = D, so NM <= M+=+X+I+D. The same guarantee is NOT claimed for
    nM -- see the nM tests below.

    This is why the denominator is the set of columns NM charges edits against.
    Measured on PBMCs_pbio.aligned.sorted.bam -- 99.82% M-op records with ZERO
    mixed CIGARs -- 0.032% of alignments still produced negative identity purely
    from indel-heavy alignments, because I and D were excluded from a denominator
    while NM counted them. Fixing the mixed-CIGAR case alone did not address it.
    """
    got = alignment_per_id(_read(header, cigar, nm))
    assert got is not None
    assert 0.0 <= got <= 100.0, f"{cigar} NM:i:{nm} gave {got}"


def test_every_aligned_base_mismatched_is_exactly_zero(header):
    """The lower bound is reachable, so the clamp is not hiding a broken formula."""
    assert alignment_per_id(_read(header, "10X", 10)) == pytest.approx(0.0)


def test_rdna_mask_identity_agrees_with_shared_definition(header):
    """The one duplicate that cannot be removed must not be allowed to drift.

    RdnaMask parses raw minimap2 SAM before any pysam object exists, and
    Util_funcs imports RdnaMask, so it cannot call back without a circular
    import. It therefore reimplements the formula, and this test is the only
    thing keeping the two in agreement. It previously summed MDN=X -- counting
    introns, omitting insertions -- while its docstring claimed to match.
    """
    from pylib.RdnaMask import _sam_hit_identity

    for cigar, nm in (
        ("100M", 1),
        ("95=5X", 5),
        ("627=1M25=", 1),
        ("50=1000N50=", 1),      # N must not enter either denominator
        ("10M5I10M", 15),
        ("10M5D10M", 15),
        ("7=1X36=1X23=1X73=1D31=4D627=1M604=", 9),
    ):
        read = _read(header, cigar, nm)
        fields = ["q", "0", "chrT", "101", "60", cigar, "*", "0", "0",
                  read.query_sequence, "*", f"NM:i:{nm}"]
        assert _sam_hit_identity(fields) == pytest.approx(
            alignment_per_id(read), abs=1e-9
        ), f"drifted on {cigar}"


def test_rdna_mask_intron_does_not_inflate_identity():
    """A spliced hit must not score ~100% purely because its intron is long.

    The old MDN=X denominator added intron length, so this record scored
    99.98% and cleared any floor. It has 10 mismatches over 100 aligned
    columns and is a 90% identity alignment.
    """
    from pylib.RdnaMask import _sam_hit_identity

    cigar = "50=50X"
    long_intron = ["q", "0", "chrT", "101", "60", "50=50000N50X", "*", "0", "0",
                   "A" * 100, "*", "NM:i:50"]
    no_intron = ["q", "0", "chrT", "101", "60", cigar, "*", "0", "0",
                 "A" * 100, "*", "NM:i:50"]
    assert _sam_hit_identity(long_intron) == pytest.approx(50.0)
    assert _sam_hit_identity(no_intron) == pytest.approx(50.0)


def test_nm_lowercase_yields_the_same_gap_aware_identity(header):
    """STAR's nM counts mismatches only, so indels must not enter its denominator.

    nM is not NM. Charging a mismatch-only numerator against an indel-inclusive
    column count divides by more than the tag measures and overstates identity by
    the indel fraction. Here 2 mismatches over 20 aligned bases is 90%; counting
    the 5 inserted bases too would report 92%.
    """
    aln = _read(header, "10M5I10M", None)
    aln.set_tag("nM", 2)
    # 2 substitutions + 5 inserted = 7 edits over 25 alignment columns.
    assert alignment_per_id(aln) == pytest.approx(72.0)


def test_nm_and_nm_lowercase_reach_the_same_identity(header):
    """The metric must not depend on which aligner wrote the bam.

    NM:i:7 and nM:i:2 describe the SAME alignment of `10M5I10M` -- 2 substituted
    bases plus 5 inserted. Gap-aware identity is 7 edits over 25 columns for
    both. Reading nM raw over aligned bases alone would give 90.0, a
    gap-EXCLUDED identity under the same name, selected silently by producer.
    """
    with_nm = _read(header, "10M5I10M", 7)
    with_lower = _read(header, "10M5I10M", None)
    with_lower.set_tag("nM", 2)
    assert alignment_per_id(with_nm) == pytest.approx(72.0)
    assert alignment_per_id(with_lower) == pytest.approx(72.0)


def test_nm_uppercase_preferred_when_both_present(header):
    """NM is the standard tag and wins, so a bam carrying both is read as the
    standard defines it rather than by tag order."""
    aln = _read(header, "10M5I10M", 7)
    aln.set_tag("nM", 99)
    assert alignment_per_id(aln) == pytest.approx(72.0)


@pytest.mark.parametrize(
    "cigar,nm_lower",
    [("100M", 1), ("95=5X", 5), ("627=1M25=", 1), ("10M5I10M", 2), ("20X", 20)],
)
def test_unpaired_nm_lowercase_is_bounded(header, cigar, nm_lower):
    """For an unpaired record nM cannot exceed its own aligned bases, so it is bounded.

    Stated narrowly on purpose. STAR documents nM as mismatches per PAIRED read,
    summed across both mates, while a CIGAR describes one mate -- so for a proper
    pair the numerator can exceed this record's aligned bases and the result can
    go negative. That is a property of the tag's definition, not something one
    record can resolve, and it is documented at the call site rather than clamped.
    """
    aln = _read(header, cigar, None)
    aln.set_tag("nM", nm_lower)
    assert aln.is_paired is False
    got = alignment_per_id(aln)
    assert 0.0 <= got <= 100.0, f"{cigar} nM:i:{nm_lower} gave {got}"


def test_substitution_count_does_not_subtract_indels_from_nm(header):
    """nM already excludes indels, so subtracting them destroys real substitutions.

    `10M5I10M` with nM:i:2 has two substituted bases. NM-style subtraction gives
    max(0, 2 - 5 - 0) = 0 and scores a wrong alignment as perfect. Under NM the
    same integer legitimately means 2 edits of which 5 are insertions, so zero
    substitutions is the right answer there -- the tag has to pick the arithmetic.
    """
    from pylib.Util_funcs import substitution_count

    with_nm_lower = _read(header, "10M5I10M", None)
    with_nm_lower.set_tag("nM", 2)
    assert substitution_count(with_nm_lower) == 2

    with_nm = _read(header, "10M5I10M", 2)
    assert substitution_count(with_nm) == 0

    assert substitution_count(_read(header, "100M", None)) is None


def test_rescue_scoring_uses_nm_substitutions_correctly(header):
    """The rescue scorers must see 2 substitutions for nM:i:2, not 0.

    These decide whether a read is reassigned to a transcriptome alignment, so a
    read whose substitutions vanish scores as a perfect match and can win a
    comparison it should lose.
    """
    import sys
    from pathlib import Path as _P

    sys.path.insert(0, str(_P(__file__).resolve().parents[1] / "pylib"))
    import IsoformReadRescue as R

    read = _read(header, "10M5I10M", None)
    read.set_tag("nM", 2)
    # 20 aligned bases, 2 substituted -> 18 explained; span 25 -> 0.72
    assert R._explained_read_bases(read) == 18
    assert R._gap_aware_identity(read) == pytest.approx(18 / 25)

    # Under NM the same integer means 2 edits, 5 of them insertions -> 0 subs.
    nm_read = _read(header, "10M5I10M", 2)
    assert R._explained_read_bases(nm_read) == 20
    assert R._gap_aware_identity(nm_read) == pytest.approx(20 / 25)


def test_alignment_score_charges_indels_under_either_tag(header):
    """Equivalent NM and nM records must score identically.

    A score should be reduced by an indel, which is the opposite of what the
    identity metrics need. Reading the raw tag made the score blind to indels
    under nM, so two candidates differing only in indel content tied -- the exact
    comparison this decides. `10M5I10M` with nM:i:2 is the same alignment as
    NM:i:7 (2 substitutions + 5 inserted bases), so both must score 20 - 7 = 13.
    """
    import sys
    from pathlib import Path as _P

    sys.path.insert(0, str(_P(__file__).resolve().parents[1] / "pylib"))
    import IsoformReadRescue as R

    with_nm = _read(header, "10M5I10M", 7)
    with_nm_lower = _read(header, "10M5I10M", None)
    with_nm_lower.set_tag("nM", 2)
    assert R._alignment_score(with_nm) == 13
    assert R._alignment_score(with_nm_lower) == 13, "nM must be charged for indels too"

    # A candidate with the same substitutions but no indels must score HIGHER,
    # which is the discrimination the raw-tag version lost.
    clean = _read(header, "20M", None)
    clean.set_tag("nM", 2)
    assert R._alignment_score(clean) == 18
    assert R._alignment_score(clean) > R._alignment_score(with_nm_lower)

    with_as = _read(header, "10M5I10M", 7)
    with_as.set_tag("AS", 123)
    assert R._alignment_score(with_as) == 123

    from pylib.Util_funcs import alignment_edit_count

    assert alignment_edit_count(with_nm) == 7
    assert alignment_edit_count(with_nm_lower) == 7
    assert alignment_edit_count(_read(header, "100M", None)) is None


def test_inconsistent_nm_and_cigar_resolves_permissively(header):
    """NM < I+D means the tag and CIGAR disagree; that must not fabricate a count.

    NM counts every inserted and deleted base, so NM >= I+D holds for any
    self-consistent record. Below that the record is internally inconsistent, and
    the choice here is the permissive reading -- zero substitutions -- matching
    how a missing tag is treated. Pinned so the behaviour is a decision rather
    than an accident of clamping.
    """
    from pylib.Util_funcs import substitution_count

    # 5 inserted bases but NM:i:1 -- impossible for a consistent record.
    bad = _read(header, "10M5I10M", 1)
    assert substitution_count(bad) == 0
    # Identity stays defined and bounded rather than going negative.
    got = alignment_per_id(bad)
    assert got is not None and 0.0 <= got <= 100.0

    # The consistent case is unaffected: NM:i:7 with 5 insertions -> 2 subs.
    assert substitution_count(_read(header, "10M5I10M", 7)) == 2


def test_recorded_measurements_match_the_measured_files(header):
    """The worked example in the comment must stay true to the corpus it cites.

    `627=1M25=` with NM:i:1 is the record that exposed the mixed-CIGAR defect on
    XP132160.ucsc.bam. If this drifts, the comment's numbers are describing
    something the code no longer does.
    """
    read = _read(header, "627=1M25=", 1)
    # 653 aligned bases, no indels, so aligned columns == aligned bases.
    from pylib.Util_funcs import aligned_base_count

    assert aligned_base_count(read) == 653
    assert alignment_per_id(read) == pytest.approx(100 - (1 / 653) * 100)
    assert alignment_per_id(read) > 99.8


@pytest.mark.parametrize(
    "cigar,kwargs",
    [
        ("100M", {"nm": 1}),
        ("95=5X", {"nm": 5}),
        ("627=1M25=", {"nm": 1}),
        ("10M5I10M", {"nm": 7}),
        ("10M5I10M", {"nm_lower": 2}),
        ("10M5D10M", {"nm": 7}),
        ("10M5D10M", {"nm_lower": 2}),
        ("50=1000N50=", {"nm": 1}),
        ("5=100D5=", {"nm": 100}),
        ("20X", {"nm": 20}),
    ],
)
def test_matches_the_gap_aware_definition_and_the_other_implementation(header, cigar, kwargs):
    """Pins the metric AND ties two modules' definitions together.

    The intended metric is gap-aware percent identity:

        100 - (mismatches + inserted + deleted) / (M + = + X + I + D) * 100

    IsoformReadRescue._gap_aware_identity() computes the same quantity by the
    algebraically equivalent route of matched read bases over span, so the two
    must agree exactly on every record. They are in different modules with
    different call sites, which is how they drifted before; this equality is what
    would catch it happening again.
    """
    import sys
    from pathlib import Path as _P

    sys.path.insert(0, str(_P(__file__).resolve().parents[1] / "pylib"))
    import IsoformReadRescue as R
    from pylib.Util_funcs import alignment_edit_count

    read = _read(header, cigar, kwargs.get("nm"))
    if "nm_lower" in kwargs:
        read.set_tag("nM", kwargs["nm_lower"])

    stats = read.get_cigar_stats()[0]
    length = stats[0] + stats[7] + stats[8] + stats[1] + stats[2]
    expected = 100 - (alignment_edit_count(read) / length) * 100

    assert alignment_per_id(read) == pytest.approx(expected, abs=1e-9)
    assert alignment_per_id(read) == pytest.approx(
        R._gap_aware_identity(read) * 100, abs=1e-9
    ), "Util_funcs and IsoformReadRescue have drifted apart"
