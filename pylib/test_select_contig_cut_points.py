#!/usr/bin/env python3

"""Tests for util/misc/select_contig_cut_points.py.

The selector places one cut per target, where a target sits at a multiple of the
segment span. A position is legal only if it ends a depth window and lies outside
every annotated locus; among legal positions the one the fewest retained primary
alignments span wins, ties going to proximity then to the lower coordinate. These
tests hold it to that, to the joint constraints that couple neighbouring targets,
and to the accounting rule that nothing is dropped without being counted and named.
"""

import importlib.util
import os
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

import pysam
import pytest


def _load(name, relative):
    path = Path(__file__).resolve().parents[1] / relative
    loader = SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


selector = _load("select_contig_cut_points_under_test", "util/misc/select_contig_cut_points.py")
extractor = selector.extractor


class Fixture:
    """Synthetic contig with an indexed FASTA, BAM and GTF."""

    def __init__(self, tmp_path, contig="chrT", length=30000):
        self.tmp_path = tmp_path
        self.contig = contig
        self.length = length
        self.reads = []
        self.transcripts = []

    def add_read(self, name, strand, blocks, flag=None):
        self.reads.append((name, strand, [tuple(b) for b in blocks], flag))
        return self

    def spanning_read(self, name, position, strand="+", reach=40):
        """A read that straddles ``position``, so a cut there severs it."""
        return self.add_read(name, strand, [(position - reach, position + reach)])

    def add_transcript(self, gene_id, transcript_id, strand, exons):
        self.transcripts.append(
            (gene_id, transcript_id, strand, [tuple(e) for e in exons])
        )
        return self

    def build(self):
        self.fasta = str(self.tmp_path / "genome.fa")
        sequence = ("ACGT" * (self.length // 4 + 1))[: self.length]
        with open(self.fasta, "wt") as ofh:
            print(">{}".format(self.contig), file=ofh)
            for i in range(0, self.length, 60):
                print(sequence[i : i + 60], file=ofh)
        pysam.faidx(self.fasta)

        self.bam = str(self.tmp_path / "reads.bam")
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": self.contig, "LN": self.length}],
        }
        records = [self._alignment(*read) for read in self.reads]
        records.sort(key=lambda aln: aln.reference_start)
        with pysam.AlignmentFile(self.bam, "wb", header=header) as ofh:
            for record in records:
                ofh.write(record)
        pysam.index(self.bam)

        self.gtf = str(self.tmp_path / "annot.gtf") if self.transcripts else None
        if self.gtf:
            with open(self.gtf, "wt") as ofh:
                for gene_id, transcript_id, strand, exons in self.transcripts:
                    lend = min(e[0] for e in exons)
                    rend = max(e[1] for e in exons)
                    for feature, attrs in (
                        ("gene", 'gene_id "{}";'.format(gene_id)),
                        (
                            "transcript",
                            'gene_id "{}"; transcript_id "{}";'.format(
                                gene_id, transcript_id
                            ),
                        ),
                    ):
                        print(
                            "\t".join(
                                (
                                    self.contig,
                                    "test",
                                    feature,
                                    str(lend),
                                    str(rend),
                                    ".",
                                    strand,
                                    ".",
                                    attrs,
                                )
                            ),
                            file=ofh,
                        )
                    tx_attrs = 'gene_id "{}"; transcript_id "{}";'.format(
                        gene_id, transcript_id
                    )
                    for number, (exon_lend, exon_rend) in enumerate(sorted(exons), 1):
                        print(
                            "\t".join(
                                (
                                    self.contig,
                                    "test",
                                    "exon",
                                    str(exon_lend),
                                    str(exon_rend),
                                    ".",
                                    strand,
                                    ".",
                                    tx_attrs + " exon_number {};".format(number),
                                )
                            ),
                            file=ofh,
                        )
        return self

    def _alignment(self, name, strand, blocks, flag):
        aln = pysam.AlignedSegment()
        aln.query_name = name
        aln.flag = (16 if strand == "-" else 0) if flag is None else flag
        aln.reference_id = 0
        aln.reference_start = blocks[0][0] - 1
        aln.mapping_quality = 60
        cigar = []
        for i, (lend, rend) in enumerate(blocks):
            if i:
                cigar.append((3, lend - blocks[i - 1][1] - 1))
            cigar.append((0, rend - lend + 1))
        aln.cigar = cigar
        query_length = sum(rend - lend + 1 for lend, rend in blocks)
        aln.query_sequence = "A" * query_length
        aln.query_qualities = pysam.qualitystring_to_array("I" * query_length)
        aln.set_tag("NM", 0)
        return aln


def _tag_mismatches(bam_path, mismatches_by_name):
    """Rewrite a built fixture bam with per-read NM tags.

    Fixture alignments are all NM 0, so nothing in them can fail an identity
    floor; a test that emission stays stricter than the tally needs one record
    that does.
    """
    with pysam.AlignmentFile(bam_path, "rb") as reader:
        header = reader.header
        records = list(reader.fetch(until_eof=True))
    for aln in records:
        if aln.query_name in mismatches_by_name:
            aln.set_tag("NM", mismatches_by_name[aln.query_name])
    with pysam.AlignmentFile(bam_path, "wb", header=header) as writer:
        for aln in records:
            writer.write(aln)
    pysam.index(bam_path)


def _select(fixture, **kwargs):
    kwargs.setdefault("strand", "+")
    kwargs.setdefault("depth_window", 100)
    kwargs.setdefault("margin", 0)
    return selector.select_cut_points(
        bam_filename=fixture.bam,
        chrom=fixture.contig,
        contig_length=fixture.length,
        gtf=fixture.gtf,
        **kwargs
    )


# -- grid alignment -------------------------------------------------------------


def test_chosen_positions_end_a_depth_window_relative_to_the_grid_origin():
    """The condition is ``(b - origin) % depth_window == 0``, not ``b % 100 == 0``.

    Asserted at a grid origin that is NOT itself a multiple of the depth window,
    and the test additionally requires the chosen positions to FAIL the
    origin-free form. Otherwise the assertion would hold for the wrong reason
    whenever the origin happened to be grid-aligned, and a selector that ignored
    the origin entirely would pass.
    """

    positions = selector.grid_positions(1, 1000, 100, 37)
    assert positions
    assert all((p - 37) % 100 == 0 for p in positions)
    assert not any(p % 100 == 0 for p in positions)


def test_cuts_land_on_the_grid_at_a_non_aligned_origin(tmp_path):
    fixture = Fixture(tmp_path, length=6000).build()
    selection = _select(
        fixture, segment_span=2000, wiggle=1000, grid_origin=37, minimum_span=500
    )

    assert selection.cuts
    for cut in selection.cuts:
        assert (cut.position - 37) % 100 == 0, cut
        assert cut.position % 100 != 0, cut


def test_a_cut_is_the_last_base_of_a_window_so_no_window_spans_it(tmp_path):
    """``b`` ends a window and ``b + 1`` begins the next.

    The invariant the whole scheme protects: a window drawing bases from two
    chunks would give a different median, a different acceptance probability, and
    different reads kept. Checked here in the 0-based frame the normalizer's grid
    actually lives in.
    """

    fixture = Fixture(tmp_path, length=6000).build()
    selection = _select(fixture, segment_span=2000, wiggle=1000, minimum_span=500)

    depth_window = selection.depth_window
    for cut in selection.cuts:
        # left chunk ends at 1-based cut.position == 0-based cut.position - 1;
        # the right chunk begins at 1-based cut.position + 1 == 0-based position
        boundary_zero_based = cut.position
        assert boundary_zero_based % depth_window == 0
        # the window containing the left chunk's last base ends exactly there
        last_base = cut.position - 1  # 0-based
        assert (last_base // depth_window + 1) * depth_window == boundary_zero_based


# -- the annotation hard constraint --------------------------------------------


def test_no_cut_falls_inside_an_annotated_isoform(tmp_path):
    fixture = Fixture(tmp_path, length=9000)
    # loci blanketing most of each wiggle window, leaving narrow legal gaps
    fixture.add_transcript("g1", "t1", "+", [(2600, 2950), (3050, 3400)])
    fixture.add_transcript("g2", "t2", "+", [(5600, 5950), (6050, 6400)])
    fixture.build()

    selection = _select(fixture, segment_span=3000, wiggle=1000, minimum_span=1000)
    annotation = extractor.load_gtf(fixture.gtf, fixture.contig, "+")

    assert selection.cuts
    assert not selection.unplaced
    for cut in selection.cuts:
        for transcript in annotation.transcripts.values():
            assert not (
                transcript.lend <= cut.position < transcript.rend
            ), "cut {} is inside {}".format(cut.position, transcript.transcript_id)
        # and the stronger form the extractor actually enforces: not inside the
        # GENE span either, since a gene is emitted whole or not at all
        for gene in annotation.genes.values():
            assert not (gene.lend <= cut.position < gene.rend)


def test_gene_span_blocks_even_where_no_isoform_does(tmp_path):
    """A position between two isoforms of one gene is still refused.

    ``genes_contained`` emits a gene whole or not at all, so cutting a gene
    between its isoforms would drop the gene from BOTH neighbours -- the silent
    loss the hard constraint exists to prevent. The gene span is therefore the
    operative unit, strictly stronger than the per-isoform rule.
    """

    fixture = Fixture(tmp_path, length=9000)
    # one gene, two isoforms, with a gap at 3000-3400 that no isoform covers
    fixture.add_transcript("g1", "t1a", "+", [(2500, 2999)])
    fixture.add_transcript("g1", "t1b", "+", [(3401, 3900)])
    fixture.build()

    annotation = extractor.load_gtf(fixture.gtf, fixture.contig, "+")
    islands = extractor.find_islands(annotation, fixture.contig, "+", margin=0)
    zones = extractor.cut_zones(islands, 1, fixture.length)

    # 3100 is inside no isoform, but inside the gene, so it is not a legal cut
    assert not any(lo <= 3100 <= hi for lo, hi in zones)
    assert all(
        not (t.lend <= 3100 < t.rend) for t in annotation.transcripts.values()
    )

    selection = _select(fixture, segment_span=3000, wiggle=1000, minimum_span=1000)
    for cut in selection.cuts:
        assert not (2500 <= cut.position < 3900)


def test_a_target_with_no_compliant_position_is_reported_not_skipped(tmp_path):
    fixture = Fixture(tmp_path, length=9000)
    # blanket the entire first wiggle window [2500, 3500] with one locus
    fixture.add_transcript("gBlanket", "tBlanket", "+", [(2400, 3600)])
    fixture.build()

    selection = _select(fixture, segment_span=3000, wiggle=1000, minimum_span=1000)

    assert [item.target for item in selection.unplaced] == [3000]
    unplaced = selection.unplaced[0]
    assert unplaced.grid_positions > 0
    assert unplaced.annotation_blocked == unplaced.grid_positions
    assert "inside an annotated locus" in unplaced.reason
    assert 3000 not in [cut.target for cut in selection.cuts]

    payload = selector.selection_to_dict(selection)
    assert payload["counts"]["targets_unplaced"] == 1
    assert payload["unplaced_targets"][0]["target"] == 3000
    assert "UNPLACED" in selector.format_report(selection)


# -- the soft objective --------------------------------------------------------


def test_the_minimum_spanning_compliant_position_wins(tmp_path):
    """A farther position with fewer spanning reads beats a nearer crowded one."""

    # 5000 bp at a 3000 bp span puts exactly one target, at 3000
    fixture = Fixture(tmp_path, length=5000)
    # target 3000 is on-grid and nearest, but four reads span it
    for i in range(4):
        fixture.spanning_read("crowd{}".format(i), 3000)
    # 2900 and 3100 carry one read each; everything else in the window is free
    fixture.spanning_read("one_left", 2900)
    fixture.spanning_read("one_right", 3100)
    fixture.build()

    selection = _select(fixture, segment_span=3000, wiggle=1000, minimum_span=1000)

    assert len(selection.cuts) == 1
    cut = selection.cuts[0]
    assert cut.spanning_dropped == 0
    assert cut.position not in (2900, 3000, 3100)
    # nearest zero-cost position on either side, lower coordinate on a tie
    assert cut.position == 2800


def test_ties_break_by_proximity_then_by_lower_coordinate(tmp_path):
    """No reads anywhere, so every compliant position costs zero."""

    fixture = Fixture(tmp_path, length=5000).build()

    # target 2900 with a 250 bp window admits 2900 exactly -> distance 0
    selection = _select(fixture, segment_span=2900, wiggle=250, minimum_span=1000)
    assert [cut.position for cut in selection.cuts] == [2900]
    assert selection.cuts[0].offset == 0

    # target 2950 with a 250 bp window admits 2900 and 3000, equidistant at 50:
    # the lower coordinate must win, deterministically
    selection = _select(fixture, segment_span=2950, wiggle=250, minimum_span=1000)
    assert [cut.position for cut in selection.cuts] == [2900]
    assert selection.cuts[0].offset == -50


def test_the_window_is_not_widened_to_reach_a_zero_cost_position(tmp_path):
    """Inside the window the minimum is taken, even when it is not zero."""

    fixture = Fixture(tmp_path, length=5000)
    # every grid position in the window [2900, 3100] is spanned: 3000 by one
    # read, 2900 and 3100 by three each. Zero-cost positions exist outside it.
    for position in (2900, 3100):
        for i in range(3):
            fixture.spanning_read("crowd_{}_{}".format(position, i), position)
    fixture.spanning_read("single", 3000)
    fixture.build()

    selection = _select(fixture, segment_span=3000, wiggle=200, minimum_span=1000)

    assert len(selection.cuts) == 1
    cut = selection.cuts[0]
    assert cut.position == 3000
    assert cut.spanning_dropped == 1
    assert cut.window_lend == 2900 and cut.window_rend == 3100
    # not widened: a free position sits just outside the window, and was refused
    with pysam.AlignmentFile(fixture.bam, "rb") as bam:
        outside = selector.spanning_counts(
            bam, fixture.contig, "+", [cut.window_lend - 100, cut.window_rend + 100]
        )
    assert outside == [0, 0]
    assert cut.spanning_dropped > min(outside)


def test_a_compromise_forced_by_the_annotation_is_flagged(tmp_path):
    """``compromised`` means the hard constraint cost reads, and is quantified."""

    fixture = Fixture(tmp_path, length=5000)
    # the window is exactly [2900, 3100]. 3000 is free of reads but sits inside a
    # locus; 2900 and 3100 are legal but each carry a read, so obeying the
    # annotation costs one dropped read that ignoring it would not have cost.
    fixture.add_transcript("gAt3000", "tAt3000", "+", [(2950, 3050)])
    fixture.spanning_read("at2900", 2900)
    fixture.spanning_read("at3100", 3100)
    fixture.build()

    selection = _select(fixture, segment_span=3000, wiggle=200, minimum_span=1000)

    assert len(selection.cuts) == 1
    cut = selection.cuts[0]
    assert cut.unconstrained_best_spanning == 0
    assert cut.spanning_dropped == 1
    assert cut.compromised is True
    assert cut.annotation_blocked == 1
    assert cut.grid_positions == 3


def test_no_compromise_is_flagged_when_the_annotation_costs_nothing(tmp_path):
    fixture = Fixture(tmp_path, length=5000)
    fixture.spanning_read("at3000", 3000)
    fixture.build()

    selection = _select(fixture, segment_span=3000, wiggle=400, minimum_span=1000)
    assert len(selection.cuts) == 1
    assert selection.cuts[0].spanning_dropped == 0
    assert selection.cuts[0].compromised is False


# -- joint selection -----------------------------------------------------------


def test_joint_selection_respects_the_minimum_span_across_overlapping_windows(
    tmp_path,
):
    """The independently cheapest pair is too close together, so it is refused.

    Windows [500, 3500] and [2500, 5500] overlap. Only 2600 and 3400 are free of
    spanning reads, and per-target argmin would take 2600 for the target at 2000
    and 3400 for the target at 4000 -- 800 apart, a sliver well under the 1500
    floor. The joint program has to give up a read somewhere instead, which is
    precisely the coupling a per-target choice cannot see.
    """

    fixture = Fixture(tmp_path, length=6000)
    minimum_span = 1500
    free = (2600, 3400)
    for position in range(500, 5501, 100):
        if position in free:
            continue
        fixture.spanning_read("blk{}".format(position), position, reach=40)
    fixture.build()

    selection = _select(
        fixture,
        segment_span=2000,
        wiggle=3000,  # +/- 1500
        minimum_span=minimum_span,
    )

    positions = [cut.position for cut in selection.cuts]
    assert not selection.unplaced
    assert len(positions) == 2
    assert positions[0] < positions[1]
    assert positions[1] - positions[0] >= minimum_span
    # the pair each target would have chosen alone is exactly what got refused
    assert positions != list(free)
    # and the cost of obeying the floor is one read, counted rather than hidden
    assert selection.total_dropped == 1

    for segment in selection.segments:
        assert segment.span >= minimum_span, segment


@pytest.mark.parametrize(
    "length,segment_span,wiggle,minimum_span",
    [
        (6000, 2000, 3000, 1500),
        (6000, 1000, 1800, 800),
        (9000, 1500, 2800, 700),
        (12000, 2000, 4000, 1000),
    ],
)
def test_positions_stay_strictly_increasing_however_far_windows_overlap(
    tmp_path, length, segment_span, wiggle, minimum_span
):
    """Structural invariant, over parameter sets whose windows heavily overlap.

    Costs are made uneven so the objective genuinely pulls positions around,
    rather than every candidate tying and the order falling out of the scan.
    """

    fixture = Fixture(tmp_path, length=length)
    for index, position in enumerate(range(200, length - 200, 100)):
        for _ in range(index % 4):
            fixture.spanning_read(
                "r{}_{}".format(position, _), position, reach=40
            )
    fixture.build()

    selection = _select(
        fixture,
        segment_span=segment_span,
        wiggle=wiggle,
        minimum_span=minimum_span,
    )

    positions = [cut.position for cut in selection.cuts]
    assert positions == sorted(set(positions))
    for left, right in zip(positions, positions[1:]):
        assert right - left >= minimum_span
    for segment in selection.segments:
        assert segment.span >= minimum_span, segment
    assert selection.segments[0].lend == 1
    assert selection.segments[-1].rend == length
    # placed + unplaced accounts for every target, with no silent loss
    assert len(selection.cuts) + len(selection.unplaced) == len(selection.targets)


def test_every_realised_segment_meets_the_minimum_span(tmp_path):
    fixture = Fixture(tmp_path, length=9500).build()
    selection = _select(fixture, segment_span=2000, wiggle=1000)

    assert selection.minimum_span == selector.minimum_span_for(2000, 100)
    for segment in selection.segments:
        assert segment.span >= selection.minimum_span, segment
    # segments tile the contig with no gap and no overlap
    assert selection.segments[0].lend == 1
    assert selection.segments[-1].rend == fixture.length
    for left, right in zip(selection.segments, selection.segments[1:]):
        assert right.lend == left.rend + 1


def test_a_trailing_sliver_is_merged_into_its_predecessor_and_reported(tmp_path):
    """The anti-sliver floor applies at the contig end too, and is not silent."""

    # 5200 with a 2000 bp span puts a target at 4000 whose residual is 1200,
    # below the 1000 floor? no -- pick a length whose tail is genuinely short
    fixture = Fixture(tmp_path, length=5500).build()
    selection = _select(fixture, segment_span=2000, wiggle=400, minimum_span=1000)

    # targets would be 2000 and 4000; 5500 - 4000 = 1500 >= 1000, so nothing merges
    assert selection.tail_merged == []

    fixture = Fixture(tmp_path, length=4600).build()
    selection = _select(fixture, segment_span=2000, wiggle=400, minimum_span=1000)
    # targets 2000 and 4000; 4600 - 4000 = 600 < 1000, so 4000 merges away
    assert selection.tail_merged == [4000]
    assert [cut.target for cut in selection.cuts] == [2000]
    assert len(selection.segments) == 2
    report = selector.format_report(selection)
    assert "tail-merged" in report
    assert selector.selection_to_dict(selection)["counts"]["targets_tail_merged"] == 1


def test_chr20_scale_target_layout_is_span_based_with_no_cap():
    """Sizing follows contig length alone: chr20 -> 6 chunks, chr1 -> 25."""

    span = 10 * selector.MB
    floor = selector.minimum_span_for(span, 100)

    chr20_targets, chr20_merged = selector.cut_targets(64444167, span, floor)
    assert len(chr20_targets) == 5
    assert chr20_merged == [60000000]  # tail of 4.44 Mb is below the 5 Mb floor

    chr1_targets, _ = selector.cut_targets(248956422, span, floor)
    assert len(chr1_targets) + 1 == 25


# -- accounting ----------------------------------------------------------------


def test_dropped_read_names_are_emitted_and_match_the_per_cut_tally(tmp_path):
    fixture = Fixture(tmp_path, length=5000)
    # force the cut to 3000 by blanketing every other grid position in the window
    for position in (2900, 3100):
        for i in range(5):
            fixture.spanning_read("blk_{}_{}".format(position, i), position)
    fixture.spanning_read("severed_a", 3000)
    fixture.spanning_read("severed_b", 3000)
    fixture.build()

    selection = _select(fixture, segment_span=3000, wiggle=200, minimum_span=1000)

    assert len(selection.cuts) == 1
    cut = selection.cuts[0]
    assert cut.position == 3000
    assert cut.spanning_dropped == 2
    assert set(selection.dropped_read_names) == {"severed_a", "severed_b"}
    assert selection.total_dropped == cut.spanning_dropped

    names_path = tmp_path / "cuts.dropped_reads.txt"
    selector.write_dropped_read_names(selection, str(names_path))
    emitted = names_path.read_text().split()
    assert emitted == sorted({"severed_a", "severed_b"})
    assert len(emitted) == cut.spanning_dropped

    detail_path = tmp_path / "cuts.dropped_reads.tsv"
    selector.write_dropped_read_detail(selection, str(detail_path))
    rows = [
        line.split("\t")
        for line in detail_path.read_text().splitlines()[1:]
        if line
    ]
    assert {row[0] for row in rows} == {"severed_a", "severed_b"}
    assert all(row[3] == "3000" for row in rows)


@pytest.mark.parametrize("strand", [None, ""], ids=["none", "empty"])
def test_the_severed_bam_holds_what_the_tally_reports_without_a_strand(
    tmp_path, strand
):
    """The count a cut is scored on and the records emitted for it must agree.

    They did not.  Cost scoring uses retained_for_extraction, which reads a falsy
    strand as "any orientation"; emission uses quant_discard_reason, which
    required an exact match and so returned wrong_strand for every record.  The
    chunked pipeline omits --strand deliberately -- its input is already
    orientation-pure -- and the CLI turned argparse's None into "", so every run
    reported reads dropped at its cuts and wrote an empty severed bam beside the
    report.  A zero there reads as "no read was severed", which is the one
    conclusion the artifact exists to support.

    Both falsy spellings are covered because the CLI produced one and the
    function signature documents the other.
    """

    fixture = Fixture(tmp_path, length=5000)
    for position in (2900, 3100):
        for i in range(5):
            fixture.spanning_read("blk_{}_{}".format(position, i), position)
    fixture.spanning_read("severed_a", 3000)
    fixture.spanning_read("severed_b", 3000)
    fixture.build()

    sink = []
    selection = _select(
        fixture,
        strand=strand,
        segment_span=3000,
        wiggle=200,
        minimum_span=1000,
        severed_sink=sink,
        min_per_id=0,
        min_mapping_quality=0,
    )

    assert selection.total_dropped == 2
    assert len(sink) == selection.total_dropped
    assert {aln.query_name for aln in sink} == {"severed_a", "severed_b"}


def test_emission_stays_stricter_than_the_tally_it_accompanies(tmp_path):
    """Falsy strand may not become a free pass: the thresholds still apply.

    The fix widens what an unset strand means, so the test that it did not widen
    anything else is the other half of it.  Scoring counts both reads; emission
    keeps only the one quantification would use.
    """

    fixture = Fixture(tmp_path, length=5000)
    for position in (2900, 3100):
        for i in range(5):
            fixture.spanning_read("blk_{}_{}".format(position, i), position)
    fixture.spanning_read("kept", 3000)
    fixture.spanning_read("too_noisy", 3000)
    fixture.build()

    _tag_mismatches(fixture.bam, {"kept": 0, "too_noisy": 40})

    sink = []
    selection = _select(
        fixture,
        strand=None,
        segment_span=3000,
        wiggle=200,
        minimum_span=1000,
        severed_sink=sink,
        min_per_id=97,
        min_mapping_quality=0,
    )

    assert selection.total_dropped == 2
    assert {aln.query_name for aln in sink} == {"kept"}


def test_every_count_carries_its_denominator(tmp_path):
    fixture = Fixture(tmp_path, length=9000)
    fixture.spanning_read("severed", 3000)
    for i in range(3):
        fixture.add_read("elsewhere{}".format(i), "+", [(500 + i * 100, 560 + i * 100)])
    fixture.add_read("other_strand", "-", [(2980, 3020)])
    fixture.build()

    selection = _select(fixture, segment_span=3000, wiggle=200, minimum_span=1000)
    payload = selector.selection_to_dict(selection)

    # denominator counts only this orientation's retained primary alignments
    assert selection.total_retained_primary == 4
    assert payload["counts"]["retained_primary_alignments"] == 4
    assert payload["counts"]["alignments_dropped_at_cuts"] <= 4
    assert payload["counts"]["targets"] == len(selection.targets)
    assert (
        payload["counts"]["cuts_placed"] + payload["counts"]["targets_unplaced"]
        == payload["counts"]["targets"]
    )
    assert payload["counts"]["segments"] == payload["counts"]["cuts_placed"] + 1
    assert "of 4 retained primary alignment(s)" in selector.format_report(selection)


def test_manifest_reports_realised_spans_and_offsets_from_target(tmp_path):
    fixture = Fixture(tmp_path, length=9000)
    fixture.spanning_read("at3000", 3000)
    fixture.build()

    selection = _select(fixture, segment_span=3000, wiggle=600, minimum_span=1000)
    payload = selector.selection_to_dict(selection)

    for entry, cut in zip(payload["cuts"], selection.cuts):
        assert entry["offset_from_target"] == cut.position - cut.target
        assert entry["position"] == cut.position
        assert entry["window"] == [cut.window_lend, cut.window_rend]
    spans = [entry["span"] for entry in payload["segments"]]
    assert sum(spans) == fixture.length
    assert all(span > 0 for span in spans)


def test_window_origin_for_each_segment_is_its_rebase_offset(tmp_path):
    """What the driver hands normalize_bam_by_strand.py --window_origin.

    Absolute 0-based == chunk-local 0-based + offset, and every boundary is
    grid-aligned, so every offset is a whole number of depth windows and the
    chunk grid coincides with the absolute grid.
    """

    fixture = Fixture(tmp_path, length=9000).build()
    selection = _select(fixture, segment_span=3000, wiggle=600, minimum_span=1000)
    payload = selector.selection_to_dict(selection)

    for entry, segment in zip(payload["segments"], selection.segments):
        assert entry["window_origin"] == segment.lend - 1
        assert entry["window_origin"] % selection.depth_window == 0


# -- the extractor honours what the selector chose -----------------------------


def test_a_selected_cut_is_accepted_by_the_extractor(tmp_path):
    """The shared hard constraint, checked end to end.

    The selector screens positions with ``find_islands``/``cut_zones`` and the
    extractor refuses a boundary with ``admissibility_offenders``. They must be
    the same rule, or the selector would hand out positions the extractor rejects.
    """

    fixture = Fixture(tmp_path, length=9000)
    fixture.add_transcript("gA", "tA", "+", [(2600, 2900)])
    fixture.add_transcript("gB", "tB", "+", [(3100, 3300)])
    fixture.add_transcript("gC", "tC", "+", [(5900, 6200)])
    fixture.spanning_read("severed", 3000)
    fixture.add_read("insideA", "+", [(2650, 2850)])
    fixture.add_read("insideC", "+", [(5950, 6150)])
    fixture.build()

    selection = _select(fixture, segment_span=3000, wiggle=1000, minimum_span=1000)
    assert selection.cuts and not selection.unplaced

    emitted = {}
    dropped = set()
    for segment in selection.segments:
        manifest = extractor.extract_partition(
            genome_fa=fixture.fasta,
            bam=fixture.bam,
            region="{}+:{}-{}".format(fixture.contig, segment.lend, segment.rend),
            output_prefix=str(tmp_path / "chunk{}".format(segment.index)),
            gtf=fixture.gtf,
            margin=0,
        )
        assert manifest["offset"] == segment.lend - 1
        assert manifest["window_origin"] == segment.lend - 1
        dropped.update(manifest["dropped_read_names"])
        with pysam.AlignmentFile(manifest["files"]["bam"], "rb") as bam:
            for aln in bam:
                assert aln.query_name not in emitted
                emitted[aln.query_name] = segment.index
                # in range on the mini contig, by construction
                assert aln.reference_start >= 0
                assert aln.reference_end <= manifest["mini_length"]

    # the selector predicted exactly the set the extractor dropped
    assert dropped == set(selection.dropped_read_names)
    # and every retained read is either emitted once or dropped, never neither
    source = {name for name, strand, _, _ in fixture.reads if strand == "+"}
    assert set(emitted) | dropped == source
    assert not (set(emitted) & dropped)


def test_selector_cli_writes_report_manifest_and_dropped_names(tmp_path, capsys):
    fixture = Fixture(tmp_path, length=5000)
    # blanket the neighbouring grid positions so the cut lands on 3000, where the
    # severed read is, making the dropped-name output non-empty and predictable
    for position in (2900, 3100):
        for i in range(3):
            fixture.spanning_read("blk_{}_{}".format(position, i), position)
    fixture.spanning_read("severed", 3000)
    # a locus well clear of the window, so --gtf has a real file to read
    fixture.add_transcript("gFar", "tFar", "+", [(600, 900)])
    fixture.build()

    prefix = str(tmp_path / "sel")
    rc = selector.main(
        [
            "--bam",
            fixture.bam,
            "--gtf",
            fixture.gtf,
            "--contig",
            fixture.contig,
            "--strand",
            "+",
            "--approx_MB_per_cut",
            "0.003",
            "--approx_MB_per_cut_wiggle_window",
            "0.0002",
            "--minimum_span",
            "1000",
            "--margin",
            "0",
            "--output_prefix",
            prefix,
        ]
    )
    assert rc == 0

    payload = json.loads(Path(prefix + ".cuts.json").read_text())
    assert len(payload) == 1
    assert payload[0]["chrom"] == fixture.contig
    assert payload[0]["counts"]["cuts_placed"] == 1
    assert payload[0]["cuts"][0]["position"] == 3000

    assert Path(prefix + ".cuts.txt").exists()
    assert Path(prefix + ".dropped_reads.txt").read_text().split() == ["severed"]
    assert payload[0]["counts"]["alignments_dropped_at_cuts"] == 1
    detail = Path(prefix + ".dropped_reads.tsv").read_text().splitlines()
    assert detail[0].startswith("read_name")
    assert detail[1].split("\t")[0] == "severed"


def test_a_read_severed_by_a_cut_is_reported_by_both_neighbours_but_named_once(
    tmp_path,
):
    """Dedup is a real requirement, not tidiness.

    A read straddling a boundary is contained by neither chunk, so per-chunk
    accounting names it twice. A filter set built from those names must contain it
    once, or the pruned baseline would be built from a malformed list.
    """

    fixture = Fixture(tmp_path, length=5000)
    for position in (2900, 3100):
        for i in range(3):
            fixture.spanning_read("blk_{}_{}".format(position, i), position)
    fixture.spanning_read("severed", 3000, reach=60)
    fixture.build()

    selection = _select(fixture, segment_span=3000, wiggle=200, minimum_span=1000)
    assert selection.cuts[0].position == 3000

    per_chunk = []
    for segment in selection.segments:
        manifest = extractor.extract_partition(
            genome_fa=fixture.fasta,
            bam=fixture.bam,
            region="{}+:{}-{}".format(fixture.contig, segment.lend, segment.rend),
            output_prefix=str(tmp_path / "c{}".format(segment.index)),
            gtf=fixture.gtf,
            margin=0,
        )
        per_chunk.extend(manifest["dropped_read_names"])

    assert per_chunk == ["severed", "severed"]

    path = tmp_path / "union.txt"
    selector.write_dropped_read_names(selection, str(path))
    assert path.read_text().split() == ["severed"]


def test_long_intron_alignments_do_not_influence_selection(tmp_path):
    """The pipeline discards them upstream, so counting them would misdirect cuts.

    A 200 kb intron spans a fifth of a default wiggle window; an alignment the
    strand split has already thrown away must not push a cut off the position it
    would otherwise take.
    """

    fixture = Fixture(tmp_path, length=9000)
    # spans 3000 via a single enormous intron
    fixture.add_read("long_intron", "+", [(2000, 2100), (8000, 8100)])
    fixture.build()

    tight = _select(
        fixture,
        segment_span=3000,
        wiggle=200,
        minimum_span=1000,
        max_intron_length=1000,
    )
    assert tight.cuts[0].spanning_dropped == 0
    assert tight.total_retained_primary == 0
    assert not tight.dropped_read_names

    loose = _select(
        fixture, segment_span=3000, wiggle=200, minimum_span=1000, max_intron_length=0
    )
    assert loose.total_retained_primary == 1
    assert loose.cuts[0].spanning_dropped == 1


# -- severed alignment emission -------------------------------------------------


def test_the_severed_reads_bam_preserves_junction_structure(tmp_path):
    """Spans cannot answer compatibility; exon blocks and junctions decide it.

    A consumer asking which severed reads linked two genes needs the block
    structure, not the extent: two reads with identical start and end can be
    compatible with different genes.  So the emission has to be the records.
    """

    fixture = Fixture(tmp_path, length=5000)
    # one spliced read across the only target, its intron away from the boundary
    fixture.add_read("spliced", "+", [(2600, 2900), (3300, 3600)])
    fixture.build()

    out_bam = tmp_path / "severed.bam"
    sink = []
    # wiggle 0 leaves the target itself as the only candidate, so the read cannot
    # be dodged and the cut is forced through it
    selection = _select(
        fixture,
        segment_span=3000,
        wiggle=0,
        minimum_span=1000,
        severed_sink=sink,
    )

    assert "spliced" in selection.dropped_read_names, (
        "fixture must actually sever the spliced read for this to test anything"
    )

    selector.write_severed_alignments_bam(
        pysam.AlignmentFile(fixture.bam, "rb").header, sink, str(out_bam)
    )

    with pysam.AlignmentFile(str(out_bam), "rb") as bam:
        emitted = list(bam.fetch())

    assert [aln.query_name for aln in emitted] == ["spliced"]
    # the junction survives: two blocks with an N between them, not one span
    assert emitted[0].get_blocks() == [(2599, 2900), (3299, 3600)]
    assert 3 in [op for op, _length in emitted[0].cigartuples]


def test_a_read_severed_by_two_cuts_is_emitted_once(tmp_path):
    """One record per alignment, or a consumer double-counts it."""

    fixture = Fixture(tmp_path, length=9000)
    # long read reaching across two adjacent targets
    fixture.add_read("long", "+", [(2500, 6500)])
    fixture.build()

    sink = []
    selection = _select(
        fixture,
        segment_span=3000,
        wiggle=500,
        minimum_span=1000,
        severed_sink=sink,
    )

    assert len(selection.dropped_read_names["long"]) == 2, (
        "fixture must sever the same read at two cuts"
    )
    assert [aln.query_name for aln in sink] == ["long"]


def test_no_sink_means_no_collection(tmp_path):
    """The default path must not pay for records nobody asked for."""

    fixture = Fixture(tmp_path, length=5000)
    fixture.spanning_read("severed", 3000)
    fixture.build()

    selection = _select(fixture, segment_span=3000, wiggle=0, minimum_span=1000)
    assert selection.dropped_read_names, "expected a severed read to be named"


def test_byte_identical_rows_are_both_emitted(tmp_path):
    """Multiplicity is data: two identical BAM rows are two severed alignments.

    Any dedupe keyed on record content collapses them -- however complete the key
    -- and the emitted BAM then stops representing every alignment the cuts sever.
    Repeats are collapsed by cut ownership instead, which is a property of
    coordinates and so cannot confuse two rows with one row seen twice.
    """

    bam_path = tmp_path / "identical.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrT", "LN": 5000}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        for _ in range(2):
            aln = pysam.AlignedSegment()
            aln.query_name = "shared_name"
            aln.flag = 0
            aln.reference_id = 0
            aln.reference_start = 2599
            aln.mapping_quality = 60
            aln.cigar = [(0, 1000)]
            aln.query_sequence = "A" * 1000
            aln.query_qualities = pysam.qualitystring_to_array("I" * 1000)
            aln.set_tag("NM", 0)
            out.write(aln)
    pysam.index(str(bam_path))

    sink = []
    selector.select_cut_points(
        bam_filename=str(bam_path),
        chrom="chrT",
        contig_length=5000,
        strand="+",
        segment_span=3000,
        wiggle=0,
        depth_window=100,
        margin=0,
        minimum_span=1000,
        severed_sink=sink,
    )

    assert len(sink) == 2, "two identical rows are two alignments, not one"


def test_an_alignment_spanning_two_cuts_is_emitted_once(tmp_path):
    """The repeat that ownership exists to collapse, in the same run as the above."""

    bam_path = tmp_path / "wide.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrT", "LN": 9000}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        aln = pysam.AlignedSegment()
        aln.query_name = "wide"
        aln.flag = 0
        aln.reference_id = 0
        aln.reference_start = 2499
        aln.mapping_quality = 60
        aln.cigar = [(0, 4000)]
        aln.query_sequence = "A" * 4000
        aln.query_qualities = pysam.qualitystring_to_array("I" * 4000)
        aln.set_tag("NM", 0)
        out.write(aln)
    pysam.index(str(bam_path))

    sink = []
    selection = selector.select_cut_points(
        bam_filename=str(bam_path),
        chrom="chrT",
        contig_length=9000,
        strand="+",
        segment_span=3000,
        wiggle=0,
        depth_window=100,
        margin=0,
        minimum_span=1000,
        severed_sink=sink,
    )

    assert len(selection.dropped_read_names["wide"]) == 2, "must span two cuts"
    assert len(sink) == 1, "one alignment, emitted at the first cut it spans"


def _one_read_bam(tmp_path, name, nm, mapq=60, length=1000, start=2599):
    """A single primary alignment whose NM tag sets its percent identity."""

    bam_path = tmp_path / "{}.bam".format(name)
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrT", "LN": 5000}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        aln = pysam.AlignedSegment()
        aln.query_name = name
        aln.flag = 0
        aln.reference_id = 0
        aln.reference_start = start
        aln.mapping_quality = mapq
        aln.cigar = [(0, length)]
        aln.query_sequence = "A" * length
        aln.query_qualities = pysam.qualitystring_to_array("I" * length)
        aln.set_tag("NM", nm)
        out.write(aln)
    pysam.index(str(bam_path))
    return bam_path


def _severed(bam_path, **kwargs):
    sink = []
    selection = selector.select_cut_points(
        bam_filename=str(bam_path),
        chrom="chrT",
        contig_length=5000,
        strand="+",
        segment_span=3000,
        wiggle=0,
        depth_window=100,
        margin=0,
        minimum_span=1000,
        severed_sink=sink,
        **kwargs
    )
    return selection, sink


def test_emission_applies_the_quant_percent_identity_threshold(tmp_path):
    """A read the quant step rejects must not be offered as severed evidence.

    Cut scoring deliberately uses a superset predicate omitting percent identity --
    over-counting a cut's cost only biases toward safer positions.  Emission cannot
    use it: a read rejected on identity never reaches a component, so offering it
    as evidence a boundary dissolved one is a false positive, not a conservative
    one.  NM=100 over 1000 aligned bases is 90% identity: kept at the default 80,
    rejected at the 97 that --HiFi sets.
    """

    bam_path = _one_read_bam(tmp_path, "pid90", nm=100)

    selection, sink = _severed(bam_path, min_per_id=80.0)
    assert "pid90" in selection.dropped_read_names
    assert [a.query_name for a in sink] == ["pid90"]

    selection, sink = _severed(bam_path, min_per_id=97.0)
    assert "pid90" in selection.dropped_read_names, (
        "cost accounting still counts it: the superset is what scored the cut"
    )
    assert sink == [], "quant would reject it at 97, so it must not be emitted"


def test_emission_applies_the_quant_mapping_quality_threshold(tmp_path):
    bam_path = _one_read_bam(tmp_path, "mq5", nm=0, mapq=5)

    _, sink = _severed(bam_path, min_mapping_quality=0)
    assert [a.query_name for a in sink] == ["mq5"]

    _, sink = _severed(bam_path, min_mapping_quality=30)
    assert sink == [], "quant would reject mapq 5, so it must not be emitted"


def test_the_emitted_bam_records_its_scope(tmp_path):
    """The artifact states what it is, because consumers have over-read it."""

    bam_path = _one_read_bam(tmp_path, "kept", nm=0)
    _, sink = _severed(bam_path, min_per_id=80.0)
    out = tmp_path / "severed.bam"
    with pysam.AlignmentFile(str(bam_path), "rb") as src:
        selector.write_severed_alignments_bam(
            src.header, sink, str(out), min_per_id=80.0, min_mapping_quality=0
        )

    with pysam.AlignmentFile(str(out), "rb") as bam:
        notes = "\n".join(bam.header.to_dict().get("CO", []))

    assert "PRE-normalization" in notes
    assert "min_per_id=80.0" in notes
    assert "NOT the quant-visible set" in notes
    assert "refuse rather than approximate" in notes



# -- strandless selection -------------------------------------------------------
#
# One set of cuts over the RAW bam, serving both strands, so the strand split can
# move inside the chunk. Nothing in the arithmetic is new -- a falsy strand
# already unions the islands and costs both orientations -- so these tests exist
# to CHECK that rather than assume it, and each one contrasts the strandless
# answer with the per-strand answer over the same fixture. A test that only
# asserted "strandless produces cuts" would pass against a selector that quietly
# ignored the other strand entirely.


def _union_fixture(tmp_path):
    """A locus per strand, each blocking a different position in one window.

    Window 2900-3100 at depth_window 100, so exactly three candidates. The minus
    locus blocks 2900 and 3000, the plus locus blocks 3200-3300 (outside the
    window, present only so the plus-filtered annotation is not empty).
    """

    fixture = Fixture(tmp_path, length=5000)
    fixture.add_transcript("gMinus", "tMinus", "-", [(2900, 3100)])
    fixture.add_transcript("gPlus", "tPlus", "+", [(3200, 3300)])
    fixture.add_read("fwd", "+", [(1000, 1200)])
    fixture.add_read("rev", "-", [(1400, 1600)])
    return fixture.build()


def test_strandless_islands_block_a_locus_on_either_strand(tmp_path):
    """Union blocking, contrasted against the per-strand answer.

    The minus locus covers the target itself. A plus-filtered selection never
    sees it and cuts at 3000; a strandless selection must clear it, because the
    chunk it bounds is the sole container of BOTH strands' loci over that
    interval and a locus it splits is lost from both neighbours exactly as a
    same-strand one would be.
    """

    fixture = _union_fixture(tmp_path)
    common = dict(segment_span=3000, wiggle=200, minimum_span=1000)

    plus = _select(fixture, strand="+", **common)
    assert [cut.position for cut in plus.cuts] == [3000]

    strandless = _select(fixture, strand="", strandless=True, **common)
    assert [cut.position for cut in strandless.cuts] == [3100]
    # and it is the annotation that moved it: the window held three candidates
    # and the minus locus removed two of them
    assert strandless.cuts[0].grid_positions == 3
    assert strandless.cuts[0].annotation_blocked == 2
    assert plus.cuts[0].annotation_blocked == 0


def test_strandless_cut_cost_counts_alignments_on_both_strands(tmp_path):
    """A shared boundary is priced by everything it severs, either orientation.

    Three minus-strand reads straddle the target and one plus-strand read
    straddles each neighbour. Filtered to plus, the target looks free and wins.
    Strandless, the target costs 3 against a neighbour's 1, so the selector pays
    the smaller price -- which is the correct objective, because the two strand
    chunks share this boundary and therefore share its cost.
    """

    fixture = Fixture(tmp_path, length=5000)
    fixture.spanning_read("plus_at_2900", 2900, strand="+")
    fixture.spanning_read("plus_at_3100", 3100, strand="+")
    for i in range(3):
        fixture.spanning_read("minus_at_3000_{}".format(i), 3000, strand="-")
    fixture.build()

    common = dict(segment_span=3000, wiggle=200, minimum_span=1000)

    plus = _select(fixture, strand="+", **common)
    assert plus.cuts[0].position == 3000
    assert plus.cuts[0].spanning_dropped == 0, "the minus reads are invisible to it"

    strandless = _select(fixture, strand="", strandless=True, **common)
    assert strandless.cuts[0].position == 2900
    assert strandless.cuts[0].spanning_dropped == 1
    # the position the plus run chose would have severed all three minus reads
    with pysam.AlignmentFile(fixture.bam, "rb") as bam:
        assert selector.spanning_counts(bam, fixture.contig, "", [3000]) == [3]


def test_strandless_selection_labels_its_segments_without_a_strand(tmp_path):
    """The region strings are consumed by the extractor, so they must parse.

    Before strandless selection existed the CLI threaded argparse's None into the
    format string and emitted "chrTNone:1-3000" -- a region no consumer could
    parse, unnoticed because the driver rebuilt the region from the strand it
    already knew. A strandless run has no such fallback: these strings ARE the
    interface.
    """

    fixture = _union_fixture(tmp_path)
    selection = _select(
        fixture,
        strand="",
        strandless=True,
        segment_span=3000,
        wiggle=200,
        minimum_span=1000,
    )
    payload = selector.selection_to_dict(selection)

    assert payload["strand"] is None
    assert payload["strandless"] is True
    regions = [segment["region"] for segment in payload["segments"]]
    assert regions == ["chrT:1-3100", "chrT:3101-5000"]
    for region in regions:
        parsed = extractor.parse_region(region)
        assert parsed.strand == ""
        assert parsed.chrom == fixture.contig
    assert "None" not in selector.format_report(selection)


def test_an_unset_strand_also_stops_emitting_the_word_None(tmp_path):
    """The same regression, on the path the chunked driver has always used.

    Omitting --strand is still "already strand-split, count everything"; it is
    not the strandless declaration. What it must not do is name a region after
    Python's None.
    """

    fixture = _union_fixture(tmp_path)
    selection = _select(
        fixture, strand=None, segment_span=3000, wiggle=200, minimum_span=1000
    )
    payload = selector.selection_to_dict(selection)

    assert payload["strandless"] is False
    assert all("None" not in s["region"] for s in payload["segments"])
    assert "None" not in selector.format_report(selection)


def test_strandless_reports_the_denominator_per_orientation(tmp_path):
    """Both orientations counted, and they sum to the denominator.

    This is the evidence that the input really was raw. A strandless selection
    over a strand-separated bam would place perfectly reasonable cuts and serve
    one strand, and the only thing that distinguishes the two cases is this
    split.
    """

    fixture = Fixture(tmp_path, length=5000)
    for i in range(3):
        fixture.add_read("fwd{}".format(i), "+", [(500 + i * 200, 600 + i * 200)])
    for i in range(2):
        fixture.add_read("rev{}".format(i), "-", [(1500 + i * 200, 1600 + i * 200)])
    fixture.build()

    selection = _select(
        fixture,
        strand="",
        strandless=True,
        segment_span=3000,
        wiggle=200,
        minimum_span=1000,
    )
    counts = selector.selection_to_dict(selection)["counts"]

    assert counts["retained_primary_forward"] == 3
    assert counts["retained_primary_reverse"] == 2
    assert counts["retained_primary_alignments"] == 5
    assert "3 forward, 2 reverse" in selector.format_report(selection)


def test_strandless_over_a_strand_separated_bam_is_reported_not_absorbed(
    tmp_path, capsys
):
    """One orientation means the caller passed a split bam. Say so.

    The cuts are still valid coordinates, so this is a report rather than a
    refusal -- but silence would look exactly like success, and the chunks would
    come out with an empty partner after the in-chunk split.
    """

    fixture = Fixture(tmp_path, length=5000)
    for i in range(3):
        fixture.add_read("fwd{}".format(i), "+", [(500 + i * 200, 600 + i * 200)])
    fixture.build()

    selection = _select(
        fixture,
        strand="",
        strandless=True,
        segment_span=3000,
        wiggle=200,
        minimum_span=1000,
    )

    assert selection.cuts, "the cuts are still placed"
    stderr = capsys.readouterr().err
    assert "only one orientation" in stderr
    assert "already strand-separated" in stderr


def test_strandless_and_strand_are_incompatible(tmp_path):
    """Costing a shared cut against one orientation prices half the reads."""

    fixture = _union_fixture(tmp_path)
    with pytest.raises(selector.SelectionError) as excinfo:
        _select(
            fixture,
            strand="+",
            strandless=True,
            segment_span=3000,
            wiggle=200,
            minimum_span=1000,
        )
    assert "cannot also filter by orientation" in str(excinfo.value)

    with pytest.raises(SystemExit):
        selector.main(
            [
                "--bam",
                fixture.bam,
                "--contig",
                fixture.contig,
                "--strandless",
                "--strand",
                "+",
                "--output_prefix",
                str(tmp_path / "never"),
            ]
        )


def test_strandless_cli_writes_a_manifest_the_extractor_can_consume(tmp_path):
    """End to end: the CLI's own output drives a strandless extraction."""

    fixture = _union_fixture(tmp_path)
    prefix = str(tmp_path / "sel")
    rc = selector.main(
        [
            "--bam",
            fixture.bam,
            "--gtf",
            fixture.gtf,
            "--contig",
            fixture.contig,
            "--strandless",
            "--approx_MB_per_cut",
            "0.003",
            "--approx_MB_per_cut_wiggle_window",
            "0.0002",
            "--minimum_span",
            "1000",
            "--margin",
            "0",
            "--output_prefix",
            prefix,
        ]
    )
    assert rc == 0

    payload = json.loads(Path(prefix + ".cuts.json").read_text())[0]
    assert payload["strandless"] is True
    assert payload["strand"] is None
    assert payload["cuts"][0]["position"] == 3100
    assert payload["counts"]["retained_primary_forward"] == 1
    assert payload["counts"]["retained_primary_reverse"] == 1

    for index, segment in enumerate(payload["segments"]):
        manifest = extractor.extract_partition(
            genome_fa=fixture.fasta,
            bam=fixture.bam,
            region=segment["region"],
            output_prefix=str(tmp_path / "chunk{}".format(index)),
            gtf=fixture.gtf,
            margin=0,
        )
        assert manifest["strand"] is None
        assert manifest["strand_split_required"] is True


# -- severing is a COST, never a veto ------------------------------------------
#
# An earlier revision made zero severed a HARD constraint in discovery: severing
# positions were struck from the candidate set and a target whose window held none
# was DECLINED. That contract was REJECTED, because at depth every base is covered
# by some read, so the rule would decline every cut and silently disable chunking
# exactly where it pays most. The tests below encode what replaced it: a severing
# position IS selectable, the run proceeds, and what it severs is counted, split by
# structure and named.
#
# The alignments here are SYNTHETIC on purpose, and the mono/spliced mix is chosen
# per test. On real corpora it cannot be chosen: measured at 2 Mb spacing on HG002
# PacBio Kinnex, 98-100% of severed alignments are spliced (chr21 20 kb: 14 mono
# against 926 spliced; chr1 200 kb: 2 against 85), because spanning probability
# scales with genomic span and a monoexonic read is a kb or two long. Real data
# therefore cannot exercise the weighting crossover at all, in either direction.


def _every_position_severed(tmp_path):
    """A window in which no grid position is free of spanning alignments.

    Target 3000 with a 200 bp window admits exactly 2900, 3000 and 3100; 3000
    carries one read and the two edges carry three each. Zero-cost positions exist
    OUTSIDE the maximum window, so this is the case where the selector has to take
    a severing position and price it.
    """

    fixture = Fixture(tmp_path, length=5000)
    for position in (2900, 3100):
        for i in range(3):
            fixture.spanning_read("crowd_{}_{}".format(position, i), position)
    fixture.spanning_read("single", 3000)
    return fixture.build()


class _RecordingBam:
    """A pysam.AlignmentFile that records the intervals fetched through it.

    Annulus-only fetching is invisible in the ANSWER -- re-scanning from the target
    at every rung chooses the same position -- so the only way to hold the code to
    it is to observe the fetches.
    """

    def __init__(self, path, mode, log):
        self._bam = pysam.AlignmentFile(path, mode)
        self._log = log

    def fetch(self, chrom, start=None, end=None, **kwargs):
        self._log.append((start, end))
        return self._bam.fetch(chrom, start, end, **kwargs)

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self._bam.close()
        return False

    def __getattr__(self, name):
        return getattr(self._bam, name)


def _record_fetches(monkeypatch):
    log = []
    monkeypatch.setattr(
        selector, "_open_bam", lambda path, mode: _RecordingBam(path, mode, log)
    )
    return log


def _search_fetches(log, depth_window=100):
    """The scoring fetches, without the two-base per-cut name lookups."""

    return [(start, end) for start, end in log if end - start > depth_window]


def test_a_severing_position_is_selected_when_nothing_better_exists(tmp_path):
    """The contract the rejected hard rule got wrong: price it, do not refuse it.

    Every position in the window severs something. The cut is placed at the
    cheapest one, the read it severs is dropped, counted and NAMED, and the run
    goes on to produce two chunks. Under the rejected rule this target came back
    declined with one chunk covering the whole contig.
    """

    fixture = _every_position_severed(tmp_path)

    selection = _select(fixture, segment_span=3000, wiggle=200, minimum_span=1000)

    assert [cut.position for cut in selection.cuts] == [3000]
    assert selection.cuts[0].spanning_dropped == 1
    assert selection.total_dropped == 1
    assert set(selection.dropped_read_names) == {"single"}
    assert not selection.unplaced
    # the cut happened, so the contig really is partitioned
    assert [(s.lend, s.rend) for s in selection.segments] == [(1, 3000), (3001, 5000)]


def test_a_window_covered_end_to_end_still_places_its_cut(tmp_path):
    """The regime the hard rule would have broken: no clean position anywhere.

    One alignment covers the entire maximum window, which is what deep data looks
    like everywhere. The selector must still cut, take the collateral, and say what
    it took.
    """

    fixture = Fixture(tmp_path, length=20000)
    fixture.add_read("blanket", "+", [(8500, 11500)])
    fixture.build()

    selection = _select(fixture, segment_span=10000, wiggle=2000, minimum_span=5000)

    assert [cut.position for cut in selection.cuts] == [10000]
    assert selection.cuts[0].spanning_dropped == 1
    assert selection.cuts[0].severed_monoexonic == 1
    assert selection.cuts[0].severed_multiexon == 0
    assert not selection.unplaced
    assert selection.dropped_read_names == {"blanket": [10000]}


# -- the multi-exon weighting ---------------------------------------------------


def _mono_against_multiexon(tmp_path, monoexonic_count):
    """Target 3000 severing one SPLICED read; 2900 severing N monoexonic reads.

    The other three grid positions in the window carry 20 monoexonic reads each,
    so neither of the two positions under test can win by default. Nothing here is
    clean, which is the point: the choice between dirty positions is the decision
    the weighting exists to make.
    """

    fixture = Fixture(tmp_path, length=5000)
    # the spliced read's INTRON crosses 3000; its blocks do not reach it
    fixture.add_read("spliced_at_3000", "+", [(2950, 2980), (3020, 3050)])
    for i in range(monoexonic_count):
        fixture.spanning_read("mono_2900_{}".format(i), 2900)
    for position in (2800, 3100, 3200):
        for i in range(20):
            fixture.spanning_read("filler_{}_{}".format(position, i), position)
    return fixture.build()


def test_one_severed_multiexon_read_loses_to_three_monoexonic_ones(tmp_path):
    """K=10: a junction-bearing read is worth more than three reads of depth."""

    fixture = _mono_against_multiexon(tmp_path, 3)

    selection = _select(
        fixture,
        segment_span=3000,
        wiggle=400,
        minimum_span=1000,
        multiexon_weight=10,
    )

    cut = selection.cuts[0]
    assert cut.position == 2900, "3000 severs one spliced read, which costs 10"
    assert (cut.severed_monoexonic, cut.severed_multiexon) == (3, 0)
    assert cut.spanning_dropped == 3
    assert cut.severed_weight == 3


@pytest.mark.parametrize(
    "weight,monoexonic_count,expected",
    (
        # unweighted: one alignment beats three, so the spliced read is chosen
        (1, 3, 3000),
        # K=10 moves the crossover to "more than ten monoexonic reads"
        (10, 3, 2900),
        (10, 9, 2900),
        (10, 11, 3000),
        # K=4 moves it again, and in the direction K predicts
        (4, 3, 2900),
        (4, 5, 3000),
    ),
)
def test_the_weighting_crossover_moves_with_K(
    tmp_path, weight, monoexonic_count, expected
):
    """The crossover sits at K monoexonic reads, and K is what moves it.

    Asserting the position rather than the cost, because the position is what the
    pipeline consumes. A weight that were merely recorded and not applied would
    pass a cost assertion and fail this one.
    """

    fixture = _mono_against_multiexon(tmp_path, monoexonic_count)

    selection = _select(
        fixture,
        segment_span=3000,
        wiggle=400,
        minimum_span=1000,
        multiexon_weight=weight,
    )

    assert [cut.position for cut in selection.cuts] == [expected]


def test_a_severed_spliced_read_is_counted_by_its_reference_span(tmp_path):
    """The invariant the whole cost rests on: the SPAN severs, not the blocks.

    An alignment's reference span contains its introns, and the splice graph's
    edges come from exactly those N ops. A cost that counted aligned BLOCKS would
    let an intron cross a cut for free while the alignment's two halves stayed
    connected through the edge -- so the cut would sever a locus and be scored as
    clean. This read has no aligned base anywhere in the window; only its intron is
    there. It must still be severed, and it must be classified MULTI-EXON.
    """

    fixture = Fixture(tmp_path, length=6000)
    # blocks 2000-2100 and 3900-4000, intron 2101-3899 across the whole window
    fixture.add_read("spliced", "+", [(2000, 2100), (3900, 4000)])
    fixture.build()

    window = (2900, 3100)
    for block in ((2000, 2100), (3900, 4000)):
        assert block[1] < window[0] or block[0] > window[1], "no aligned base here"

    selection = _select(
        fixture,
        segment_span=3000,
        wiggle=200,
        minimum_span=1000,
        multiexon_weight=10,
    )

    cut = selection.cuts[0]
    assert cut.position == 3000
    assert cut.spanning_dropped == 1
    assert (cut.severed_monoexonic, cut.severed_multiexon) == (0, 1)
    assert cut.severed_weight == 10
    assert selection.dropped_read_names == {"spliced": [3000]}


def test_the_weight_must_not_make_a_spliced_read_cheaper(tmp_path):
    """0 would make severing junctions free, which is the one thing K must not do."""

    fixture = Fixture(tmp_path, length=5000).build()
    with pytest.raises(selector.SelectionError):
        _select(fixture, segment_span=3000, minimum_span=1000, multiexon_weight=0)


# -- progressive expansion ------------------------------------------------------


def _zero_only_beyond_2kb(tmp_path):
    """One read covering target 20000 +/- 2 kb, in a 4 kb maximum radius.

    Every position within 2 kb of the target severs it, so the search cannot stop
    before the third rung; positions beyond 2 kb are clean. The nearest clean
    positions are 17900 and 22100, tying on distance, and the lower coordinate wins.
    """

    fixture = Fixture(tmp_path, length=40000)
    fixture.add_read("blanket", "+", [(17950, 22050)])
    return fixture.build()


def test_only_the_new_annulus_is_fetched_as_the_search_widens(
    tmp_path, monkeypatch
):
    """The measured reason expansion is safe at depth, asserted on the FETCHES.

    Re-scanning from the target at every rung returns the same answer, so a test on
    the chosen position cannot tell the two apart. It is the I/O that differs, and
    it differs in the wrong direction: the ladder's radii sum to more than the
    window, so naive re-scanning costs MORE than never expanding at all. Here the
    three rungs would cost 2002 + 4002 + 8002 = 14006 bases re-scanned against 8002
    for one flat scan of the whole window; the annuli cost 7610.
    """

    fixture = _zero_only_beyond_2kb(tmp_path)
    log = _record_fetches(monkeypatch)

    selection = _select(
        fixture,
        segment_span=20000,
        wiggle=8000,
        minimum_span=10000,
        count_denominator=False,
        expansion_radii=(1000, 2000, 4000),
    )

    assert [cut.position for cut in selection.cuts] == [17900]
    assert selection.cuts[0].search_radius == 4000

    # rung 1 is one interval around the target; rungs 2 and 3 are two annuli each,
    # and NONE of them re-reads the middle
    assert _search_fetches(log) == [
        (18999, 21001),  # radius 1000
        (17999, 18901),  # radius 2000, left annulus
        (21099, 22001),  # radius 2000, right annulus
        (15999, 17901),  # radius 4000, left annulus
        (22099, 24001),  # radius 4000, right annulus
    ]
    fetched = sum(end - start for start, end in _search_fetches(log))
    assert fetched == 7610
    # the flat cost of scanning the maximum window once, which annulus fetching
    # must not exceed and naive re-scanning does exceed by 1.75x
    assert fetched < 8002


def test_the_search_stops_at_the_first_rung_that_severs_nothing(
    tmp_path, monkeypatch
):
    """Early termination: one fetch, not five, when the target itself is clean."""

    fixture = Fixture(tmp_path, length=40000)
    # far from the window, so nothing in it is severed
    fixture.spanning_read("elsewhere", 5000)
    fixture.build()
    log = _record_fetches(monkeypatch)

    selection = _select(
        fixture,
        segment_span=20000,
        wiggle=8000,
        minimum_span=10000,
        count_denominator=False,
        expansion_radii=(1000, 2000, 4000),
    )

    assert [cut.position for cut in selection.cuts] == [20000]
    assert selection.cuts[0].spanning_dropped == 0
    assert selection.cuts[0].search_radius == 1000
    assert _search_fetches(log) == [(18999, 21001)]


def test_reaching_the_maximum_window_takes_the_best_available(tmp_path):
    """No clean position at any radius: take the cheapest and own the collateral.

    Distinct costs, so "best" means something: every position severs the blanket
    read, and the target additionally severs four more. The answer is the cheapest
    position nearest the target, at the maximum radius, with the severed reads
    counted and named rather than the target declined.
    """

    fixture = Fixture(tmp_path, length=40000)
    fixture.add_read("blanket", "+", [(15900, 24100)])
    for i in range(4):
        fixture.spanning_read("extra_{}".format(i), 20000)
    fixture.build()

    selection = _select(
        fixture,
        segment_span=20000,
        wiggle=8000,
        minimum_span=10000,
        expansion_radii=(1000, 2000, 4000),
    )

    cut = selection.cuts[0]
    assert cut.position == 19900, "20000 severs five, its neighbours one"
    assert cut.spanning_dropped == 1
    assert cut.search_radius == 4000, "the maximum was reached, not exceeded"
    assert not selection.unplaced
    assert set(selection.dropped_read_names) == {"blanket"}


def test_expansion_reaches_the_position_a_flat_scan_would(tmp_path):
    """Progressive search is an I/O strategy, not a different objective.

    A single rung at the maximum radius IS the flat scan the module used to do.
    Whatever the ladder, the position must match it.
    """

    fixture = _zero_only_beyond_2kb(tmp_path)

    common = dict(segment_span=20000, wiggle=8000, minimum_span=10000)
    laddered = _select(fixture, expansion_radii=(1000, 2000, 4000), **common)
    flat = _select(fixture, expansion_radii=(4000,), **common)

    assert [c.position for c in laddered.cuts] == [c.position for c in flat.cuts]
    assert laddered.total_dropped == flat.total_dropped


def test_the_rung_ladder_never_exceeds_the_maximum_and_always_reaches_it():
    """The window width is a promise: the last rung is the maximum, exactly."""

    assert selector.expansion_rungs(500000) == [5000, 25000, 100000, 250000, 500000]
    # clipped, and the ladder stops as soon as it reaches the maximum
    assert selector.expansion_rungs(30000) == [5000, 25000, 30000]
    assert selector.expansion_rungs(3000) == [3000]
    assert selector.expansion_rungs(0) == [0]
    # beyond the last default radius the maximum is still reached in one more step
    assert selector.expansion_rungs(700000)[-1] == 700000


def test_an_annulus_never_re_covers_what_the_inner_radius_already_did():
    """The arithmetic annulus fetching rests on, checked without a bam."""

    assert selector.annulus_intervals(1000, None, 100, 1, 5000) == [(900, 1100)]
    assert selector.annulus_intervals(1000, 100, 300, 1, 5000) == [
        (700, 899),
        (1101, 1300),
    ]
    # clipped to the window, and an interval the window already excludes is dropped
    assert selector.annulus_intervals(1000, 100, 300, 950, 1150) == [(1101, 1150)]
    assert selector.annulus_intervals(1000, 300, 300, 1, 5000) == []


# -- the maximum window is ABSOLUTE --------------------------------------------


@pytest.mark.parametrize("spacing_MB", (0.2, 2, 10, 20, 100))
def test_the_maximum_window_does_not_track_the_spacing(tmp_path, spacing_MB):
    """The default window is 1 Mb at every spacing, and deliberately so.

    A proportional rule looks tidy because 10% of the shipped 10 Mb spacing is
    exactly the 1 Mb default, and it is wrong: at 2 Mb spacing it would give
    200 kb, which still severs 743 alignments on chr21 where 1 Mb severs none,
    and chr1 and chr21 disagree 8.5-fold at identical parameters (87 against 743).
    The distance to a read-free position is a property of the genome and the
    library in BASES. This test exists to fail if anyone couples the two.
    """

    fixture = Fixture(tmp_path, length=1000).build()

    selection = _select(
        fixture, segment_span=int(spacing_MB * selector.MB), minimum_span=1000
    )

    assert selection.wiggle == 1 * selector.MB


def test_an_explicit_maximum_window_is_honoured_exactly(tmp_path):
    """Testing a finer spacing means passing a smaller window, and getting it."""

    fixture = Fixture(tmp_path, length=1000).build()

    selection = _select(
        fixture, segment_span=2 * selector.MB, wiggle=20000, minimum_span=1000
    )

    assert selection.wiggle == 20000


# -- the annotation is the only remaining decline reason ------------------------


def test_a_fully_blocked_window_is_declined_for_the_annotation(tmp_path):
    """A read can be dropped and counted. A locus cannot, so this stays hard.

    ``genes_contained`` emits a gene whole or not at all, so a locus straddling a
    boundary is emitted by NEITHER neighbouring chunk. A target whose whole window
    is inside an annotated locus is therefore declined, its chunks stay joined, and
    that is accepted behaviour rather than a defect.
    """

    fixture = Fixture(tmp_path, length=5000)
    fixture.add_transcript("g1", "t1", "+", [(2500, 3500)])
    fixture.spanning_read("noise", 3000)
    fixture.build()

    selection = _select(fixture, segment_span=3000, wiggle=200, minimum_span=1000)

    assert selection.cuts == []
    assert len(selection.unplaced) == 1
    declined = selection.unplaced[0]
    assert declined.target == 3000
    assert declined.declined_annotation is True
    assert declined.best_spanning is None, "no admissible position was costed"
    assert selection.total_dropped == 0
    # the two chunks it would have separated stay joined
    assert len(selection.segments) == 1
    assert selection.segments[0].span == fixture.length

    payload = selector.selection_to_dict(selection)
    assert payload["counts"]["targets_declined_annotation"] == 1
    assert payload["unplaced_targets"][0]["declined_annotation"] is True
    report = selector.format_report(selection)
    assert "1 declined for annotation" in report
    assert "DECLINED for the annotation" in report


def test_a_severing_cut_is_never_reported_as_declined(tmp_path):
    """The rejected contract's fingerprint, asserted absent."""

    fixture = _every_position_severed(tmp_path)
    selection = _select(fixture, segment_span=3000, wiggle=200, minimum_span=1000)

    payload = selector.selection_to_dict(selection)
    assert payload["counts"]["cuts_placed"] == 1
    assert payload["counts"]["targets_unplaced"] == 0
    assert payload["counts"]["targets_declined_annotation"] == 0
    assert payload["counts"]["alignments_dropped_at_cuts"] == 1


# -- the accounting a reader sees ----------------------------------------------


def test_the_manifest_and_report_split_severing_by_structure(tmp_path):
    """Severing is expected now, so the report is the only place its cost shows.

    A bare total cannot distinguish a cut that dropped three reads of depth from
    one that dropped three junctions, and those are not the same event.
    """

    fixture = Fixture(tmp_path, length=6000)
    fixture.add_read("spliced_a", "+", [(2900, 2960), (3040, 3100)])
    fixture.spanning_read("mono_a", 3000)
    fixture.spanning_read("mono_b", 3000)
    fixture.build()

    selection = _select(
        fixture,
        segment_span=3000,
        wiggle=0,
        minimum_span=1000,
        multiexon_weight=10,
    )

    cut = selection.cuts[0]
    assert cut.position == 3000
    assert (cut.severed_monoexonic, cut.severed_multiexon) == (2, 1)
    assert cut.spanning_dropped == 3
    assert cut.severed_weight == 12

    payload = selector.selection_to_dict(selection)
    assert payload["params"]["severed_multiexon_weight"] == 10
    assert payload["cuts"][0]["severed_monoexonic"] == 2
    assert payload["cuts"][0]["severed_multiexon"] == 1
    assert payload["cuts"][0]["severed_weighted_cost"] == 12
    assert payload["counts"]["alignments_dropped_monoexonic"] == 2
    assert payload["counts"]["alignments_dropped_multiexon"] == 1
    assert payload["counts"]["severed_weighted_cost_at_cuts"] == 12

    report = selector.format_report(selection)
    assert "severed multi-exon weight 10" in report
    assert "severing is a COST, never a veto" in report
    assert "2 monoexonic and 1 multi-exon, weighted cost 12" in report


def test_the_cli_flag_reaches_the_selection(tmp_path):
    """--severed_multiexon_weight is the pipeline's lever on the cost; it must work."""

    fixture = _mono_against_multiexon(tmp_path, 3)

    def run(weight):
        prefix = str(tmp_path / "w{}".format(weight))
        assert (
            selector.main(
                [
                    "--bam",
                    fixture.bam,
                    "--genome_fa",
                    fixture.fasta,
                    "--strand",
                    "+",
                    "--approx_MB_per_cut",
                    "0.003",
                    "--approx_MB_per_cut_wiggle_window",
                    "0.0004",
                    "--depth_window",
                    "100",
                    "--margin",
                    "0",
                    "--minimum_span",
                    "1000",
                    "--severed_multiexon_weight",
                    str(weight),
                    "--output_prefix",
                    prefix,
                ]
            )
            == 0
        )
        with open(prefix + ".cuts.json", "rt") as fh:
            return json.load(fh)[0]

    weighted = run(10)
    assert weighted["params"]["severed_multiexon_weight"] == 10
    assert weighted["cuts"][0]["position"] == 2900
    assert weighted["counts"]["alignments_dropped_monoexonic"] == 3

    unweighted = run(1)
    assert unweighted["cuts"][0]["position"] == 3000
    assert unweighted["counts"]["alignments_dropped_multiexon"] == 1


def test_a_reference_contig_absent_from_the_bam_is_skipped_not_fatal(tmp_path, capsys):
    """A reference carrying sequences the alignment does not is the NORMAL case.

    The contig list comes from the genome fasta when one is given, but every candidate is
    priced by FETCHING it from the bam, and pysam raises ValueError on a name the bam
    header lacks. Decoys, EBV and alt scaffolds routinely sit in a reference a given bam
    never used: GRCh38_no_alt.fa has 195 sequences against a 10x SBX PBMC bam's 194, and
    chrEBV -- zero reads, zero gencode records -- killed a whole-genome chunked run 246 s
    in, because one failed make-chunks unit makes the pipeline refuse to extract or merge.

    Such a contig cannot be partitioned anyway: no alignments means no depth to place a
    cut against. So it is skipped and COUNTED, not refused.
    """

    fixture = Fixture(tmp_path, length=9000)
    for i in range(12):
        fixture.add_read("r{}".format(i), "+", [(400 + 200 * i, 460 + 200 * i)])
    fixture.add_transcript("g1", "t1", "+", [(300, 900)])
    fixture.build()

    # Append a contig to the FASTA that the bam header has never heard of. The bam is
    # left exactly as built, which is the whole point: reference and alignment disagree.
    with open(fixture.fasta, "at") as ofh:
        print(">cAbsentFromBam", file=ofh)
        print("ACGT" * 30, file=ofh)
    os.remove(fixture.fasta + ".fai")
    pysam.faidx(fixture.fasta)

    argv = [
        "--bam", fixture.bam,
        "--genome_fa", fixture.fasta,
        "--gtf", fixture.gtf,
        "--approx_MB_per_cut", "0.002",
        "--strandless",
        "--output_prefix", str(tmp_path / "cuts"),
    ]
    rc = selector.main(argv)

    assert rc == 0, "a reference contig the bam lacks must not fail the selection"
    err = capsys.readouterr().err
    assert "absent from the bam header" in err, "the skip must be reported, not silent"
    assert "cAbsentFromBam" in err, "and must name what was skipped"
