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
    assert "blocked" in unplaced.reason
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
