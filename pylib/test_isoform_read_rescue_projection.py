#!/usr/bin/env python3

import importlib.machinery
import importlib.util
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import networkx as nx
import pysam

import LRAA_Globals
from GenomeFeature import Exon, Intron, PolyAsite, TSS
from IsoformReadRescue import (
    _build_transcript_models,
    _explained_read_bases,
    _collect_read_sequences,
    _merge_contiguous_genomic_segments,
    _normalize_read_identifier,
    _parse_rescue_alignments,
)
from LRAA import LRAA
from Splice_graph import Splice_graph
from Transcript import Transcript


def test_merge_contiguous_genomic_segments_collapses_split_exons():
    segments = [
        (100, 110),
        (111, 125),
        (200, 210),
        (211, 211),
        (300, 350),
    ]

    assert _merge_contiguous_genomic_segments(segments) == [
        (100, 125),
        (200, 211),
        (300, 350),
    ]


def test_merge_contiguous_genomic_segments_preserves_true_introns():
    segments = [
        (100, 150),
        (250, 300),
    ]

    assert _merge_contiguous_genomic_segments(segments) == segments


def _build_minus_strand_boundary_fixture():
    Exon.reset_counter()
    splice_graph = Splice_graph()
    splice_graph._contig_acc = "chr1"
    splice_graph._contig_strand = "-"

    polya = PolyAsite("chr1", 100, 100, "-", 10)
    exon_low = Exon("chr1", 100, 149, "-", 10)
    intron = Intron("chr1", 150, 199, "-", 10)
    exon_high = Exon("chr1", 200, 249, "-", 10)
    tss = TSS("chr1", 249, 249, "-", 10)

    splice_graph._splice_graph = nx.DiGraph(
        [
            (polya, exon_low),
            (exon_low, intron),
            (intron, exon_high),
            (exon_high, tss),
        ]
    )
    splice_graph._intron_objs = {"150:199": intron}
    splice_graph._finalize_splice_graph()

    transcript = Transcript("chr1", [[100, 149], [200, 249]], "-")
    transcript.set_gene_id("g1")
    transcript.set_transcript_id("t1")
    transcript.set_simple_path(
        [
            polya.get_id(),
            exon_low.get_id(),
            intron.get_id(),
            exon_high.get_id(),
            tss.get_id(),
        ]
    )
    transcript_models = _build_transcript_models(splice_graph, [transcript], "A" * 300)

    return splice_graph, transcript_models, polya, exon_low, intron, exon_high, tss


def test_minus_strand_rescue_soft_clips_follow_transcript_ends(tmp_path, monkeypatch):
    monkeypatch.setitem(LRAA_Globals.config, "max_soft_clip_at_TSS", 0)
    monkeypatch.setitem(LRAA_Globals.config, "max_soft_clip_at_PolyA", 0)
    (
        splice_graph,
        transcript_models,
        polya,
        exon_low,
        intron,
        exon_high,
        tss,
    ) = _build_minus_strand_boundary_fixture()

    rescue_sam = tmp_path / "minus_rescue.sam"
    sequence = "C" * 100
    rescue_sam.write_text(
        "@HD\tVN:1.6\tSO:unknown\n"
        "@SQ\tSN:g1^t1\tLN:100\n"
        f"clip_tx_3prime\t0\tg1^t1\t1\t60\t95M5S\t*\t0\t0\t{sequence}\t*\tAS:i:95\tNM:i:0\n"
        f"clip_tx_5prime\t0\tg1^t1\t6\t60\t5S95M\t*\t0\t0\t{sequence}\t*\tAS:i:95\tNM:i:0\n"
    )

    _, read_details = _parse_rescue_alignments(
        str(rescue_sam),
        splice_graph,
        transcript_models,
        read_path_mapper=LRAA(splice_graph)._map_read_to_graph,
    )
    paths_by_read = {
        read_name: read_details["read_name_to_multipaths"][
            _normalize_read_identifier(read_name)
        ][0].get_simple_path()
        for read_name in ("clip_tx_3prime", "clip_tx_5prime")
    }

    assert paths_by_read["clip_tx_3prime"] == [
        exon_low.get_id(),
        intron.get_id(),
        exon_high.get_id(),
        tss.get_id(),
    ]
    assert paths_by_read["clip_tx_5prime"] == [
        polya.get_id(),
        exon_low.get_id(),
        intron.get_id(),
        exon_high.get_id(),
    ]


def _supplementary_then_primary_bam(path):
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}

    def _alignment(flag, start, cigar, seq):
        aln = pysam.AlignedSegment()
        aln.query_name = "r1"
        aln.flag = flag
        aln.reference_id = 0
        aln.reference_start = start
        aln.mapping_quality = 60
        aln.cigar = cigar
        aln.query_sequence = seq
        aln.set_tag("NM", 0)
        return aln

    with pysam.AlignmentFile(str(path), "wb", header=header) as bam_writer:
        bam_writer.write(_alignment(2048, 100, [(5, 10), (0, 40)], "G" * 40))
        bam_writer.write(_alignment(0, 300, [(0, 50)], "T" * 50))
    pysam.index(str(path))


def test_rescue_sequences_come_from_the_primary_record(tmp_path):
    bam_path = tmp_path / "supplementary.bam"
    _supplementary_then_primary_bam(bam_path)

    (
        read_name_to_seq,
        read_name_to_genome_explained,
        read_name_to_genome_gap_id,
        read_name_to_allowed_target_ids,
    ) = _collect_read_sequences(str(bam_path), "chr1", None, None, {"r1"}, "+")

    assert read_name_to_seq == {"r1": "T" * 50}
    # the primary is a clean 50-base match, so it explains the whole read and agrees
    # with the genome perfectly over the span it covers
    assert read_name_to_genome_explained == {"r1": 50}
    assert read_name_to_genome_gap_id == {"r1": 1.0}
    # No exon overlap index passed, so no locality restriction is computed and
    # _parse_rescue_alignments imposes none. Only rescue supplies one.
    assert read_name_to_allowed_target_ids is None


def test_rescue_declined_when_it_explains_less_of_the_read_than_the_genome():
    # A target that disagrees with part of a read can soft-clip the disagreement away
    # and still score 100% identity over what remains, which is how a read carrying a
    # novel terminal exon gets reshaped onto a reference model. Rescue is additive, so
    # such an alignment cannot correct anything and only adds a competing path. It is
    # declined by comparing explained read bases, not identity over the aligned part.
    import pysam

    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "g1^t1", "LN": 100}]}

    def _aln(cigar, nm):
        segment = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(header))
        segment.query_name = "r1"
        segment.reference_id = 0
        segment.reference_start = 0
        segment.cigarstring = cigar
        segment.query_sequence = "A" * 100
        segment.set_tag("NM", nm)
        return segment

    # whole read matched: explains every base
    assert _explained_read_bases(_aln("100M", 0)) == 100
    # thirty bases clipped away: those bases are not explained, though identity over
    # the aligned portion is still a perfect 100%
    assert _explained_read_bases(_aln("30S70M", 0)) == 70
    # substitutions reduce the count, insertions and deletions do not double-count
    assert _explained_read_bases(_aln("100M", 4)) == 96
    assert _explained_read_bases(_aln("96M4I", 4)) == 96


def _gap_baseline_rescue_sam(tmp_path):
    """One rescue alignment that ties the explained-bases baseline and loses on gaps.

    40M20I40M over a 100-base target: 80 read bases are matched, so explained bases are
    80 either way, which is exactly the tie the explained-bases test cannot break. The 20
    inserted bases enter the gap-aware span, so the alignment agrees with its target over
    100 positions rather than 80 -- 0.80 where a clean genome alignment scores 1.00.
    """
    sam = tmp_path / "gap_baseline.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:unknown\n"
        "@SQ\tSN:g1^t1\tLN:100\n"
        "r1\t0\tg1^t1\t1\t60\t40M20I40M\t*\t0\t0\t{}\t*\tAS:i:80\tNM:i:20\n".format(
            "A" * 100
        )
    )
    return sam


@pytest.mark.parametrize(
    "genome_gap_id, expect_declined",
    [
        # The genome agrees with the genome perfectly, the rescue agrees with its target
        # over a longer span at 0.80: the rescue is the worse account and is declined.
        (1.0, True),
        # Non-vacuity control for the case above. Same alignment, same explained-bases
        # tie, same locality -- only the baseline moves, and now the rescue is accepted.
        # Without this, a fixture that never reached the rule would look like a pass.
        (0.5, False),
    ],
)
def test_rescue_declined_when_it_agrees_worse_with_its_target_than_the_genome(
    tmp_path, monkeypatch, genome_gap_id, expect_declined
):
    """The gap-aware baseline is enforced, not merely recorded, and is unconditional.

    Both genome baselines used to be emptied for callers that passed genome target
    gating, which was how the whole-genome-versus-whole-transcriptome modes ran; those
    modes and that flag are gone, so there is no longer any way to reach rescue with the
    baselines switched off. Measured on chr21 HiFi, these baselines decline 668 candidate
    rescues against 16 for locality, so they are the larger of the two rules and the one
    whose silent loss would matter most.

    The two content thresholds are pinned rather than inherited: at their defaults the
    20-base insertion is refused as a long indel and 0.80 identity sits exactly on the
    LowFi floor, and either would decline the alignment before the rule under test ran.
    """
    monkeypatch.setitem(LRAA_Globals.config, "rescue_unassigned_max_indel_length", 0)
    monkeypatch.setitem(LRAA_Globals.config, "rescue_unassigned_min_per_id", 50)
    splice_graph, transcript_models = _build_minus_strand_boundary_fixture()[:2]

    rescued_mps, details = _parse_rescue_alignments(
        str(_gap_baseline_rescue_sam(tmp_path)),
        splice_graph,
        transcript_models,
        read_path_mapper=LRAA(splice_graph)._map_read_to_graph,
        read_name_to_allowed_target_ids={"r1": {"g1^t1"}},
        # Explained bases tie at 80, so this rule must pass the alignment through.
        read_name_to_genome_explained={"r1": 80},
        read_name_to_genome_gap_id={"r1": genome_gap_id},
    )

    rejections = details["alignment_rejections"]
    assert (
        "explains_less_than_genome" not in rejections
    ), "the explained-bases rule must tie here, or this fixture tests the wrong rule"
    if expect_declined:
        assert rejections == {"agrees_worse_than_genome": 1}, rejections
        assert rescued_mps == []
    else:
        assert "agrees_worse_than_genome" not in rejections, rejections
        assert len(rescued_mps) == 1, rejections


def test_rescue_raises_rather_than_skipping_when_minimap2_is_absent(monkeypatch):
    """Returning empty here would silently drop every read rescue was to recover."""

    import shutil as rescue_shutil

    import IsoformReadRescue

    monkeypatch.setattr(
        IsoformReadRescue.shutil, "which", lambda name: None, raising=True
    )
    assert rescue_shutil is IsoformReadRescue.shutil

    with pytest.raises(RuntimeError, match="minimap2"):
        IsoformReadRescue.rescue_unassigned_reads_to_transcriptome(
            None, [], "", "reads.bam", "chr1", 1, 1000, set()
        )


def test_lraa_exits_before_doing_work_when_rescue_has_no_minimap2(monkeypatch):
    """The CLI guard must fire when rescue is enabled, and only then."""

    lraa_cli = _load_lraa_cli()
    monkeypatch.setattr(lraa_cli.shutil, "which", lambda name: None, raising=True)

    with pytest.raises(SystemExit) as exit_info:
        lraa_cli._require_minimap2_for_transcriptome_alignment(True)
    assert "minimap2" in str(exit_info.value)
    assert "transcriptome read rescue" in str(exit_info.value)
    assert "--no_rescue_unassigned_reads_via_transcriptome_alignment" in str(
        exit_info.value
    ), "and must name the way to run without it"

    # rescue turned off never aligns to transcripts, so it must still run
    lraa_cli._require_minimap2_for_transcriptome_alignment(False)


def _load_lraa_cli():
    """Import the top-level LRAA script, which is not importable by name."""

    loader = importlib.machinery.SourceFileLoader("lraa_cli", str(REPO_ROOT / "LRAA"))
    spec = importlib.util.spec_from_loader("lraa_cli", loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module
