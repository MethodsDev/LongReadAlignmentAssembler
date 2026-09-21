#!/usr/bin/env python3

import os
import subprocess
import sys

import pytest

from Ascii_isoform_illustrator import (
    Track,
    View,
    compare_structures,
    exons_from_simple_path,
    illustrate,
    intron_chain,
    track_from_simple_path,
    transcripts_from_gtf,
)
from GenomeFeature import Exon, Intron, PolyAsite, TSS
from LRAA_Globals import SPACER
from Splice_graph import Splice_graph
from Transcript import Transcript


def _drawing_rows(text):
    """Only the rows that carry a drawing (they end in a track name)."""
    return [line for line in text.split("\n") if "=" in line or "#" in line]


def _body(line, width, label_width=21):
    """Just the drawing columns of a rendered row: no coordinate label, no name."""
    return line[label_width + 2 : label_width + 2 + width]


def test_intron_chain_matches_transcript_convention():
    # the same identity SQANTI_like_annotator keys on; if these ever diverge,
    # a "shared junction" here would not be a shared junction there
    exons = [(100, 200), (300, 400), (500, 600)]
    transcript = Transcript("contig", [list(e) for e in exons], "+")
    assert intron_chain(exons) == transcript.get_introns()


def test_identical_and_subchain_and_novel_relationships():
    ref = [(100, 200), (300, 400), (500, 600), (700, 800)]

    same = compare_structures(list(ref), ref)
    assert same.relationship == "identical_splice_pattern"
    assert same.query_only_introns == []

    # a 5'-truncated model: contiguous run of the reference's introns
    ism = compare_structures([(320, 400), (500, 600), (700, 800)], ref)
    assert ism.relationship == "contained_subchain"
    assert ism.query_only_introns == []
    assert len(ism.ref_only_introns) == 1

    # exon skipping using only reference splice sites, in a new arrangement
    recombined = compare_structures([(100, 200), (300, 400), (700, 800)], ref)
    assert recombined.relationship == "novel_junctions"
    assert recombined.query_only_introns == [(401, 699)]

    # a donor the reference never uses
    novel = compare_structures([(100, 250), (300, 400), (500, 600), (700, 800)], ref)
    assert novel.relationship == "novel_junctions"
    assert novel.query_only_introns == [(251, 299)]
    assert len(novel.shared_introns) == 2


def test_retained_intron_reuses_known_junctions_without_creating_one():
    # intron 2 retained: exons 300-400 and 500-600 fused.  Both surviving introns
    # exist in the reference, but no longer as a contiguous run of its chain, so
    # this is neither novel nor a truncation and must not be reported as either.
    ref = [(100, 200), (300, 400), (500, 600), (700, 800)]
    retained = compare_structures([(100, 200), (300, 600), (700, 800)], ref)
    assert retained.relationship == "known_junctions_recombined"
    assert retained.query_only_introns == []
    assert retained.shared_introns == [(201, 299), (601, 699)]

    # exon skipping, by contrast, welds two exons across a junction the reference
    # does not contain, while the rest of the chain still matches
    skipping = compare_structures([(100, 200), (500, 600), (700, 800)], ref)
    assert skipping.relationship == "novel_junctions"
    assert skipping.query_only_introns == [(201, 499)]
    assert skipping.shared_introns == [(601, 699)]


def test_query_extending_past_the_reference_is_a_superchain():
    short_ref = [(100, 200), (300, 400), (500, 600)]
    longer = compare_structures([(100, 200), (300, 400), (500, 600), (700, 800)], short_ref)
    assert longer.relationship == "contains_reference_subchain"
    assert longer.query_only_introns == [(601, 699)]


def test_disjoint_chains_are_reported_as_sharing_nothing():
    # every disjoint query chain also has query-only introns, so this must be
    # tested before the novel-junction case or the label is never emitted; a
    # locus-level mismatch and a one-junction difference are not the same finding
    result = compare_structures(
        [(5000, 5100), (5300, 5400)], [(100, 200), (300, 400), (500, 600)]
    )
    assert result.relationship == "no_shared_junctions"
    assert result.shared_introns == []
    assert result.query_only_introns == [(5101, 5299)]
    assert result.overlap_bp == 0


def test_monoexonic_reports_overlap_rather_than_junctions():
    result = compare_structures([(150, 250)], [(100, 200), (300, 400)])
    assert result.relationship == "monoexonic_query"
    assert result.overlap_bp == 51
    assert "overlap=51bp" in result.summary()


def test_arrowhead_marks_the_three_prime_end_and_vanishes_when_clipped():
    plus = View("contig", 1, 1000, width=40, mode="proportional")
    plus.add_transcript("fwd", [(100, 200), (400, 500)], "+")
    plus.add_transcript("rev", [(100, 200), (400, 500)], "-")
    rows = _drawing_rows(plus.render())
    fwd, rev = _body(rows[0], 40), _body(rows[1], 40)
    assert fwd.rstrip()[-1] == ">"
    assert rev.strip()[0] == "<"

    # the same '+' model seen through a window that ends before it does: drawing
    # an arrowhead at the edge would assert a 3' end that is not there
    clipped = View("contig", 1, 300, width=40, mode="proportional")
    clipped.add_transcript("fwd", [(100, 200), (400, 500)], "+")
    assert ">" not in _body(_drawing_rows(clipped.render())[0], 40)


def test_compressed_mode_keeps_a_small_exon_visible_in_a_large_locus():
    exons = [(1000, 1050), (99000, 99050)]
    true_scale = View("contig", 1, 100000, width=60, mode="proportional")
    true_scale.add_transcript("t", exons, "+")
    compressed = View("contig", 1, 100000, width=60, mode="compressed")
    compressed.add_transcript("t", exons, "+")

    # at true scale a 51 bp exon in a 100 kb window is a rounding error; the whole
    # point of compressed mode is that it is not
    assert _body(_drawing_rows(true_scale.render())[0], 60).count("=") <= 2
    assert _body(_drawing_rows(compressed.render())[0], 60).count("=") >= 20


def test_highlighted_introns_are_drawn_with_the_highlight_char():
    view = View("contig", 1, 1000, width=60, mode="proportional")
    view.add_transcript(
        "t", [(100, 200), (400, 500), (700, 800)], "+", highlight=[(201, 399)]
    )
    body = _body(_drawing_rows(view.render())[0], 60)
    assert "*" in body
    assert "-" in body  # the non-highlighted intron keeps the default glyph


def test_illustrate_annotates_every_row_against_the_reference():
    text = illustrate(
        [
            ("ref", [(100, 200), (300, 400), (500, 600)], "+"),
            ("same", [(100, 200), (300, 400), (500, 600)], "+"),
            ("novel", [(100, 250), (300, 400), (500, 600)], "+"),
        ],
        reference="ref",
        contig_acc="contig",
        width=60,
    )
    assert "(reference)" in text
    assert "identical_splice_pattern" in text
    assert "novel_junctions" in text
    assert "*" in text  # the novel intron is marked


def _write_gtf(path, rows):
    with open(path, "w") as ofh:
        for contig, feature, lend, rend, strand, attrs in rows:
            ofh.write(
                "\t".join(
                    [contig, "test", feature, str(lend), str(rend), ".", strand, ".", attrs]
                )
                + "\n"
            )
    return path


def test_region_selects_overlapping_transcripts_not_only_contained_ones(tmp_path):
    gtf = _write_gtf(
        os.path.join(tmp_path, "models.gtf"),
        [
            ("chr1", "exon", 100, 200, "+", 'gene_id "g1"; transcript_id "inside";'),
            ("chr1", "exon", 300, 400, "+", 'gene_id "g1"; transcript_id "inside";'),
            ("chr1", "exon", 100, 200, "+", 'gene_id "g1"; transcript_id "spanning";'),
            ("chr1", "exon", 9000, 9100, "+", 'gene_id "g1"; transcript_id "spanning";'),
            ("chr1", "exon", 20000, 20100, "+", 'gene_id "g2"; transcript_id "elsewhere";'),
        ],
    )
    names = {
        t.get_transcript_id()
        for t in transcripts_from_gtf(gtf, contig_acc="chr1", lend=50, rend=500)
    }
    # a model running off the right edge is exactly what you are looking at when
    # you ask for a window; the parser's containment rule would hide it
    assert names == {"inside", "spanning"}


def test_gene_and_transcript_filters(tmp_path):
    gtf = _write_gtf(
        os.path.join(tmp_path, "models.gtf"),
        [
            ("chr1", "exon", 100, 200, "+", 'gene_id "g1"; gene_name "ALPHA"; transcript_id "t1";'),
            ("chr1", "exon", 100, 200, "+", 'gene_id "g2"; gene_name "BETA"; transcript_id "t2";'),
        ],
    )
    assert [t.get_transcript_id() for t in transcripts_from_gtf(gtf, gene="ALPHA")] == ["t1"]
    assert [t.get_transcript_id() for t in transcripts_from_gtf(gtf, transcript_ids=["t2"])] == [
        "t2"
    ]


def test_unstranded_records_need_an_explicit_strand(tmp_path):
    # LRAA's multipath debug dumps write '?' because the dump has no strand; the
    # parser must keep refusing that by default and accept it only when asked
    gtf = _write_gtf(
        os.path.join(tmp_path, "mpgns.gtf"),
        [
            ("contig", "exon", 100, 200, "?", 'gene_id "g"; transcript_id "mp1";'),
            ("contig", "exon", 300, 400, "?", 'gene_id "g"; transcript_id "mp1";'),
        ],
    )
    with pytest.raises(AssertionError):
        transcripts_from_gtf(gtf)

    transcripts = transcripts_from_gtf(gtf, unstranded_as="+")
    assert [t.get_strand() for t in transcripts] == ["+"]


def _simple_path_splice_graph():
    sg = Splice_graph()
    Exon.reset_counter()
    Intron.intron_id_counter = 0
    TSS.TSS_id_counter = 0
    PolyAsite.polyA_counter = 0

    def register(node):
        sg._node_id_to_node[node.get_id()] = node
        return node.get_id()

    ids = {
        "tss": register(TSS("contig", 100, 100, "+", 5)),
        # two abutting exon nodes: one exon split at an internal splice site
        "e1a": register(Exon("contig", 100, 150, "+", 10)),
        "e1b": register(Exon("contig", 151, 200, "+", 10)),
        "i1": register(Intron("contig", 201, 299, "+", 8)),
        "e2": register(Exon("contig", 300, 400, "+", 10)),
        "e3": register(Exon("contig", 700, 800, "+", 10)),
        "polya": register(PolyAsite("contig", 800, 800, "+", 4)),
    }
    return sg, ids


def test_simple_path_fuses_adjacent_exon_nodes_and_keeps_introns():
    sg, ids = _simple_path_splice_graph()
    path = [ids["tss"], ids["e1a"], ids["e1b"], ids["i1"], ids["e2"], ids["polya"]]
    exons, spacers, site_marks = exons_from_simple_path(sg, path)
    # E:1 and E:2 abut, so they are one drawn exon, not two with a phantom intron
    assert exons == [(100, 200), (300, 400)]
    assert spacers == []
    assert site_marks == [(100, "T"), (800, "A")]


def test_spacer_in_a_simple_path_is_drawn_as_unresolved():
    sg, ids = _simple_path_splice_graph()
    path = [ids["e1a"], ids["e1b"], SPACER, ids["e3"]]
    exons, spacers, _ = exons_from_simple_path(sg, path)
    assert exons == [(100, 200), (700, 800)]
    assert spacers == [(201, 699)]

    track = track_from_simple_path(sg, path, "read_path")
    assert track.junc_chars == {(201, 699): "?"}

    view = View("contig", 1, 1000, width=60, mode="proportional")
    view.add(track)
    body = view.render().split("\n")[-1]
    assert "?" in body  # the gap is drawn as unknown, not as a confident intron


def test_track_rejects_a_multi_character_glyph():
    with pytest.raises(ValueError):
        Track("t", [(1, 10)], "+", glyph="==")


def test_mark_legend_prints_the_coordinate_once():
    view = View("contig", 1, 1000, width=40, mode="proportional")
    view.add_transcript("t", [(100, 200), (400, 500)], "+")
    view.mark(450)
    legend = view.render().split("\n")[-1]
    assert legend.split() == ["450", "^"]

    view_with_label = View("contig", 1, 1000, width=40, mode="proportional")
    view_with_label.add_transcript("t", [(100, 200), (400, 500)], "+")
    view_with_label.mark(450, "^", "novel acceptor")
    assert view_with_label.render().split("\n")[-1].split() == [
        "450",
        "^",
        "novel",
        "acceptor",
    ]


def test_illustrate_refuses_to_overlay_two_contigs():
    chr1_model = Transcript("chr1", [[100, 200], [300, 400]], "+")
    chr2_model = Transcript("chr2", [[100, 200], [300, 400]], "+")
    with pytest.raises(ValueError, match="more than one contig"):
        illustrate([chr1_model, chr2_model])


def test_cli_refuses_a_selection_spanning_contigs(tmp_path):
    gtf = _write_gtf(
        os.path.join(tmp_path, "two_contigs.gtf"),
        [
            ("chr1", "exon", 100, 200, "+", 'gene_id "g"; transcript_id "t1";'),
            ("chr1", "exon", 300, 400, "+", 'gene_id "g"; transcript_id "t1";'),
            ("chr2", "exon", 100, 200, "+", 'gene_id "g"; transcript_id "t2";'),
            ("chr2", "exon", 300, 400, "+", 'gene_id "g"; transcript_id "t2";'),
        ],
    )
    cli = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "..", "util", "ascii_isoform_view.py"
    )
    spanning = subprocess.run(
        [sys.executable, cli, "-g", gtf], capture_output=True, text=True
    )
    assert spanning.returncode != 0
    assert "spans 2 contigs" in spanning.stderr

    restricted = subprocess.run(
        [sys.executable, cli, "-g", gtf, "-r", "chr2:1-1000"],
        capture_output=True,
        text=True,
    )
    assert restricted.returncode == 0, restricted.stderr
    assert "t2" in restricted.stdout and "t1" not in restricted.stdout
