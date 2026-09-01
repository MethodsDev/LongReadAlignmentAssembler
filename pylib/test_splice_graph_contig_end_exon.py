#!/usr/bin/env python3

"""An exon reaching the contig's last base must not be reported one base past it.

``_contig_base_cov`` is sized ``contig_len + 1`` so that position p lives at index p
(Splice_graph.py:493). The scan that turns coverage into exon segments visits
i = 1..contig_len and every close inside the loop uses i or i - 1, correctly. Only
the tail case -- a segment still open when the scan runs out of array -- used the
array's LENGTH as a coordinate, producing an exon node at contig_len + 1 and a GTF
record whose end lies outside the contig.

Latent for as long as it existed. With
``fracture_splice_graph_at_input_transcript_bounds`` ON, which is the default, an
input transcript's rend is a registered boundary and the segment is closed there
before the scan finishes. ``merge_LRAA_GTFs.py`` turns that off for every merge
(:218), so a model reaching the contig's last base fell through to the tail case:
MEASURED, an input model spanning 1-5000 on a 5,000 bp contig merged to 1-5001,
while an interior model was unaffected under either setting. That is an invalid GTF
coordinate, and tools that validate or extract sequence are entitled to reject it.
"""

import LRAA_Globals
import pytest
from Splice_graph import Splice_graph
from Transcript import Transcript

CONTIG = "chrX"
CONTIG_LEN = 5000


def _seq(n):
    return "A" * n


def _transcript(tid, lend, rend):
    t = Transcript(CONTIG, [[lend, rend]], "+")
    t.set_transcript_id(tid)
    t.set_gene_id("g_" + tid)
    return t


def _nodes(spans, fracture):
    """Exon node coordinates for models at `spans`, with the fracture knob set."""

    previous = LRAA_Globals.config["fracture_splice_graph_at_input_transcript_bounds"]
    LRAA_Globals.config["fracture_splice_graph_at_input_transcript_bounds"] = fracture
    try:
        sg = Splice_graph()
        sg.build_splice_graph_for_contig(
            CONTIG,
            "+",
            _seq(CONTIG_LEN),
            alignments_bam_file=None,
            region_lend=None,
            region_rend=None,
            input_transcripts=[
                _transcript("t{}".format(i), lend, rend)
                for i, (lend, rend) in enumerate(spans)
            ],
            quant_mode=False,
        )
        return sorted(n.get_coords() for n in sg._splice_graph.nodes())
    finally:
        LRAA_Globals.config[
            "fracture_splice_graph_at_input_transcript_bounds"
        ] = previous


@pytest.mark.parametrize("fracture", [True, False])
def test_an_exon_reaching_the_contig_end_stops_at_the_contig_end(fracture):
    """Whichever way the fracture knob is set, the node ends at the last base.

    Parametrized on the knob because the bug was reachable through exactly one of
    its settings, and the assertion is the same either way: a coordinate outside the
    contig is wrong regardless of how the graph was built.
    """

    assert _nodes([(1, CONTIG_LEN), (1, CONTIG_LEN)], fracture) == [(1, CONTIG_LEN)]


@pytest.mark.parametrize("fracture", [True, False])
def test_an_exon_ending_short_of_the_contig_end_is_unchanged(fracture):
    """The control. The tail case is specific to running out of array, so an
    interior segment must be untouched by the fix -- if this moved, the change was
    to the scan's arithmetic rather than to its boundary condition."""

    assert _nodes([(1000, 2000), (1000, 2000)], fracture) == [(1000, 2000)]


def test_no_exon_node_can_exceed_the_contig_length():
    """Stated as the invariant rather than as a coordinate.

    A model abutting the last base, one ending a base short, and one well inside,
    all in the same graph: nothing may be reported beyond CONTIG_LEN.
    """

    spans = [
        (1, CONTIG_LEN),
        (100, CONTIG_LEN - 1),
        (2000, 3000),
    ]
    for fracture in (True, False):
        for lend, rend in _nodes(spans, fracture):
            assert rend <= CONTIG_LEN, (fracture, lend, rend)
            assert lend >= 1, (fracture, lend, rend)
