#!/usr/bin/env python3

"""Guide TSS/PolyA boundary seeding is gated to the modes that own the input GTF (v0.43.0).

_integrate_input_transcript_structures turns an input transcript's annotated TSS/PolyA
boundaries into graph vertices (path #2, distinct from the read-derived veto reprieve).
As of v0.43.0 that seeding runs ONLY where the input GTF is authoritative:

  * quant-only (quant_mode=True): the GTF is the fixed model set reads are assigned to.
  * merge (LRAA_Globals.LRAA_MODE == "MERGE"): the merge reconciles constituent TSS/PolyA.

In DISCOVERY (quant_mode=False and not MERGE) the input GTF is a STRUCTURAL GUIDE only --
its termini are not cleavage-validated, and seeding them re-introduced A-rich internal-
priming ends as PolyA vertices (the scg internal-priming excess). This exercises the gate
end to end by building a graph from a PolyA+TSS-annotated multi-exon guide with NO bam
(the merge's own call shape) and asserting which modes seed a boundary vertex.

Fracturing/introns/coverage are applied regardless of the gate; only boundary vertex
seeding is gated, so discovery keeps the structural guidance either way.
"""

import pytest

import LRAA_Globals
import Splice_graph
from Transcript import GTF_contig_to_transcripts


# multi-exon guide (one intron), matching the residual multi-exonic scg case; the 3' end
# (PolyA) is 1000 and the 5' end (TSS) is 500 on '+'. Clean genomic context, so nothing
# here depends on the veto -- seeding is veto-exempt (from_reads=False) by construction.
_GUIDE_GTF = (
    'chr1\tLRAA\ttranscript\t500\t1000\t.\t+\t.\t'
    'gene_id "g1"; transcript_id "t1"; TSS "True"; PolyA "True";\n'
    'chr1\tLRAA\texon\t500\t700\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    'chr1\tLRAA\texon\t801\t1000\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
)


@pytest.fixture(autouse=True)
def _restore_mode():
    saved = LRAA_Globals.LRAA_MODE
    yield
    LRAA_Globals.LRAA_MODE = saved


def _seeded_boundaries(tmp_path, lraa_mode, quant_mode):
    gtf = tmp_path / "g.gtf"
    gtf.write_text(_GUIDE_GTF)
    guides = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(str(gtf))["chr1"]

    LRAA_Globals.LRAA_MODE = lraa_mode
    sg = Splice_graph.Splice_graph()
    sg.build_splice_graph_for_contig(
        "chr1", "+", "C" * 2000,
        None,            # no bam: input_transcripts only, as the merge builds
        None, None,
        guides,
        quant_mode=quant_mode,
    )
    polyA = sorted(o.get_coords()[0] for o in sg._PolyA_objs)
    tss = sorted(o.get_coords()[0] for o in sg._TSS_objs)
    return polyA, tss


def test_discovery_does_not_seed_guide_boundaries(tmp_path):
    """quant_mode=False, non-merge: the guide is a structural hint; no TSS/PolyA vertex."""
    polyA, tss = _seeded_boundaries(tmp_path, "ID-init_norm_reads", quant_mode=False)
    assert polyA == []
    assert tss == []


def test_quant_only_seeds_guide_boundaries(tmp_path):
    """quant_mode=True: the GTF is the model set, so its boundaries become vertices."""
    polyA, tss = _seeded_boundaries(tmp_path, "ID-init_quant", quant_mode=True)
    assert polyA == [1000]
    assert tss == [500]


def test_merge_seeds_guide_boundaries(tmp_path):
    """LRAA_MODE=='MERGE' (quant_mode=False): the merge must still read TSS+PolyA off the
    input GTFs -- the explicit requirement this gate must not break."""
    polyA, tss = _seeded_boundaries(tmp_path, "MERGE", quant_mode=False)
    assert polyA == [1000]
    assert tss == [500]
