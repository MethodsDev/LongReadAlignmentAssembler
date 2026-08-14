#!/usr/bin/env python3

"""The guide annotation is an input to graph construction, so it must be in the cache key.

It was not. The cache root is per-output-prefix only, so two ref-guided runs sharing a
prefix but supplying different GTFs could reuse each other's graphs. That matters more
since `spare_polyA_veto_at_known_3prime`, because a reference 3' coordinate now decides
whether an A-rich PolyA candidate becomes a vertex at all -- so the same BAM and config
can legitimately produce different graphs from different guides.
"""

import importlib.machinery
import importlib.util
import os

import pytest

from Transcript import Transcript


def _load_lraa():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    loader = importlib.machinery.SourceFileLoader("lraa_cachekey_test", os.path.join(here, "LRAA"))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


lraa = _load_lraa()


def _guide(transcript_id, exons, strand="+"):
    transcript = Transcript("chr1", exons, strand)
    transcript.set_gene_id("g")
    transcript.set_transcript_id(transcript_id)
    return transcript


def _key(transcripts, reference_transcripts=None, quant_mode=False):
    cache_key, _signature = lraa._compute_splice_graph_cache_entry(
        "chr1",
        "+",
        "/tmp/nonexistent.bam",
        1,
        100000,
        transcripts,
        quant_mode,
        reference_transcripts=reference_transcripts,
    )
    return cache_key


def test_quant_mode_changes_the_key():
    """The assembly build and the final-quant build must not share a cache entry.

    quant_mode gates pruning during construction (Splice_graph.py:327-328, :332-333), so
    the two stages build different graphs from the same inputs.  While it was absent from
    the signature the keys collided and the quant stage loaded the discovery graph: every
    read then mapped to a SPACER path, was recorded as having no usable genome path, and
    was sent to transcriptome rescue for the wrong reason.  Guide digests happened to
    separate the two stages for ref-guided runs, which is why this went unnoticed -- so
    hold the invariant directly, with an identical transcript set on both sides.
    """
    guides = [_guide("t1", [[500, 1000], [1500, 2000]])]
    assert _key(guides, quant_mode=False) != _key(guides, quant_mode=True)


def test_changing_a_guide_three_prime_end_changes_the_key():
    """The coordinate the veto exemption reads. Same BAM, same config, different guide
    terminus -- the graphs genuinely differ, so the keys must."""
    a = _key([_guide("t1", [[500, 1000]])])
    b = _key([_guide("t1", [[500, 1200]])])
    assert a != b


def test_changing_internal_structure_changes_the_key():
    """Termini alone are what the veto reads, but intron chains shape the rest of the
    graph, so the digest covers exon segments too."""
    a = _key([_guide("t1", [[500, 700], [900, 1000]])])
    b = _key([_guide("t1", [[500, 750], [900, 1000]])])
    assert a != b


def test_adding_a_guide_transcript_changes_the_key():
    one = _key([_guide("t1", [[500, 1000]])])
    two = _key([_guide("t1", [[500, 1000]]), _guide("t2", [[2000, 2500]])])
    assert one != two


def test_the_same_guide_gives_a_stable_key():
    """Otherwise nothing would ever cache-hit."""
    assert _key([_guide("t1", [[500, 1000]])]) == _key([_guide("t1", [[500, 1000]])])


def test_guide_order_does_not_change_the_key():
    """Parse order is not meaningful; a key that depended on it would miss constantly."""
    forward = _key([_guide("t1", [[500, 1000]]), _guide("t2", [[2000, 2500]])])
    reverse = _key([_guide("t2", [[2000, 2500]]), _guide("t1", [[500, 1000]])])
    assert forward == reverse


def test_a_different_strand_guide_does_not_affect_this_strand_key():
    """The signature is per contig/strand, and only same-strand structures contribute."""
    plus_only = _key([_guide("t1", [[500, 1000]], "+")])
    plus_and_minus = _key(
        [_guide("t1", [[500, 1000]], "+"), _guide("t2", [[3000, 3500]], "-")]
    )
    assert plus_only == plus_and_minus


def test_ref_free_and_ref_guided_keys_differ():
    assert _key(None) != _key([_guide("t1", [[500, 1000]])])


def test_the_transcript_id_alone_does_not_change_the_key():
    """Structure decides the graph, not names; keying on ids would miss on a re-annotated
    but structurally identical guide."""
    assert _key([_guide("nameA", [[500, 1000]])]) == _key([_guide("nameB", [[500, 1000]])])


# ------------------------------------------------ provenance, not just structure


def test_the_same_structures_with_a_different_provenance_split_change_the_key():
    """The decisive case for spare_polyA_veto_at_known_3prime.

    The integrated set decides graph structure; the genuine reference subset decides
    which A-rich PolyA candidates keep their vertex. Identical combined structures can
    arise with different provenance -- one run where a model is reconstructed, another
    where the same coordinates come from the annotation -- and those build DIFFERENT
    graphs. A single digest over the mixture cannot tell them apart.
    """
    combined = [_guide("t1", [[500, 1000]]), _guide("t2", [[2000, 2500]])]
    reference_is_t1 = _key(combined, reference_transcripts=[combined[0]])
    reference_is_t2 = _key(combined, reference_transcripts=[combined[1]])
    assert reference_is_t1 != reference_is_t2


def test_no_reference_differs_from_a_reference_over_the_same_structures():
    combined = [_guide("t1", [[500, 1000]])]
    assert _key(combined, reference_transcripts=None) != _key(
        combined, reference_transcripts=combined
    )


# --------------------------------- boundary metadata also steers graph construction


def _parse_guides(tmp_path, attributes, name="g.gtf"):
    """Build guides the way production does -- through the GTF parser -- so the imported
    TSS/PolyA flags, their support counts and TPM are genuinely set. Constructing a bare
    Transcript and calling a setter does not reproduce that state: has_PolyA() reflects
    the imported annotation, not the read-count field.
    """
    from Transcript import GTF_contig_to_transcripts

    path = tmp_path / name
    path.write_text(
        'chr1\tLRAA\ttranscript\t500\t1000\t.\t+\t.\t'
        'gene_id "g1"; transcript_id "t1"; ' + attributes + '\n'
        'chr1\tLRAA\texon\t500\t1000\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    )
    return GTF_contig_to_transcripts.parse_GTF_to_Transcripts(str(path))["chr1"]


def test_an_annotated_polyA_flag_changes_the_key_at_identical_exon_segments(tmp_path):
    """Same single exon 500-1000 in both; only the annotated boundary differs, which
    _integrate_input_transcript_structures turns into PolyA evidence."""
    without = _parse_guides(tmp_path, 'TSS "False"; PolyA "False";', "a.gtf")
    with_polya = _parse_guides(tmp_path, 'TSS "False"; PolyA "True";', "b.gtf")
    assert [t.get_exon_segments() for t in without] == [
        t.get_exon_segments() for t in with_polya
    ], "fixtures must differ only in boundary metadata"
    assert _key(without) != _key(with_polya)


def test_annotated_boundary_support_changes_the_key(tmp_path):
    low = _parse_guides(tmp_path, 'TSS "True"; PolyA "True"; PolyA_read_count "3";', "c.gtf")
    high = _parse_guides(tmp_path, 'TSS "True"; PolyA "True"; PolyA_read_count "99";', "d.gtf")
    assert _key(low) != _key(high)


def _remap(transcripts):
    """Reproduce what LRAA.assign_transcripts_paths_in_graph does to these objects.

    It sets a simple path and calls refresh_boundary_annotations_from_simple_path, which
    clears _imported_has_TSS/_imported_has_POLYA.  build_ME_transcripts and
    build_SE_transcripts do this to the reference objects IN PLACE, before
    run_transcript_assembly computes the final graph's cache entry over those same
    objects -- so this, not the freshly parsed state, is what the digest usually sees.
    """
    for transcript in transcripts:
        transcript.set_simple_path(["E:1"])
        transcript.refresh_boundary_annotations_from_simple_path()
    return transcripts


def test_annotated_boundary_metadata_still_changes_the_key_after_remapping(tmp_path):
    without = _remap(_parse_guides(tmp_path, 'TSS "False"; PolyA "False";', "e.gtf"))
    with_polya = _remap(_parse_guides(tmp_path, 'TSS "False"; PolyA "True";', "f.gtf"))
    assert _key(without) != _key(with_polya)


def test_graph_derived_boundary_state_does_not_enter_the_key():
    """has_TSS/has_PolyA read the simple path, which is an OUTPUT of a splice graph.

    Keying on it makes the key for a graph depend on a graph, and on a transcript with no
    simple path those accessors assert -- the digest then records the same "ERR" for
    every record and the field distinguishes nothing.
    """
    plain = _guide("t1", [[500, 1000]])
    with_path = _guide("t1", [[500, 1000]])
    with_path.set_simple_path(["POLYA:chr1:1000:+", "E:1"])
    assert _key([plain]) == _key([with_path])
