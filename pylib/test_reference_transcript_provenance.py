#!/usr/bin/env python3

import pytest
import networkx as nx

import LRAA_Globals
import TranscriptFiltering
from GenomeFeature import Exon, Intron
from LRAA import LRAA
from Splice_graph import Splice_graph
from Transcript import Transcript


def _splice_graph(strand="+"):
    splice_graph = Splice_graph()
    splice_graph._contig_acc = "chr1"
    splice_graph._contig_strand = strand
    return splice_graph


def _transcript(transcript_id, path, exons=None, strand="+"):
    if exons is None:
        exons = [[100, 150], [200, 250]]
    transcript = Transcript("chr1", exons, strand)
    transcript.set_transcript_id(transcript_id)
    transcript.set_simple_path(path)
    return transcript


@pytest.mark.parametrize(
    ("strand", "draft_path", "reference_path"),
    [
        (
            "+",
            ["TSS:draft", "E:upstream", "E:1", "I:1", "E:2", "POLYA:draft"],
            ["TSS:reference", "E:1", "I:1", "E:2", "POLYA:reference"],
        ),
        (
            "-",
            ["POLYA:draft", "E:2", "I:1", "E:1", "E:upstream", "TSS:draft"],
            ["POLYA:reference", "E:2", "I:1", "E:1", "TSS:reference"],
        ),
    ],
)
def test_reference_provenance_ignores_reference_terminal_boundaries(
    strand, draft_path, reference_path
):
    lraa = LRAA(_splice_graph(strand))
    draft = _transcript("draft", draft_path, strand=strand)
    reference = _transcript("reference", reference_path, strand=strand)

    lraa.assign_reference_transcript_provenance([draft], [reference])

    assert draft.includes_reference_transcript()
    assert draft.get_source_reference_transcript_ids() == {"reference"}


def test_reference_provenance_requires_contiguous_containment():
    lraa = LRAA(_splice_graph())
    draft = _transcript("draft", ["E:1", "I:1", "E:2"])
    reference = _transcript(
        "reference", ["TSS:reference", "E:1", "I:other", "E:2", "POLYA:reference"]
    )

    lraa.assign_reference_transcript_provenance([draft], [reference])

    assert not draft.includes_reference_transcript()
    assert draft.get_source_reference_transcript_ids() == set()


def test_reference_provenance_rejects_internal_boundary_nodes():
    lraa = LRAA(_splice_graph())
    draft = _transcript("draft", ["E:1", "I:1", "E:2"])
    reference = _transcript(
        "reference",
        ["TSS:terminal", "E:1", "TSS:internal", "I:1", "E:2", "POLYA:terminal"],
    )

    lraa.assign_reference_transcript_provenance([draft], [reference])

    assert not draft.includes_reference_transcript()


def test_failed_remap_clears_stale_path_and_excludes_rescue_target(monkeypatch):
    monkeypatch.setitem(LRAA_Globals.config, "show_progress_assign_transcripts", False)
    monkeypatch.setitem(LRAA_Globals.config, "infer_TSS", False)
    monkeypatch.setitem(LRAA_Globals.config, "infer_PolyA", False)

    splice_graph = _splice_graph()
    exon_left = Exon("chr1", 100, 150, "+", 1)
    intron = Intron("chr1", 151, 199, "+", 1)
    exon_right = Exon("chr1", 200, 250, "+", 1)
    splice_graph._splice_graph = nx.DiGraph([(exon_left, intron), (intron, exon_right)])
    splice_graph._intron_objs = {"151:199": intron}
    splice_graph._finalize_splice_graph()

    mapped_draft = _transcript("mapped", ["E:old"])
    unmapped_draft = _transcript(
        "unmapped", ["E:stale"], exons=[[300, 350], [400, 450]]
    )
    lraa = LRAA(splice_graph)

    mapped = lraa.assign_transcripts_paths_in_graph([mapped_draft, unmapped_draft])
    mapped_object_ids = {id(transcript) for transcript in mapped}
    # Rescue targets mirror the post-assembly quant step: the draft models that
    # mapped into the final splice graph, which is also the set being quantified.
    rescue_targets = [
        transcript
        for transcript in (mapped_draft, unmapped_draft)
        if id(transcript) in mapped_object_ids
    ]

    assert mapped == [mapped_draft]
    assert rescue_targets == [mapped_draft]
    with pytest.raises(AssertionError):
        unmapped_draft.get_simple_path()


def test_reference_provenance_survives_quant_reset_and_drives_filter(monkeypatch):
    monkeypatch.setitem(LRAA_Globals.config, "HiFi", False)
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 1)
    monkeypatch.setitem(
        LRAA_Globals.config, "ref_trans_filter_mode", "retain_expressed"
    )
    transcript = _transcript("draft", ["E:1"], exons=[[100, 150]])
    transcript.set_source_reference_transcript_ids({"reference"})

    transcript.init_quant_info()
    transcript.set_read_counts_assigned(1)
    retained = TranscriptFiltering.filter_monoexonic_isoforms_by_TPM_threshold(
        [transcript], min_TPM=2e6
    )

    assert retained == [transcript]
    assert transcript.get_source_reference_transcript_ids() == {"reference"}
