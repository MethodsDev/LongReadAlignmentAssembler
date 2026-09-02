import pytest

from Transcript import Transcript
import LRAA_Globals


def make_mock_transcript(idx, lend, rend):
    # Minimal Transcript construction: emulate a monoexonic transcript
    t = Transcript("chrTEST", [(lend, rend)], "+")
    t.set_simple_path(["E:shared"])
    t.set_gene_id(f"g{idx}")
    t.set_transcript_id(f"t{idx}")
    return t


def test_large_cluster_skips_leiden_and_uses_dsu(monkeypatch):
    # Force a large cluster triggering DSU fallback.
    #
    # monkeypatch, NOT direct assignment: LRAA_Globals.config is process-global, and a
    # threshold of 10 left behind sends every later test in the session down the DSU
    # path. That is not hypothetical -- assigning these directly made
    # test_recluster_extent_invariance fail whenever it ran after this file, because
    # its fixture cluster exceeds 10 and its recorded Leiden call never happened.
    monkeypatch.setitem(LRAA_Globals.config, "use_community_clustering", True)
    monkeypatch.setitem(
        LRAA_Globals.config, "max_transcripts_for_community_clustering", 10
    )

    transcripts = []
    # Create 25 transcripts with overlapping spans so they end up in one initial cluster
    for i in range(25):
        transcripts.append(make_mock_transcript(i, 100 + i * 5, 200 + i * 5))

    reclustered = Transcript.recluster_transcripts_to_genes(
        transcripts, contig_acc="chrTEST", contig_strand="+"
    )

    # Expect DSU path used: a single gene cluster mapping all transcripts
    gene_ids = {t.get_gene_id() for t in reclustered}
    assert len(gene_ids) == 1, "All transcripts should be reclustered into a single gene"
    assert len(reclustered) == 25
