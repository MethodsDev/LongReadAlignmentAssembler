import pytest

from Transcript import Transcript
import LRAA_Globals


def make_mock_transcript(idx, lend, rend):
    # Minimal Transcript construction: emulate a monoexonic transcript
    t = Transcript(
        contig_acc="chrTEST",
        lend=lend,
        rend=rend,
        orient="+",
        exon_segments=[(lend, rend)],
        simple_path=[idx],
        scored_path_obj=None,
    )
    t.set_gene_id(f"g{idx}")
    t.set_transcript_id(f"t{idx}")
    return t


def test_large_cluster_skips_leiden_and_uses_dsu():
    # Force a large cluster triggering DSU fallback.
    LRAA_Globals.config["use_community_clustering"] = True
    LRAA_Globals.config["max_transcripts_for_community_clustering"] = 10  # very small threshold for test

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
