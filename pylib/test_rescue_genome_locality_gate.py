#!/usr/bin/env python3

"""The genomic locality criterion on transcriptome read rescue.

Rescue realigns an unassigned read against local transcript sequences and credits it to
whichever target it scores best against. Nothing in that scoring knows where the read
is: on chr21 HG002 a read whose only primary alignment sat at 6,500,980 scored best
against a model at 43,109,302, 36.6 Mb away, and across that run 22 of the 33
(read, transcript) credits rescue made were to models the read does not overlap.

The read's genome alignment is trusted for LOCALITY but not for optimality: rescue
exists precisely because that alignment is a poor fit, so it refines placement within
the locus the read already occupies and never moves it to another one. That is why
plain coordinate overlap of the read's aligned blocks against the target's exons is
the whole criterion -- no minimum overlap, no fraction of the read, no intron-chain or
path-compatibility reasoning, all of which would be assertions about structure that
the suboptimal genome alignment cannot support.

These tests go through rescue_unassigned_reads_to_transcriptome, the entry point the
pipeline calls (pylib/LRAA.py:324 and the quant path in LRAA). A test that built its
own invocation instead would validate a flow no user reaches, and would keep passing
with the criterion absent from the one that ships.

The fixture is synthetic because the defect needs a read whose genome alignment and
whose best transcript match disagree about locality by megabases, and no small real
corpus is guaranteed to contain a clean one.
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import networkx as nx
import pysam

import IsoformReadRescue
from GenomeFeature import Exon, Intron
from IsoformReadRescue import (
    _normalize_read_identifier,
    rescue_unassigned_reads_to_transcriptome,
)
from LRAA import LRAA
from Splice_graph import Splice_graph
from Transcript import Transcript


CONTIG = "chr_locality"
CONTIG_LEN = 501000

# Gene A, where both reads' genome alignments are.
A_EXON_1 = (1001, 1300)
A_INTRON = (1301, 1600)
A_EXON_2 = (1601, 1900)

# Gene B, 0.5 Mb away, which one read's sequence matches exactly. Far enough that no
# alignment block of a read at gene A can reach it, which is the whole point.
B_EXON_1 = (500001, 500300)
B_INTRON = (500301, 500600)
B_EXON_2 = (500601, 500900)

LOCAL_READ = "read_matching_the_local_model"
DISTANT_READ = "read_matching_the_distant_model"


def _contig_seq():
    """A deterministic contig whose two genes share no sequence.

    Random rather than repetitive so that minimap2 has exactly one place to put each
    read: if the two transcripts resembled each other the test would be measuring
    aligner ambiguity instead of the gate.
    """
    import random

    rng = random.Random(20260817)
    return "".join(rng.choice("ACGT") for _ in range(CONTIG_LEN))


def _build_two_gene_fixture(contig_seq):
    Exon.reset_counter()
    splice_graph = Splice_graph()
    splice_graph._contig_acc = CONTIG
    splice_graph._contig_strand = "+"

    exon_a1 = Exon(CONTIG, *A_EXON_1, "+", 10)
    intron_a = Intron(CONTIG, *A_INTRON, "+", 10)
    exon_a2 = Exon(CONTIG, *A_EXON_2, "+", 10)
    exon_b1 = Exon(CONTIG, *B_EXON_1, "+", 10)
    intron_b = Intron(CONTIG, *B_INTRON, "+", 10)
    exon_b2 = Exon(CONTIG, *B_EXON_2, "+", 10)

    splice_graph._splice_graph = nx.DiGraph(
        [
            (exon_a1, intron_a),
            (intron_a, exon_a2),
            (exon_b1, intron_b),
            (intron_b, exon_b2),
        ]
    )
    splice_graph._intron_objs = {
        "{}:{}".format(*A_INTRON): intron_a,
        "{}:{}".format(*B_INTRON): intron_b,
    }
    splice_graph._finalize_splice_graph()

    def _transcript(gene_id, transcript_id, exons, nodes):
        transcript = Transcript(CONTIG, [list(coords) for coords in exons], "+")
        transcript.set_gene_id(gene_id)
        transcript.set_transcript_id(transcript_id)
        transcript.set_simple_path([node.get_id() for node in nodes])
        return transcript

    transcripts = [
        _transcript("geneA", "txA", [A_EXON_1, A_EXON_2], [exon_a1, intron_a, exon_a2]),
        _transcript("geneB", "txB", [B_EXON_1, B_EXON_2], [exon_b1, intron_b, exon_b2]),
    ]

    paths = {
        "A": [exon_a1.get_id(), intron_a.get_id(), exon_a2.get_id()],
        "B": [exon_b1.get_id(), intron_b.get_id(), exon_b2.get_id()],
    }
    return splice_graph, transcripts, paths


def _transcript_seq(contig_seq, exons):
    return "".join(contig_seq[lend - 1 : rend] for lend, rend in exons)


def _write_genome_bam(path, contig_seq):
    """Both reads spliced across gene A's exons; only their sequences differ.

    Identical coordinates on purpose: the two reads are indistinguishable to every
    part of rescue except the locality check and the aligner, so if collection dropped
    one of them it dropped both and the positive control below fails rather than the
    locality assertion passing for the wrong reason.

    NM:i:0 is deliberate. It makes the genome alignment explain the whole read at
    perfect gap-aware identity, which is the most permissive the pre-gate
    explained-bases and identity baseline can be -- so a rescue to gene B clears that
    baseline and locality is the only thing left that can decline it.
    """
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}]}
    spliced_cigar = [
        (0, A_EXON_1[1] - A_EXON_1[0] + 1),
        (3, A_INTRON[1] - A_INTRON[0] + 1),
        (0, A_EXON_2[1] - A_EXON_2[0] + 1),
    ]

    def _record(read_name, seq):
        aln = pysam.AlignedSegment()
        aln.query_name = read_name
        aln.flag = 0
        aln.reference_id = 0
        aln.reference_start = A_EXON_1[0] - 1
        aln.mapping_quality = 60
        aln.cigar = spliced_cigar
        aln.query_sequence = seq
        aln.set_tag("NM", 0)
        return aln

    with pysam.AlignmentFile(str(path), "wb", header=header) as bam_writer:
        bam_writer.write(
            _record(LOCAL_READ, _transcript_seq(contig_seq, [A_EXON_1, A_EXON_2]))
        )
        bam_writer.write(
            _record(DISTANT_READ, _transcript_seq(contig_seq, [B_EXON_1, B_EXON_2]))
        )
    pysam.index(str(path))


def _paths_by_read(rescued_mps):
    """Keyed by read name via the compact id MultiPath actually stores.

    get_read_names() only returns names when LRAA_Globals.READ_NAME_STORE is populated,
    which it is not in a unit test, and falls back to stringified ids there.
    """
    id_to_name = {
        _normalize_read_identifier(read_name): read_name
        for read_name in (LOCAL_READ, DISTANT_READ)
    }
    paths = {}
    for rescued_mp in rescued_mps:
        for read_id in rescued_mp.get_read_ids():
            read_name = id_to_name.get(read_id, read_id)
            paths.setdefault(read_name, []).append(rescued_mp.get_simple_path())
    return paths


@pytest.fixture
def two_gene_rescue(tmp_path):
    contig_seq = _contig_seq()
    splice_graph, transcripts, paths = _build_two_gene_fixture(contig_seq)
    bam_path = tmp_path / "genome.bam"
    _write_genome_bam(bam_path, contig_seq)
    return {
        "splice_graph": splice_graph,
        "transcripts": transcripts,
        "node_paths": paths,
        "contig_seq": contig_seq,
        "bam": str(bam_path),
        "read_path_mapper": LRAA(splice_graph)._map_read_to_graph,
    }


def test_rescue_declines_a_target_the_read_does_not_overlap(two_gene_rescue):
    """The pipeline's rescue must not move a read 0.5 Mb to a better-scoring model.

    Fails when the locality criterion is absent from
    rescue_unassigned_reads_to_transcriptome: the distant read matches gene B exactly
    and gene A not at all, so without it rescue credits the read to B.
    """
    rescued_mps = rescue_unassigned_reads_to_transcriptome(
        two_gene_rescue["splice_graph"],
        two_gene_rescue["transcripts"],
        two_gene_rescue["contig_seq"],
        two_gene_rescue["bam"],
        CONTIG,
        None,
        None,
        {LOCAL_READ, DISTANT_READ},
        read_path_mapper=two_gene_rescue["read_path_mapper"],
    )
    paths = _paths_by_read(rescued_mps)

    # Positive control first. An empty result would satisfy the locality assertion for
    # every reason other than the gate -- no minimap2 hit, an empty transcript fasta, a
    # projection that failed -- so the run has to be shown capable of rescuing at all.
    assert paths.get(LOCAL_READ) == [
        two_gene_rescue["node_paths"]["A"]
    ], "rescue produced nothing for the read that does overlap its target; the locality assertion below would then be vacuous"

    distant_b_nodes = set(two_gene_rescue["node_paths"]["B"])
    for path in paths.get(DISTANT_READ, []):
        assert not distant_b_nodes.intersection(
            path
        ), f"read at {A_EXON_1[0]}-{A_EXON_2[1]} was credited to the model at {B_EXON_1[0]}-{B_EXON_2[1]}: {path}"


def test_the_fixture_misassigns_when_the_locality_predicate_is_neutralised(
    two_gene_rescue, monkeypatch
):
    """Proof the assertion above is not vacuous: neutralise locality and B wins.

    An ablation of exactly one mechanism, through the same production entry point:
    _get_alignment_overlapping_targets is replaced by one that reports every target in
    the index, so every read is allowed everywhere and nothing else changes. If the
    misassignment does not reappear, the fixture cannot fail for locality reasons and
    the test above proves nothing, so this failing is as informative as that one.
    """

    def _permit_every_target(read, exon_overlap_index):
        return {interval.data: 1 for interval in exon_overlap_index}

    monkeypatch.setattr(
        IsoformReadRescue,
        "_get_alignment_overlapping_targets",
        _permit_every_target,
    )

    rescued_mps = rescue_unassigned_reads_to_transcriptome(
        two_gene_rescue["splice_graph"],
        two_gene_rescue["transcripts"],
        two_gene_rescue["contig_seq"],
        two_gene_rescue["bam"],
        CONTIG,
        None,
        None,
        {LOCAL_READ, DISTANT_READ},
        read_path_mapper=two_gene_rescue["read_path_mapper"],
    )
    paths = _paths_by_read(rescued_mps)

    assert paths.get(DISTANT_READ) == [
        two_gene_rescue["node_paths"]["B"]
    ], "neutralising locality did not reproduce the misassignment, so the test above cannot fail and proves nothing"


def test_locality_leaves_the_existing_genome_baseline_untouched(tmp_path):
    """One added criterion, nothing else disturbed.

    Locality was added to the acceptance path rescue already used rather than by switching
    rescue onto the genome-gated target flow, because that flow also emptied the
    explained-bases and gap-aware-identity baseline and applied read filters of its own.
    Measured on chr21 HiFi, that baseline declines 668 rescues where locality declines 16,
    so losing it silently would have been the larger change and not the one asked for.
    That gated flow has since been deleted along with the whole-genome-versus-whole-
    transcriptome modes it served, which is what makes these baselines unconditional.
    Asserts the same collector still reports both baselines alongside the new
    allowed-target sets, over the same offered reads.
    """
    from IsoformReadRescue import (
        _build_exon_overlap_index,
        _build_transcript_models,
        _collect_read_sequences,
    )

    contig_seq = _contig_seq()
    splice_graph, transcripts, _ = _build_two_gene_fixture(contig_seq)
    transcript_models = _build_transcript_models(splice_graph, transcripts, contig_seq)
    bam_path = tmp_path / "genome.bam"
    _write_genome_bam(bam_path, contig_seq)

    def _collect(index):
        return _collect_read_sequences(
            str(bam_path),
            CONTIG,
            None,
            None,
            {LOCAL_READ, DISTANT_READ},
            "+",
            exon_overlap_index=index,
        )

    without = _collect(None)
    withidx = _collect(_build_exon_overlap_index(transcript_models))

    # Both reads are spliced across gene A's two 300-base exons with NM:i:0, so the
    # genome alignment explains all 600 read bases at perfect gap-aware identity.
    assert withidx[1] == {LOCAL_READ: 600, DISTANT_READ: 600}
    assert withidx[2] == {LOCAL_READ: 1.0, DISTANT_READ: 1.0}
    # Same offered reads and same baselines with the criterion as without it.
    assert withidx[0] == without[0]
    assert withidx[1] == without[1]
    assert withidx[2] == without[2]
    # Only the allowed-target sets are new, and gene B is absent from both reads even
    # though one of them matches gene B's sequence exactly.
    assert without[3] is None
    assert withidx[3] == {LOCAL_READ: {"geneA^txA"}, DISTANT_READ: {"geneA^txA"}}


def test_locality_is_overlap_only_and_admits_a_single_shared_base(tmp_path):
    """The gate asks where the read is, not how much of the target it covers.

    A read reaching one base into a target's exon is admitted. Guards against the check
    being tightened into a bp threshold or a fraction-of-read minimum, which is a
    different acceptance policy from the one asked for.
    """
    from IsoformReadRescue import (
        _build_exon_overlap_index,
        _get_alignment_overlapping_targets,
    )

    contig_seq = _contig_seq()
    splice_graph, transcripts, _ = _build_two_gene_fixture(contig_seq)
    from IsoformReadRescue import _build_transcript_models

    transcript_models = _build_transcript_models(splice_graph, transcripts, contig_seq)
    index = _build_exon_overlap_index(transcript_models)

    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}]}
    aln = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(header))
    aln.query_name = "one_base_into_geneA"
    aln.flag = 0
    aln.reference_id = 0
    # Ends on A_EXON_1's first base: 100 bases of upstream sequence plus that one.
    aln.reference_start = A_EXON_1[0] - 101
    aln.mapping_quality = 60
    aln.cigar = [(0, 101)]
    aln.query_sequence = "A" * 101
    aln.set_tag("NM", 0)

    overlaps = _get_alignment_overlapping_targets(aln, index)
    assert overlaps == {"geneA^txA": 1}
