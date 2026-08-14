#!/usr/bin/env python3

import logging
import importlib.util
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path
import itertools
from collections import defaultdict
from contextlib import contextmanager

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA_Globals
from GenomeFeature import Exon
from IsoformReadRescue import _normalize_read_identifier
from MultiPath import MultiPath
from MultiPathCounter import MultiPathCounter
from Quantify import Quantify
from Splice_graph import Splice_graph


def _load_lraa_cli_module():
    loader = SourceFileLoader("lraa_cli_test_module", str(REPO_ROOT / "LRAA"))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def _build_splice_graph():
    Exon.reset_counter()
    splice_graph = Splice_graph()
    e1 = Exon("chr1", 100, 150, "+", 1)
    e2 = Exon("chr1", 200, 250, "+", 1)
    e3 = Exon("chr1", 300, 350, "+", 1)
    for exon in (e1, e2, e3):
        splice_graph._node_id_to_node[exon.get_id()] = exon
    return splice_graph, e1, e2, e3


def _make_mp_counter(splice_graph, simple_path, read_name):
    mp_counter = MultiPathCounter()
    multipath = MultiPath(
        splice_graph,
        [simple_path],
        read_types={"PacBio"},
        read_names={read_name},
        read_count=1,
    )
    mp_counter.add(multipath)
    return mp_counter


def _make_tx_multipath(splice_graph, simple_path, read_name):
    return MultiPath(
        splice_graph,
        [simple_path],
        read_types={"PacBio"},
        read_names={read_name},
        read_count=1,
    )


def test_read_overlap_threshold_resolves_config_at_assignment_time():
    original_fraction = LRAA_Globals.config["fraction_read_align_overlap"]
    original_aggressive_assignment = LRAA_Globals.config["aggressively_assign_reads"]
    original_show_progress = LRAA_Globals.config["show_progress_quant_assign"]

    def build_assignment_fixture():
        splice_graph, e1, e2, _ = _build_splice_graph()
        splice_graph._contig_acc = "chr1"
        splice_graph._contig_strand = "+"

        transcript = MultiPath(splice_graph, [[e1.get_id()]]).toTranscript()
        transcript.set_gene_id("g1")
        transcript.set_transcript_id("t1")

        read_multipath = MultiPath(
            splice_graph,
            [[e1.get_id(), e2.get_id()]],
            read_types={"PacBio"},
            read_names={1},
            read_count=1,
        )
        mp_counter = MultiPathCounter()
        mp_counter.add(read_multipath)

        quantify = Quantify(False, 1)
        quantify._assign_path_nodes_to_gene([transcript])
        return quantify, splice_graph, mp_counter, read_multipath, transcript

    try:
        LRAA_Globals.config["aggressively_assign_reads"] = False
        LRAA_Globals.config["show_progress_quant_assign"] = False
        LRAA_Globals.config["fraction_read_align_overlap"] = 0.4

        quantify, splice_graph, mp_counter, read_multipath, transcript = (
            build_assignment_fixture()
        )
        quantify._assign_reads_to_transcripts(splice_graph, mp_counter)
        assert quantify._mp_to_transcripts[read_multipath] == [transcript]

        quantify, splice_graph, mp_counter, read_multipath, _ = (
            build_assignment_fixture()
        )
        quantify._assign_reads_to_transcripts(
            splice_graph, mp_counter, fraction_read_align_overlap=0.6
        )
        assert read_multipath not in quantify._mp_to_transcripts
    finally:
        LRAA_Globals.config["fraction_read_align_overlap"] = original_fraction
        LRAA_Globals.config["aggressively_assign_reads"] = (
            original_aggressive_assignment
        )
        LRAA_Globals.config["show_progress_quant_assign"] = original_show_progress


def test_explode_mp_counter_normalizes_genome_read_keys():
    lraa_cli = _load_lraa_cli_module()
    splice_graph, e1, e2, _ = _build_splice_graph()
    raw_read_name = "LR44258.383_AT5G48000.3"
    genome_mp_counter = _make_mp_counter(
        splice_graph, [e1.get_id(), e2.get_id()], raw_read_name
    )

    read_name_to_multipaths = lraa_cli._explode_mp_counter_to_read_multipaths(
        splice_graph, genome_mp_counter
    )

    assert set(read_name_to_multipaths.keys()) == {
        _normalize_read_identifier(raw_read_name)
    }


def test_genome_tx_arb_tied_per_id_keeps_genome_for_same_normalized_read():
    lraa_cli = _load_lraa_cli_module()
    splice_graph, e1, e2, e3 = _build_splice_graph()
    raw_read_name = "LR27826.152_AT3G09440.4"
    read_key = _normalize_read_identifier(raw_read_name)

    genome_mp_counter = _make_mp_counter(
        splice_graph, [e1.get_id(), e2.get_id()], raw_read_name
    )
    tx_multipath = _make_tx_multipath(
        splice_graph, [e1.get_id(), e3.get_id()], raw_read_name
    )

    chosen_mp_counter, stats = lraa_cli._arbitrate_genome_vs_transcriptome_read_paths(
        splice_graph,
        genome_mp_counter,
        {
            "read_name_to_multipaths": {read_key: [tx_multipath]},
            "read_name_to_best_per_id": {read_key: 1.0},
            "read_name_to_primary_per_id": {read_key: 1.0},
        },
    )

    assert stats["reads_total"] == 1
    assert stats["reads_kept_genome"] == 1
    assert stats["reads_selected_tx_total"] == 0
    assert stats["reads_selected_tx_missing_genome"] == 0
    assert stats["reads_tx_present_but_kept_genome"] == 1
    assert stats["reads_tx_tied_per_id_kept_genome"] == 1
    chosen_multipaths = lraa_cli._explode_mp_counter_to_read_multipaths(
        splice_graph, chosen_mp_counter
    )
    assert list(chosen_multipaths.keys()) == [read_key]
    assert len(chosen_multipaths[read_key]) == 1
    assert chosen_multipaths[read_key][0].get_simple_path() == [
        e1.get_id(),
        e2.get_id(),
    ]


def test_genome_tx_arb_normalizes_failed_genome_read_names():
    lraa_cli = _load_lraa_cli_module()
    splice_graph, e1, e2, e3 = _build_splice_graph()
    raw_read_name = "LR1617.1_AT1G17360.1"
    read_key = _normalize_read_identifier(raw_read_name)

    genome_mp_counter = _make_mp_counter(
        splice_graph, [e1.get_id(), e2.get_id()], raw_read_name
    )
    tx_multipath = _make_tx_multipath(
        splice_graph, [e1.get_id(), e3.get_id()], raw_read_name
    )

    chosen_mp_counter, stats = lraa_cli._arbitrate_genome_vs_transcriptome_read_paths(
        splice_graph,
        genome_mp_counter,
        {
            "read_name_to_multipaths": {read_key: [tx_multipath]},
            "read_name_to_best_per_id": {read_key: 0.95},
            "read_name_to_primary_per_id": {read_key: 0.99},
        },
        failed_genome_read_names={raw_read_name},
    )

    assert stats["reads_total"] == 1
    assert stats["reads_selected_tx_total"] == 1
    assert stats["reads_selected_tx_failed_genome"] == 1
    assert stats["reads_selected_tx_missing_genome"] == 0
    assert stats["reads_kept_genome"] == 0
    chosen_multipaths = lraa_cli._explode_mp_counter_to_read_multipaths(
        splice_graph, chosen_mp_counter
    )
    assert list(chosen_multipaths.keys()) == [read_key]
    assert chosen_multipaths[read_key][0].get_simple_path() == [
        e1.get_id(),
        e3.get_id(),
    ]


def test_genome_tx_arb_higher_tx_per_id_is_not_miscounted_as_missing_genome():
    lraa_cli = _load_lraa_cli_module()
    splice_graph, e1, e2, e3 = _build_splice_graph()
    raw_read_name = "LR571.1000_AT1G06640.3"
    read_key = _normalize_read_identifier(raw_read_name)

    genome_mp_counter = _make_mp_counter(
        splice_graph, [e1.get_id(), e2.get_id()], raw_read_name
    )
    tx_multipath = _make_tx_multipath(
        splice_graph, [e1.get_id(), e3.get_id()], raw_read_name
    )

    chosen_mp_counter, stats = lraa_cli._arbitrate_genome_vs_transcriptome_read_paths(
        splice_graph,
        genome_mp_counter,
        {
            "read_name_to_multipaths": {read_key: [tx_multipath]},
            "read_name_to_best_per_id": {read_key: 0.99},
            "read_name_to_primary_per_id": {read_key: 0.95},
        },
    )

    assert stats["reads_total"] == 1
    assert stats["reads_selected_tx_total"] == 1
    assert stats["reads_selected_tx_higher_per_id"] == 1
    assert stats["reads_selected_tx_missing_genome"] == 0
    assert stats["reads_kept_genome"] == 0
    chosen_multipaths = lraa_cli._explode_mp_counter_to_read_multipaths(
        splice_graph, chosen_mp_counter
    )
    assert list(chosen_multipaths.keys()) == [read_key]
    assert chosen_multipaths[read_key][0].get_simple_path() == [
        e1.get_id(),
        e3.get_id(),
    ]


def _build_single_transcript_fixture():
    """One transcript on the first exon; the third exon belongs to no transcript."""
    splice_graph, e1, _, e3 = _build_splice_graph()
    splice_graph._contig_acc = "chr1"
    splice_graph._contig_strand = "+"

    transcript = MultiPath(splice_graph, [[e1.get_id()]]).toTranscript()
    transcript.set_gene_id("g1")
    transcript.set_transcript_id("t1")

    quantify = Quantify(False, 1)
    quantify._assign_path_nodes_to_gene([transcript])
    return quantify, splice_graph, e1, e3


def _assign_read_on_exon(quantify, splice_graph, exon, read_count):
    read_multipath = MultiPath(
        splice_graph,
        [[exon.get_id()]],
        read_types={"PacBio"},
        read_names={1},
        read_count=read_count,
    )
    mp_counter = MultiPathCounter()
    mp_counter.add(read_multipath)
    quantify._assign_reads_to_transcripts(splice_graph, mp_counter)
    return read_multipath


def _warnings_from(caplog):
    return [
        record.getMessage()
        for record in caplog.records
        if record.levelno >= logging.WARNING
    ]


def test_unanchored_read_paths_report_the_reads_they_drop(caplog):
    original_show_progress = LRAA_Globals.config["show_progress_quant_assign"]
    try:
        LRAA_Globals.config["show_progress_quant_assign"] = False

        quantify, splice_graph, _, unused_exon = _build_single_transcript_fixture()
        with caplog.at_level(logging.WARNING):
            unanchored = _assign_read_on_exon(quantify, splice_graph, unused_exon, 4)

        assert unanchored not in quantify._mp_to_transcripts
        # the reads carried by the dropped path are reported, not just the path count
        assert any("carrying 4 reads" in message for message in _warnings_from(caplog))
    finally:
        LRAA_Globals.config["show_progress_quant_assign"] = original_show_progress


def test_anchored_read_paths_do_not_report_dropped_reads(caplog):
    original_show_progress = LRAA_Globals.config["show_progress_quant_assign"]
    try:
        LRAA_Globals.config["show_progress_quant_assign"] = False

        quantify, splice_graph, transcript_exon, _ = _build_single_transcript_fixture()
        with caplog.at_level(logging.WARNING):
            anchored = _assign_read_on_exon(quantify, splice_graph, transcript_exon, 4)

        assert quantify._mp_to_transcripts[anchored]
        assert not any("matched no gene" in m for m in _warnings_from(caplog))
    finally:
        LRAA_Globals.config["show_progress_quant_assign"] = original_show_progress


###############################################################
## EM over read-sharing components of genes, not over one gene
###############################################################


@contextmanager
def _progress_silenced():
    original_show_progress = LRAA_Globals.config["show_progress_quant_assign"]
    LRAA_Globals.config["show_progress_quant_assign"] = False
    try:
        yield
    finally:
        LRAA_Globals.config["show_progress_quant_assign"] = original_show_progress


_EXON_COORDS = ((100, 150), (200, 250), (300, 350), (400, 450), (500, 550), (600, 650))


def _build_exon_splice_graph():
    """Splice graph of six non-overlapping exons, addressed by 1-based index."""
    Exon.reset_counter()
    splice_graph = Splice_graph()
    splice_graph._contig_acc = "chr1"
    splice_graph._contig_strand = "+"
    exons = dict()
    for i, (lend, rend) in enumerate(_EXON_COORDS, start=1):
        exon = Exon("chr1", lend, rend, "+", 1)
        splice_graph._node_id_to_node[exon.get_id()] = exon
        exons[i] = exon
    return splice_graph, exons


def _make_transcript(splice_graph, exons, exon_indices, gene_id, transcript_id):
    node_ids = [exons[i].get_id() for i in exon_indices]
    transcript = MultiPath(splice_graph, [node_ids]).toTranscript()
    transcript.set_gene_id(gene_id)
    transcript.set_transcript_id(transcript_id)
    return transcript


def _make_read_multipath(splice_graph, exons, exon_indices, read_id, read_count=1):
    """A multipath standing for `read_count` reads, with that many distinct read ids.

    The ids matter, not just the count. Support is measured with get_read_weight(), which
    sums the per-read normalization weight and falls back to the number of distinct read
    ids when nothing was normalized. Handing one id a count of six would make weight and
    count disagree -- a state the real pipeline does not produce (verified: zero diverging
    multipaths over a real unnormalized chr20 region) but which a synthetic shortcut can
    invent, and which would then quietly change what EM is given here and nowhere else.

    Ids are spread by `read_id * 1000 + k` so two multipaths asking for overlapping ranges
    cannot collide and silently merge each other's reads.
    """
    node_ids = [exons[i].get_id() for i in exon_indices]
    return MultiPath(
        splice_graph,
        [node_ids],
        read_types={"PacBio"},
        read_names={read_id * 1000 + k for k in range(read_count)},
        read_count=read_count,
    )


def _counter_of(multipaths):
    mp_counter = MultiPathCounter()
    for multipath in multipaths:
        mp_counter.add(multipath)
    return mp_counter


def _gene_to_transcripts(transcripts):
    gene_to_transcripts = defaultdict(list)
    for transcript in transcripts:
        gene_to_transcripts[transcript.get_gene_id()].append(transcript)
    return gene_to_transcripts


def _build_overlapping_gene_fixture(include_shared_read):
    """A readthrough-style overlap: geneA and geneB both carry exon 3.

    With include_shared_read, one read lands on exon 3 alone and is therefore
    compatible with transcripts of both genes -- the population the component EM
    exists for.  Without it the two genes share no read at all.
    """
    splice_graph, exons = _build_exon_splice_graph()

    tA1 = _make_transcript(splice_graph, exons, (1, 2, 3), "geneA", "tA1")
    tA2 = _make_transcript(splice_graph, exons, (1, 3), "geneA", "tA2")
    tB1 = _make_transcript(splice_graph, exons, (3, 4), "geneB", "tB1")

    read_multipaths = [
        _make_read_multipath(splice_graph, exons, (1, 2, 3), 2, read_count=6),
        _make_read_multipath(splice_graph, exons, (1, 3), 3, read_count=3),
        _make_read_multipath(splice_graph, exons, (3, 4), 4, read_count=1),
    ]
    shared_multipath = None
    if include_shared_read:
        shared_multipath = _make_read_multipath(splice_graph, exons, (3,), 1)
        read_multipaths.insert(0, shared_multipath)

    return (
        splice_graph,
        [tA1, tA2, tB1],
        _counter_of(read_multipaths),
        shared_multipath,
    )


def _quantified_by_transcript_id(transcripts):
    return {
        transcript.get_transcript_id(): (
            transcript.get_read_counts_assigned(),
            transcript.get_isoform_fraction(),
        )
        for transcript in transcripts
    }


def test_genes_sharing_a_compatible_read_are_one_component_quantified_jointly():
    with _progress_silenced():
        splice_graph, transcripts, mp_counter, shared_multipath = (
            _build_overlapping_gene_fixture(include_shared_read=True)
        )
        quantify = Quantify(True, 1000)
        quantify.quantify(splice_graph, transcripts, mp_counter)

        # one multipath compatible with transcripts of both genes is the edge
        assert quantify._mp_to_transcripts[shared_multipath] == [
            transcripts[1],
            transcripts[0],
            transcripts[2],
        ]
        assert quantify._build_read_sharing_gene_components(
            _gene_to_transcripts(transcripts)
        ) == [["geneA", "geneB"]]

        # the joint EM apportions the shared read by relative expression, so the
        # far weaker geneB takes a sliver of it rather than the whole read that
        # an independent per-gene EM would have handed it
        assert 1.0 < transcripts[2].get_read_counts_assigned() < 1.2


def test_genes_sharing_no_read_stay_separate_and_match_independent_em():
    with _progress_silenced():
        splice_graph, transcripts, mp_counter, _ = _build_overlapping_gene_fixture(
            include_shared_read=False
        )
        quantify = Quantify(True, 1000)
        quantify.quantify(splice_graph, transcripts, mp_counter)

        assert quantify._build_read_sharing_gene_components(
            _gene_to_transcripts(transcripts)
        ) == [["geneA"], ["geneB"]]
        jointly_run = _quantified_by_transcript_id(transcripts)

        # the same fixture one gene at a time: exactly what per-gene EM did
        independently_run = dict()
        for gene_id in ("geneA", "geneB"):
            splice_graph, transcripts, mp_counter, _ = _build_overlapping_gene_fixture(
                include_shared_read=False
            )
            gene_transcripts = [
                transcript
                for transcript in transcripts
                if transcript.get_gene_id() == gene_id
            ]
            Quantify(True, 1000).quantify(splice_graph, gene_transcripts, mp_counter)
            independently_run.update(_quantified_by_transcript_id(gene_transcripts))

        # exact equality, not approximate: singleton components must not drift
        assert jointly_run == independently_run


def test_shared_read_assignment_sums_to_one_across_both_genes():
    with _progress_silenced():
        splice_graph, transcripts, mp_counter, shared_multipath = (
            _build_overlapping_gene_fixture(include_shared_read=True)
        )
        frac_assignments = Quantify(True, 1000).quantify(
            splice_graph, transcripts, mp_counter
        )

        genes_taking_a_share = {
            transcript.get_gene_id()
            for transcript in transcripts
            if frac_assignments[transcript.get_transcript_id()].get(shared_multipath)
        }
        assert genes_taking_a_share == {"geneA", "geneB"}

        assert sum(
            frac_assignments[transcript.get_transcript_id()][shared_multipath]
            for transcript in transcripts
        ) == pytest.approx(1.0, abs=1e-9)


def test_isoform_fraction_is_scoped_to_its_own_gene_not_the_component():
    with _progress_silenced():
        splice_graph, transcripts, mp_counter, _ = _build_overlapping_gene_fixture(
            include_shared_read=True
        )
        Quantify(True, 1000).quantify(splice_graph, transcripts, mp_counter)

        tA1, tA2, tB1 = transcripts
        component_reads = sum(t.get_read_counts_assigned() for t in transcripts)
        gene_a_reads = tA1.get_read_counts_assigned() + tA2.get_read_counts_assigned()

        for transcript in (tA1, tA2):
            assert transcript.get_isoform_fraction() == pytest.approx(
                transcript.get_read_counts_assigned() / gene_a_reads
            )
            # the component-scoped number is a visibly different number
            assert transcript.get_isoform_fraction() != pytest.approx(
                transcript.get_read_counts_assigned() / component_reads, abs=1e-3
            )

        # geneB's sole isoform owns all of geneB whatever else the component holds
        assert tB1.get_isoform_fraction() == pytest.approx(1.0)
        assert tB1.get_read_counts_assigned() / component_reads < 0.2


def _build_four_gene_fixture():
    """geneA and geneB share exon 3; geneZ and geneM share nothing with anyone.

    Gene ids are deliberately out of step with coordinate order -- geneZ is
    leftmost and geneM rightmost -- so a component list ordered by name would not
    match one ordered by coordinate.
    """
    splice_graph, exons = _build_exon_splice_graph()
    transcripts = [
        _make_transcript(splice_graph, exons, (1,), "geneZ", "tZ"),
        _make_transcript(splice_graph, exons, (2, 3), "geneA", "tA"),
        _make_transcript(splice_graph, exons, (3, 4), "geneB", "tB"),
        _make_transcript(splice_graph, exons, (5, 6), "geneM", "tM"),
    ]
    mp_counter = _counter_of(
        [
            _make_read_multipath(splice_graph, exons, (1,), 1, read_count=2),
            _make_read_multipath(splice_graph, exons, (2, 3), 2, read_count=4),
            _make_read_multipath(splice_graph, exons, (3,), 3),
            _make_read_multipath(splice_graph, exons, (3, 4), 4, read_count=3),
            _make_read_multipath(splice_graph, exons, (5, 6), 5, read_count=2),
        ]
    )
    return splice_graph, transcripts, mp_counter


def test_component_construction_is_deterministic_under_shuffled_transcript_order():
    expected = [["geneZ"], ["geneA", "geneB"], ["geneM"]]

    with _progress_silenced():
        for permuted in itertools.permutations(range(4)):
            splice_graph, transcripts, mp_counter = _build_four_gene_fixture()
            shuffled = [transcripts[i] for i in permuted]

            quantify = Quantify(True, 1000)
            quantify._assign_path_nodes_to_gene(shuffled)
            quantify._assign_reads_to_transcripts(splice_graph, mp_counter)

            assert (
                quantify._build_read_sharing_gene_components(
                    _gene_to_transcripts(shuffled)
                )
                == expected
            )


def _quant_log_messages(caplog, include_shared_read):
    caplog.clear()
    with _progress_silenced():
        splice_graph, transcripts, mp_counter, _ = _build_overlapping_gene_fixture(
            include_shared_read=include_shared_read
        )
        with caplog.at_level(logging.INFO, logger="Quantify"):
            Quantify(True, 1000).quantify(splice_graph, transcripts, mp_counter)
    return [record.getMessage() for record in caplog.records]


def test_joint_em_is_named_in_the_log_and_single_gene_lines_are_unchanged(caplog):
    # a joint EM is exactly the situation someone chasing an unexpected number
    # needs to see, so it names its member genes and its combined span
    joint = _quant_log_messages(caplog, include_shared_read=True)
    assert any(
        "quant estimates for isoforms of read-sharing component of 2 genes "
        "{geneA,geneB} chr1+:100-450" in message
        for message in joint
    )

    # and an ordinary run reads exactly as it did before
    singletons = _quant_log_messages(caplog, include_shared_read=False)
    assert any(
        "quant estimates for isoforms of geneA chr1+:100-350" in message
        for message in singletons
    )
    assert any(
        "quant estimates for isoforms of geneB chr1+:300-450" in message
        for message in singletons
    )
    assert not any("read-sharing component" in message for message in singletons)