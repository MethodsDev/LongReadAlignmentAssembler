#!/usr/bin/env python3

import importlib.util
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

from GenomeFeature import Exon
from IsoformReadRescue import _normalize_read_identifier
from MultiPath import MultiPath
from MultiPathCounter import MultiPathCounter
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
