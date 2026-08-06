#!/usr/bin/env python3

import os
import sys
from pathlib import Path

import pysam

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA_Globals

# Splice_graph must be imported first: it and Pretty_alignment_manager import each other.
from Splice_graph import Splice_graph
from Pretty_alignment_manager import Pretty_alignment_manager


def _alignment(read_name, flag, start, cigar, seq, mismatches):
    aln = pysam.AlignedSegment()
    aln.query_name = read_name
    aln.flag = flag
    aln.reference_id = 0
    aln.reference_start = start
    aln.mapping_quality = 60
    aln.cigar = cigar
    aln.query_sequence = seq
    aln.set_tag("NM", mismatches)
    return aln


def _write_bam(path, alignments):
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam_writer:
        for alignment in alignments:
            bam_writer.write(alignment)
    pysam.index(str(path))


def _retrieve(bam_path, cache_dir):
    manager = Pretty_alignment_manager(
        Splice_graph(), alignment_cache_dir=str(cache_dir)
    )
    alignments = manager.retrieve_pretty_alignments(
        "chr1", "+", None, str(bam_path), use_cache=True
    )
    return manager, alignments


def _low_per_id_bam(tmp_path, monkeypatch):
    monkeypatch.setitem(LRAA_Globals.config, "min_per_id", 80)
    bam_path = tmp_path / "reads.bam"
    _write_bam(
        bam_path,
        [
            _alignment("kept", 0, 100, [(0, 50)], "A" * 50, 0),
            _alignment("low_identity", 0, 200, [(0, 50)], "C" * 50, 25),
        ],
    )
    return bam_path


def test_discarded_read_names_survive_a_cache_hit(tmp_path, monkeypatch):
    bam_path = _low_per_id_bam(tmp_path, monkeypatch)
    cache_dir = tmp_path / "cache"

    cold_manager, cold_alignments = _retrieve(bam_path, cache_dir)
    warm_manager, warm_alignments = _retrieve(bam_path, cache_dir)

    assert len(cold_alignments) == 1
    assert len(warm_alignments) == len(cold_alignments)

    cold_discards = cold_manager.get_last_discarded_read_names_by_reason()
    warm_discards = warm_manager.get_last_discarded_read_names_by_reason()
    assert cold_discards.get("low_perID") == {"low_identity"}
    assert warm_discards == cold_discards


def test_cache_without_discard_provenance_is_reextracted(tmp_path, monkeypatch):
    bam_path = _low_per_id_bam(tmp_path, monkeypatch)
    cache_dir = tmp_path / "cache"

    _retrieve(bam_path, cache_dir)

    discard_caches = sorted(cache_dir.glob("*.discards.pkl*"))
    assert discard_caches
    for path in discard_caches:
        os.remove(path)

    manager, alignments = _retrieve(bam_path, cache_dir)
    assert len(alignments) == 1
    assert manager.get_last_discarded_read_names_by_reason().get("low_perID") == {
        "low_identity"
    }
