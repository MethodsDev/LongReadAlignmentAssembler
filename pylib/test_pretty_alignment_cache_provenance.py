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


def test_two_bams_sharing_a_basename_do_not_share_a_cache(tmp_path, monkeypatch):
    """The cache is keyed on file identity, not on `os.path.basename(bam_file)`.

    sampleA/aligned.bam and sampleB/aligned.bam are the ordinary shape of a cohort,
    and while the key carried only the basename the second run read the first one's
    extraction. One alignment each, at different positions, so a collision returns
    the wrong read rather than merely the wrong count.
    """
    cache_dir = tmp_path / "cache"
    first = tmp_path / "sampleA" / "aligned.bam"
    second = tmp_path / "sampleB" / "aligned.bam"
    first.parent.mkdir()
    second.parent.mkdir()
    _write_bam(first, [_alignment("from_A", 0, 100, [(0, 50)], "A" * 50, 0)])
    _write_bam(second, [_alignment("from_B", 0, 400, [(0, 50)], "A" * 50, 0)])

    _, first_alignments = _retrieve(first, cache_dir)
    _, second_alignments = _retrieve(second, cache_dir)

    assert [a.get_read_name() for a in first_alignments] == ["from_A"]
    assert [a.get_read_name() for a in second_alignments] == ["from_B"]


def test_a_cache_from_an_earlier_extraction_version_is_not_reused(tmp_path, monkeypatch):
    """A code change to what extraction returns leaves every other key field alone.

    Without a method version in the token, a pickle written before such a change is
    a valid hit, so the change is silently undone on any resumed or --no_cleanup run
    -- which is exactly how the region off-by-one would have survived its own fix.
    Here the cache is written, then the version is bumped, and the stale entry must
    be ignored rather than served.
    """
    import Pretty_alignment_manager as pam

    cache_dir = tmp_path / "cache"
    bam_path = tmp_path / "reads.bam"
    _write_bam(bam_path, [_alignment("only", 0, 100, [(0, 50)], "A" * 50, 0)])

    _retrieve(bam_path, cache_dir)
    written = sorted(p.name for p in cache_dir.iterdir())
    assert any(f".v{pam.EXTRACTION_METHOD_VERSION}." in name for name in written)

    monkeypatch.setattr(
        pam, "EXTRACTION_METHOD_VERSION", pam.EXTRACTION_METHOD_VERSION + 1
    )
    _, alignments = _retrieve(bam_path, cache_dir)

    assert [a.get_read_name() for a in alignments] == ["only"]
    after = sorted(p.name for p in cache_dir.iterdir())
    assert len(after) > len(written), "bumping the version must not reuse the old stem"


def _one_alignment_bam(tmp_path):
    bam_path = tmp_path / "reads.bam"
    _write_bam(bam_path, [_alignment("only", 0, 100, [(0, 50)], "A" * 50, 0)])
    return bam_path


def _cache_stems(cache_dir):
    return {p.name.split(".restrict-")[0] for p in cache_dir.iterdir() if p.suffix == ".pkl"}


def test_changing_min_per_id_does_not_reuse_the_cache(tmp_path, monkeypatch):
    """--HiFi raises min_per_id from 80 to 97, so it changes which alignments survive.

    Extraction reads it at Bam_alignment_extractor.py:164. While it was absent from the
    cache token, a run with --HiFi and one without shared a stem, and whichever ran
    second was served the first one's alignments.
    """
    bam_path = _one_alignment_bam(tmp_path)
    cache_dir = tmp_path / "cache"

    monkeypatch.setitem(LRAA_Globals.config, "min_per_id", 80)
    _retrieve(bam_path, cache_dir)
    lenient = _cache_stems(cache_dir)

    monkeypatch.setitem(LRAA_Globals.config, "min_per_id", 97)
    _retrieve(bam_path, cache_dir)

    assert _cache_stems(cache_dir) - lenient, "raising min_per_id must not hit the cache"


def test_changing_max_intron_length_does_not_reuse_the_cache(tmp_path, monkeypatch):
    """Extraction discards alignments carrying an intron over the cap, at :58-59.

    Sweeping the cap is the obvious way to ask what it costs, and while it was absent
    from the token every value after the first answered with the first one's pickle.
    """
    bam_path = _one_alignment_bam(tmp_path)
    cache_dir = tmp_path / "cache"

    monkeypatch.setitem(LRAA_Globals.config, "max_intron_length", 200000)
    _retrieve(bam_path, cache_dir)
    default_cap = _cache_stems(cache_dir)

    monkeypatch.setitem(LRAA_Globals.config, "max_intron_length", 10000)
    _retrieve(bam_path, cache_dir)

    assert _cache_stems(cache_dir) - default_cap, "changing the cap must not hit the cache"



def _retrieve_correcting(bam_path, cache_dir, splice_graph=None, contig_seq=None, **kwargs):
    manager = Pretty_alignment_manager(
        Splice_graph() if splice_graph is None else splice_graph,
        alignment_cache_dir=str(cache_dir),
    )
    alignments = manager.retrieve_pretty_alignments(
        "chr1", "+", contig_seq, str(bam_path), use_cache=True, **kwargs
    )
    return manager, alignments


def test_the_token_is_built_only_after_correction_is_finally_decided(tmp_path, monkeypatch):
    """Oversimplify force-disables correction, and the token must already know that.

    While the override ran after the token had been assembled, an oversimplify contig
    wrote UNCORRECTED alignments under a stem that said corrected: a later ordinary
    corrected run found a valid hit and was served content that had never been through
    the corrector. That entry is poisoned rather than stale -- wrong, not merely old,
    and nothing about it looks wrong.

    Asserted on stem equality rather than on any particular field, so the invariant
    survives a change of token format: an oversimplify run must land where a run with
    correction switched off lands, and away from a run that really does correct.
    """
    bam_path = _one_alignment_bam(tmp_path)

    def stems_written(cache_name, oversimplify, try_correct_alignments):
        monkeypatch.setitem(LRAA_Globals.config, "oversimplify_enabled", oversimplify)
        monkeypatch.setitem(
            LRAA_Globals.config, "oversimplify_contigs", ["chr1"] if oversimplify else []
        )
        cache_dir = tmp_path / cache_name
        _retrieve_correcting(
            bam_path, cache_dir, try_correct_alignments=try_correct_alignments
        )
        return _cache_stems(cache_dir)

    oversimplified = stems_written("oversimplified", True, True)
    uncorrected = stems_written("uncorrected", False, False)
    corrected = stems_written("corrected", False, True)

    assert oversimplified, "the oversimplify run must have written a cache to compare"
    assert oversimplified == uncorrected, (
        "oversimplify stores uncorrected alignments, so it must share the stem of a run "
        f"that asked for none: {sorted(oversimplified)} vs {sorted(uncorrected)}"
    )
    assert oversimplified.isdisjoint(corrected), (
        "a genuinely corrected run must not read the oversimplify run's pickle: "
        f"{sorted(oversimplified)} vs {sorted(corrected)}"
    )