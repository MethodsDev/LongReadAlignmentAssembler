#!/usr/bin/env python3

"""The internal-priming veto waives itself ONLY where a trusted polyA reference calls a
cleavage site.

A read-derived PolyA candidate in A-rich genomic context is normally rejected before it
can become a graph vertex. As of v0.43.0 the only thing that waives that veto is an
explicit --polyA_known trusted cleavage set (Splice_graph._reference_three_prime_ends,
populated from LRAA_Globals.known_polyA_three_prime_ends). The veto is NO LONGER waived by
the --gtf structural annotation's transcript termini: a structural guide (GENCODE, or the
cluster-guided de-novo init GTF) is not cleavage-validated, and trusting its ends let
A-rich internal-priming termini within +/-25 nt be re-blessed as PolyA sites (the "halo").

"At that position" means within max_dist_between_alt_polyA_sites / 2, i.e. +/-25 nt, not an
exact coordinate match; see Splice_graph._reference_endorses_polyA_site.

Inert without --polyA_known, so the default vetoes every A-rich read candidate.
"""

import pytest

import LRAA_Globals
import Splice_graph
from Transcript import Transcript


def _guide(transcript_id, exons, strand):
    transcript = Transcript("chr1", exons, strand)
    transcript.set_gene_id("g")
    transcript.set_transcript_id(transcript_id)
    return transcript


def _polyA_sites(contig_seq_str, counter, strand="+", known=None, spare=True):
    # `known` mimics the --polyA_known set that populates _reference_three_prime_ends.
    sg = Splice_graph.Splice_graph()
    sg._contig_seq_str = contig_seq_str
    LRAA_Globals.config["spare_polyA_veto_at_known_3prime"] = spare
    sg._reference_three_prime_ends = sorted(set(known)) if (known and spare) else []
    sg._incorporate_PolyA_objects("chr1", strand, counter, from_reads=True)
    return sorted(obj.get_coords()[0] for obj in sg._PolyA_objs)


A_RICH_AT_1000 = "C" * 1000 + "A" * 20 + "C" * 1980


@pytest.fixture(autouse=True)
def _restore_config():
    keys = ("spare_polyA_veto_at_known_3prime", "polyA_known",
            "max_dist_between_alt_polyA_sites")
    saved = {k: LRAA_Globals.config.get(k) for k in keys}
    yield
    for k, v in saved.items():
        LRAA_Globals.config[k] = v
    LRAA_Globals._KNOWN_POLYA_ENDS_CACHE = None
    LRAA_Globals._KNOWN_POLYA_ENDS_CACHE_PATH = None


# --- the veto reprieve, now driven by the --polyA_known set --------------------------


def test_without_a_known_set_the_veto_still_fires():
    """Default (no --polyA_known): every A-rich read candidate is vetoed."""
    assert _polyA_sites(A_RICH_AT_1000, {1000: 40}, known=None) == []


def test_a_known_polyA_at_the_same_position_waives_the_veto():
    assert _polyA_sites(A_RICH_AT_1000, {1000: 40}, known=[1000]) == [1000]


def test_a_known_site_elsewhere_does_not_waive_it():
    assert _polyA_sites(A_RICH_AT_1000, {1000: 40}, known=[2500]) == []


def test_the_exemption_can_be_switched_off():
    assert _polyA_sites(A_RICH_AT_1000, {1000: 40}, known=[1000], spare=False) == []


def test_the_tolerance_is_the_alt_polyA_window_half():
    """+/-25 (max_dist_between_alt_polyA_sites / 2), inclusive, both sides."""
    tol = int(LRAA_Globals.config["max_dist_between_alt_polyA_sites"] / 2)
    for offset in (tol, -tol):
        assert _polyA_sites(A_RICH_AT_1000, {1000: 40}, known=[1000 + offset]) == [1000], offset
    for offset in (tol + 1, -(tol + 1)):
        assert _polyA_sites(A_RICH_AT_1000, {1000: 40}, known=[1000 + offset]) == [], offset


def test_the_tolerance_tracks_the_configured_window(monkeypatch):
    monkeypatch.setitem(LRAA_Globals.config, "max_dist_between_alt_polyA_sites", 4)
    assert _polyA_sites(A_RICH_AT_1000, {1000: 40}, known=[1002]) == [1000]
    assert _polyA_sites(A_RICH_AT_1000, {1000: 40}, known=[1003]) == []


def test_a_clean_context_candidate_is_unaffected_either_way():
    """The reprieve only ever loosens the veto, never tightens anything."""
    clean = "C" * 3000
    assert _polyA_sites(clean, {1000: 40}, known=None) == [1000]
    assert _polyA_sites(clean, {1000: 40}, known=[1000]) == [1000]


# --- the reprieve SOURCE is --polyA_known only; --gtf guides no longer endorse --------


def test_gtf_guides_do_not_endorse_without_polyA_known():
    """v0.43.0: a structural --gtf guide's transcript termini no longer waive the veto."""
    sg = Splice_graph.Splice_graph()
    sg._contig_acc = "chr1"
    LRAA_Globals.config["spare_polyA_veto_at_known_3prime"] = True
    LRAA_Globals.config["polyA_known"] = None
    LRAA_Globals._KNOWN_POLYA_ENDS_CACHE = None
    LRAA_Globals._KNOWN_POLYA_ENDS_CACHE_PATH = None
    assert sg._collect_reference_three_prime_ends([_guide("g", [[500, 1000]], "+")], "+") == []


def test_polyA_known_supplies_the_reprieve_set(tmp_path):
    bed = tmp_path / "known.bed"
    bed.write_text("chr1\t999\t1000\tk:1000:+\t.\t+\nchr1\t1999\t2000\tk:2000:-\t.\t-\n")
    LRAA_Globals.config["polyA_known"] = str(bed)
    LRAA_Globals._KNOWN_POLYA_ENDS_CACHE = None
    LRAA_Globals._KNOWN_POLYA_ENDS_CACHE_PATH = None
    sg = Splice_graph.Splice_graph()
    sg._contig_acc = "chr1"
    LRAA_Globals.config["spare_polyA_veto_at_known_3prime"] = True
    assert sg._collect_reference_three_prime_ends(None, "+") == [1000]
    assert sg._collect_reference_three_prime_ends(None, "-") == [2000]  # strand-specific
    # guides passed alongside are ignored; only the known file counts
    assert sg._collect_reference_three_prime_ends([_guide("g", [[100, 300]], "+")], "+") == [1000]


def test_exemption_off_returns_empty_even_with_known(tmp_path):
    bed = tmp_path / "known.bed"
    bed.write_text("chr1\t999\t1000\tk\t.\t+\n")
    LRAA_Globals.config["polyA_known"] = str(bed)
    LRAA_Globals._KNOWN_POLYA_ENDS_CACHE = None
    LRAA_Globals._KNOWN_POLYA_ENDS_CACHE_PATH = None
    sg = Splice_graph.Splice_graph()
    sg._contig_acc = "chr1"
    LRAA_Globals.config["spare_polyA_veto_at_known_3prime"] = False
    assert sg._collect_reference_three_prime_ends(None, "+") == []
