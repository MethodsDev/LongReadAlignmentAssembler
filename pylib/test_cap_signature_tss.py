#!/usr/bin/env python3

"""prune_uncapped_alt_TSS: the cap signature at alternative TSSs.

A read starting at a capped RNA end carries one more untemplated G than the library's
uncapped reads. Pretty_alignment records the 5' G run before the strip zeroes the clip;
Splice_graph calibrates the library's baseline from read starts away from TSS peaks and
prunes a TSS downstream of a capped site when its own reads are rarely capped.
"""

import sys
from collections import defaultdict
from pathlib import Path

import networkx as nx
import pysam
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA_Globals
from GenomeFeature import TSS, Exon, Intron
from Pretty_alignment import Pretty_alignment
from Splice_graph import Splice_graph


@pytest.fixture
def restore_config():
    keep = dict(LRAA_Globals.config)
    yield
    LRAA_Globals.config.clear()
    LRAA_Globals.config.update(keep)


def _alignment(cigar, seq, is_reverse=False):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100000}]}
    )
    a = pysam.AlignedSegment(header)
    a.query_name = "r"
    a.reference_id = 0
    a.reference_start = 1000
    a.mapping_quality = 60
    a.flag = 16 if is_reverse else 0
    a.cigarstring = cigar
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    return a


def _G_run(alignment):
    return Pretty_alignment(alignment, [[1001, 1100]]).get_five_prime_G_run()


# --------------------------------------------------------------------------
# the per-read G run
# --------------------------------------------------------------------------


def test_the_G_run_is_recorded_before_the_strip(restore_config):
    """The strip zeroes the clip; the G run must survive it."""
    LRAA_Globals.config["max_untemplated_G_at_TSS"] = 3
    assert _G_run(_alignment("100M", "A" * 100)) == 0
    assert _G_run(_alignment("1S100M", "G" + "A" * 100)) == 1
    assert _G_run(_alignment("2S100M", "GG" + "A" * 100)) == 2


def test_only_the_run_touching_the_alignment_counts(restore_config):
    """ONT: adapter + primer + G run; the run next to the alignment is the signature."""
    assert _G_run(_alignment("8S100M", "ACGTGGGG" + "A" * 100)) == 4
    assert _G_run(_alignment("3S100M", "GGA" + "A" * 100)) == 0


def test_reverse_reads_count_the_complement_at_the_right_end(restore_config):
    assert _G_run(_alignment("100M2S", "A" * 100 + "CC", is_reverse=True)) == 2
    assert _G_run(_alignment("100M2S", "A" * 100 + "GG", is_reverse=True)) == 0


# --------------------------------------------------------------------------
# calibration
# --------------------------------------------------------------------------


def _bare_graph():
    Splice_graph._TSS_cap_calibration.clear()
    sg = Splice_graph.__new__(Splice_graph)
    sg._contig_acc = "chr1"
    sg._contig_strand = "+"
    return sg


def test_the_baseline_is_the_common_run_away_from_peaks(restore_config):
    """Bulk PacBio-like: uncapped starts carry one G, capped two."""
    counter = defaultdict(lambda: defaultdict(float))
    counter[1000][(2, True)] = 70
    counter[1000][(1, True)] = 30
    for pos in range(5000, 5300):
        counter[pos][(1, True)] = 0.9
        counter[pos][(2, True)] = 0.1
    baseline, active, stats = _bare_graph()._calibrate_TSS_cap_signature(
        counter, [1000], 50
    )
    assert baseline == 1
    assert active
    assert stats["capped_peak"] == pytest.approx(0.7)


def test_no_signature_leaves_the_rule_off(restore_config):
    """Simulated or trimmed reads: peaks and body look alike."""
    counter = defaultdict(lambda: defaultdict(float))
    counter[1000][(0, True)] = 100
    for pos in range(5000, 5300):
        counter[pos][(0, True)] = 1
    baseline, active, _ = _bare_graph()._calibrate_TSS_cap_signature(
        counter, [1000], 50
    )
    assert baseline == 0
    assert not active


def test_too_few_calibration_reads_leave_the_rule_off(restore_config):
    counter = defaultdict(lambda: defaultdict(float))
    counter[1000][(1, True)] = 100
    counter[5000][(0, True)] = 10
    _, active, _ = _bare_graph()._calibrate_TSS_cap_signature(counter, [1000], 50)
    assert not active


def _strand_with_a_deep_uncapped_pile():
    """Twelve promoters, 40% capped, 20 reads each; one 5' truncation pile of 2000
    uncapped reads; body starts 5% capped (FTD chr10+ around VIM)."""
    counter = defaultdict(lambda: defaultdict(float))
    peaks = []
    for i in range(12):
        p = 10000 + 1000 * i
        counter[p][(1, True)] = 8
        counter[p][(0, True)] = 12
        peaks.append(p)
    counter[50000][(0, True)] = 2000
    peaks.append(50000)
    for pos in range(60000, 60400):
        counter[pos][(0, True)] = 0.95
        counter[pos][(1, True)] = 0.05
    return counter, peaks


def test_one_deep_uncapped_pile_does_not_switch_the_strand_off(restore_config):
    """Pooled, the pile makes the peaks 7% capped against 5% elsewhere; per peak, the
    upper quartile is 40%."""
    counter, peaks = _strand_with_a_deep_uncapped_pile()
    baseline, active, stats = _bare_graph()._calibrate_TSS_cap_signature(
        counter, peaks, 50
    )
    assert stats["capped_peak"] < 0.10
    assert stats["peaks_judged"] == 13
    assert (baseline, active) == (0, True)


def test_per_peak_contrast_stays_off_without_a_signature(restore_config):
    """Twenty peaks at the body rate, with binomial scatter: no signature."""
    counter = defaultdict(lambda: defaultdict(float))
    peaks = []
    for i in range(20):
        p = 10000 + 1000 * i
        counter[p][(1, True)] = i % 3  # 0-2 of 20 marked
        counter[p][(0, True)] = 20 - i % 3
        peaks.append(p)
    for pos in range(60000, 60400):
        counter[pos][(0, True)] = 0.95
        counter[pos][(1, True)] = 0.05
    _, active, stats = _bare_graph()._calibrate_TSS_cap_signature(counter, peaks, 50)
    assert stats["peaks_judged"] == 20
    assert not active


def test_TSS_cap_signature_on_needs_only_a_baseline(restore_config):
    """Cluster runs: too few reads to calibrate, the library is known to be capped."""
    LRAA_Globals.config["TSS_cap_signature"] = "on"
    counter = defaultdict(lambda: defaultdict(float))
    counter[1000][(1, True)] = 10
    for pos in range(5000, 5030):
        counter[pos][(0, True)] = 1
    baseline, active, _ = _bare_graph()._calibrate_TSS_cap_signature(
        counter, [1000], 50
    )
    assert (baseline, active) == (0, True)

    LRAA_Globals.config["TSS_cap_signature"] = "off"
    _, active, _ = _bare_graph()._calibrate_TSS_cap_signature(counter, [1000], 50)
    assert not active


# --------------------------------------------------------------------------
# pruning
# --------------------------------------------------------------------------


def _tss(pos, informative, marked):
    t = TSS("chr1", pos, pos, "+", informative)
    t.set_cap_counts(informative, marked)
    return t


def _graph_with(TSS_specs):
    """+ strand: exon1 -> intron -> exon2 -> exon3 (unbranched segment); a TSS is
    attached to the exon named in its spec."""
    sg = _bare_graph()
    g = nx.DiGraph()
    e1 = Exon("chr1", 100, 200, "+", 10)
    i1 = Intron("chr1", 201, 999, "+", 10)
    e2 = Exon("chr1", 1000, 1100, "+", 10)
    e3 = Exon("chr1", 1101, 1200, "+", 10)
    g.add_edge(e1, i1)
    g.add_edge(i1, e2)
    g.add_edge(e2, e3)
    exons = {"e1": e1, "e2": e2, "e3": e3}
    tss_objs = []
    for name, (exon, informative, marked) in TSS_specs.items():
        t = _tss(exons[exon].get_coords()[0], informative, marked)
        g.add_edge(t, exons[exon])
        tss_objs.append(t)
    sg._splice_graph = g
    sg._TSS_objs = list(tss_objs)
    return sg, dict(zip(TSS_specs, tss_objs))


def test_an_uncapped_site_downstream_of_a_capped_one_is_pruned(restore_config):
    """Past a splice junction too: truncation sites are not confined to the first exon."""
    LRAA_Globals.config["min_TSS_cap_informative_reads"] = 5
    LRAA_Globals.config["min_TSS_cap_frac"] = 0.10
    sg, t = _graph_with(
        {"U": ("e1", 100, 50), "D2": ("e2", 300, 3), "D3": ("e3", 50, 25)}
    )
    kept = sg._prune_uncapped_alt_TSS(list(t.values()))
    assert t["D2"] not in kept and t["D2"] not in sg._splice_graph
    assert t["U"] in kept and t["D3"] in kept


def test_the_upstream_most_site_is_never_judged(restore_config):
    """Uncapped spike-ins: no capped site upstream, nothing pruned."""
    sg, t = _graph_with({"U": ("e1", 100, 1), "D2": ("e2", 300, 3)})
    kept = sg._prune_uncapped_alt_TSS(list(t.values()))
    assert set(kept) == set(t.values())


def test_a_site_without_enough_informative_reads_is_not_judged(restore_config):
    LRAA_Globals.config["min_TSS_cap_informative_reads"] = 5
    sg, t = _graph_with({"U": ("e1", 100, 50), "D2": ("e2", 4, 0)})
    kept = sg._prune_uncapped_alt_TSS(list(t.values()))
    assert t["D2"] in kept


def test_the_defaults():
    """Pruning on everywhere (it switches itself off without a cap signature); the
    new-TSS gate only where the single-cell workflow turns it on."""
    assert LRAA_Globals.config["prune_uncapped_alt_TSS"] is True
    assert LRAA_Globals.config["require_cap_for_new_TSS"] is False


# --------------------------------------------------------------------------
# new TSSs in cluster-guided runs
# --------------------------------------------------------------------------


def test_a_new_TSS_needs_the_cap_signature(restore_config):
    """require_cap_for_new_TSS: sites the guide lacks need >= 2 capped reads and
    min_TSS_cap_frac; sites at a guide 5' end are left alone."""
    LRAA_Globals.config.update(
        {
            "require_cap_for_new_TSS": True,
            "min_TSS_cap_informative_reads": 5,
            "min_new_TSS_capped_reads": 2,
            "min_TSS_cap_frac": 0.10,
        }
    )
    sg = _bare_graph()
    sg._guide_TSS_positions = [1000]
    at_guide = _tss(1010, 50, 0)
    capped = _tss(3000, 10, 3)
    one_capped = _tss(5000, 5, 1)
    sparse = _tss(7000, 4, 4)
    sg._TSS_objs = [at_guide, capped, one_capped, sparse]
    sg._reject_uncapped_new_TSS(list(sg._TSS_objs), 50)
    assert sg._TSS_objs == [at_guide, capped]


def test_calibration_ignores_unspliced_reads_and_reuses_the_last_one(restore_config):
    """Unspliced truncation piles dilute the contrast; a build without spliced reads
    (the SE graph) takes the contig strand's earlier calibration."""
    sg = _bare_graph()
    counter = defaultdict(lambda: defaultdict(float))
    counter[1000][(1, True)] = 60
    counter[1000][(0, True)] = 40
    counter[1000][(0, False)] = 5000  # an unspliced pile at the peak
    for pos in range(5000, 5300):
        counter[pos][(0, True)] = 0.95
        counter[pos][(1, True)] = 0.05
    baseline, active, _ = sg._calibrate_TSS_cap_signature(counter, [1000], 50)
    assert (baseline, active) == (0, True)

    unspliced_only = defaultdict(lambda: defaultdict(float))
    unspliced_only[1000][(0, False)] = 100
    baseline, active, stats = sg._calibrate_TSS_cap_signature(
        unspliced_only, [1000], 50
    )
    assert (baseline, active) == (0, True) and "reused" in stats
