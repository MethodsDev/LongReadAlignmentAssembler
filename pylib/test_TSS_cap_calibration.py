#!/usr/bin/env python3

"""TSS_cap_calibration: the library-level cap-signature decision (pre-pass).

The pre-pass reads 5' clips itself instead of building Pretty_alignments, so its clip
rule must stay Pretty_alignment's; the decision must call two separated modes "on", one
mode "off", and too little data "undetermined"; a saved scan table must decide as the
scan did. testing/TSS_cap_calibration runs the standalone command end to end.
"""

import random
import sys
from pathlib import Path

import pysam
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA_Globals
import TSS_cap_calibration as T
from Pretty_alignment import Pretty_alignment


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


CLIPS = ["", "G", "GG", "GGG", "GGGG", "A", "GA", "AG", "ACGTGGG", "ACGTGG", "ACGTGGGG", "TTGGGC"]


@pytest.mark.parametrize("clip", CLIPS)
def test_the_clip_rule_is_Pretty_alignments(clip, restore_config):
    """Same G run, same verdict on whether the read can define a TSS, both strands."""
    body = "A" * 100
    for reverse in (False, True):
        if reverse:
            comp = clip.translate(str.maketrans("ACGT", "TGCA"))[::-1]
            seq = body + comp
            cigar = "100M{}S".format(len(clip)) if clip else "100M"
            clip_seq = comp
        else:
            seq = clip + body
            cigar = "{}S100M".format(len(clip)) if clip else "100M"
            clip_seq = clip
        pa = Pretty_alignment(_alignment(cigar, seq, reverse), [[1001, 1100]])
        residual_pa = pa.right_soft_clipping if reverse else pa.left_soft_clipping
        G_run, residual = T.five_prime_clip(clip_seq, not reverse)
        assert G_run == pa.get_five_prime_G_run(), (clip, reverse)
        assert (residual == 0) == (residual_pa == 0), (clip, reverse)


def _per_contig(peak_fracs, n_reads=40, off=None, seed=1):
    """A scan result with peaks at the given capped fractions (baseline 0)."""
    rnd = random.Random(seed)
    peaks = []
    for i, f in enumerate(peak_fracs):
        m = sum(rnd.random() < f for _ in range(n_reads))
        peaks.append((1000 * (i + 1), {0: n_reads - m, 1: m} if m else {0: n_reads}))
    return {"chr1": {"+": {"off_peak": off or {0: 950, 1: 50}, "peaks": peaks, "n_starts": None}}}


def test_two_separated_modes_are_on(restore_config):
    fracs = [0.02] * 300 + [0.5] * 300
    baseline, active, st = T.decide(_per_contig(fracs))
    assert (baseline, active, st["decision"]) == (0, True, "on")
    assert st["capped_mode_median"] > 0.4 and st["uncapped_mode_median"] < 0.1
    assert sum(st["histogram"]) == st["peaks_judged"] == 600


def test_one_mode_at_zero_is_off(restore_config):
    _, active, st = T.decide(_per_contig([0.01] * 600))
    assert (active, st["decision"]) == (False, "off")


def test_a_sparse_capped_tail_is_off(restore_config):
    """MORF-like: 2% of peaks capped-looking is not a library signature."""
    _, active, st = T.decide(_per_contig([0.0] * 2400 + [0.4] * 50))
    assert st["decision"] == "off"
    assert "of peaks" in st["reason"]


def test_too_few_peaks_is_undetermined(restore_config):
    _, active, st = T.decide(_per_contig([0.02] * 20 + [0.5] * 18))
    assert (active, st["decision"]) == (False, "undetermined")
    overrides, _ = T.overrides_from(0, active, st)
    assert overrides == {}


def test_a_thin_capped_mode_is_undetermined(restore_config):
    """A clear signature with fewer capped-mode peaks than the floor: thin data, not absence."""
    _, active, st = T.decide(_per_contig([0.02] * 100 + [0.5] * 30))
    assert st["decision"] == "undetermined"


def test_decisions_become_config_keys(restore_config):
    b, a, st = T.decide(_per_contig([0.02] * 300 + [0.5] * 300))
    assert T.overrides_from(b, a, st)[0] == {"TSS_cap_signature": "on", "TSS_cap_G_baseline": 0}
    b, a, st = T.decide(_per_contig([0.01] * 600))
    assert T.overrides_from(b, a, st)[0] == {"TSS_cap_signature": "off", "TSS_cap_G_baseline": 0}


def test_the_scan_table_decides_as_the_scan(tmp_path, restore_config):
    pc = _per_contig([0.02] * 300 + [0.5] * 300 + [0.3] * 5, n_reads=25)
    pc["chr1"]["+"]["peaks"].append((999999, {0: 3}))  # below the floor
    f = str(tmp_path / "x.scan.tsv.gz")
    T.write_scan_table(f, pc, min_peak_reads=20)
    back = T.read_scan_table(f)
    a, b = T.decide(pc)[2], T.decide(back)[2]
    for k in ("decision", "baseline_G_run", "peaks_judged", "histogram", "capped_mode_median"):
        assert a[k] == b[k], k
    assert all(sum(r.values()) >= 20 for _, r in back["chr1"]["+"]["peaks"])


def test_the_report_files(tmp_path, restore_config):
    b, a, st = T.decide(_per_contig([0.02] * 300 + [0.5] * 300))
    ov, st = T.overrides_from(b, a, st)
    out = T.write_report(str(tmp_path / "lib"), ov, st)
    assert set(out) == {"json", "tsv", "svg"}
    import xml.dom.minidom

    xml.dom.minidom.parse(out["svg"])
