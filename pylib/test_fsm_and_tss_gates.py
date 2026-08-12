#!/usr/bin/env python3

"""Two evidence quantities the filters can be pointed at, and the boundaries of each.

`min_FSM_reads_gate` swaps an absolute count of reads matching an isoform's exact
intron chain for the unique-read fraction, which is relative to gene depth.
`max_untemplated_G_at_TSS` strips the non-genomic G run reverse transcriptase adds
opposite the cap, which the shipped `max_soft_clip_at_TSS` of 0 otherwise treats as
disqualifying -- rejecting the reads that evidence a real transcript start.

Both default off. These tests cover what turns on when they are not.
"""

import sys
from pathlib import Path

import pysam
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA_Globals
import TranscriptFiltering
from Pretty_alignment import Pretty_alignment


# --------------------------------------------------------------------------
# full-splice-match counting
# --------------------------------------------------------------------------


class _FakeMultiPath:
    def __init__(self, path, read_count):
        self._path = list(path)
        self._read_count = read_count

    def get_simple_path(self):
        return self._path

    def get_simple_path_tuple(self):
        return tuple(self._path)

    def get_read_count(self):
        return self._read_count


class _FakeTranscript:
    def __init__(self, transcript_id, path):
        self._transcript_id = transcript_id
        self._path = list(path)

    def get_transcript_id(self):
        return self._transcript_id

    def get_simple_path(self):
        return self._path


SPLICED = ["TSS:1", "E:1", "I:1", "E:2", "I:2", "E:3", "POLYA:1"]


def test_only_an_identical_intron_chain_counts():
    """A read must match the whole chain -- not a prefix, subset, or superset."""
    tx = _FakeTranscript("t1", SPLICED)
    assigned = {
        "t1": {
            _FakeMultiPath(["E:1", "I:1", "E:2", "I:2", "E:3"], 5): 1.0,  # exact
            _FakeMultiPath(["E:1", "I:1", "E:2"], 90): 1.0,  # prefix, ISM
            _FakeMultiPath(["E:2", "I:2", "E:3", "I:9", "E:4"], 90): 1.0,  # extends
            _FakeMultiPath(["E:1", "I:7", "E:3", "I:2", "E:9"], 90): 1.0,  # differs
        }
    }
    assert TranscriptFiltering.get_isoform_FSM_read_count(tx, assigned) == 5


def test_terminal_variation_does_not_break_the_match():
    """Differing TSS/PolyA nodes are not part of the splice pattern."""
    tx = _FakeTranscript("t1", SPLICED)
    assigned = {
        "t1": {
            _FakeMultiPath(["POLYA:9", "E:1", "I:1", "E:2", "I:2", "E:3", "TSS:9"], 4): 1.0
        }
    }
    assert TranscriptFiltering.get_isoform_FSM_read_count(tx, assigned) == 4


def test_each_read_counts_once_however_it_is_assigned():
    """The count is evidence for a structure, not a share of it.

    A multipath reachable from several carriers must not be counted per carrier,
    and the fraction EM assigned is deliberately ignored -- it is recomputed every
    filtering round, so gating on it would protect a model once and drop it next.
    """
    tx = _FakeTranscript("t1", SPLICED)
    mp = _FakeMultiPath(["E:1", "I:1", "E:2", "I:2", "E:3"], 7)
    assert TranscriptFiltering.get_isoform_FSM_read_count(tx, {"t1": {mp: 0.0}}) == 7


def test_a_monoexonic_model_has_no_splice_pattern():
    """Zero, not an empty-chain match against every other monoexon."""
    tx = _FakeTranscript("t1", ["TSS:1", "E:1", "POLYA:1"])
    mono = _FakeMultiPath(["E:1"], 50)
    assert TranscriptFiltering.get_splice_pattern(tx) == ()
    assert TranscriptFiltering.get_isoform_FSM_read_count(tx, {"t1": {mono: 1.0}}) == 0


def test_an_unassigned_isoform_counts_nothing():
    tx = _FakeTranscript("absent", SPLICED)
    assert TranscriptFiltering.get_isoform_FSM_read_count(tx, {}) == 0


# --------------------------------------------------------------------------
# untemplated G stripping at the TSS
# --------------------------------------------------------------------------


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


def _clips(alignment):
    pa = Pretty_alignment(alignment, [[1001, 1100]])
    return pa.left_soft_clipping, pa.right_soft_clipping


def test_g_run_is_kept_by_default(restore_config):
    """Off unless asked for: the shipped behaviour is unchanged."""
    LRAA_Globals.config["max_untemplated_G_at_TSS"] = 0
    assert _clips(_alignment("2S100M", "GG" + "A" * 100))[0] == 2


def test_a_leading_g_run_is_stripped_when_enabled(restore_config):
    LRAA_Globals.config["max_untemplated_G_at_TSS"] = 3
    assert _clips(_alignment("2S100M", "GG" + "A" * 100))[0] == 0
    assert _clips(_alignment("1S100M", "G" + "A" * 100))[0] == 0


def test_reverse_strand_looks_for_the_complement(restore_config):
    """A reverse alignment stores the reverse complement, so the 5' G run is a
    C run at the right end. Getting this backwards would silently do nothing."""
    LRAA_Globals.config["max_untemplated_G_at_TSS"] = 3
    assert _clips(_alignment("100M2S", "A" * 100 + "CC", is_reverse=True))[1] == 0
    # G's at the right end of a reverse read are not the cap signature
    assert _clips(_alignment("100M2S", "A" * 100 + "GG", is_reverse=True))[1] == 2


def test_the_wrong_end_is_left_alone(restore_config):
    """Only the TSS side is eligible; the polyA side has its own rule."""
    LRAA_Globals.config["max_untemplated_G_at_TSS"] = 3
    assert _clips(_alignment("100M2S", "A" * 100 + "GG"))[1] == 2


def test_a_mixed_clip_is_not_the_cap_signature(restore_config):
    """This is the specificity the rule buys over tolerating N clipped bases:
    GG is the library artifact, AT is a misalignment."""
    LRAA_Globals.config["max_untemplated_G_at_TSS"] = 3
    assert _clips(_alignment("2S100M", "AT" + "A" * 100))[0] == 2
    assert _clips(_alignment("2S100M", "AG" + "A" * 100))[0] == 2


def test_a_run_longer_than_allowed_is_not_stripped(restore_config):
    LRAA_Globals.config["max_untemplated_G_at_TSS"] = 2
    assert _clips(_alignment("3S100M", "GGG" + "A" * 100))[0] == 3
    LRAA_Globals.config["max_untemplated_G_at_TSS"] = 3
    assert _clips(_alignment("3S100M", "GGG" + "A" * 100))[0] == 0
