#!/usr/bin/env python3

"""The unique-read gate must count READS, not the support they stand for.

`min_unique_reads_novel_isoform` is the companion to `min_reads_novel_isoform`. Both
default to 2 and `--min_reads_novel` sets both, but they compare different quantities:

  min_reads_novel_isoform         Transcript.get_read_counts_assigned(), EM's fractional
                                  assignment total -- follows a weight
  min_unique_reads_novel_isoform  the sum of mp.get_read_count() over the multipaths whose
                                  fractional assignment clears unique_read_filter_min_frac
                                  (TranscriptFiltering.get_isoform_unique_assigned_read_count)
                                  -- a tally of literal reads

The second must stay a read tally. Coverage normalization retains one read in place of
several and records the ratio in XW, so if the tally summed get_read_weight() instead, a
single read kept at p = 1/15 would arrive carrying 15 and clear a threshold of 2 on its
own. A structure witnessed exactly once at a deep locus would then be reported as though
fifteen reads had independently attested to it -- which is the opposite of what a
"unique read" requirement is for. A confidence statement cannot be satisfied by
re-weighting a single observation.

Three things already prevent weights from reaching this gate in a real run: the pass that
feeds it reads `--bam` rather than the normalized bam, `--bam` is required to carry no XW
tag, and that pass is pinned `weight_reads=False`. This file tests the remaining link --
that the gate itself is counting reads -- because those three are properties of the
caller and this is a property of the gate. It is exercised by driving
`filter_isoforms_by_min_isoform_fraction`, the function the nested helper lives in, with a
registry that makes one read weigh far more than it counts.
"""

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

import LRAA_Globals
import TranscriptFiltering
from GenomeFeature import Exon
from MultiPath import MultiPath
from Quantify import Quantify
from Splice_graph import Splice_graph

CONTIG = "chr1"
# Two spliced structures sharing no node, so each read is unique to one isoform and the
# gate's fractional-assignment precondition is trivially met for every read.
EXON_COORDS = [(100, 150), (200, 250), (400, 450), (500, 550)]

# The dominant isoform's reads each weigh 1; the scarce one's single read weighs this.
# Well above the threshold, so a gate summing weights would wave it straight through.
SOLO_READ_WEIGHT = 15.0
N_DOMINANT_READS = 50


@pytest.fixture(autouse=True)
def _isolate_the_unique_gate():
    """Neutralize every OTHER filter, so a verdict can only come from the unique gate.

    `filtered_low_unique_or_isoform_fraction` is tested BEFORE the unique gate and would
    otherwise claim the scarce isoform first -- its unique-read fraction of the gene is
    tiny either way -- and the verdict would say nothing about which quantity the unique
    gate summed. The FSM retention exemption is disabled for the same reason.
    """
    keys = (
        "min_frac_gene_unique_reads",
        "min_FSM_reads_gate",
        "min_FSM_reads_retain_isoform",
        "min_unique_reads_novel_isoform",
        "unique_read_filter_min_frac",
        "aggressively_assign_reads",
        "show_progress_quant_assign",
        "fraction_read_align_overlap",
    )
    saved = {k: LRAA_Globals.config[k] for k in keys}
    LRAA_Globals.config["min_frac_gene_unique_reads"] = 0.0
    LRAA_Globals.config["min_FSM_reads_gate"] = 0
    LRAA_Globals.config["min_FSM_reads_retain_isoform"] = 0
    LRAA_Globals.config["min_unique_reads_novel_isoform"] = 2
    LRAA_Globals.config["aggressively_assign_reads"] = False
    LRAA_Globals.config["show_progress_quant_assign"] = False
    LRAA_Globals.config["fraction_read_align_overlap"] = 0.4
    LRAA_Globals.reset_read_weight_registry()
    try:
        yield
    finally:
        LRAA_Globals.config.update(saved)
        LRAA_Globals.reset_read_weight_registry()


def _graph():
    Exon.reset_counter()
    sg = Splice_graph()
    exons = [Exon(CONTIG, lend, rend, "+", 1) for lend, rend in EXON_COORDS]
    for exon in exons:
        sg._node_id_to_node[exon.get_id()] = exon
    sg._contig_acc = CONTIG
    sg._contig_strand = "+"
    return sg, [e.get_id() for e in exons]


def _transcript(sg, node_ids, exon_idxs, transcript_id):
    t = MultiPath(sg, [[node_ids[i] for i in exon_idxs]]).toTranscript()
    t.set_gene_id("g1")
    t.set_transcript_id(transcript_id)
    return t


def _assign(sg, node_ids, transcript, reads):
    """Attach one single-read multipath per read name directly to a transcript.

    _estimate_isoform_read_support reads multipaths off the transcript objects, so
    assigning them here is what a completed compatibility pass would have done and lets
    this test control exactly which read belongs to which isoform.
    """
    mps = []
    for read_name, exon_idxs in reads:
        mps.append(
            MultiPath(
                sg,
                [[node_ids[i] for i in exon_idxs]],
                read_types={"PacBio"},
                read_names={read_name},
                read_count=1,
            )
        )
    transcript.add_multipaths_evidence_assigned(mps)
    return mps


def _fixture():
    """A dominant isoform with 50 plain reads and a scarce one with a single heavy read."""
    sg, node_ids = _graph()

    dominant = _transcript(sg, node_ids, [0, 1], "t_dominant")
    scarce = _transcript(sg, node_ids, [2, 3], "t_scarce")

    _assign(
        sg, node_ids, dominant,
        [(f"dom{i}", [0, 1]) for i in range(N_DOMINANT_READS)],
    )

    solo_id = MultiPath._coerce_read_identifier(None, "solo")
    LRAA_Globals.register_read_weight(solo_id, SOLO_READ_WEIGHT)
    solo_mps = _assign(sg, node_ids, scarce, [("solo", [2, 3])])

    return sg, dominant, scarce, solo_mps[0]


def test_the_fixture_makes_one_read_weigh_far_more_than_it_counts():
    """Vacuity guard: without this gap the gate cannot be caught summing the wrong one."""
    _sg, _dominant, _scarce, solo_mp = _fixture()
    assert solo_mp.get_read_count() == 1
    assert solo_mp.get_read_weight() == pytest.approx(SOLO_READ_WEIGHT)
    assert SOLO_READ_WEIGHT >= 2 * LRAA_Globals.config["min_unique_reads_novel_isoform"], (
        "the solo read must weigh enough to clear the threshold on its own, or a gate "
        "summing weights would be filtered here too and the test would prove nothing"
    )


def test_a_single_heavy_read_does_not_satisfy_the_two_unique_read_requirement():
    """The gate under test, driven through the function its helper lives in.

    The scarce isoform holds exactly ONE read, which weighs 15. With
    min_unique_reads_novel_isoform at 2, a gate counting reads sees 1 and removes it; a
    gate summing weights sees 15 and keeps it. Switch
    get_isoform_unique_assigned_read_count from mp.get_read_count() to
    mp.get_read_weight() and this is the assertion that fails.
    """
    sg, dominant, scarce, _solo_mp = _fixture()

    retained = TranscriptFiltering.filter_isoforms_by_min_isoform_fraction(
        [dominant, scarce],
        min_isoform_fraction=0.0,
        run_EM=False,
        max_EM_iterations=1,
    )
    retained_ids = {t.get_transcript_id() for t in retained}

    assert "t_dominant" in retained_ids, (
        "the 50-read isoform must survive, or this fixture is filtering on something "
        "other than the gate under test"
    )
    assert "t_scarce" not in retained_ids, (
        "a single read weighing 15 satisfied a 2-unique-read requirement: the gate is "
        "summing get_read_weight() where it must sum get_read_count()"
    )


def test_two_literal_reads_do_satisfy_it():
    """The converse, so the test above is not passing because everything is filtered.

    Same scarce isoform, now witnessed by two actual reads of weight 1 each. It must
    survive, which is what makes the removal above attributable to the read TALLY rather
    than to the isoform being scarce, novel, or lightly weighted.
    """
    sg, node_ids = _graph()
    dominant = _transcript(sg, node_ids, [0, 1], "t_dominant")
    scarce = _transcript(sg, node_ids, [2, 3], "t_scarce")

    _assign(
        sg, node_ids, dominant,
        [(f"dom{i}", [0, 1]) for i in range(N_DOMINANT_READS)],
    )
    _assign(sg, node_ids, scarce, [("s0", [2, 3]), ("s1", [2, 3])])

    retained = TranscriptFiltering.filter_isoforms_by_min_isoform_fraction(
        [dominant, scarce],
        min_isoform_fraction=0.0,
        run_EM=False,
        max_EM_iterations=1,
    )
    retained_ids = {t.get_transcript_id() for t in retained}

    assert "t_scarce" in retained_ids, (
        "two unweighted reads must clear a two-unique-read requirement"
    )
