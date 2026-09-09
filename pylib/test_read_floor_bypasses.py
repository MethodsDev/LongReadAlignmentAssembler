#!/usr/bin/env python3

"""The two places "every model carries a read" could still be false.

min_reads_retain_isoform is applied at the top of run_transcript_assembly, which is
not everywhere a model can be written:

  1. OVERSIMPLIFY. A ref-guided run on an oversimplified contig skips discovery
     entirely and writes the provided annotation straight to the output GTF
     (LRAA, the `oversimplify_enabled and not QUANT_ONLY` branch). chrM is
     oversimplified by default in the single-cell workflow, so this is not a corner
     case: without a floor there, every annotated chrM model is reported whatever the
     reads did.

  2. RE-EM. filter_isoforms_by_min_isoform_fraction reruns EM, so the counts the floor
     checked are not the counts it produced -- dropping a competitor redistributes its
     mass and a model that cleared one read on the old abundances can hold less than
     one on the new ones.

Both are asserted on read counts, not on TPM: see test_reference_retention_read_floor
for why those are not the same number.
"""

import importlib.util
from importlib.machinery import SourceFileLoader
from io import StringIO
from pathlib import Path

import pytest

import LRAA_Globals
import TranscriptFiltering
from Transcript import Transcript


def _load_lraa_module():
    lraa_path = Path(__file__).resolve().parents[1] / "LRAA"
    loader = SourceFileLoader("lraa_read_floor_bypass_test", str(lraa_path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


class _FakeAlignment:
    """One alignment over the first transcript's exon."""

    def __init__(self, read_name):
        self._read_name = read_name

    def get_pretty_alignment_segments(self):
        return [[1, 100]]

    def get_read_name(self):
        return self._read_name


def _manager_yielding(*read_names):
    """A Pretty_alignment_manager stand-in serving one alignment per name given.

    Parametrised by read count because the spool only engages above a floor of 1, so a
    model must be able to carry more than one read to be RETAINED on that path.
    """

    class _Manager:
        def __init__(self, splice_graph=None):
            pass

        def retrieve_pretty_alignments(self, *args, **kwargs):
            return [_FakeAlignment(name) for name in read_names]

    return _Manager


_FakePrettyAlignmentManager = _manager_yielding("cell-1^umi-1^read-1")


# -- 1. oversimplify ----------------------------------------------------------


@pytest.fixture
def oversimplify_models():
    supported = Transcript("chrM", [[1, 100]], "+")
    supported.set_gene_id("g_supported")
    supported.set_transcript_id("t_supported")

    # No read overlaps this one, so best-overlap assigns it nothing.
    starved = Transcript("chrM", [[400, 500]], "+")
    starved.set_gene_id("g_starved")
    starved.set_transcript_id("t_starved")

    return supported, starved


def _run_oversimplify(lraa, models, min_reads=0, quant=None, track=None):
    return lraa._run_oversimplify_best_overlap(
        "chrM", "+", "A" * 600, "unused.bam", None, None,
        list(models),
        quant if quant is not None else StringIO(),
        track if track is not None else StringIO(),
        min_reads_retain_isoform=min_reads,
    )


def _quant_transcript_ids(quant_output):
    return {
        line.split("\t")[1]
        for line in quant_output.getvalue().splitlines()
        if line.strip()
    }


def _tracking_transcript_ids(tracking_output):
    return {
        line.split("\t")[1]
        for line in tracking_output.getvalue().splitlines()
        if line.strip()
    }


def test_oversimplify_assigns_nothing_to_an_unread_model(
    monkeypatch, oversimplify_models
):
    """The precondition: best-overlap leaves the starved model on zero reads.

    Without this the assertions below would pass for the wrong reason.
    """
    lraa = _load_lraa_module()
    monkeypatch.setattr(lraa, "Pretty_alignment_manager", _FakePrettyAlignmentManager)
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 1)

    supported, starved = oversimplify_models
    _run_oversimplify(lraa, (supported, starved))

    assert supported.get_assigned_read_count() == 1.0
    assert starved.get_assigned_read_count() == 0.0


def test_the_floor_removes_a_zero_read_model_from_quant_and_the_returned_set(
    monkeypatch, oversimplify_models
):
    """chrM is oversimplified by default, so this is the default path for every
    mitochondrial model. A missing floor here reports them all regardless of support."""
    lraa = _load_lraa_module()
    monkeypatch.setattr(lraa, "Pretty_alignment_manager", _FakePrettyAlignmentManager)
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 1)

    supported, starved = oversimplify_models
    quant = StringIO()
    reported = _run_oversimplify(lraa, (supported, starved), min_reads=1.0, quant=quant)

    assert [t.get_transcript_id() for t in reported] == ["t_supported"]
    assert _quant_transcript_ids(quant) == {"t_supported"}


def test_the_three_outputs_name_the_same_models(monkeypatch, oversimplify_models):
    """The invariant a caller-side filter cannot give you.

    Filtering only the GTF emission left the dropped models in quant.expr, so the two
    outputs disagreed about which models exist. The GTF is now written from the return
    value and quant is written from the same list, so the sets are equal by
    construction -- this asserts they stay that way, tracking included.
    """
    lraa = _load_lraa_module()
    monkeypatch.setattr(lraa, "Pretty_alignment_manager", _FakePrettyAlignmentManager)
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 1)

    supported, starved = oversimplify_models
    quant, track = StringIO(), StringIO()
    reported = _run_oversimplify(
        lraa, (supported, starved), min_reads=1.0, quant=quant, track=track
    )

    gtf_ids = {t.get_transcript_id() for t in reported}
    assert gtf_ids == _quant_transcript_ids(quant) == {"t_supported"}
    assert _tracking_transcript_ids(track) == {"t_supported"}


def test_a_floor_above_one_withholds_the_tracking_of_a_dropped_model(monkeypatch):
    """The spool path: a model can be dropped while holding tracking rows.

    At a floor of 1 or less this cannot happen -- a dropped model has zero reads and so
    wrote nothing -- which is why that path still streams. At 2 a one-read model is
    dropped from the GTF and quant while having produced a row, and streaming would
    leave it in tracking alone. Exercises the disk spool, which the default floor never
    touches.
    """
    lraa = _load_lraa_module()
    monkeypatch.setattr(lraa, "Pretty_alignment_manager", _FakePrettyAlignmentManager)
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 1)

    one_read = Transcript("chrM", [[1, 100]], "+")
    one_read.set_gene_id("g_one")
    one_read.set_transcript_id("t_one_read")

    quant, track = StringIO(), StringIO()
    reported = _run_oversimplify(
        lraa, (one_read,), min_reads=2.0, quant=quant, track=track
    )

    assert one_read.get_assigned_read_count() == 1.0, "fixture must give it one read"
    assert reported == []
    assert _quant_transcript_ids(quant) == set()
    assert _tracking_transcript_ids(track) == set(), (
        "a model absent from the GTF and quant must not survive in tracking"
    )


def test_a_floor_above_one_replays_the_tracking_of_a_retained_model(monkeypatch):
    """The spool must REPLAY, not merely discard.

    Two reads and a floor of 2, so the model is retained on the path where its rows
    went to the spool rather than straight out. Without this a spool that dropped
    everything would satisfy the withholding test above and lose every tracking row on
    any run configured above the default floor.
    """
    lraa = _load_lraa_module()
    monkeypatch.setattr(
        lraa, "Pretty_alignment_manager", _manager_yielding("read-1", "read-2")
    )
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 2)

    kept = Transcript("chrM", [[1, 100]], "+")
    kept.set_gene_id("g_kept")
    kept.set_transcript_id("t_kept")

    quant, track = StringIO(), StringIO()
    reported = _run_oversimplify(lraa, (kept,), min_reads=2.0, quant=quant, track=track)

    assert kept.get_assigned_read_count() == 2.0, "fixture must give it two reads"
    assert [t.get_transcript_id() for t in reported] == ["t_kept"]
    assert _quant_transcript_ids(quant) == {"t_kept"}

    tracking_rows = [l for l in track.getvalue().splitlines() if l.strip()]
    assert len(tracking_rows) == 2, "both spooled rows must be replayed"
    assert {r.split("\t")[5] for r in tracking_rows} == {"read-1", "read-2"}
    assert {r.split("\t")[1] for r in tracking_rows} == {"t_kept"}


def test_quant_only_semantics_report_every_transcript_asked_about(
    monkeypatch, oversimplify_models
):
    """The default of 0 is quant-only's contract, and it must stay reachable.

    A quantification asked about N transcripts answers about N, including the ones that
    got nothing. That is a different question from discovery's "does this model exist".
    """
    lraa = _load_lraa_module()
    monkeypatch.setattr(lraa, "Pretty_alignment_manager", _FakePrettyAlignmentManager)
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 1)

    supported, starved = oversimplify_models
    quant = StringIO()
    reported = _run_oversimplify(lraa, (supported, starved), quant=quant)

    assert {t.get_transcript_id() for t in reported} == {"t_supported", "t_starved"}
    assert _quant_transcript_ids(quant) == {"t_supported", "t_starved"}


def test_oversimplify_counts_ignore_an_imported_tpm(monkeypatch):
    """The accumulator reads assigned counts, not the GTF's TPM attribute.

    _run_oversimplify_best_overlap runs on transcripts parsed straight from the input
    GTF, so this is the path where an imported TPM is actually present. Reading the
    running total through get_read_counts_assigned() would start the count at the TPM
    and hand every mitochondrial model a support figure it never earned.
    """
    lraa = _load_lraa_module()
    monkeypatch.setattr(lraa, "Pretty_alignment_manager", _FakePrettyAlignmentManager)
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 1)

    transcript = Transcript("chrM", [[1, 100]], "+")
    transcript.set_gene_id("g1")
    transcript.set_transcript_id("t1")
    transcript._imported_TPM_val = 5000.0

    _run_oversimplify(lraa, (transcript,))

    assert transcript.get_assigned_read_count() == 1.0


# -- 2. re-EM -----------------------------------------------------------------


def _model(read_counts, transcript_id):
    transcript = Transcript("chr1", [[100, 150], [200, 250]], "+")
    transcript.set_transcript_id(transcript_id)
    transcript.set_read_counts_assigned(read_counts)
    return transcript


def test_a_survivor_demoted_by_re_em_is_removed_by_the_reapplied_floor():
    """The gap the single up-front application leaves.

    Both models clear the floor on the counts filtering starts from. The
    isoform-fraction pass then drops one and EM moves its mass, leaving the other
    below a read. Applying the floor only once would report it.
    """
    survivor = _model(1.5, "survivor")
    competitor = _model(40.0, "competitor")

    first_pass = TranscriptFiltering.filter_isoforms_by_min_assigned_reads(
        [survivor, competitor], 1.0
    )
    assert first_pass == [survivor, competitor], "both clear the floor initially"

    # What the EM inside filter_isoforms_by_min_isoform_fraction does: the competitor
    # is filtered, and re-estimating without it moves mass off the survivor.
    survivor.set_read_counts_assigned(0.6)
    after_isoform_fraction_em = [survivor]

    second_pass = TranscriptFiltering.filter_isoforms_by_min_assigned_reads(
        after_isoform_fraction_em, 1.0
    )
    assert second_pass == []


def test_re_em_promotion_is_not_undone_by_the_reapplied_floor():
    """The other direction: EM can also move mass ONTO a model.

    A model below the floor at the start is already gone, but one that gains mass must
    not be removed by the second application for having been low at the first.
    """
    gainer = _model(1.2, "gainer")

    assert TranscriptFiltering.filter_isoforms_by_min_assigned_reads([gainer], 1.0) == [
        gainer
    ]
    gainer.set_read_counts_assigned(9.0)
    assert TranscriptFiltering.filter_isoforms_by_min_assigned_reads([gainer], 1.0) == [
        gainer
    ]


def test_the_floor_is_wired_into_both_emitting_paths():
    """Source-level: where the floor is called, which behaviour tests cannot show.

    Two paths emit models. run_transcript_assembly must apply the floor before the
    filters AND again after the isoform-fraction EM, because that EM moves counts. The
    oversimplify branch must hand its floor to _run_oversimplify_best_overlap rather
    than filtering afterwards, because filtering afterwards is what let quant.expr and
    the GTF disagree.
    """
    lines = (Path(__file__).resolve().parents[1] / "LRAA").read_text().splitlines()

    floor_calls = [
        i for i, l in enumerate(lines) if "filter_isoforms_by_min_assigned_reads(" in l
    ]
    isoform_fraction_call = next(
        i for i, l in enumerate(lines) if "filter_isoforms_by_min_isoform_fraction(" in l
    )

    assert len(floor_calls) == 2, (
        "expected the floor before the filters and again after the isoform-fraction "
        "EM; found {}".format(len(floor_calls))
    )
    assert any(i < isoform_fraction_call for i in floor_calls), "no up-front floor"
    assert any(i > isoform_fraction_call for i in floor_calls), (
        "no floor after the isoform-fraction EM: a model demoted by that EM would be "
        "reported below the floor"
    )

    oversimplify_calls = [
        i
        for i, l in enumerate(lines)
        if "_run_oversimplify_best_overlap(" in l and not l.lstrip().startswith("def ")
    ]
    # the quant-only caller and the ref-guided caller
    assert len(oversimplify_calls) == 2, oversimplify_calls
    passes_floor = [
        i
        for i in oversimplify_calls
        if "min_reads_retain_isoform=" in "\n".join(lines[i : i + 20])
    ]
    assert len(passes_floor) == 1, (
        "exactly the ref-guided caller passes a floor; quant-only must keep its "
        "row-per-requested-transcript contract. found {}".format(passes_floor)
    )
