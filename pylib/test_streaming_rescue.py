#!/usr/bin/env python3

"""Streaming transcriptome rescue: which reads it is offered, and what it refuses to load.

Two things here are load bearing and neither is visible in a run's output.

The candidate branch.  Rescue only means anything if it is offered the same reads the
batch path collects, and the biggest of those categories -- reads the extractor discarded
for low percent identity, 10,119 of 14,455 on ONT chr20 -- is one the streaming loop drops
before anything can consider it.  So the discard test has to branch.  The danger in
branching there is that Util_funcs.quant_discard_reason is the SINGLE retention policy
shared by quantification, coverage normalization and the chunked pipeline: a branch that
changed what it effectively retains would move depth measurement and XW sampling weights
with it, silently.  These tests pin both halves -- that low_perID and only low_perID is
offered, and that the retention counters are byte-identical with the branch active.

The aligner assertions.  Every failure mode mappy has here is silent.  A failed index load
is falsy rather than an exception, an index built without sequences answers every query
with zero hits and no error, and a multi-part index exposes only its first part with no
warning.  Each of those would present as "rescue found nothing", which on ONT is also what
success looks like.  So each has to raise.
"""

import os
import shutil
import subprocess
import sys

import pysam
import pytest

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import IsoformReadRescue
import LRAA_Globals
import StreamingQuant
from LRAA_Globals import SPACER

mappy = pytest.importorskip("mappy")
import StreamingRescue  # noqa: E402  (needs mappy present)

CONTIG = "chr1"
CONTIG_LEN = 5000

LOW_PER_ID = IsoformReadRescue.RESCUE_CANDIDATE_LOW_PER_ID
SPACER_PATH = IsoformReadRescue.RESCUE_CANDIDATE_SPACER_PATH
NO_GRAPH_PATH = IsoformReadRescue.RESCUE_CANDIDATE_NO_GRAPH_PATH
UNASSIGNED = IsoformReadRescue.RESCUE_CANDIDATE_UNASSIGNED_TO_TARGETS

# Each read is placed at its own start coordinate, because _map_read_to_graph is handed
# alignment segments and never a read name -- the coordinate is the only handle a fake
# mapper has for deciding what to return.
START_ASSIGNED = 100
START_SPACER = 200
START_NO_PATH = 300
START_UNASSIGNABLE = 400
START_LOW_PER_ID = 500
START_DUPLICATE = 600
START_SUPPLEMENTARY = 700


class _FakeSpliceGraph:
    def canonical_simple_path(self, simple_path):
        return "canon:" + ",".join(simple_path)

    def get_contig_acc(self):
        return CONTIG

    def get_contig_strand(self):
        return "+"


class _FakeLRAA:
    """_map_read_to_graph keyed on where the read starts."""

    def __init__(self):
        self._splice_graph = _FakeSpliceGraph()

    def _map_read_to_graph(self, segments, **kwargs):
        start = segments[0][0]
        if start == START_ASSIGNED + 1:
            return ["E:1"]
        if start == START_SPACER + 1:
            return ["E:2", SPACER, "E:3"]
        if start == START_NO_PATH + 1:
            return None
        if start == START_UNASSIGNABLE + 1:
            return ["E:9"]
        if start == START_LOW_PER_ID + 1:
            return ["E:1"]
        return ["E:1"]


class _Census:
    """The rescuer interface, recording only. Same object the parity measurement uses."""

    def __init__(self):
        self.offered = StreamingRescue.RescueCandidateCensus()

    def offer(self, read, category):
        return self.offered.offer(read, category)

    @property
    def by_category(self):
        return {c: set(n) for c, n in self.offered.candidate_read_names.items()}


class _Rescues:
    """Accepts every read it is offered, onto one fixed path."""

    def __init__(self, path=("E:1",)):
        self._path = tuple(path)
        self.offered = []

    def offer(self, read, category):
        self.offered.append((read.query_name, category))
        return (self._path,)


def _write_bam(path):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"},
         "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}]}
    )
    specs = [
        # (name, start, NM, flags)
        ("assigned", START_ASSIGNED, 1, {}),
        ("spacer", START_SPACER, 1, {}),
        ("nopath", START_NO_PATH, 1, {}),
        ("unassignable", START_UNASSIGNABLE, 1, {}),
        # NM 30 over 100 aligned bases is 70% identity, under the default min_per_id 80
        ("lowperid", START_LOW_PER_ID, 30, {}),
        ("duplicate", START_DUPLICATE, 1, {"is_duplicate": True}),
        ("supplementary", START_SUPPLEMENTARY, 1, {"is_supplementary": True}),
    ]
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for name, start, nm, flags in specs:
            a = pysam.AlignedSegment(header)
            a.query_name = name
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigarstring = "100M"
            a.query_sequence = "ACGT" * 25
            a.query_qualities = pysam.qualitystring_to_array("I" * 100)
            a.set_tag("NM", nm, value_type="i")
            for attr, value in flags.items():
                setattr(a, attr, value)
            fh.write(a)
    pysam.index(str(path))
    return path


def _table():
    table = StreamingQuant.AssignmentTable()
    table.add("canon:E:1", [("g1", "t1", "h1", 2, "mp1", 1.0, 1.0, 0)])
    # Seen and matched nothing, which is a real answer the table must carry as one.
    table.add("canon:E:9", [])
    return table


def _run(bam, tmp_path, rescuer=None, resolver=None):
    tracking = tmp_path / "tracking.tsv"
    with open(tracking, "wt") as fh:
        totals = StreamingQuant.stream_assign(
            str(bam), CONTIG, "+", "A" * CONTIG_LEN, _FakeLRAA(), _table(), fh,
            resolver=resolver if resolver is not None else (lambda path: []),
            try_correct_alignments=False,
            rescuer=rescuer,
        )
    return totals, tracking.read_text().splitlines()


@pytest.fixture
def bam(tmp_path):
    return _write_bam(tmp_path / "reads.bam")


def test_only_low_per_id_is_rescued_among_discarded_reads(bam, tmp_path):
    """A read rejected ONLY on percent identity is a candidate; other reasons still discard.

    The duplicate and supplementary reads are the control. Both are discarded by the same
    predicate as the low-identity read, and offering them would mean the branch keyed on
    "was discarded" rather than on "was discarded for this one reason".
    """
    census = _Census()
    _run(bam, tmp_path, rescuer=census)

    assert census.by_category[LOW_PER_ID] == {"lowperid"}
    offered = set().union(*census.by_category.values())
    assert "duplicate" not in offered
    assert "supplementary" not in offered


def test_retention_is_unchanged_by_the_branch(bam, tmp_path):
    """The hard boundary: rescue may change who is offered rescue and nothing else.

    quant_discard_reason is shared with coverage normalization and the chunked pipeline, so
    a branch that quietly retained one more read would move depth measurement with it.
    Compared against a run with no rescuer at all, which is the pre-change behaviour.
    """
    baseline_totals, baseline_rows = _run(bam, tmp_path)
    census_totals, census_rows = _run(bam, tmp_path, rescuer=_Census())

    for counter in (
        "reads_streamed", "reads_filtered", "reads_no_path", "reads_spacer_path",
        "reads_unassignable", "reads_assigned", "paths_resolved_in_stream",
        "reads_on_stream_resolved_paths",
    ):
        assert getattr(census_totals, counter) == getattr(baseline_totals, counter), counter
    assert census_rows == baseline_rows
    assert dict(census_totals.frac_sum) == dict(baseline_totals.frac_sum)


def test_a_rescued_read_does_not_disturb_the_reads_that_stream(bam, tmp_path):
    """With rescue accepting everything, every non-candidate read is accounted for as before."""
    baseline_totals, baseline_rows = _run(bam, tmp_path)
    rescued_totals, rescued_rows = _run(bam, tmp_path, rescuer=_Rescues())

    for counter in ("reads_streamed", "reads_filtered", "reads_assigned"):
        assert getattr(rescued_totals, counter) == getattr(baseline_totals, counter), counter
    # Every baseline row still present and unchanged; rescue only adds.
    assert set(baseline_rows) <= set(rescued_rows)
    assert rescued_totals.rescue_rows_written == len(rescued_rows) - len(baseline_rows)


def test_spacer_and_no_path_reads_land_in_distinct_categories(bam, tmp_path):
    """Streaming's branch is exclusive where the batch path's is not.

    The batch path adds a spacer read at the spacer site AND again at the no-path site,
    because its paths_list_no_spacers comes out empty for it. Streaming records it once. So
    the comparable quantity across the two is the UNION -- measured 4,336 on ONT chr20,
    identical read for read -- and per-site counts are not.
    """
    census = _Census()
    _run(bam, tmp_path, rescuer=census)

    assert census.by_category[SPACER_PATH] == {"spacer"}
    assert census.by_category[NO_GRAPH_PATH] == {"nopath"}


def test_unassigned_to_targets_is_off_unless_its_own_flag_is_set(bam, tmp_path, monkeypatch):
    """The fourth category is an extension, so it must not ride along by default.

    Under --stream_reads the batch path derives this category from a first pass over the
    coverage-normalized bam while the stream reads the full one, and separately the batch
    path omits every read whose path anchors to no gene (Quantify.py, where top_genes is
    None continues without recording). Measured on ONT chr20, batch reports 3,442 where
    streaming finds 11,195 -- so enabling it targets reads the batch path does not.
    """
    census = _Census()
    _run(bam, tmp_path, rescuer=census)
    assert UNASSIGNED not in census.by_category

    monkeypatch.setitem(
        LRAA_Globals.config, "stream_reads_rescue_unassigned_to_targets", True
    )
    census_on = _Census()
    _run(bam, tmp_path, rescuer=census_on)
    assert census_on.by_category[UNASSIGNED] == {"unassignable"}


# --------------------------------------------------------------------------------------
# Aligner construction
# --------------------------------------------------------------------------------------


class _Transcript:
    def __init__(self, tid, gene_id, exons, strand="+"):
        self._tid, self._gene_id, self._exons, self._strand = tid, gene_id, exons, strand

    def has_introns(self):
        return len(self._exons) > 1

    def get_simple_path(self):
        return [f"E:{lend}-{rend}" for lend, rend in self._exons]

    def get_strand(self):
        return self._strand

    def get_gene_id(self):
        return self._gene_id

    def get_transcript_id(self):
        return self._tid

    def get_exon_segments(self):
        return list(self._exons)


class _Node:
    def __init__(self, coords):
        self._coords = coords

    def get_coords(self):
        return self._coords


class _ModelSpliceGraph(_FakeSpliceGraph):
    def get_node_obj_via_id(self, node_id):
        lend, rend = node_id[2:].split("-")
        return _Node((int(lend), int(rend)))


def _transcripts(n=2):
    # Distinct sequences per transcript, so the index is not degenerate.
    return [
        _Transcript(f"t{i}", f"g{i}", [(1 + 400 * i, 200 + 400 * i), (250 + 400 * i, 380 + 400 * i)])
        for i in range(n)
    ]


CONTIG_SEQ = "".join(
    "ACGTTGCAAGGCTTACGATCCGTTAAGCTCGGATCCAAGTTGCACGTAT"[(i * 7) % 49] for i in range(4000)
)


def _rescuer(tmp_path, transcripts=None, **kwargs):
    return StreamingRescue.StreamingRescuer(
        _ModelSpliceGraph(),
        _transcripts() if transcripts is None else transcripts,
        CONTIG_SEQ,
        None,
        tmp_dir=str(tmp_path),
        **kwargs,
    )


def test_a_working_index_reports_every_target(tmp_path):
    """The baseline the failure cases are measured against."""
    rescuer = _rescuer(tmp_path)
    try:
        assert rescuer.stats()["n_targets"] == 2
    finally:
        rescuer.close()


def test_no_transcript_models_declines_rather_than_raises(tmp_path):
    """A contig-strand with no multi-exonic target is ordinary, and the batch path returns
    early on the same condition. Rescue then accepts nothing, which is what having nothing
    to rescue against means."""
    rescuer = _rescuer(tmp_path, transcripts=[])
    try:
        assert rescuer.stats()["n_targets"] == 0
    finally:
        rescuer.close()


def test_an_index_that_fails_to_load_raises_rather_than_reading_as_empty(tmp_path, monkeypatch):
    """mappy.Aligner.__bool__ reports whether the index pointer is non-NULL, so a failed
    load raises nothing and behaves as an aligner that finds no hits."""
    real_aligner = mappy.Aligner

    def failing_aligner(*args, **kwargs):
        # A real failed load: mappy returns a falsy Aligner for an unreadable index.
        return real_aligner(fn_idx_in=str(tmp_path / "does-not-exist.mmi"))

    monkeypatch.setattr(mappy, "Aligner", failing_aligner)
    with pytest.raises(RuntimeError, match="could not load a transcript index"):
        _rescuer(tmp_path)


def test_a_truncated_index_raises_rather_than_offering_fewer_targets(tmp_path, monkeypatch):
    """mappy reads only the FIRST part of a multi-part index, with no warning, so a read
    would be offered a subset of the targets it is entitled to and simply fail to rescue.
    Faked by reporting a short n_seq, since provoking a real second part needs a reference
    larger than the default -I 8G batch size."""
    real_aligner = mappy.Aligner

    class Truncated:
        def __init__(self, *args, **kwargs):
            self._inner = real_aligner(*args, **kwargs)

        def __bool__(self):
            return bool(self._inner)

        def __getattr__(self, name):
            return getattr(self._inner, name)

        @property
        def n_seq(self):
            return self._inner.n_seq - 1

    monkeypatch.setattr(mappy, "Aligner", Truncated)
    with pytest.raises(RuntimeError, match="index is truncated"):
        _rescuer(tmp_path)


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="needs the minimap2 binary")
def test_a_sequence_less_index_raises_rather_than_finding_nothing(tmp_path, monkeypatch):
    """The quietest failure of the three, and the reason it is checked with a REAL index
    rather than a fake: mappy always requests base-level alignment, and map() returns
    immediately when the index holds no sequences. bool(aligner) is True, n_seq is right,
    and every read maps nowhere."""
    models = IsoformReadRescue._build_transcript_models(
        _ModelSpliceGraph(), _transcripts(), CONTIG_SEQ
    )
    fasta = tmp_path / "tx.fa"
    IsoformReadRescue._write_transcript_fasta(str(fasta), models)
    mmi = tmp_path / "noseq.mmi"
    subprocess.run(
        ["minimap2", "--idx-no-seq", "-d", str(mmi), str(fasta)],
        check=True, capture_output=True,
    )

    real_aligner = mappy.Aligner

    def from_seqless_index(*args, **kwargs):
        kwargs.pop("fn_idx_in", None)
        return real_aligner(fn_idx_in=str(mmi), **kwargs)

    monkeypatch.setattr(mappy, "Aligner", from_seqless_index)
    with pytest.raises(RuntimeError, match="retains no sequences"):
        _rescuer(tmp_path)


def test_hits_are_labelled_as_minimap2_would_label_them(tmp_path):
    """minimap2 calls the first hit whose id equals its parent the SAM primary and every
    later such hit supplementary (mm_set_sam_pri); mappy exposes only is_primary, so that
    split is rebuilt from hit order. It matters because batch acceptance drops
    supplementary hits and keeps secondary ones -- collapsing the two would admit hits the
    batch path rejects.
    """
    rescuer = _rescuer(tmp_path)
    try:
        models = rescuer._transcript_models
        target = next(iter(models.values()))
        records = rescuer._pysam_records_for("q1", target["sequence"])
        assert records, "a transcript must align to its own index"

        first = records[0]
        assert not first.is_secondary and not first.is_supplementary
        assert first.reference_name in models
        assert first.has_tag("NM")
        # No SEQ is carried, so read length has to come off the CIGAR.
        assert first.infer_read_length() == len(target["sequence"])
        # Soft clips only; a reference skip would mean the preset spliced the cDNA.
        assert all(op != 3 for op, _ in first.cigartuples)

        primaries = [r for r in records if not (r.is_secondary or r.is_supplementary)]
        assert len(primaries) == 1, "exactly one SAM primary, as mm_set_sam_pri gives"
    finally:
        rescuer.close()


def test_reverse_strand_hits_swap_their_soft_clips(tmp_path):
    """q_st/q_en stay on the original query strand (PAF convention) while a SAM record
    stores the reverse complement, so the two clip lengths exchange. Getting this backwards
    shifts the projected genomic interval by the clip difference and stays plausible."""
    rescuer = _rescuer(tmp_path)
    try:
        target = next(iter(rescuer._transcript_models.values()))
        forward = target["sequence"]
        # A prefix, so the forward hit is clipped on one side only and the swap is visible.
        query = mappy.revcomp(forward[: len(forward) - 40])
        records = rescuer._pysam_records_for("q_rev", query)
        reverse_hits = [r for r in records if r.is_reverse]
        if not reverse_hits:
            pytest.skip("no reverse hit produced for this fixture")
        rec = reverse_hits[0]
        left = rec.cigartuples[0][1] if rec.cigartuples[0][0] == 4 else 0
        right = rec.cigartuples[-1][1] if rec.cigartuples[-1][0] == 4 else 0
        assert left + right + rec.query_alignment_length == len(query)
        assert rec.infer_read_length() == len(query)
    finally:
        rescuer.close()


def test_candidate_names_are_not_recorded_unless_asked(tmp_path):
    """At a billion reads a candidate name set is exactly the unbounded per-read state the
    streaming pass exists not to hold, so the parity measurement's recording must be opt in."""
    rescuer = _rescuer(tmp_path)
    off_default = rescuer._record_candidate_names
    rescuer.close()
    recording = _rescuer(tmp_path, record_candidate_names=True)
    try:
        assert off_default is False
        assert recording._record_candidate_names is True
    finally:
        recording.close()


# --------------------------------------------------------------------------------------
# Locality
# --------------------------------------------------------------------------------------


def _distinct_contig_seq(length=1200, seed=20260817):
    """Non-repetitive sequence, so two transcripts drawn from it are distinguishable.

    CONTIG_SEQ above is deliberately simple and is in fact periodic, which is fine where
    only "does anything align" matters. It is useless here: under a short repeat every
    target matches every read and locality cannot be told from chance.
    """
    out = []
    state = seed
    for _ in range(length):
        state = (state * 1103515245 + 12345) & 0x7FFFFFFF
        out.append("ACGT"[(state >> 16) & 3])
    return "".join(out)


LOCAL_SEQ = _distinct_contig_seq()

# Two transcripts at disjoint loci. A read placed on NEAR must not be credited to FAR.
NEAR_EXONS = [(1, 200), (251, 400)]
FAR_EXONS = [(701, 900), (951, 1100)]


def _locality_rescuer(tmp_path):
    transcripts = [
        _Transcript("near", "gNear", NEAR_EXONS),
        _Transcript("far", "gFar", FAR_EXONS),
    ]
    return StreamingRescue.StreamingRescuer(
        _ModelSpliceGraph(), transcripts, LOCAL_SEQ, None, tmp_dir=str(tmp_path)
    )


def _genome_read_over_near(sequence):
    """A primary genome record whose aligned blocks sit inside NEAR's first exon."""
    header = pysam.AlignmentHeader.from_references([CONTIG], [len(LOCAL_SEQ)])
    read = pysam.AlignedSegment(header)
    read.query_name = "displaced"
    read.reference_id = 0
    read.reference_start = 0
    read.mapping_quality = 60
    read.cigarstring = f"{len(sequence)}M"
    read.query_sequence = sequence
    read.set_tag("NM", 0, value_type="i")
    return read


def test_the_index_sees_only_the_targets_a_read_overlaps(tmp_path):
    """The input streaming supplies to the shared locality rule, checked on its own.

    Everything downstream is the batch path's; what streaming adds is this per-read set,
    computed from the read's own genome blocks against the same exon overlap index.
    """
    rescuer = _locality_rescuer(tmp_path)
    try:
        read = _genome_read_over_near(LOCAL_SEQ[:150])
        overlapping = set(
            IsoformReadRescue._get_alignment_overlapping_targets(
                read, rescuer._exon_overlap_index
            )
        )
        assert overlapping == {"gNear^near"}, overlapping
    finally:
        rescuer.close()


def test_a_hit_to_a_target_the_read_does_not_overlap_is_declined(tmp_path):
    """The rejection itself, end to end through the shared acceptance function.

    The read's genome alignment sits on NEAR while its sequence is FAR's, so its only
    good transcriptome hit is at a locus the read is not at -- the case locality exists
    to refuse. Worth a test of its own because on ONT chr20 the gate fires zero times:
    the aligned-fraction threshold leaves six reads whose hits are all local, so an
    end-to-end run there cannot tell this gate working from it being absent.
    """
    rescuer = _locality_rescuer(tmp_path)
    try:
        far_sequence = "".join(LOCAL_SEQ[lend - 1 : rend] for lend, rend in FAR_EXONS)
        read = _genome_read_over_near(far_sequence)
        rescued = rescuer.offer(read, LOW_PER_ID)

        assert rescued == (), "a non-local target must not be credited"
        assert rescuer.alignment_rejections["locality"] >= 1
        assert rescuer.reads_locality_declined == 1
        # Nothing local survived either, so this read's placement really was displaced.
        assert rescuer.reads_locality_displaced == 1
        assert rescuer.stats()["alignment_rejections"]["locality"] >= 1
    finally:
        rescuer.close()


# --------------------------------------------------------------------------------------
# The refusals
# --------------------------------------------------------------------------------------

LRAA_SCRIPT = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "LRAA")


def _guard_inputs(tmp_path):
    gtf = tmp_path / "a.gtf"
    gtf.write_text(
        f'{CONTIG}\ttest\texon\t101\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    )
    genome = tmp_path / "g.fa"
    genome.write_text(f">{CONTIG}\n" + "A" * CONTIG_LEN + "\n")
    return _write_bam(tmp_path / "guard.bam"), gtf, genome


def _guard_run(tmp_path, *extra):
    bam, gtf, genome = _guard_inputs(tmp_path)
    cmd = [
        sys.executable, LRAA_SCRIPT, "--quant_only", "--bam", str(bam),
        "--gtf", str(gtf), "--genome", str(genome),
        "--output_prefix", str(tmp_path / "out"),
        # Explicit on both counts: these guards are about --stream_reads's
        # interaction with transcriptome rescue specifically, and --chunk now
        # also defaults on, which would otherwise dispatch every one of these
        # invocations into the chunked pipeline before it ever reaches the
        # guards below -- a confound that did not exist when this harness was
        # written.
        "--stream_reads", "--no_chunk",
    ] + list(extra)
    return subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))


def test_streaming_rescues_unassigned_reads_by_default(tmp_path):
    """The pre-existing refusal only fired when the flag was left unstated, and
    --stream_reads_rescue_unassigned's default no longer leaves it unstated: unset, it
    now resolves to whichever way transcriptome rescue itself resolved, on by default --
    so this exact combination (streaming, rescue on, the flag not given) runs rescue
    inside the stream instead of refusing. The refusal itself is still live; it now
    takes an explicit --no_stream_reads_rescue_unassigned to reach it (see
    test_streaming_rescue_needs_rescue_left_enabled and
    test_the_fourth_category_flag_is_refused_on_its_own for the guards that still fire).
    """
    r = _guard_run(tmp_path)
    combined = r.stdout + r.stderr
    assert r.returncode == 0, combined[-3000:]
    assert "cannot reproduce transcriptome rescue" not in combined
    assert list(tmp_path.glob("out*quant.expr")), "must emit results"


def test_the_retired_assignment_mode_flag_is_rejected_rather_than_ignored(tmp_path):
    """--quant_read_assignment_mode is gone. A script still passing `genome` was asking
    for rescue OFF, and rescue is on by default -- so silently accepting and ignoring the
    flag would turn rescue on for that run with nothing to say so. argparse must refuse
    the unknown argument instead."""
    r = _guard_run(tmp_path, "--quant_read_assignment_mode", "genome")
    assert r.returncode != 0
    out = r.stdout + r.stderr
    assert "unrecognized arguments" in out
    assert "--quant_read_assignment_mode" in out
    assert not list(tmp_path.glob("out*quant.expr")), "and must not emit results"


def test_the_fourth_category_flag_is_refused_on_its_own(tmp_path):
    """It extends streaming rescue and does nothing without it. Ignoring it would leave a
    run that asked for the extension, got the base behaviour, and said so nowhere."""
    r = _guard_run(
        tmp_path,
        "--no_rescue_unassigned_reads_via_transcriptome_alignment",
        "--stream_reads_rescue_unassigned_to_targets",
    )
    assert r.returncode != 0
    assert "does nothing without it" in (r.stdout + r.stderr)


def test_streaming_rescue_needs_rescue_left_enabled(tmp_path):
    """Asking to rescue inside the stream while turning rescue off asks for opposite
    things. This used to be caught because --no_rescue... downgraded the assignment mode;
    with the modes gone the contradiction is named directly."""
    r = _guard_run(
        tmp_path,
        "--stream_reads_rescue_unassigned",
        "--no_rescue_unassigned_reads_via_transcriptome_alignment",
    )
    assert r.returncode != 0
    assert "transcriptome rescue left enabled" in (r.stdout + r.stderr)
