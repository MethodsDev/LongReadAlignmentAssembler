#!/usr/bin/env python3

"""The oversimplify paths stream their reads rather than holding the contig's library.

chrM is oversimplified by default, and both paths read the full unthinned --bam: PBMC
Kinnex chrM is 5.5M records on one strand. They took those reads from
Pretty_alignment_manager.retrieve_pretty_alignments, which returns a list, so the whole
strand was resident at once although each read is handled on its own.

  1. PARITY. LRAA._oversimplify_alignments yields the reads that call returned, in its
     order and with its segments.
  2. STREAMING. Each path writes a read's tracking row before it asks for the next
     read, so nothing upstream has to hold more than one.
  3. Bam_alignment_extractor.get_read_alignments is iter_read_alignments collected.
"""

import importlib.util
from importlib.machinery import SourceFileLoader
from io import StringIO
from pathlib import Path

import pysam

import LRAA_Globals
from Bam_alignment_extractor import Bam_alignment_extractor
from Transcript import Transcript


def _load_lraa_module():
    lraa_path = Path(__file__).resolve().parents[1] / "LRAA"
    loader = SourceFileLoader("lraa_oversimplify_streams_test", str(lraa_path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def _write_bam(bam_path):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chrM", "LN": 5000}],
    }
    # (name, flag, start, cigar, mapq, NM)
    records = [
        ("unspliced", 0, 10, [(0, 80)], 60, 0),
        ("spliced", 0, 20, [(0, 50), (3, 100), (0, 50)], 60, 0),
        # intron 1,500 > max_intron_length below: discarded by the retention policy
        ("long_terminal", 0, 30, [(0, 40), (3, 1500), (0, 60), (3, 100), (0, 40)], 60, 0),
        ("secondary", 256, 40, [(0, 80)], 60, 0),
        ("low_identity", 0, 50, [(0, 80)], 60, 40),
        ("softclipped", 0, 60, [(4, 20), (0, 60)], 60, 0),
    ]
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam_writer:
        for name, flag, start, cigar, mapq, nm in records:
            read_len = sum(n for op, n in cigar if op in (0, 4))
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = "A" * read_len
            read.flag = flag
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = mapq
            read.cigartuples = cigar
            read.query_qualities = pysam.qualitystring_to_array("I" * read_len)
            read.set_tag("NM", nm)
            bam_writer.write(read)
    pysam.index(str(bam_path))
    return str(bam_path)


def _oversimplify_config(monkeypatch, tmp_path):
    monkeypatch.setenv("LRAA_TMP_DIR", str(tmp_path))
    monkeypatch.setitem(LRAA_Globals.config, "oversimplify_enabled", True)
    monkeypatch.setitem(LRAA_Globals.config, "oversimplify_contigs", ["chrM"])
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 1)
    monkeypatch.setitem(LRAA_Globals.config, "max_intron_length", 1000)
    monkeypatch.setitem(LRAA_Globals.config, "min_per_id", 80)


def _rendered(alignments):
    return [
        (pa.get_read_name(), [list(s) for s in pa.get_pretty_alignment_segments()])
        for pa in alignments
    ]


def test_streamed_alignments_match_retrieve_pretty_alignments(tmp_path, monkeypatch):
    lraa = _load_lraa_module()
    _oversimplify_config(monkeypatch, tmp_path)
    bam_path = _write_bam(tmp_path / "reads.bam")

    listed = lraa.Pretty_alignment_manager(
        splice_graph=None, alignment_cache_dir=str(tmp_path / "cache")
    ).retrieve_pretty_alignments(
        "chrM",
        "+",
        "A" * 5000,
        bam_path,
        region_lend=None,
        region_rend=None,
        use_cache=False,
        restrict_splice_type=None,
        try_correct_alignments=False,
        SE_read_encapsulation_mask=None,
        per_id_QC_raise_error=False,
    )
    streamed = list(lraa._oversimplify_alignments("chrM", "+", bam_path, None, None))

    assert _rendered(streamed) == _rendered(listed)
    # The precondition that makes the comparison mean something: both the kept and
    # the discarded sets are non-trivial.
    names = [name for name, _ in _rendered(streamed)]
    assert names == ["unspliced", "spliced", "softclipped"]


def _streaming_source(tracking):
    """Yields three reads, recording how many tracking rows existed at each request."""

    rows_seen_at_request = []

    class _Alignment:
        def __init__(self, name):
            self._name = name

        def get_pretty_alignment_segments(self):
            return [[1, 100]]

        def get_read_name(self):
            return self._name

    def _alignments(*args, **kwargs):
        for name in ("r1", "r2", "r3"):
            rows_seen_at_request.append(len(tracking.getvalue().splitlines()))
            yield _Alignment(name)

    return _alignments, rows_seen_at_request


def test_best_overlap_writes_each_read_before_reading_the_next(monkeypatch):
    lraa = _load_lraa_module()
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 1)
    tracking = StringIO()
    source, rows_seen_at_request = _streaming_source(tracking)
    monkeypatch.setattr(lraa, "_oversimplify_alignments", source)

    transcript = Transcript("chrM", [[1, 100]], "+")
    transcript.set_gene_id("g1")
    transcript.set_transcript_id("t1")
    lraa._run_oversimplify_best_overlap(
        "chrM", "+", "A" * 200, "unused.bam", None, None, [transcript], StringIO(), tracking
    )

    assert rows_seen_at_request == [0, 1, 2]
    assert transcript.get_assigned_read_count() == 3.0


def test_aggregate_writes_each_read_before_reading_the_next(monkeypatch):
    lraa = _load_lraa_module()
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 1)
    tracking = StringIO()
    quant = StringIO()
    source, rows_seen_at_request = _streaming_source(tracking)
    monkeypatch.setattr(lraa, "_oversimplify_alignments", source)

    lraa._run_oversimplify_aggregate(
        "chrM", "+", "A" * 200, "unused.bam", None, None, StringIO(), quant, tracking
    )

    assert rows_seen_at_request == [0, 1, 2]
    # all_reads is counted in the same pass that wrote the rows
    assert quant.getvalue().rstrip().split("\t")[2:4] == ["3", "3.0"]


def test_get_read_alignments_is_the_iterator_collected(tmp_path, monkeypatch):
    _oversimplify_config(monkeypatch, tmp_path)
    bam_path = _write_bam(tmp_path / "reads.bam")

    for pretty in (True, False):
        listed = Bam_alignment_extractor(bam_path).get_read_alignments(
            "chrM", "+", pretty=pretty
        )
        streamed = list(
            Bam_alignment_extractor(bam_path).iter_read_alignments(
                "chrM", "+", pretty=pretty
            )
        )
        name = (lambda x: x.get_read_name()) if pretty else (lambda x: x.query_name)
        assert [name(x) for x in streamed] == [name(x) for x in listed]

    extractor = Bam_alignment_extractor(bam_path)
    list(extractor.iter_read_alignments("chrM", "+", pretty=True))
    assert extractor.get_last_discarded_read_names_by_reason()["low_perID"] == {
        "low_identity"
    }
