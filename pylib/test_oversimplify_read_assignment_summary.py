#!/usr/bin/env python3

"""An oversimplified work unit reports the reads it quantified, and counting them costs no memory.

The two oversimplify paths wrote no read-assignment summary, so the _contig_job_runner
backstop wrote one for them: reads counted in a second pass over the FIRST-PASS bam
(the coverage-normalized one under --stream_reads) with zero assigned, while the paths
quantify the full bam_file_for_quant. On a deep chrM that is a far smaller population
than the one in quant.expr. The backstop also collected the reads' names to count them,
about 0.9 GB on PBMC Kinnex chrM.

  1. Each path writes its own row: reads_total is the reads it streamed,
     reads_kept_genome the ones it credited to a model.
  2. Through the chunked route, where the unit's first-pass bam is a normalized slice,
     the merged total is the full library, not the thinned one.
  3. Bam_alignment_extractor.count_retained_reads counts what the name set counted:
     unpaired records one each, a proper pair once, low_perID only when asked.
"""

import csv
import importlib.util
import json
import os
from importlib.machinery import SourceFileLoader
from io import StringIO
from pathlib import Path

import pysam

import ChunkedRun
import LRAA_Globals
from Bam_alignment_extractor import Bam_alignment_extractor
from Transcript import Transcript


def _load_lraa_module():
    lraa_path = Path(__file__).resolve().parents[1] / "LRAA"
    loader = SourceFileLoader("lraa_oversimplify_summary_test", str(lraa_path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def _worker_row(path):
    with open(path) as fh:
        rows = [
            r for r in csv.DictReader(fh, delimiter="\t") if r["row_type"] == "worker"
        ]
    assert len(rows) == 1, rows
    return rows[0]


def _source(segments_by_name):
    class _Alignment:
        def __init__(self, name, segs):
            self._name = name
            self._segs = segs

        def get_pretty_alignment_segments(self):
            return self._segs

        def get_read_name(self):
            return self._name

    def _alignments(*args, **kwargs):
        for name, segs in segments_by_name:
            yield _Alignment(name, segs)

    return _alignments


def test_best_overlap_reports_reads_seen_and_assigned(tmp_path, monkeypatch):
    lraa = _load_lraa_module()
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 1)
    # r3 overlaps no model's exons: seen, not assigned.
    monkeypatch.setattr(
        lraa,
        "_oversimplify_alignments",
        _source([("r1", [[1, 100]]), ("r2", [[20, 90]]), ("r3", [[500, 600]])]),
    )
    transcript = Transcript("chrM", [[1, 100]], "+")
    transcript.set_gene_id("g1")
    transcript.set_transcript_id("t1")

    summary = tmp_path / "unit.read_assignment.summary.tsv"
    lraa._run_oversimplify_best_overlap(
        "chrM",
        "+",
        "A" * 1000,
        "unused.bam",
        None,
        None,
        [transcript],
        StringIO(),
        StringIO(),
        read_assignment_summary_path=str(summary),
    )

    row = _worker_row(summary)
    assert (int(row["reads_total"]), int(row["reads_kept_genome"])) == (3, 2)
    assert int(row["reads_rescue_requested"]) == 0


def test_aggregate_reports_every_read_assigned(tmp_path, monkeypatch):
    lraa = _load_lraa_module()
    monkeypatch.setitem(LRAA_Globals.config, "num_total_reads", 1)
    monkeypatch.setattr(
        lraa,
        "_oversimplify_alignments",
        _source([("r1", [[1, 100]]), ("r2", [[20, 90]]), ("r3", [[500, 600]])]),
    )

    summary = tmp_path / "unit.read_assignment.summary.tsv"
    lraa._run_oversimplify_aggregate(
        "chrM",
        "+",
        "A" * 1000,
        "unused.bam",
        None,
        None,
        StringIO(),
        StringIO(),
        StringIO(),
        read_assignment_summary_path=str(summary),
    )

    row = _worker_row(summary)
    assert (int(row["reads_total"]), int(row["reads_kept_genome"])) == (3, 3)


def _write_chrM(root, records):
    seq = ("ACGT" * 2000)[:8000]
    fasta = root / "g.fa"
    fasta.write_text(">chrM\n{}\n".format(seq))
    pysam.faidx(str(fasta))
    bam = root / "r.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chrM", "LN": 8000}],
    }
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for name, start, span, flag in sorted(records, key=lambda r: r[1]):
            aln = pysam.AlignedSegment()
            aln.query_name = name
            aln.flag = flag
            aln.reference_id = 0
            aln.reference_start = start
            aln.mapping_quality = 60
            aln.cigarstring = "{}M".format(span)
            aln.query_sequence = seq[start : start + span]
            aln.query_qualities = pysam.qualitystring_to_array("I" * span)
            out.write(aln)
    pysam.index(str(bam))
    return str(fasta), str(bam)


def test_chunked_oversimplified_total_is_the_full_library_not_the_normalized_one(
    tmp_path,
):
    """chrM stacked 60 deep against a normalization cap of 5.

    The unit's first-pass bam is the normalized slice, which the backstop counted
    over: well under 60. The paths quantify the full slice, and their row says so.
    """

    records = [("m_{}".format(i), 1000, 1000, 0) for i in range(60)]
    records += [("mr_{}".format(i), 4000, 800, 16) for i in range(20)]
    fasta, bam = _write_chrM(tmp_path, records)

    outdir = tmp_path / "work"
    ChunkedRun.run(
        ChunkedRun.default_args(
            bam=bam,
            genome_fa=fasta,
            output_dir=str(outdir),
            discovery=True,
            num_total_reads=len(records),
            no_reuse_source_bam=True,
            stop_after_make_chunks=True,
            cpu_budget=2,
            normalize_max_cov_level=5,
        )
    )
    with open(os.path.join(str(outdir), "chunk_plan.json")) as fh:
        plan = json.load(fh)
    units = []
    for entry in plan["chunks"]:
        ChunkedRun.run(
            ChunkedRun.default_args(
                output_dir=str(outdir),
                only_chunk=entry["chunk_id"],
                cpu_budget=1,
                # Stage 4 normalizes inside the per-chunk run, so the cap goes here.
                normalize_max_cov_level=5,
            )
        )
        with open(
            os.path.join(str(outdir), "chunks", entry["chunk_id"], "units.json")
        ) as fh:
            units.extend(json.load(fh)["units"])

    rows = {}
    for unit in units:
        row = _worker_row("{}.read_assignment.summary.tsv".format(unit["quant_prefix"]))
        rows[unit["strand"]] = (int(row["reads_total"]), int(row["reads_kept_genome"]))

    assert rows == {"+": (60, 60), "-": (20, 20)}


def _write_pair_bam(root):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chrM", "LN": 5000}],
    }
    PAIR = 0x1 | 0x2
    # (name, flag, start, NM)
    records = [
        ("u1", 0, 10, 0),
        ("u2", 0, 20, 0),
        ("u_low", 0, 30, 40),
        ("u_secondary", 256, 40, 0),
        ("p1", PAIR | 0x40, 50, 0),
        ("p1", PAIR | 0x80, 60, 0),
        # one mate kept, one low identity: one read, counted as kept
        ("p2", PAIR | 0x40, 70, 0),
        ("p2", PAIR | 0x80, 80, 40),
        ("p_low", PAIR | 0x40, 90, 40),
        ("p_low", PAIR | 0x80, 100, 40),
    ]
    with pysam.AlignmentFile(str(root / "pairs.bam"), "wb", header=header) as out:
        for name, flag, start, nm in records:
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = "A" * 80
            read.flag = flag
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigartuples = [(0, 80)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 80)
            read.set_tag("NM", nm)
            out.write(read)
    pysam.index(str(root / "pairs.bam"))
    return str(root / "pairs.bam")


def test_count_retained_reads_counts_what_the_name_set_counted(tmp_path, monkeypatch):
    monkeypatch.setitem(LRAA_Globals.config, "min_per_id", 80)
    bam = _write_pair_bam(tmp_path)

    # The distinct names the replaced method collected, from the reads the
    # extractor keeps and the ones it discards for low_perID.
    extractor = Bam_alignment_extractor(bam)
    kept_names = {r.query_name for r in extractor.iter_read_alignments("chrM", "+")}
    low_names = extractor.get_last_discarded_read_names_by_reason()["low_perID"]
    assert kept_names == {"u1", "u2", "p1", "p2"}

    counter = Bam_alignment_extractor(bam)
    assert counter.count_retained_reads("chrM", "+") == len(kept_names)
    assert counter.count_retained_reads("chrM", "+", include_low_perID=True) == len(
        kept_names | low_names
    )
