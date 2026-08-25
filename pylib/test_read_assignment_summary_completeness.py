#!/usr/bin/env python3

"""The read-assignment summary must account for EVERY read of a run.

``<output_prefix>.read_assignment.summary.tsv`` is the run's statement of where
its reads went. ``run_quant_only`` has three early returns that fire before its
summary writer is reached -- no input transcripts for this contig/strand, an empty
splice graph, and no input transcript left after mapping into the graph -- and a
work unit taking any of them used to contribute NOTHING to the merged table.

That is not a rounding error. A whole-genome chunked run has ~600 quant units, and
the third case in particular can still HOLD reads: an intergenic chunk, or one
whose every annotated transcript was filtered out, has a bam full of reads and
nothing to assign them to. Those reads existed and were simply missing from
``reads_total``, with nothing in the file saying a subset was being reported.

The fix is COUNT-THEN-REPORT: each early return counts the reads it actually saw
and writes a summary reporting them with zero assigned. Not all-zeros -- a
fabricated zero ``reads_total`` for a shard that held reads would make the merged
total claim completeness while still undercounting, which is strictly worse than
the absent row it replaces.

The headline test here is the parity one: merged chunked ``reads_total`` must EQUAL
unchunked whole-contig ``reads_total`` on a fixture that deliberately contains a
region with reads and no quantifiable transcript. It fails the moment any early
return stops writing its summary. The rest pin the three returns individually, the
scope of the count, and stage 6's refusal of a unit that produced no summary at
all.
"""

import csv
import importlib.util
import json
import os
import random
import subprocess
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

import networkx as nx
import pysam
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "pylib"))

import ChunkedRun  # noqa: E402
import LRAA as LRAA_module  # noqa: E402
import LRAA_Globals  # noqa: E402
from GenomeFeature import Exon  # noqa: E402
from Splice_graph import Splice_graph  # noqa: E402
from Transcript import Transcript  # noqa: E402


def _load_lraa_driver():
    """The LRAA driver script, as a module.

    ``run_quant_only`` and its summary writers live in the driver, not in a
    package, and the three early returns are branches inside one function -- so
    they are exercised by CALLING it. Driving them through a whole-pipeline
    subprocess instead would reach each branch only by accident of fixture
    geometry, and could not distinguish "the branch wrote nothing" from "the
    branch was never taken".
    """

    loader = SourceFileLoader("lraa_driver_read_assignment_test", str(REPO / "LRAA"))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


lraa = _load_lraa_driver()


CONTIG = "chrT"
CONTIG_LEN = 24000

# Three loci, deliberately far apart so cut selection has read-free ground to cut
# in and no alignment straddles a chunk boundary. Only locus A is annotated:
#
#   A  1001-2500   10 spliced forward reads, the ONE annotated transcript
#   B  9001-10000   8 unspliced reverse reads, and no reverse annotation at all
#   C 17001-18000   6 unspliced forward reads, unannotated
#
# So a chunk holding C has reads and nothing to quantify -- the case the whole
# change exists for -- and the reverse strand takes the same early return
# everywhere. 24 reads total, well under normalize_max_cov_level (1000), so
# coverage normalization thins nothing and the read population the first
# assignment pass sees is the whole library on both routes.
LOCUS_A_READS = 10
LOCUS_B_READS = 8
LOCUS_C_READS = 6
TOTAL_READS = LOCUS_A_READS + LOCUS_B_READS + LOCUS_C_READS
FORWARD_READS = LOCUS_A_READS + LOCUS_C_READS
REVERSE_READS = LOCUS_B_READS


def _write_inputs(root):
    """A one-contig genome, an indexed bam over the three loci, and a GTF naming
    only locus A."""

    rng = random.Random(7)
    seq = [rng.choice("ACGT") for _ in range(CONTIG_LEN)]
    seq[1500:1502] = ["G", "T"]  # donor of locus A's intron [1501, 2000]
    seq[1998:2000] = ["A", "G"]  # acceptor
    seq = "".join(seq)

    fasta = root / "genome.fa"
    fasta.write_text(">{}\n{}\n".format(CONTIG, seq))
    pysam.faidx(str(fasta))

    records = []
    for i in range(LOCUS_A_READS):
        records.append(
            ("A_plus_{}".format(i), 1000, "500M500N500M",
             seq[1000:1500] + seq[2000:2500], 0)
        )
    for i in range(LOCUS_B_READS):
        records.append(
            ("B_minus_{}".format(i), 9000, "1000M", seq[9000:10000], 16)
        )
    for i in range(LOCUS_C_READS):
        records.append(
            ("C_plus_{}".format(i), 17000, "1000M", seq[17000:18000], 0)
        )
    records.sort(key=lambda r: r[1])

    bam = root / "reads.bam"
    _write_bam(bam, seq, records)

    # Only locus A. Nothing annotates locus C or the reverse strand, which is what
    # makes those work units take an early return while still holding reads.
    gtf = root / "annot.gtf"
    attrs = 'gene_id "gA"; transcript_id "tA";'
    with open(gtf, "wt") as fh:
        for feature, start, end in (
            ("transcript", 1001, 2500),
            ("exon", 1001, 1500),
            ("exon", 2001, 2500),
        ):
            print(
                "\t".join(
                    (CONTIG, "test", feature, str(start), str(end), ".", "+", ".",
                     attrs)
                ),
                file=fh,
            )

    return {"fasta": str(fasta), "bam": str(bam), "gtf": str(gtf), "seq": seq}


def _write_bam(path, seq, records):
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, start, cigar, query, flag in records:
            aln = pysam.AlignedSegment()
            aln.query_name = name
            aln.flag = flag
            aln.reference_id = 0
            aln.reference_start = start
            aln.mapping_quality = 60
            aln.cigarstring = cigar
            aln.query_sequence = query
            aln.query_qualities = pysam.qualitystring_to_array("I" * len(query))
            out.write(aln)
    pysam.index(str(path))


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    root = tmp_path_factory.mktemp("read_assignment_summary")
    doc = _write_inputs(root)
    doc["root"] = root
    return doc


def _summary_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _total_reads(path):
    totals = [r for r in _summary_rows(path) if r["row_type"] == "TOTAL"]
    assert len(totals) == 1, path
    return int(totals[0]["reads_total"])


def _run_quant_only(
    summary_path,
    tmpdir,
    contig_strand="+",
    input_transcripts=None,
    prebuilt_splice_graph=None,
    bam_for_sg=None,
    bam_for_quant=None,
    region=(None, None),
    contig_seq=None,
):
    """One ``run_quant_only`` call, with only what an early return reads.

    The quant and tracking handles are opened because the signature demands them;
    every branch under test returns before writing to either.
    """

    os.makedirs(tmpdir, exist_ok=True)
    with open(os.path.join(tmpdir, "quant.out"), "wt") as ofh_quant, open(
        os.path.join(tmpdir, "tracking.out"), "wt"
    ) as ofh_track:
        return lraa.run_quant_only(
            CONTIG,
            contig_strand,
            contig_seq,
            bam_for_sg,
            bam_for_quant,
            region[0],
            region[1],
            input_transcripts,
            ofh_quant,
            ofh_track,
            1,
            True,
            tmpdir,
            prebuilt_splice_graph=prebuilt_splice_graph,
            read_assignment_summary_path=summary_path,
        )


def _worker_row(path):
    rows = [r for r in _summary_rows(path) if r["row_type"] == "worker"]
    assert len(rows) == 1, rows
    return rows[0]


def _non_empty_splice_graph():
    """A splice graph with nodes, built directly.

    ``sg.is_empty()`` is ``len(self._splice_graph) == 0`` and nothing between it
    and the branch under test reads any other member, so the two exons here are
    the whole precondition: a graph that exists.
    """

    sg = Splice_graph()
    sg._contig_acc = CONTIG
    sg._contig_strand = "+"
    sg._splice_graph = nx.DiGraph(
        [
            (
                Exon(CONTIG, 1001, 1500, "+", 10),
                Exon(CONTIG, 2001, 2500, "+", 10),
            )
        ]
    )
    return sg


# -- the three early returns, individually


def test_no_input_transcripts_reports_the_reads_it_saw(inputs, tmp_path):
    """Early return 1: nothing to quantify on this contig/strand.

    The forward strand holds 16 reads. With no input transcripts the unit assigns
    none of them, and the honest report is 16 seen and 0 assigned -- not an absent
    file, and not a zero total.
    """

    summary = tmp_path / "case1.tsv"
    result = _run_quant_only(
        str(summary),
        str(tmp_path / "work"),
        input_transcripts=None,
        bam_for_sg=inputs["bam"],
        bam_for_quant=inputs["bam"],
        contig_seq=inputs["seq"],
    )

    assert result == (None, None)
    row = _worker_row(str(summary))
    assert int(row["reads_total"]) == FORWARD_READS
    assert int(row["reads_kept_genome"]) == 0
    assert int(row["reads_selected_tx_total"]) == 0


def test_an_empty_splice_graph_reports_the_reads_it_saw(inputs, tmp_path):
    """Early return 2: a graph was attempted and came out empty.

    Handed an empty prebuilt graph so the branch is reached deterministically
    rather than by finding inputs that happen to build nothing. The reads are
    still in the bam either way, which is the point.
    """

    empty = Splice_graph()
    empty._contig_acc = CONTIG
    empty._contig_strand = "+"
    empty._splice_graph = nx.DiGraph()

    summary = tmp_path / "case2.tsv"
    result = _run_quant_only(
        str(summary),
        str(tmp_path / "work"),
        input_transcripts=[Transcript(CONTIG, [[1001, 1500], [2001, 2500]], "+")],
        prebuilt_splice_graph=empty,
        bam_for_sg=inputs["bam"],
        bam_for_quant=inputs["bam"],
        contig_seq=inputs["seq"],
    )

    assert result == (None, None)
    row = _worker_row(str(summary))
    assert int(row["reads_total"]) == FORWARD_READS
    assert int(row["reads_kept_genome"]) == 0


def test_transcripts_filtered_out_of_the_graph_report_the_reads_they_saw(
    inputs, tmp_path, monkeypatch
):
    """Early return 3: the graph exists, and no input transcript survived mapping
    into it.

    The one that mattered most. This shard built a graph, so its bam held enough
    coverage to build one, and its annotation was still emptied -- an intergenic
    chunk, or one where every annotated model was filtered out. It holds reads and
    quantifies none of them.

    ``assign_transcripts_paths_in_graph`` is stubbed to return an empty list, with
    a real transcript going in. That models the case the branch is for -- a
    non-empty annotation emptied by mapping -- without depending on which
    coordinates a full graph build would reject; the branch's own precondition is
    the empty list it receives, not how the list came to be empty.
    """

    monkeypatch.setattr(
        LRAA_module.LRAA,
        "assign_transcripts_paths_in_graph",
        lambda self, transcripts: [],
    )

    summary = tmp_path / "case3.tsv"
    result = _run_quant_only(
        str(summary),
        str(tmp_path / "work"),
        input_transcripts=[Transcript(CONTIG, [[1001, 1500], [2001, 2500]], "+")],
        prebuilt_splice_graph=_non_empty_splice_graph(),
        bam_for_sg=inputs["bam"],
        bam_for_quant=inputs["bam"],
        contig_seq=inputs["seq"],
    )

    assert result == (None, None)
    row = _worker_row(str(summary))
    assert int(row["reads_total"]) == FORWARD_READS
    assert int(row["reads_kept_genome"]) == 0


# -- what the count is a count OF


@pytest.mark.parametrize(
    "contig_strand,region,expected",
    [
        ("+", (None, None), FORWARD_READS),
        ("-", (None, None), REVERSE_READS),
        ("+", (15001, 20000), LOCUS_C_READS),
        ("+", (1, 5000), LOCUS_A_READS),
    ],
    ids=["forward", "reverse", "locus_C_window", "locus_A_window"],
)
def test_the_count_is_scoped_to_this_shard(
    inputs, tmp_path, contig_strand, region, expected
):
    """reads_total is THIS unit's reads, not the bam's.

    A whole-bam count would be plausible and wrong in exactly the way the change
    exists to prevent: chunked units would then each report 24 and the merged total
    would be ten times the library. Both axes a unit is scoped on are exercised --
    orientation and window -- because a count that honoured only one of them would
    still pass a single-axis test.
    """

    summary = tmp_path / "scoped.tsv"
    _run_quant_only(
        str(summary),
        str(tmp_path / "work"),
        contig_strand=contig_strand,
        input_transcripts=None,
        bam_for_sg=inputs["bam"],
        bam_for_quant=inputs["bam"],
        region=region,
        contig_seq=inputs["seq"],
    )

    assert int(_worker_row(str(summary))["reads_total"]) == expected


def test_the_count_follows_the_bam_the_first_assignment_pass_reads(
    inputs, tmp_path, monkeypatch
):
    """Under --stream_reads the first pass reads the THINNED bam, so the count must
    too.

    ``reads_total`` on the assigning path is a count over that pass's multipaths,
    which come from ``bam_file_for_sg`` while streaming and ``bam_file_for_quant``
    otherwise. An early return that always counted one of them would disagree with
    the assigning units it is merged with -- and would do so silently, since both
    numbers look like read counts.
    """

    seq = inputs["seq"]
    thin = tmp_path / "thin.bam"
    _write_bam(
        thin,
        seq,
        [
            ("A_plus_{}".format(i), 1000, "500M500N500M",
             seq[1000:1500] + seq[2000:2500], 0)
            for i in range(4)
        ],
    )

    streaming = tmp_path / "streaming.tsv"
    _run_quant_only(
        str(streaming),
        str(tmp_path / "work_stream"),
        input_transcripts=None,
        bam_for_sg=str(thin),
        bam_for_quant=inputs["bam"],
        contig_seq=seq,
    )
    assert LRAA_Globals.config["stream_reads"] is True
    assert int(_worker_row(str(streaming))["reads_total"]) == 4

    monkeypatch.setitem(LRAA_Globals.config, "stream_reads", False)
    batched = tmp_path / "batched.tsv"
    _run_quant_only(
        str(batched),
        str(tmp_path / "work_batch"),
        input_transcripts=None,
        bam_for_sg=str(thin),
        bam_for_quant=inputs["bam"],
        contig_seq=seq,
    )
    assert int(_worker_row(str(batched))["reads_total"]) == FORWARD_READS


# -- the parity bar


@pytest.fixture(scope="module")
def unchunked_summary(inputs, tmp_path_factory):
    """One whole-contig ``LRAA --no_chunk`` run, as a subprocess.

    A subprocess because the driver writes its scratch into the working directory
    and merges its work units at exit; this is the same invocation a user makes,
    and its ``read_assignment.summary.tsv`` is the number the chunked route has to
    reproduce.
    """

    outdir = tmp_path_factory.mktemp("unchunked")
    result = subprocess.run(
        [
            sys.executable,
            str(REPO / "LRAA"),
            "--no_chunk",
            "--quant_only",
            "--gtf", inputs["gtf"],
            "--bam", inputs["bam"],
            "--genome", inputs["fasta"],
            "--output_prefix", "un",
        ],
        cwd=str(outdir),
        capture_output=True,
        text=True,
        timeout=900,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    path = outdir / "un.LRAA.quant-only.read_assignment.summary.tsv"
    assert path.exists(), sorted(p.name for p in outdir.iterdir())
    return str(path)


def _run_chunked(outdir, bam, fasta, num_total_reads, gtf=None, discovery=False,
                 approx_MB_per_cut=None):
    """make-chunks, then one leaf per chunk, then the units they produced.

    The scatter shape, in process: it is what a workflow runs, and it is the shape
    the merge is handed a unit list from.
    """

    extra = {}
    if approx_MB_per_cut is not None:
        extra["approx_MB_per_cut"] = approx_MB_per_cut
        extra["approx_MB_per_cut_wiggle_window"] = approx_MB_per_cut / 5.0
    ChunkedRun.run(
        ChunkedRun.default_args(
            bam=bam,
            genome_fa=fasta,
            gtf=gtf,
            output_dir=str(outdir),
            discovery=discovery,
            num_total_reads=num_total_reads,
            no_reuse_source_bam=True,
            stop_after_make_chunks=True,
            cpu_budget=2,
            **extra
        )
    )
    with open(os.path.join(str(outdir), "chunk_plan.json")) as fh:
        plan = json.load(fh)

    units = []
    for entry in plan["chunks"]:
        chunk_id = entry["chunk_id"]
        ChunkedRun.run(
            ChunkedRun.default_args(
                output_dir=str(outdir), only_chunk=chunk_id, cpu_budget=1
            )
        )
        with open(
            os.path.join(str(outdir), "chunks", chunk_id, "units.json")
        ) as fh:
            units.extend(json.load(fh)["units"])

    rank = {"+": 0, "-": 1}
    units.sort(key=lambda u: (rank[u["strand"]], u["order"]))
    return {"outdir": Path(str(outdir)), "plan": plan, "units": units}


@pytest.fixture(scope="module")
def chunked_units(inputs, tmp_path_factory):
    """The fixture's inputs through the chunked route.

    ``approx_MB_per_cut`` is forced down to 5 kb so a 24 kb contig is really cut
    into five chunks -- at the shipped 10 Mb it would be one chunk and the merge
    would have nothing to add up.
    """

    result = _run_chunked(
        tmp_path_factory.mktemp("chunked") / "work",
        inputs["bam"],
        inputs["fasta"],
        TOTAL_READS,
        gtf=inputs["gtf"],
        discovery=False,
        approx_MB_per_cut=0.005,
    )
    assert len(result["plan"]["chunks"]) > 1, result["plan"]["chunks"]
    return result


def test_extraction_partitions_the_library_without_dropping_a_read(chunked_units):
    """The parity test's premise, asserted separately so a failure localizes.

    Parity can only hold if the chunks are a partition: every read in exactly one
    chunk, none dropped for overhanging a boundary. The loci sit far from any cut,
    so this holds -- and if a geometry change ever breaks it, this test says so
    instead of leaving the parity failure to be misread as a summary bug.
    """

    emitted = 0
    for entry in chunked_units["plan"]["chunks"]:
        manifest_path = (
            chunked_units["outdir"] / "chunks" / entry["chunk_id"]
            / "chunk.partition.json"
        )
        with open(manifest_path) as fh:
            manifest = json.load(fh)
        assert manifest["counts"]["alignments_dropped_overhang"] == 0, entry
        emitted += manifest["counts"]["alignments_emitted"]

    assert emitted == TOTAL_READS


def test_chunked_reads_total_equals_the_unchunked_whole_contig_total(
    unchunked_summary, chunked_units, tmp_path
):
    """THE test. The same library, the same number, whichever route it took.

    Before count-then-report, every unit taking one of ``run_quant_only``'s early
    returns was simply absent from this table -- here that is the reverse strand
    everywhere plus the forward unit of the chunk holding locus C, so the chunked
    total came up 14 reads short of the unchunked one and said nothing about it.

    Fails if ANY of the three early returns stops writing its summary: the merge
    refuses a unit without one, and if it did not, the total would drop.
    """

    merged_dir = tmp_path / "merged"
    os.makedirs(merged_dir, exist_ok=True)
    merged = ChunkedRun.merge_read_assignment_summaries(
        str(merged_dir), chunked_units["units"]
    )

    assert merged is not None
    assert merged["units_absent"] == 0
    assert merged["units_merged"] == merged["units_total"] == len(
        chunked_units["units"]
    )
    assert merged["complete"] is True

    unchunked_total = _total_reads(unchunked_summary)
    assert unchunked_total == TOTAL_READS, "the fixture's own library size"
    assert _total_reads(merged["path"]) == unchunked_total

    # and the reads are there because a unit that assigned NOTHING still reported
    # them: at least one worker row carries reads with zero kept. Without this the
    # totals could agree by both being complete for the wrong reason.
    unassigned = [
        row
        for row in _summary_rows(merged["path"])
        if row["row_type"] == "worker"
        and int(row["reads_total"]) > 0
        and int(row["reads_kept_genome"]) == 0
    ]
    assert unassigned, "no unit reported seen-but-unassigned reads"
    assert sum(int(row["reads_total"]) for row in unassigned) == (
        REVERSE_READS + LOCUS_C_READS
    )


# -- stage 6 no longer tolerates a unit that reported nothing


def test_stage_six_refuses_a_unit_with_no_summary(chunked_units, tmp_path):
    """A missing summary is now a refusal, not a skip.

    It used to be tolerated on evidence -- accepted when the unit's ``quant.expr``
    held no data rows -- because that was the only way to survive the early
    returns. Those write now, so the tolerance covered nothing except the case it
    could not tell apart from a broken unit: a shard holding reads whose expr
    happens to be empty is exactly a shard whose reads would vanish from the
    total.

    The refusal has to NAME the unit, because on a 600-unit run "some unit is
    missing a file" is not actionable.
    """

    units = [dict(u) for u in chunked_units["units"]]
    victim = units[0]["unit_id"]
    units[0] = dict(units[0], quant_prefix=str(tmp_path / "nothing_here"))

    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.merge_read_assignment_summaries(str(tmp_path), units)

    message = str(err.value)
    assert victim in message
    assert "undercount" in message
    assert "read_assignment.summary.tsv" in message


def test_stage_six_refuses_even_when_no_unit_has_a_summary(tmp_path):
    """The other half of the old tolerance: returning None when NOT ONE unit had a
    summary.

    That existed so a chunk directory predating these files could still merge.
    Nothing writes such a directory any more, and the branch's effect on a real run
    was to publish quant output with no accounting at all beside it -- so it is a
    refusal too. ``units`` genuinely empty is still None, below: nothing to
    describe is not the same as everything missing.
    """

    units = [
        {"unit_id": "u{}".format(i), "quant_prefix": str(tmp_path / "u{}".format(i)),
         "offset": 0}
        for i in range(3)
    ]

    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.merge_read_assignment_summaries(str(tmp_path), units)
    assert "3 of 3" in str(err.value)

    assert ChunkedRun.merge_read_assignment_summaries(str(tmp_path), []) is None


# -- the paths that never reach run_quant_only at all


def _tiny_genome(root, contig, length, records):
    fasta = root / "g.fa"
    rng = random.Random(11)
    seq = "".join(rng.choice("ACGT") for _ in range(length))
    fasta.write_text(">{}\n{}\n".format(contig, seq))
    pysam.faidx(str(fasta))

    bam = root / "r.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": contig, "LN": length}],
    }
    records = sorted(records, key=lambda r: r[1])
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for name, start, span, flag in records:
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


def _merged_totals(units, merged_dir):
    os.makedirs(merged_dir, exist_ok=True)
    merged = ChunkedRun.merge_read_assignment_summaries(str(merged_dir), units)
    assert merged is not None
    assert merged["units_absent"] == 0
    assert merged["complete"] is True
    return _total_reads(merged["path"])


def test_a_discovery_unit_that_assembles_nothing_still_reports_its_reads(tmp_path):
    """``run_transcript_assembly`` has nine early returns of its OWN, none of which
    reaches run_quant_only.

    Nothing assembled, no draft mapped into the final splice graph, and one per
    isoform filter that emptied the set. MEASURED before this change on exactly this
    fixture -- one read per strand, too thin for any model to survive -- the unit
    wrote a quant.expr and NO summary at all, so both reads were absent from the
    run's accounting and stage 6 had no way to know.

    Fixed by counting once where a work unit FINISHES rather than at each of those
    returns, which is the only place a future early return cannot bypass.
    """

    fasta, bam = _tiny_genome(
        tmp_path,
        "chrD",
        8000,
        [("solo_f", 1000, 300, 0), ("solo_r", 5000, 300, 16)],
    )
    run = _run_chunked(tmp_path / "work", bam, fasta, 2, discovery=True)

    assert _merged_totals(run["units"], tmp_path / "merged") == 2


def test_an_oversimplified_contig_still_reports_its_reads(tmp_path):
    """chrM is oversimplified BY DEFAULT (``--oversimplify chrM,M``), and the
    oversimplify paths take no summary argument at all.

    So on any genome with a chrM, one chunk of every run produced no read-assignment
    summary -- MEASURED: both units of an oversimplified chrM chunk wrote quant
    output and no summary, and their 40 reads were missing from the total. Requiring
    a summary from every unit would have refused every such run outright, which is
    what makes this the test that keeps the refusal shippable.
    """

    records = [("m_{}".format(i), 1000 + (i % 5), 1000, 0) for i in range(30)]
    records += [("mr_{}".format(i), 4000, 800, 16) for i in range(10)]
    fasta, bam = _tiny_genome(tmp_path, "chrM", 8000, records)

    assert "chrM" in ChunkedRun.default_args().oversimplify
    run = _run_chunked(tmp_path / "work", bam, fasta, len(records), discovery=True)

    assert _merged_totals(run["units"], tmp_path / "merged") == len(records)