#!/usr/bin/env python3

"""The make-chunks phase as ONE flat queue, and the four things that can go wrong.

Stages 2 and 3 used to be serial loops. They now share one budget-sized pool over
(contig, strand) cut selections and (contig, strand, chunk) extractions, with each
contig's extractions submitted by that contig's own selection so there is no stage
barrier. Four properties carry the change, and each has a test here that fails
without it:

ENUMERATION IS NOT A CHOICE. Which contigs cut selection runs on has to stay
exactly the set the selector enumerates for itself. A smaller set -- "main
chromosomes", or "references idxstats reports reads for" -- silently drops both
the reads and the TPM denominator of what it omits, and changes the cut manifest
the parity gates diff. idxstats supplies the launch cost and nothing else, which
is tested here from both sides: a zero-mapped reference is still enumerated, and
it still gets a cost of zero.

ORDER NORMALISATION. Per-contig selections finish in whatever order the pool
finishes them in, and ``selections[key]`` feeds a monotonic ``order`` counter that
reaches the merged table's row order. So the assembled list is rebuilt in the
selector's own ``sorted(lengths)`` order, and the test shuffles completion
deliberately -- a test that let them complete in submission order would pass
against the bug. ``cut_placement_report`` is the second consumer of that list and
is checked the same way: it is a PURE function of the assembled selections, so it
needs no aggregation of its own and the shuffle must not move a line of it.

THE COLD INDEX. ``ensure_gtf_index`` has no lock: concurrent callers all miss the
same cold stamp and each build their own copy. The test asserts one build with N
callers AND, as its own positive control, N builds without the pre-warm -- because
"one build happened" proves nothing unless the unguarded case is shown to produce
more.

THE WHOLE-CONTIG SKIP. A strandless chunk spanning its whole contig would write a
mini bam that restates the source: offset 0, same contig name, nothing able to
overhang. The test holds the reused-source manifest to the EXTRACTED manifest's
own numbers rather than to a hand-written expectation, and holds the strand split
of the reused source to the record stream of a split of the extracted copy.
"""

import importlib.util
import json
import os
import random
import subprocess
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

import pysam
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "pylib"))

import ChunkedRun  # noqa: E402  (path insert must precede the import)

SEPARATE_BAM = REPO / "util" / "separate_bam_by_strand.py"


def _load(path, name):
    loader = SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


extractor = _load(
    REPO / "util" / "misc" / "extract_contig_region_inputs.py",
    "extract_contig_region_inputs_under_pcm_test",
)


# --------------------------------------------------------------------- fixtures


def write_genome(path, contigs):
    """A fasta of ``{name: length}``, indexed. Sequence is unremarkable on purpose."""

    with open(path, "wt") as ofh:
        for name, length in contigs.items():
            sequence = ("ACGT" * (length // 4 + 1))[:length]
            print(">" + name, file=ofh)
            for i in range(0, length, 60):
                print(sequence[i : i + 60], file=ofh)
    pysam.faidx(str(path))
    return str(path)


def write_reads(path, contigs, reads):
    """``reads`` is a list of (name, contig, start_1based, span, flag)."""

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": n, "LN": l} for n, l in contigs.items()],
    }
    order = list(contigs)
    records = []
    with pysam.AlignmentFile(str(path), "wb", header=header) as ofh:
        for name, contig, start, span, flag in reads:
            aln = pysam.AlignedSegment(ofh.header)
            aln.query_name = name
            aln.reference_id = order.index(contig)
            aln.reference_start = start - 1
            aln.mapping_quality = 60
            aln.cigarstring = "{}M".format(span)
            aln.query_sequence = "A" * span
            aln.query_qualities = pysam.qualitystring_to_array("I" * span)
            aln.flag = flag
            records.append(aln)
        records.sort(key=lambda a: (a.reference_id, a.reference_start))
        for aln in records:
            ofh.write(aln)
    pysam.index(str(path))
    return str(path)


def write_gtf(path, entries):
    """``entries`` is a list of (contig, gene_id, transcript_id, strand, lend, rend)."""

    with open(path, "wt") as ofh:
        for contig, gene, tx, strand, lend, rend in entries:
            for feature, attrs in (
                ("gene", 'gene_id "{}";'.format(gene)),
                (
                    "transcript",
                    'gene_id "{}"; transcript_id "{}";'.format(gene, tx),
                ),
                ("exon", 'gene_id "{}"; transcript_id "{}";'.format(gene, tx)),
            ):
                print(
                    "\t".join(
                        (
                            contig,
                            "test",
                            feature,
                            str(lend),
                            str(rend),
                            ".",
                            strand,
                            ".",
                            attrs,
                        )
                    ),
                    file=ofh,
                )
    return str(path)


def selection(chrom, contig_length, segments, retained=0):
    """A cuts.json selection as the selector writes one, reduced to what is read."""

    return {
        "chrom": chrom,
        "strand": None,
        "strandless": True,
        "contig_length": contig_length,
        "segments": [
            {
                "region": "{}:{}-{}".format(chrom, lend, rend),
                "lend": lend,
                "rend": rend,
                "span": rend - lend + 1,
                "window_origin": lend - 1,
            }
            for lend, rend in segments
        ],
        "counts": {"retained_primary_alignments": retained, "segments": len(segments)},
    }


# ------------------------------------------------------------------ enumeration


def test_enumeration_is_the_selectors_own_contig_set(tmp_path):
    """Every fasta reference, sorted, including ones holding no reads.

    The selector takes ``sorted(lengths)`` off the genome fasta. Anything narrower
    changes which contigs the run covers, and a contig that is not selected on has
    no chunk, contributes no reads and vanishes from the denominator.
    """

    contigs = {"cB": 4000, "cA": 5000, "cEmpty": 3000}
    genome = write_genome(tmp_path / "g.fa", contigs)
    bam = write_reads(
        tmp_path / "r.bam",
        contigs,
        [("r1", "cA", 100, 50, 0), ("r2", "cB", 100, 50, 16)],
    )
    args = ChunkedRun.default_args(genome_fa=genome, bam=bam)

    enumerated, lengths = ChunkedRun.enumerate_prep_contigs(args)
    assert enumerated == ["cA", "cB", "cEmpty"]
    assert lengths == contigs

    # the same set the selector would build for itself, computed the way it does
    with pysam.FastaFile(genome) as fasta:
        assert enumerated == sorted(fasta.references)


def test_idxstats_prices_the_queue_and_does_not_filter_it(tmp_path):
    """A zero-mapped reference is enumerated AND costed zero.

    Two claims that have to hold together: the reference is in the queue (so its
    empty chunk still exists, as the serial loop's did), and it sorts last (so a
    202-reference genome does not launch 92 no-ops ahead of chr1).
    """

    contigs = {"cA": 5000, "cB": 4000, "cEmpty": 3000}
    write_genome(tmp_path / "g.fa", contigs)
    bam = write_reads(
        tmp_path / "r.bam",
        contigs,
        [("r{}".format(i), "cB", 100 + 10 * i, 50, 0) for i in range(7)]
        + [("s1", "cA", 100, 50, 0)],
    )

    costs = ChunkedRun.contig_mapped_counts(bam)
    assert costs["cB"] == 7
    assert costs["cA"] == 1
    assert costs.get("cEmpty", 0) == 0

    # cB is the SHORTEST of the two that carry reads and still launches first:
    # the proxy is reads, not bases. That is the chrM shape -- 16.5 kb holding
    # 6.6 % of a library -- and launching it last is what leaves it in the tail.
    import CpuBudget

    units = [
        CpuBudget.WorkUnit(contig_acc=c, contig_strand="strandless", cost=costs.get(c, 0))
        for c in sorted(contigs)
    ]
    assert [u.contig_acc for u in CpuBudget.order_longest_first(units)] == [
        "cB",
        "cA",
        "cEmpty",
    ]


def test_an_absent_contig_restriction_is_refused_by_name(tmp_path):
    contigs = {"cA": 5000}
    genome = write_genome(tmp_path / "g.fa", contigs)
    bam = write_reads(tmp_path / "r.bam", contigs, [("r1", "cA", 100, 50, 0)])
    args = ChunkedRun.default_args(genome_fa=genome, bam=bam, contig="cZ")

    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.enumerate_prep_contigs(args)
    assert "cZ" in str(err.value)
    assert "cA" in str(err.value)


def test_a_contig_restriction_reduces_the_enumeration_to_one(tmp_path):
    contigs = {"cA": 5000, "cB": 4000}
    genome = write_genome(tmp_path / "g.fa", contigs)
    bam = write_reads(tmp_path / "r.bam", contigs, [("r1", "cA", 100, 50, 0)])
    args = ChunkedRun.default_args(genome_fa=genome, bam=bam, contig="cB")

    assert ChunkedRun.enumerate_prep_contigs(args)[0] == ["cB"]


def test_an_empty_genome_is_refused_rather_than_producing_an_empty_pool(tmp_path):
    empty = tmp_path / "empty.fa"
    empty.write_text("")
    (tmp_path / "empty.fa.fai").write_text("")
    args = ChunkedRun.default_args(genome_fa=str(empty), bam="unused.bam")

    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.enumerate_prep_contigs(args)
    assert "no contigs" in str(err.value)


# ----------------------------------------------------------- order normalisation


def _sources():
    return [("", ChunkedRun.STRANDLESS_TAG, "/dev/null", "tok")]


def test_selection_order_is_the_contig_order_not_the_completion_order():
    """The assembled list is what a serial run's would be, whatever order finishes.

    This is the one correctness requirement of the whole change: ``selections[key]``
    -> ``planned_chunks`` -> the ``order`` counter -> ``chunk_quant_units`` ->
    ``ordered_units`` -> ``merge_and_translate``. Reassembling in completion order
    moves the merged table's rows.
    """

    contigs = ["cA", "cB", "cC"]
    by_contig = {
        ("", "cA"): selection("cA", 900, [(1, 400), (401, 900)]),
        ("", "cB"): selection("cB", 500, [(1, 500)]),
        ("", "cC"): selection("cC", 700, [(1, 300), (301, 700)]),
    }
    serial = {"": [by_contig[("", c)] for c in contigs]}
    reference = [
        (p["chunk_id"], p["order"]) for p in ChunkedRun.planned_chunks(_sources(), serial)
    ]
    assert [cid for cid, _ in reference] == [
        "cA_00",
        "cA_01",
        "cB_00",
        "cC_00",
        "cC_01",
    ]
    assert [order for _, order in reference] == [0, 1, 2, 3, 4]

    rng = random.Random(20260817)
    for _ in range(24):
        completion = list(by_contig)
        rng.shuffle(completion)
        # what the pool holds: an unordered mapping filled in completion order
        pooled = {}
        for key in completion:
            pooled[key] = by_contig[key]
        # the normalisation under test
        assembled = {"": [pooled[("", c)] for c in contigs]}
        assert [
            (p["chunk_id"], p["order"])
            for p in ChunkedRun.planned_chunks(_sources(), assembled)
        ] == reference

    # POSITIVE CONTROL: completion order really can differ from contig order, so
    # the assertion above is not vacuous. Assembling in a shuffled order moves the
    # counters, which is exactly the regression.
    shuffled = {"": [by_contig[("", c)] for c in ["cB", "cC", "cA"]]}
    assert [
        (p["chunk_id"], p["order"])
        for p in ChunkedRun.planned_chunks(_sources(), shuffled)
    ] != reference


def test_quant_unit_order_survives_shuffled_completion():
    """``ordered_units`` sorts on ``order``, so a moved counter moves the merge."""

    contigs = ["cA", "cB"]
    by_contig = {
        ("", "cA"): selection("cA", 900, [(1, 400), (401, 900)]),
        ("", "cB"): selection("cB", 500, [(1, 500)]),
    }

    def units_for(order):
        selections = {"": [by_contig[("", c)] for c in order]}
        chunks = []
        for planned in ChunkedRun.planned_chunks(_sources(), selections):
            chunks.append(
                {
                    "units": ChunkedRun.chunk_quant_units(
                        planned["chunk_id"],
                        "/tmp/nowhere",
                        "/tmp/nowhere/chunk",
                        planned["key"],
                        0,
                        planned["order"],
                    )
                }
            )
        return [u["unit_id"] for u in ChunkedRun.ordered_units(chunks)]

    assert units_for(contigs) == [
        "cA_00_plus",
        "cA_01_plus",
        "cB_00_plus",
        "cA_00_minus",
        "cA_01_minus",
        "cB_00_minus",
    ]
    assert units_for(["cB", "cA"]) != units_for(contigs)


def _report_selection(chrom, contig_length, segments, cuts, declined):
    """A selection carrying everything ``cut_placement_report`` reads.

    ``selection`` above is reduced to what the chunk planner reads. The placement
    report reads the target accounting beside it and the per-cut severing under
    that. ``cuts`` is a list of ``(position, target, monoexonic, multi-exon)``.
    """

    sel = selection(chrom, contig_length, segments)
    sel["counts"].update(
        {
            "targets": len(segments) - 1 + len(declined),
            "cuts_placed": len(cuts),
            "targets_tail_merged": 0,
            "segments": len(segments),
            "alignments_dropped_at_cuts": sum(mono + multi for _p, _t, mono, multi in cuts),
            "alignments_dropped_monoexonic": sum(mono for _p, _t, mono, _x in cuts),
            "alignments_dropped_multiexon": sum(multi for _p, _t, _m, multi in cuts),
        }
    )
    sel["cuts"] = [
        {
            "position": position,
            "target": target,
            "offset_from_target": position - target,
            "spanning_alignments_dropped": mono + multi,
            "severed_monoexonic": mono,
            "severed_multiexon": multi,
            "search_radius": 500,
        }
        for position, target, mono, multi in cuts
    ]
    sel["unplaced_targets"] = [
        {
            "target": target,
            "reason": "the annotation left this window with no admissible position",
            "declined_annotation": True,
            "best_spanning_in_window": 7,
        }
        for target in declined
    ]
    return sel


def test_the_placement_report_is_the_serial_one_whatever_order_completes():
    """The report is derived from the ASSEMBLED selections, so it inherits their order.

    ``cut_placement_report`` is a PURE function of ``(selections, discovery)`` and
    needs no aggregation of its own: the pool contributes a mapping from
    ``(key, contig)`` to that contig's selection, and the contig-order reassembly
    is what fixes the report. Completion order is shuffled deliberately -- a
    report accumulated as each unit finished would be a second aggregation, and it
    would order its per-contig lines and its ``per_selection`` list by whichever
    contig happened to finish first.
    """

    contigs = ["cA", "cB", "cC"]
    by_contig = {
        ("", "cA"): _report_selection(
            "cA", 900, [(1, 400), (401, 900)], [(400, 402, 3, 1)], []
        ),
        ("", "cB"): _report_selection("cB", 500, [(1, 500)], [], [250]),
        ("", "cC"): _report_selection(
            "cC", 700, [(1, 300), (301, 700)], [(300, 300, 0, 2)], []
        ),
    }
    reference = ChunkedRun.cut_placement_report(
        {"": [by_contig[("", c)] for c in contigs]}, True
    )
    # Not a vacuous comparison: the report saw cuts, severing and a decline, so
    # there is something in it that an order could move.
    assert reference[1]["cuts_placed"] == 2
    assert reference[1]["alignments_severed"] == 6
    assert reference[1]["alignments_severed_multiexon"] == 3
    assert reference[1]["cuts_declined_annotation"] == 1
    assert [s["label"] for s in reference[1]["per_selection"]] == ["cA", "cB", "cC"]

    rng = random.Random(20260818)
    for _ in range(24):
        completion = list(by_contig)
        rng.shuffle(completion)
        # what the pool holds: an unordered mapping filled in completion order
        pooled = {}
        for key in completion:
            pooled[key] = by_contig[key]
        assembled = {"": [pooled[("", c)] for c in contigs]}
        assert ChunkedRun.cut_placement_report(assembled, True) == reference

    # POSITIVE CONTROL: order really does reach the report, so the equality above
    # is a statement about the reassembly rather than about a report that cannot
    # tell its inputs apart.
    shuffled = {"": [by_contig[("", c)] for c in ["cB", "cC", "cA"]]}
    assert ChunkedRun.cut_placement_report(shuffled, True) != reference


def test_a_whole_contig_segment_is_recognised_and_a_partial_one_is_not():
    """``spans_whole_contig`` is what routes a chunk past extraction."""

    whole = selection("cM", 16569, [(1, 16569)])
    split = selection("cA", 900, [(1, 400), (401, 900)])
    flags = {
        p["chunk_id"]: p["spans_whole_contig"]
        for p in ChunkedRun.planned_chunks(
            _sources(), {"": [split, whole]}
        )
    }
    assert flags == {"cA_00": False, "cA_01": False, "cM_00": True}


def test_a_selection_without_a_contig_length_never_claims_whole_contig():
    """Absent evidence is not a whole contig. An old cuts.json has no length."""

    stale = selection("cA", 900, [(1, 900)])
    del stale["contig_length"]
    planned = list(ChunkedRun.planned_chunks(_sources(), {"": [stale]}))
    assert planned[0]["spans_whole_contig"] is False


# ------------------------------------------------- whole-invocation aggregation
#
# Five of the selector's artifacts are per-INVOCATION rather than per-contig, so
# N invocations produce N of each and something has to put them back together.
# `.cuts.json` is reassembled in contig order (above). `.dropped_reads.txt` is
# unioned into a SET, which is the one aggregation a naive parallelisation gets
# wrong quietly: the severed-read accounting compares selection NAMED against
# extraction DROPPED, so a lost or duplicated file changes a set comparison into
# a failure or, worse, a control that subtracts the wrong records.


def test_severed_names_union_across_per_contig_files(tmp_path):
    """N per-contig files must give the same SET one whole-invocation file gave."""

    serial_dir = tmp_path / "serial"
    serial_dir.mkdir()
    (serial_dir / "strandless.dropped_reads.txt").write_text(
        "readA\nreadB\nreadC\nreadD\n"
    )

    parallel_dir = tmp_path / "parallel"
    parallel_dir.mkdir()
    # per contig, out of order, one contig severing nothing, and readC severed on
    # two contigs -- a read name is not unique to a contig and the union has to
    # collapse it exactly as the one-file form did
    (parallel_dir / "cB.strandless.dropped_reads.txt").write_text("readC\nreadA\n")
    (parallel_dir / "cA.strandless.dropped_reads.txt").write_text("readB\nreadC\n")
    (parallel_dir / "cC.strandless.dropped_reads.txt").write_text("")
    (parallel_dir / "cD.strandless.dropped_reads.txt").write_text("readD\n")

    assert ChunkedRun.severed_read_names(str(serial_dir)) == ChunkedRun.severed_read_names(
        str(parallel_dir)
    )
    assert len(ChunkedRun.severed_read_names(str(parallel_dir))) == 4


def test_the_severed_accounting_gate_still_closes_over_per_contig_files(tmp_path):
    """The accounting reads N per-contig files, in every direction it can go.

    This is the correctness gate the whole change is checked against, so it is
    tested in the shape the change puts it in rather than only in the shape it had.
    Since d6ae82f the two directions of a mismatch are not the same defect: reads
    NAMED by selection that no chunk dropped are still fatal, and reads DROPPED by
    extraction that selection never named are tolerated, persisted and reported.
    Both directions are exercised here over the multi-file shape, because the union
    is what the pool changed and either direction could be read off the wrong set
    of files.

    The accounting itself is ONE call over the whole chunk list, made after the pool
    drains (``run``: verify_severed_accounting(cut_dir, chunks)), which is what makes
    a single ``EXTRACTION_ONLY_DROPS`` file correct: a per-contig invocation would
    overwrite it and the baseline would quietly stop subtracting the earlier
    contigs' tolerated reads.
    """

    cut_dir = tmp_path / "cuts"
    cut_dir.mkdir()
    (cut_dir / "cA.strandless.dropped_reads.txt").write_text("readA\n")
    (cut_dir / "cB.strandless.dropped_reads.txt").write_text("readB\nreadC\n")
    tolerated = cut_dir / ChunkedRun.EXTRACTION_ONLY_DROPS

    def chunk(cid, names):
        return {"chunk_id": cid, "manifest": {"dropped_read_names": list(names)}}

    # a severed read is dropped by BOTH neighbouring chunks, so mentions exceed
    # the distinct count -- the accounting keeps the two apart on purpose
    counts = ChunkedRun.verify_severed_accounting(
        str(cut_dir),
        [
            chunk("cA_00", ["readA"]),
            chunk("cA_01", ["readA"]),
            chunk("cB_00", ["readB"]),
            chunk("cB_01", ["readB", "readC"]),
            chunk("cB_02", ["readC"]),
        ],
    )
    assert counts["named_by_cut_selection"] == 3
    assert counts["dropped_by_extraction"] == 3
    assert counts["per_chunk_drop_mentions"] == 6
    assert counts["sets_identical"] is True
    assert counts["dropped_not_named"] == 0
    # exact accounting writes nothing: the tolerated-drops file exists only when
    # there is something to tolerate
    assert not tolerated.exists()

    # One contig's file going missing is what a lost per-contig artifact looks
    # like, and it lands in the TOLERATED direction: extraction dropped readB and
    # readC, and the surviving selection files name neither. The run continues.
    (cut_dir / "cB.strandless.dropped_reads.txt").unlink()
    absorbed = ChunkedRun.verify_severed_accounting(
        str(cut_dir),
        [chunk("cA_00", ["readA"]), chunk("cB_00", ["readB", "readC"])],
    )
    assert absorbed["dropped_not_named"] == 2
    assert absorbed["sets_identical"] is False
    assert tolerated.read_text().split() == ["readB", "readC"]

    # THE LOAD-BEARING HALF, and the reason "it did not raise" is not enough:
    # run_baseline prunes the control bam by this union, with no exclusion, so the
    # tolerated reads have to be in it or the two arms consume different records
    # while the run reports success.
    assert ChunkedRun.severed_read_names(str(cut_dir)) == {"readA", "readB", "readC"}
    # and the other half, over N files: the persisted file matches the same glob as
    # the per-contig selection artifacts, so the SELECTION question must not pick
    # it up -- otherwise a rerun reclassifies tolerated reads as predictions.
    assert ChunkedRun.severed_read_names(
        str(cut_dir), exclude=(ChunkedRun.EXTRACTION_ONLY_DROPS,)
    ) == {"readA"}

    # IDEMPOTENT over the per-contig glob: the same run again reports the same two
    # tolerated reads rather than exact accounting that nothing measured.
    again = ChunkedRun.verify_severed_accounting(
        str(cut_dir),
        [chunk("cA_00", ["readA"]), chunk("cB_00", ["readB", "readC"])],
    )
    assert again == absorbed
    assert tolerated.read_text().split() == ["readB", "readC"]

    # POSITIVE CONTROL: the gate is still a gate. A read NAMED in one contig's
    # file that no chunk dropped is the direction that stays fatal, and it has to
    # be found across the N files rather than only in the first one.
    (cut_dir / "cC.strandless.dropped_reads.txt").write_text("readD\n")
    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.verify_severed_accounting(
            str(cut_dir),
            [chunk("cA_00", ["readA"]), chunk("cB_00", ["readB", "readC"])],
        )
    assert "readD" in str(err.value)


def test_per_contig_prefixes_cannot_collide(tmp_path):
    """Every per-invocation artifact is keyed on the contig AND the orientation.

    N invocations sharing a prefix would leave only the last contig's cuts.json,
    dropped-read list and severed bam on disk -- and the run would look fine,
    because the aggregation reads whatever is there.
    """

    contigs = ["cA", "cB"]
    sources = [
        ("+", "plus", "/dev/null", "tok"),
        ("-", "minus", "/dev/null", "tok"),
    ]
    args = ChunkedRun.default_args(genome_fa="g.fa", gtf="a.gtf")
    plans = [
        ChunkedRun.cut_selection_plan(args, "/tmp/o", "/tmp/o/cuts", source, contig)
        for source in sources
        for contig in contigs
    ]
    for field in ("prefix", "log", "token", "unit_id"):
        values = [p[field] for p in plans]
        assert len(set(values)) == len(plans), field
    # and the severed bam, which is per invocation and not merged
    sev = [p["cmd"][p["cmd"].index("--severed_reads_bam") + 1] for p in plans]
    assert len(set(sev)) == len(plans)
    # every invocation is restricted to its own contig
    for plan in plans:
        assert plan["cmd"][plan["cmd"].index("--contig") + 1] == plan["contig"]


# ---------------------------------------------------------------- gtf pre-warm


CALLER = '''
import os, sys, time
from importlib.machinery import SourceFileLoader
import importlib.util
loader = SourceFileLoader("ex", sys.argv[1])
spec = importlib.util.spec_from_loader(loader.name, loader)
ex = importlib.util.module_from_spec(spec)
loader.exec_module(ex)
gtf, cache, marker, arrivals, expected = sys.argv[2:7]
expected = int(expected)
real = ex._build_gtf_index
def _arrived():
    with open(marker, "rt") as fh:
        return len([line for line in fh if line.strip()])
def counting(gtf_filename, gz_path):
    with open(marker, "at") as fh:
        fh.write("{}\\n".format(os.getpid()))
    # SECOND barrier, inside the build. On the real input the window is wide by
    # construction -- 36.2 s to index a 1.49 GB GENCODE GTF, which every other
    # caller spends looking at a cold stamp. This fixture's GTF builds in
    # microseconds, so without the hold the first caller past the outer barrier
    # can finish and write the stamp before the rest have read it: observed 5 of
    # 6 and then 1 of 6 on a loaded box. A caller the pre-warm turned away never
    # reaches this function, so the hold cannot manufacture a duplicate build --
    # it only stops the count from measuring the scheduler.
    stop = time.time() + 60
    while _arrived() < expected and time.time() < stop:
        time.sleep(0.01)
    return real(gtf_filename, gz_path)
ex._build_gtf_index = counting
# Barrier: nobody consults the stamp until everybody has started, so a cold
# cache is cold for all of them at the moment they look. Without this the
# count would depend on process start-up jitter.
open(os.path.join(arrivals, str(os.getpid())), "wt").close()
deadline = time.time() + 60
while len(os.listdir(arrivals)) < expected and time.time() < deadline:
    time.sleep(0.01)
ex.ensure_gtf_index(gtf, cache_dir=cache)
'''


def _concurrent_callers(tmp_path, gtf, cache, n):
    """Run ``n`` real PROCESSES through ensure_gtf_index at once. Returns the pids
    that reached ``_build_gtf_index``.

    Processes rather than threads, and not for tidiness: ``_build_gtf_index``
    scopes its temporary files by PID precisely so that concurrent builders do not
    clobber each other, so N threads in one process is a collision the real
    fan-out cannot have. Threads were tried first and produced exactly that
    artefact -- four of six builders failing on each other's temps -- which would
    have made this control measure the wrong thing.
    """

    script = tmp_path / "caller.py"
    script.write_text(CALLER)
    marker = tmp_path / "builds.txt"
    marker.write_text("")
    arrivals = tmp_path / "arrivals"
    arrivals.mkdir()
    procs = [
        subprocess.Popen(
            [
                sys.executable,
                str(script),
                str(REPO / "util" / "misc" / "extract_contig_region_inputs.py"),
                gtf,
                cache,
                str(marker),
                str(arrivals),
                str(n),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        for _ in range(n)
    ]
    for proc in procs:
        out, err = proc.communicate(timeout=180)
        assert proc.returncode == 0, err.decode()[-2000:]
    return [line for line in marker.read_text().splitlines() if line.strip()]


def test_the_gtf_index_is_built_once_for_the_whole_pool(tmp_path, monkeypatch):
    """One serial pre-warm, then N concurrent units build NOTHING.

    Paired with its positive control below, because "one build" is only meaningful
    against a measured alternative.
    """

    gtf = write_gtf(tmp_path / "a.gtf", [("cA", "gA", "tA", "+", 100, 900)])
    outdir = tmp_path / "out"
    outdir.mkdir()
    args = ChunkedRun.default_args(gtf=gtf)

    builds = []
    real_build = extractor._build_gtf_index

    def counting_build(gtf_filename, gz_path):
        builds.append(gz_path)
        return real_build(gtf_filename, gz_path)

    monkeypatch.setattr(extractor, "_build_gtf_index", counting_build)
    monkeypatch.setattr(ChunkedRun, "_extractor_module", lambda: extractor)

    timing = {}
    index = ChunkedRun.warm_gtf_index(args, str(outdir), timing)
    assert index is not None
    assert len(builds) == 1
    assert timing["stages"]["gtf_index_warm"]["index"] == index

    cache = ChunkedRun.gtf_index_cache_dir(str(outdir))
    assert _concurrent_callers(tmp_path, gtf, cache, 6) == []


def test_without_the_prewarm_concurrent_callers_each_build_their_own(tmp_path):
    """The positive control. ensure_gtf_index has no lock, and this is the cost.

    Measured on a 1.49 GB GENCODE GTF with ``_build_gtf_index`` pid-marked: six
    concurrent callers against a cold cache produced SIX independent builds, at
    36.2 s and ~1.5 GB of temp each. Nothing is corrupted -- the temps are
    pid-scoped and os.replace'd into place -- so the pre-warm is a cost fix, and
    this test measures the count it fixes.

    SIX, not "more than one", and the fixture earns that. The real GTF's 36 s
    build holds the window open on its own; this one builds in microseconds, so
    ``_concurrent_callers`` holds every builder INSIDE the build until all six are
    in it. Without that hold the count measures the scheduler -- observed 5 of 6
    and then 1 of 6 on a loaded box. The hold cannot manufacture the result it
    measures: a caller the pre-warm turns away never reaches the build at all,
    which is the zero the guarded test above asserts.
    """

    gtf = write_gtf(tmp_path / "b.gtf", [("cA", "gA", "tA", "+", 100, 900)])
    cache = str(tmp_path / "cold_cache")
    os.makedirs(cache)

    pids = _concurrent_callers(tmp_path, gtf, cache, 6)
    assert len(pids) == 6
    # every build came from a DIFFERENT process, so this is concurrency rather
    # than one caller building twice
    assert len(set(pids)) == 6
    # and the result is still usable: duplicated work, never a corrupt index
    assert extractor.ensure_gtf_index(gtf, cache_dir=cache) is not None


# ---------------------------------------------------------------- memory bound


def test_the_pool_is_capped_by_memory_when_the_box_cannot_hold_the_budget():
    cap, note = ChunkedRun.prep_memory_cap(32, available_mib=2048)
    assert cap == 2048 // ChunkedRun.PREP_UNIT_PEAK_MIB
    assert cap < 32
    assert "MEMORY" in note
    assert str(ChunkedRun.PREP_UNIT_PEAK_MIB) in note
    assert "2048" in note


def test_an_ample_box_is_not_capped_at_all():
    # 62 GB against 16 workers: 16 x 300 MiB is 4.7 GiB, 13x headroom
    assert ChunkedRun.prep_memory_cap(16, available_mib=62 * 1024) == (None, None)
    # 16 GB laptop, same width, still fits
    assert ChunkedRun.prep_memory_cap(16, available_mib=16 * 1024) == (None, None)
    # 16 GB at 24-wide is 7.03 GiB, and MemAvailable on such a box is well over it
    assert ChunkedRun.prep_memory_cap(24, available_mib=15 * 1024) == (None, None)


def test_a_squeezed_box_still_runs_one_unit_rather_than_refusing():
    """One at a time is what the serial loop did. Refusing would be worse."""

    cap, note = ChunkedRun.prep_memory_cap(8, available_mib=10)
    assert cap == 1
    assert note


def test_an_unreadable_meminfo_says_so_instead_of_pretending_to_have_checked():
    cap, note = ChunkedRun.prep_memory_cap(8, available_mib=None)
    if ChunkedRun.available_memory_mib() is None:
        assert cap is None
        assert "could not be read" in note


def test_the_memory_cap_reaches_the_allocation():
    """The bound is only a bound if ``allocate`` is given it."""

    import CpuBudget

    cap, _note = ChunkedRun.prep_memory_cap(32, available_mib=1200)
    assert cap == 4
    allocation = CpuBudget.allocate(budget=32, num_units=32, max_unit_workers=cap)
    assert allocation.unit_workers == 4
    uncapped = CpuBudget.allocate(budget=32, num_units=32)
    assert uncapped.unit_workers == 32


# --------------------------------------------------- the whole-contig extraction


class WholeContigCorpus:
    """One contig, both orientations, plus records the retention filter rejects.

    The rejected records are the point: a reused source still holds the secondary,
    supplementary and duplicate records that a written mini bam would have left
    behind, so this is what proves the split of the source keeps the same set the
    extraction would have written.
    """

    contig = "cWhole"
    length = 6000

    def __init__(self, tmp_path):
        self.contigs = {self.contig: self.length, "cOther": 4000}
        self.genome = write_genome(tmp_path / "wc.fa", self.contigs)
        reads = []
        for i in range(6):
            reads.append(("keep_f{}".format(i), self.contig, 500 + 100 * i, 80, 0))
        for i in range(4):
            reads.append(("keep_r{}".format(i), self.contig, 2000 + 100 * i, 80, 16))
        # not retained: secondary, supplementary, duplicate
        reads.append(("drop_secondary", self.contig, 3000, 80, 0x100))
        reads.append(("drop_supplementary", self.contig, 3100, 80, 0x800))
        reads.append(("drop_duplicate", self.contig, 3200, 80, 0x400))
        # a different contig, which a reused source still carries and the
        # restriction has to keep out
        reads.append(("other_contig", "cOther", 500, 80, 0))
        self.bam = write_reads(tmp_path / "wc.bam", self.contigs, reads)
        self.gtf = write_gtf(
            tmp_path / "wc.gtf",
            [
                (self.contig, "gF", "tF", "+", 400, 1400),
                (self.contig, "gR", "tR", "-", 1900, 2600),
                ("cOther", "gO", "tO", "+", 400, 900),
            ],
        )
        self.retained_forward = 6
        self.retained_reverse = 4

    def region(self):
        return "{}:1-{}".format(self.contig, self.length)

    def extract(self, prefix, reuse):
        argv = [
            "--genome_fa",
            self.genome,
            "--bam",
            self.bam,
            "--gtf",
            self.gtf,
            "--region",
            self.region(),
            "--output_prefix",
            str(prefix),
            "--max_intron_length",
            "200000",
            "--secondary_alignments",
            "exclude",
        ]
        if reuse:
            argv.append("--reuse_source_bam")
        assert extractor.main(argv) == 0
        return json.loads(Path(str(prefix) + ".partition.json").read_text())


def test_a_whole_contig_chunk_reuses_the_source_and_states_the_same_numbers(tmp_path):
    """The reused manifest is held to the EXTRACTED manifest, not to a guess.

    Everything but ``files`` and the two reuse flags has to agree, because the
    reuse path runs the same loop over the same fetch with the same predicates and
    only declines to write. If the numbers were handed in from elsewhere this
    comparison would be the place it showed.
    """

    corpus = WholeContigCorpus(tmp_path)
    full_dir = tmp_path / "full"
    reuse_dir = tmp_path / "reuse"
    full_dir.mkdir()
    reuse_dir.mkdir()

    full = corpus.extract(full_dir / "chunk", reuse=False)
    reused = corpus.extract(reuse_dir / "chunk", reuse=True)

    assert full["bam_reused_from_source"] is False
    assert reused["bam_reused_from_source"] is True
    assert full["spans_whole_contig"] is True and reused["spans_whole_contig"] is True

    ignored = {"files", "bam_reused_from_source"}
    assert {k: v for k, v in reused.items() if k not in ignored} == {
        k: v for k, v in full.items() if k not in ignored
    }

    # non-vacuity: the filter really did reject records, so "same counts" is a
    # statement about a filtered set and not about a pass-through
    assert full["counts"]["alignments_emitted"] == (
        corpus.retained_forward + corpus.retained_reverse
    )
    assert full["counts"]["nonprimary_excluded"] == 2
    assert full["counts"]["alignments_dropped_overhang"] == 0
    assert full["dropped_read_names"] == []

    # no mini bam was written, and files.bam names the source
    assert not (reuse_dir / "chunk.bam").exists()
    assert reused["files"]["bam"] == os.path.abspath(corpus.bam)
    assert Path(full["files"]["bam"]).exists()

    # the cheap artifacts are still written: LRAA derives its contig list from
    # the genome fasta, so a mini fasta is not optional
    assert (reuse_dir / "chunk.fa").exists()
    assert (reuse_dir / "chunk.gtf").exists()
    assert (
        Path(reuse_dir / "chunk.fa").read_text()
        == Path(full_dir / "chunk.fa").read_text()
    )
    assert (
        Path(reuse_dir / "chunk.gtf").read_text()
        == Path(full_dir / "chunk.gtf").read_text()
    )


def test_splitting_the_reused_source_by_contig_equals_splitting_the_copy(tmp_path):
    """The downstream stage really does accept the source in place of the copy.

    ``separate_bam_by_strand`` retains exactly what ``retained_for_extraction``
    retains, so the two record streams must be equal -- with the source's extra
    contig and its rejected records both kept out. Compared as SAM text without
    the header, because a reused source keeps its own header and that is the one
    thing that legitimately differs.
    """

    corpus = WholeContigCorpus(tmp_path)
    full_dir = tmp_path / "full"
    full_dir.mkdir()
    corpus.extract(full_dir / "chunk", reuse=False)

    def split(bam, prefix, contig=None):
        cmd = [
            sys.executable,
            str(SEPARATE_BAM),
            "--bam",
            str(bam),
            "--output_prefix",
            str(prefix),
            "--max_intron_length",
            "200000",
        ]
        if contig:
            cmd += ["--contig", contig]
        subprocess.run(cmd, check=True, capture_output=True)
        return {
            s: sorted(
                subprocess.run(
                    ["samtools", "view", "{}.{}.bam".format(prefix, s)],
                    check=True,
                    capture_output=True,
                    text=True,
                ).stdout.splitlines()
            )
            for s in ("+", "-")
        }

    from_copy = split(full_dir / "chunk.bam", tmp_path / "from_copy")
    from_source = split(corpus.bam, tmp_path / "from_source", contig=corpus.contig)

    assert from_copy == from_source
    assert len(from_copy["+"]) == corpus.retained_forward
    assert len(from_copy["-"]) == corpus.retained_reverse
    assert not any("other_contig" in line for line in from_source["+"])
    assert not any("drop_" in line for line in from_source["+"] + from_source["-"])

    # POSITIVE CONTROL for the restriction: without --contig the source's other
    # contig leaks in, so the restriction is doing the work and not the fixture.
    unrestricted = split(corpus.bam, tmp_path / "unrestricted")
    assert unrestricted != from_source
    assert any("other_contig" in line for line in unrestricted["+"])


def test_a_partial_region_refuses_to_reuse_the_source(tmp_path):
    corpus = WholeContigCorpus(tmp_path)
    with pytest.raises(extractor.ExtractionError) as err:
        extractor.extract_partition(
            genome_fa=corpus.genome,
            bam=corpus.bam,
            region="{}:1-3000".format(corpus.contig),
            output_prefix=str(tmp_path / "partial"),
            gtf=corpus.gtf,
            max_intron_length=200000,
            reuse_source_bam=True,
        )
    assert "whole contig" in str(err.value)
    assert not (tmp_path / "partial.bam").exists()


def test_a_strand_suffixed_region_refuses_to_reuse_the_source(tmp_path):
    corpus = WholeContigCorpus(tmp_path)
    with pytest.raises(extractor.ExtractionError) as err:
        extractor.extract_partition(
            genome_fa=corpus.genome,
            bam=corpus.bam,
            region="{}+:1-{}".format(corpus.contig, corpus.length),
            output_prefix=str(tmp_path / "stranded"),
            gtf=corpus.gtf,
            max_intron_length=200000,
            reuse_source_bam=True,
        )
    assert "STRANDLESS" in str(err.value)


def test_an_unindexed_source_is_refused_before_a_consumer_meets_it(tmp_path):
    """A manifest naming an unindexed source is a pysam error one stage later."""

    corpus = WholeContigCorpus(tmp_path)
    with pytest.raises(extractor.ExtractionError) as err:
        extractor._require_source_index(str(tmp_path / "nonexistent.bam"))
    assert "samtools index" in str(err.value)
    # and the real one passes, so the guard is not simply always failing
    extractor._require_source_index(corpus.bam)


def test_only_strandless_whole_contig_chunks_are_planned_for_reuse(tmp_path):
    """Strand-first never reuses: its chunk bam is normalized and quantified.

    A strand-first chunk bam has consumers that cannot restrict by contig, and the
    strand bam behind it holds every other contig, so the substitution is not
    available there however whole the region is.
    """

    outdir = tmp_path / "out"
    (outdir / "logs").mkdir(parents=True)
    chunk_root = outdir / "chunks"
    chunk_root.mkdir()
    args = ChunkedRun.default_args(
        genome_fa="g.fa", bam="r.bam", gtf="a.gtf", margin=200
    )

    whole = selection("cM", 16569, [(1, 16569)])
    strandless_plans = [
        ChunkedRun.extraction_plan(args, str(outdir), str(chunk_root), p, "tok")
        for p in ChunkedRun.planned_chunks(_sources(), {"": [whole]})
    ]
    assert [p["reuse_source_bam"] for p in strandless_plans] == [True]
    assert "--reuse_source_bam" in strandless_plans[0]["cmd"]

    strand_first = [("+", "plus", "/dev/null", "tok")]
    whole_plus = selection("cM", 16569, [(1, 16569)])
    strand_plans = [
        ChunkedRun.extraction_plan(args, str(outdir), str(chunk_root), p, "tok")
        for p in ChunkedRun.planned_chunks(strand_first, {"+": [whole_plus]})
    ]
    assert [p["reuse_source_bam"] for p in strand_plans] == [False]
    assert "--reuse_source_bam" not in strand_plans[0]["cmd"]

    # the sentinel has to move with the decision, or a directory extracted one way
    # is reported reusable by a run that would have done the other
    assert strandless_plans[0]["token"] != strand_plans[0]["token"]


def test_the_split_reads_the_manifests_bam_and_restricts_it(tmp_path):
    """Stage 3b's command follows the manifest rather than assuming a mini bam."""

    cdir = tmp_path / "chunks" / "cM_00"
    cdir.mkdir(parents=True)
    args = ChunkedRun.default_args(max_intron_length=200000)
    ckpt = ChunkedRun.Checkpoints(str(tmp_path / "__ckpt"))

    def chunk_for(reused):
        return {
            "chunk_id": "cM_00",
            "chrom": "cM",
            "dir": str(cdir),
            "prefix": str(cdir / "chunk"),
            "log": str(tmp_path / "cM.log"),
            "upstream_token": "up",
            "manifest": {
                "strand": None,
                "strand_split_required": True,
                "bam_reused_from_source": reused,
                "files": {"bam": "/srv/source.bam"},
                "counts": {"alignments_emitted": 0},
            },
        }

    calls = []

    def capture(name, cmd, log, cwd, rss, append=True):
        calls.append(cmd)
        raise ChunkedRun.PipelineError("stop after the command is built")

    import unittest.mock

    with unittest.mock.patch.object(ChunkedRun, "run_step", capture):
        for reused in (True, False):
            with pytest.raises(ChunkedRun.PipelineError):
                ChunkedRun.split_chunk_by_strand(
                    args, ckpt, chunk_for(reused), 0.25
                )

    reused_cmd, copy_cmd = calls
    assert "/srv/source.bam" in reused_cmd
    assert reused_cmd[reused_cmd.index("--contig") + 1] == "cM"
    assert str(cdir / "chunk.bam") in copy_cmd
    assert "--contig" not in copy_cmd


# ----------------------------------------------------------- no stage barrier


def _prep_corpus(tmp_path, contigs):
    genome = write_genome(tmp_path / "g.fa", contigs)
    bam = write_reads(
        tmp_path / "r.bam",
        contigs,
        [("r{}".format(i), c, 100 + 50 * i, 80, 0) for i, c in enumerate(contigs)],
    )
    gtf = write_gtf(
        tmp_path / "a.gtf", [(c, "g" + c, "t" + c, "+", 200, 900) for c in contigs]
    )
    return genome, bam, gtf


def test_an_extraction_runs_while_another_contig_is_still_selecting(
    tmp_path, monkeypatch
):
    """THE design claim, checked directly rather than inferred from wall time.

    A contig's extractions are submitted by that contig's own selection, so there
    is no point at which every selection has to be finished. With a stage barrier
    the makespan floor is the largest chromosome's WHOLE prep -- chr1 at 121.6 s
    measured, an 8.3x ceiling more cores cannot move. Without one the floor is the
    longest single dependency chain, one selection plus its longest extraction.

    Made deterministic with an event rather than with sleeps: cB's selection blocks
    until cA's extraction has run, so the test can only pass if the extraction was
    reachable while a selection was outstanding, and can only hang if it was not.
    """

    import threading

    contigs = {"cA": 5000, "cB": 4000}
    genome, bam, gtf = _prep_corpus(tmp_path, contigs)
    outdir = tmp_path / "out"
    (outdir / "logs").mkdir(parents=True)
    args = ChunkedRun.default_args(
        genome_fa=genome, bam=bam, gtf=gtf, strandless_chunks=True, cpu_budget=2
    )
    ckpt = ChunkedRun.Checkpoints(str(outdir / "__ckpt"))
    monkeypatch.setattr(ChunkedRun, "warm_gtf_index", lambda *a, **k: None)

    extraction_ran = threading.Event()
    saw_barrier_free = {}

    def selecting(plan, ckpt_, outdir_, rss):
        if plan["contig"] == "cB":
            # released by cA's EXTRACTION, not by cA's selection
            assert extraction_ran.wait(timeout=60), (
                "cA's extraction never ran while cB was still selecting: the "
                "queue has a stage barrier"
            )
        return {"wall_s": 0.0}, selection(
            plan["contig"], contigs[plan["contig"]], [(1, contigs[plan["contig"]])]
        )

    def extracting(plan, ckpt_, outdir_, rss):
        saw_barrier_free[plan["chunk_id"]] = not extraction_ran.is_set()
        extraction_ran.set()
        return (
            {"wall_s": 0.0},
            {
                "strand": None,
                "strand_split_required": True,
                "offset": 0,
                "window_origin": 0,
                "counts": {"alignments_emitted": 0},
                "dropped_read_names": [],
                "files": {"bam": "/dev/null"},
            },
        )

    monkeypatch.setattr(ChunkedRun, "run_cut_selection", selecting)
    monkeypatch.setattr(ChunkedRun, "run_extraction", extracting)

    timing = {}
    selections, cut_dir, chunks, makespan, allocation = (
        ChunkedRun.run_prep_concurrently(
            args, ckpt, str(outdir), timing, _sources_for(args), 0.25
        )
    )

    assert saw_barrier_free["cA_00"] is True
    assert [c["chunk_id"] for c in chunks] == ["cA_00", "cB_00"]
    # and the assembled order is contig order even though cB finished last
    assert [s["chrom"] for s in selections[""]] == ["cA", "cB"]
    assert allocation.unit_workers == 2


def test_a_single_contig_invocation_still_pools_its_chunks(tmp_path, monkeypatch):
    """The WDL shard case: one contig, many chunks, and the pool must be wide.

    A per-contig WDL shard seeds ONE selection, so sizing the pool from the seed
    length would run its extractions one at a time -- on a chr1 shard that is
    103.9 s of measured serial extraction inside a 5-core task.
    """

    contigs = {"cA": 900000}
    genome, bam, gtf = _prep_corpus(tmp_path, contigs)
    outdir = tmp_path / "out"
    (outdir / "logs").mkdir(parents=True)
    args = ChunkedRun.default_args(
        genome_fa=genome, bam=bam, gtf=gtf, strandless_chunks=True, cpu_budget=6
    )
    ckpt = ChunkedRun.Checkpoints(str(outdir / "__ckpt"))
    monkeypatch.setattr(ChunkedRun, "warm_gtf_index", lambda *a, **k: None)

    segments = [(1 + 100000 * i, 100000 * (i + 1)) for i in range(9)]
    monkeypatch.setattr(
        ChunkedRun,
        "run_cut_selection",
        lambda plan, *a: ({"wall_s": 0.0}, selection("cA", 900000, segments)),
    )
    concurrent = {"max": 0, "now": 0}
    import threading

    seen = threading.Lock()
    # Released by the sixth concurrent unit, so the first six prove the pool is as
    # wide as the budget and units seven to nine pass straight through. A Barrier
    # would be wrong here: nine units over a barrier of six resets it and leaves
    # the last three waiting for parties that do not exist.
    saturated = threading.Event()

    def extracting(plan, *a):
        with seen:
            concurrent["now"] += 1
            concurrent["max"] = max(concurrent["max"], concurrent["now"])
            if concurrent["now"] >= 6:
                saturated.set()
        assert saturated.wait(timeout=60), (
            "only {} extraction(s) were ever concurrent at a budget of 6: the pool "
            "was sized from the seeded selection count".format(concurrent["max"])
        )
        with seen:
            concurrent["now"] -= 1
        return (
            {"wall_s": 0.0},
            {
                "strand": None,
                "strand_split_required": True,
                "offset": 0,
                "window_origin": 0,
                "counts": {"alignments_emitted": 0},
                "dropped_read_names": [],
                "files": {"bam": "/dev/null"},
            },
        )

    monkeypatch.setattr(ChunkedRun, "run_extraction", extracting)

    _sel, _cd, chunks, _mk, allocation = ChunkedRun.run_prep_concurrently(
        args, ckpt, str(outdir), {}, _sources_for(args), 0.25
    )
    assert len(chunks) == 9
    assert allocation.unit_workers == 6
    assert concurrent["max"] == 6


# ---------------------------------------------------------------- pool refusal


def test_one_failed_unit_refuses_the_whole_make_chunks_phase(tmp_path, monkeypatch):
    """Named unit, named log, nothing extracted, nothing merged.

    A partial make-chunks result is a run missing whole contigs, which merges
    perfectly happily and produces a smaller answer.
    """

    contigs = {"cA": 5000, "cB": 4000, "cC": 3000}
    genome = write_genome(tmp_path / "g.fa", contigs)
    bam = write_reads(
        tmp_path / "r.bam",
        contigs,
        [("r{}".format(i), c, 100 + 50 * i, 80, 0) for i, c in enumerate(contigs)],
    )
    gtf = write_gtf(
        tmp_path / "a.gtf",
        [(c, "g" + c, "t" + c, "+", 200, 900) for c in contigs],
    )
    outdir = tmp_path / "out"
    (outdir / "logs").mkdir(parents=True)
    args = ChunkedRun.default_args(
        genome_fa=genome,
        bam=bam,
        gtf=gtf,
        strandless_chunks=True,
        cpu_budget=3,
    )
    ckpt = ChunkedRun.Checkpoints(str(outdir / "__ckpt"))
    monkeypatch.setattr(ChunkedRun, "warm_gtf_index", lambda *a, **k: None)

    def one_contig_fails(plan, ckpt_, outdir_, rss):
        if plan["contig"] == "cB":
            raise RuntimeError("selector exploded on cB")
        return {"wall_s": 0.1}, selection(plan["contig"], contigs[plan["contig"]], [(1, contigs[plan["contig"]])])

    extracted = []

    def never(plan, ckpt_, outdir_, rss):
        extracted.append(plan["chunk_id"])
        return {"wall_s": 0.1}, {}

    monkeypatch.setattr(ChunkedRun, "run_cut_selection", one_contig_fails)
    monkeypatch.setattr(ChunkedRun, "run_extraction", never)

    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.run_prep_concurrently(
            args, ckpt, str(outdir), {}, _sources_for(args), 0.25
        )
    message = str(err.value)
    assert "cuts_cB_strandless" in message
    assert "stage2_cuts_cB_strandless.log" in message
    assert "refusing" in message
    # the failure aborts the queue rather than extracting what it can
    assert len(extracted) < 3


def _sources_for(args):
    return ChunkedRun.cut_sources(args, None, "inputs", None)


def test_a_selection_that_never_reports_is_refused_rather_than_dropped(
    tmp_path, monkeypatch
):
    """A hole in the selection list is a contig silently absent from the run."""

    contigs = {"cA": 5000, "cB": 4000}
    genome = write_genome(tmp_path / "g.fa", contigs)
    bam = write_reads(
        tmp_path / "r.bam", contigs, [("r1", "cA", 100, 80, 0), ("r2", "cB", 100, 80, 0)]
    )
    gtf = write_gtf(
        tmp_path / "a.gtf", [(c, "g" + c, "t" + c, "+", 200, 900) for c in contigs]
    )
    outdir = tmp_path / "out"
    (outdir / "logs").mkdir(parents=True)
    args = ChunkedRun.default_args(
        genome_fa=genome, bam=bam, gtf=gtf, strandless_chunks=True, cpu_budget=2
    )
    ckpt = ChunkedRun.Checkpoints(str(outdir / "__ckpt"))
    monkeypatch.setattr(ChunkedRun, "warm_gtf_index", lambda *a, **k: None)

    def silently_skips_one(plan, ckpt_, outdir_, rss):
        if plan["contig"] == "cB":
            # neither raises nor records: the shape a swallowed error would take
            raise _SilentSkip()
        return {"wall_s": 0.1}, selection("cA", 5000, [(1, 5000)])

    monkeypatch.setattr(ChunkedRun, "run_cut_selection", silently_skips_one)
    monkeypatch.setattr(
        ChunkedRun, "run_extraction", lambda *a, **k: ({"wall_s": 0.0}, {})
    )

    with pytest.raises(ChunkedRun.PipelineError) as err:
        ChunkedRun.run_prep_concurrently(
            args, ckpt, str(outdir), {}, _sources_for(args), 0.25, select_only=True
        )
    assert "cuts_cB_strandless" in str(err.value)


class _SilentSkip(Exception):
    pass
