#!/usr/bin/env python3

"""Stage 6's quant.tracking merge must sort without materialising the table.

The merge used to append every unit's rows to one Python list and sort that.
MEASURED, that list costs 588 B of Python object per 140 B of data (4.2x), which
put the serial merge at 41.6 GiB on a 63.8M-row PBMC arm -- the peak of the whole
run, in all four measured arms, independent of chunk count. It is now an external
sort: bounded in-memory runs, spilled, then ``heapq.merge``d into the writer.

The prerequisite a plain k-way merge would need DOES NOT HOLD. Per-unit tracking
files are emitted transcript-major (``Quantify.report_quant_results`` walks
transcripts by descending read count; only read names within one multipath are
sorted) and LRAA's shard merge concatenates them unreordered
(``LRAA:_append_without_first_line``). ``test_per_unit_tracking_is_not_presorted``
pins that here, on generated rows that mirror that emission order, so the sort
cannot be quietly dropped later on the strength of a wrong assumption.

Three things are proven:

1. BYTE-IDENTITY against the materialising implementation, kept in this file as
   ``_reference_merge_tracking``, on rows deliberately interleaved across units
   and carrying keys duplicated both within and across units. Both discovery
   modes, and both the spilling and the everything-fits paths.
2. The header-consistency ``PipelineError`` still fires.
3. MEMORY SHAPE, as an exact identity on the reported peak container size rather
   than as an RSS measurement: with the run cap held fixed, quadrupling the row
   count leaves the streaming peak unchanged while the materialising peak
   quadruples.
"""

import gzip
import os
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "pylib"))

import ChunkedRun  # noqa: E402
import Util_funcs  # noqa: E402

# Stage 6 REQUIRES a read-assignment summary from every quant unit
# (``ChunkedRun.merge_read_assignment_summaries``), because every LRAA path now
# writes one -- including the ``run_quant_only`` early returns that quantify
# nothing and report the reads they saw with zero assigned. These tests are about
# the TRACKING merge and care nothing for that table, but a unit directory without
# it is no longer a shape a real run produces, so the fixture would be testing the
# merge against an input it can never receive. Loaded from LRAA and driven through
# its own writer rather than hand-writing the 28-column schema, which would be
# free to drift from what LRAA emits -- the exact disagreement the merge's
# field-name check exists to catch. Same approach as
# ``test_merge_chunk_outputs.py``.
import importlib.machinery  # noqa: E402
import importlib.util  # noqa: E402


def _load_lraa_driver():
    spec = importlib.util.spec_from_loader(
        "lraa_driver_tracking_fixture",
        importlib.machinery.SourceFileLoader(
            "lraa_driver_tracking_fixture", str(REPO / "LRAA")
        ),
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

lraa = _load_lraa_driver()


def _write_read_assignment_summary(prefix, reads_total):
    """The per-unit summary a real quant unit leaves beside its quant.expr.

    Written through LRAA's worker writer and then its own merger, because that is
    how a unit's file comes to hold a ``worker`` row AND its own ``TOTAL``: stage
    6 skips the latter or it counts every unit twice.
    """

    worker = prefix + ".read_assignment.summary.worker.tsv"
    lraa._write_read_assignment_summary(
        worker,
        "chrT",
        "+",
        {"reads_total": reads_total, "reads_kept_genome": reads_total},
    )
    lraa._merge_read_assignment_summary_files(
        [worker], prefix + ".read_assignment.summary.tsv"
    )
    os.remove(worker)


EXPR_HEADER = [
    "gene_id",
    "transcript_id",
    "uniq_reads",
    "all_reads",
    "isoform_fraction",
    "unique_gene_read_fraction",
    "TPM",
    "exons",
    "introns",
    "splice_hash_code",
    "RPM_total_reads",
]
TRACKING_HEADER = [
    "gene_id",
    "transcript_id",
    "transcript_splice_hash_code",
    "num_exons",
    "mp_id",
    "read_name",
    "frac_assigned",
    "read_weight",
]

# Three models shared by every unit, as a quant-only run against one GTF gives
# them: two spliced (so the splice-hash remap has something to remap) and one
# monoexonic (whose hash IS its transcript id and must pass through untouched).
MODELS = [
    ("g1", "t1", "chrT:(+)[[100, 200], [300, 400]]", "chrT:(+)[[201, 299]]"),
    ("g1", "t2", "chrT:(+)[[500, 600], [700, 800]]", "chrT:(+)[[601, 699]]"),
    ("g2", "t3", "chrT:(+)[[900, 950]]", ""),
]


def _write_unit(root, unit_id, offset, tracking_rows):
    """One quant unit's chunk-local artifacts: expr, gzipped tracking, gtf."""

    prefix = str(root / unit_id)

    with open(prefix + ".quant.expr", "wt") as fh:
        print("\t".join(EXPR_HEADER), file=fh)
        for gene, tx, exons, introns in MODELS:
            hash_code = Util_funcs.get_hash_code(introns) if introns else tx
            print(
                "\t".join(
                    [gene, tx, "5", "5.0", "0.500", "1.000", "500000.000",
                     exons, introns, hash_code, "50.000"]
                ),
                file=fh,
            )

    # Gzipped, which is what LRAA has emitted since v0.20.0 and so what the
    # merge reads in production; the plain form is exercised by the sibling
    # whole-contig-frame test.
    with gzip.open(prefix + ".quant.tracking.gz", "wt") as fh:
        print("\t".join(TRACKING_HEADER), file=fh)
        for row in tracking_rows:
            print("\t".join(row), file=fh)

    with open(prefix + ".gtf", "wt") as fh:
        print("# LRAA version comment", file=fh)
        for gene, tx, _exons, introns in MODELS:
            attrs = 'gene_id "{}"; transcript_id "{}";'.format(gene, tx)
            end = "400" if introns else "950"
            print(
                "\t".join(
                    ("chrT", "LRAA", "transcript", "100", end, ".", "+", ".",
                     attrs)
                ),
                file=fh,
            )

    # Distinct per unit and derived from this unit's own rows, so a merged TOTAL
    # that dropped or double-counted a unit shows up as a wrong number rather
    # than as a coincidence of equal values.
    _write_read_assignment_summary(prefix, len(tracking_rows))

    return {"unit_id": unit_id, "offset": offset, "quant_prefix": prefix}


def _unit_tracking_rows(unit_index, n_units, rows_per_unit):
    """Rows shaped like a real per-unit tracking file, not like sorted input.

    Two properties the merge has to survive, both deliberate:

    INTERLEAVED -- read ``k`` lands in unit ``k % n_units``, so every unit owns
    a stride of the merged key space and no unit's rows are contiguous in the
    output. A merge that concatenated units in order would fail on any pair.

    DUPLICATE KEYS -- one read is assigned to two models within a unit, and a
    ``shared`` read carries the SAME (read_name, transcript_id, gene_id) in
    every unit. Those are the rows whose relative order is decided by sort
    stability alone, so they are what byte-identity actually rests on.

    Emission order is transcript-major and reverse-scrambled within it, which is
    what ``report_quant_results`` produces and is not the merge's sort order.
    """

    rows = []
    for j in range(rows_per_unit):
        read = "read{:06d}".format(j * n_units + unit_index)
        gene, tx, _exons, introns = MODELS[(unit_index + j) % len(MODELS)]
        hash_code = (
            Util_funcs.get_hash_code(introns) if introns else tx
        )
        rows.append(
            [gene, tx, hash_code, "2" if introns else "1",
             "MP{}".format(j), read, "1.000000", "1.000000"]
        )
        if j % 5 == 0:
            # same read, second model: a duplicate on read_name alone
            gene2, tx2, _e2, introns2 = MODELS[(unit_index + j + 1) % len(MODELS)]
            rows.append(
                [gene2, tx2,
                 Util_funcs.get_hash_code(introns2) if introns2 else tx2,
                 "2" if introns2 else "1",
                 "MP{}b".format(j), read, "0.500000", "1.000000"]
            )

    # fully duplicated key, within the unit and across every unit
    gene, tx, _exons, introns = MODELS[0]
    for copy in range(2):
        rows.append(
            [gene, tx, Util_funcs.get_hash_code(introns), "2",
             "MPdup{}".format(copy), "shared", "0.250000", "1.000000"]
        )

    # transcript-major, then reversed: the emission order, never the sort order
    rows.sort(key=lambda r: (r[1], r[5]))
    rows.reverse()
    return rows


def _make_units(root, n_units=4, rows_per_unit=10, offset_step=1000):
    os.makedirs(root, exist_ok=True)
    return [
        _write_unit(
            root,
            "c{}".format(i),
            # non-zero so the coordinate translation and therefore the splice
            # hash remap both run; a zero offset would skip them entirely
            # (``coords_already_whole_contig``) and leave the remap empty.
            (i + 1) * offset_step,
            _unit_tracking_rows(i, n_units, rows_per_unit),
        )
        for i in range(n_units)
    ]


def _reference_merge_tracking(units, hash_remap, discovery, track_out):
    """The materialising implementation, verbatim, as the byte-identity oracle.

    Signature-compatible with ``ChunkedRun._merge_tracking_streaming`` so it can
    be swapped in and the SAME surrounding ``merge_and_translate`` -- same expr
    pass, same ``hash_remap`` -- produces both files being compared. Its peak
    resident row count is the whole table, which is the number the streaming
    peak is measured against.
    """

    track_header = None
    track_rows = []
    for unit in units:
        _, header, rows = ChunkedRun.read_tsv(
            ChunkedRun.resolve_tracking(unit["quant_prefix"])
        )
        if track_header is None:
            track_header = header
        elif header != track_header:
            raise ChunkedRun.PipelineError(
                "unit {} quant.tracking header differs from the first "
                "unit's".format(unit["unit_id"])
            )
        col = {name: i for i, name in enumerate(header)}
        for row in rows:
            row = list(row)
            if discovery:
                row[col["gene_id"]] = ChunkedRun._namespace_id(
                    unit["unit_id"], row[col["gene_id"]]
                )
                row[col["transcript_id"]] = ChunkedRun._namespace_id(
                    unit["unit_id"], row[col["transcript_id"]]
                )
            old = row[col["transcript_splice_hash_code"]]
            row[col["transcript_splice_hash_code"]] = hash_remap.get(
                (unit["unit_id"], old), old
            )
            track_rows.append(row)

    tcol = {name: i for i, name in enumerate(track_header)}
    track_rows.sort(
        key=lambda r: (
            r[tcol["read_name"]],
            r[tcol["transcript_id"]],
            r[tcol["gene_id"]],
        )
    )
    with gzip.open(track_out, "wt") as ofh:
        print("\t".join(track_header), file=ofh)
        for row in track_rows:
            print("\t".join(row), file=ofh)

    return track_header, len(track_rows), len(track_rows)


def _merge(tmp_path, name, units, discovery, reference=False):
    outdir = tmp_path / name
    os.makedirs(outdir, exist_ok=True)
    if reference:
        original = ChunkedRun._merge_tracking_streaming
        ChunkedRun._merge_tracking_streaming = _reference_merge_tracking
        try:
            return ChunkedRun.merge_and_translate(
                str(outdir), units, discovery=discovery
            )
        finally:
            ChunkedRun._merge_tracking_streaming = original
    return ChunkedRun.merge_and_translate(str(outdir), units, discovery=discovery)


def _tracking_bytes(result):
    """The merged table DECOMPRESSED: gzip stores an mtime, the table does not."""

    return gzip.decompress(Path(result["quant_tracking"]).read_bytes())


def test_per_unit_tracking_is_not_presorted():
    """The k-way-merge prerequisite is absent, so the sort cannot be dropped.

    Recorded as a test rather than a comment because the whole design turns on
    it: if per-unit files were already in merge-key order the spill machinery
    would be dead weight, and someone will eventually assume they are.
    """

    rows = _unit_tracking_rows(0, 4, 10)
    keys = [(r[5], r[1], r[0]) for r in rows]
    assert keys != sorted(keys)


@pytest.mark.parametrize("discovery", [False, True])
@pytest.mark.parametrize("run_rows", [4, 7, 10_000])
def test_streaming_merge_is_byte_identical(tmp_path, discovery, run_rows):
    """Same units, both implementations, same bytes.

    ``run_rows`` spans the spilling path (4 and 7 rows per run over ~50 rows,
    with 7 leaving a partial final run) and the everything-fits path (10_000,
    where nothing spills and no temp file is created).
    """

    units = _make_units(tmp_path / "units")
    original_cap = ChunkedRun._TRACKING_SORT_RUN_ROWS
    ChunkedRun._TRACKING_SORT_RUN_ROWS = run_rows
    try:
        streamed = _merge(tmp_path, "streamed", units, discovery)
    finally:
        ChunkedRun._TRACKING_SORT_RUN_ROWS = original_cap
    reference = _merge(tmp_path, "reference", units, discovery, reference=True)

    assert _tracking_bytes(streamed) == _tracking_bytes(reference)
    assert streamed["tracking_rows"] == reference["tracking_rows"]
    assert streamed["tracking_rows"] > 0
    # the remap has to have fired, or the identity above proves nothing about it
    assert streamed["splice_hash_codes_recomputed"] > 0


def test_spilled_runs_are_cleaned_up(tmp_path):
    """No ``__tracking_sort.*`` survives, on the path that creates one."""

    units = _make_units(tmp_path / "units")
    original_cap = ChunkedRun._TRACKING_SORT_RUN_ROWS
    ChunkedRun._TRACKING_SORT_RUN_ROWS = 4
    try:
        result = _merge(tmp_path, "streamed", units, False)
    finally:
        ChunkedRun._TRACKING_SORT_RUN_ROWS = original_cap

    merged_dir = Path(result["quant_tracking"]).parent
    assert [p.name for p in merged_dir.glob("__tracking_sort.*")] == []


def test_header_mismatch_still_raises(tmp_path):
    """A unit whose tracking header differs is still refused, not merged."""

    units = _make_units(tmp_path / "units", n_units=3)
    odd = units[1]
    rows = _unit_tracking_rows(1, 3, 10)
    with gzip.open(odd["quant_prefix"] + ".quant.tracking.gz", "wt") as fh:
        print("\t".join(TRACKING_HEADER + ["extra_column"]), file=fh)
        for row in rows:
            print("\t".join(row + ["x"]), file=fh)

    original_cap = ChunkedRun._TRACKING_SORT_RUN_ROWS
    ChunkedRun._TRACKING_SORT_RUN_ROWS = 4
    try:
        with pytest.raises(ChunkedRun.PipelineError, match="quant.tracking header"):
            _merge(tmp_path, "streamed", units, False)
    finally:
        ChunkedRun._TRACKING_SORT_RUN_ROWS = original_cap

    # refused BEFORE anything was written, as the materialising version was, and
    # with no spill left behind
    merged_dir = tmp_path / "streamed" / "merged"
    assert not (merged_dir / "chunked.quant.tracking.gz").exists()
    assert [p.name for p in merged_dir.glob("__tracking_sort.*")] == []


def test_peak_resident_rows_does_not_grow_with_the_table(tmp_path):
    """The memory shape, as an identity on container size, not on RSS.

    Peak resident rows for the external sort is ``max(cap, number_of_runs)`` and
    ``number_of_runs`` is ``ceil(rows / cap)``, so the peak stays AT the cap for
    every table smaller than ``cap**2`` -- 2.5e11 rows at the shipped cap of
    500k, against the 63.8M rows of the largest arm ever measured. Held here by
    quadrupling the table under a fixed cap: the streaming peak is unchanged,
    the materialising peak quadruples with it.
    """

    small = _make_units(tmp_path / "small", n_units=4, rows_per_unit=10)
    large = _make_units(tmp_path / "large", n_units=4, rows_per_unit=40)

    original_cap = ChunkedRun._TRACKING_SORT_RUN_ROWS
    ChunkedRun._TRACKING_SORT_RUN_ROWS = 16
    try:
        small_streamed = _merge(tmp_path, "small_s", small, False)
        large_streamed = _merge(tmp_path, "large_s", large, False)
    finally:
        ChunkedRun._TRACKING_SORT_RUN_ROWS = original_cap
    small_ref = _merge(tmp_path, "small_r", small, False, reference=True)
    large_ref = _merge(tmp_path, "large_r", large, False, reference=True)

    # the tables really did grow ~4x, or the comparison below is vacuous
    assert large_streamed["tracking_rows"] > 3.5 * small_streamed["tracking_rows"]

    # materialising: peak IS the table
    assert small_ref["tracking_merge_peak_resident_rows"] == small_ref["tracking_rows"]
    assert large_ref["tracking_merge_peak_resident_rows"] == large_ref["tracking_rows"]

    # streaming: peak is the cap, on both, and far below either table
    assert small_streamed["tracking_merge_peak_resident_rows"] == 16
    assert large_streamed["tracking_merge_peak_resident_rows"] == 16
    assert 16 < small_streamed["tracking_rows"]


def test_peak_resident_rows_is_measured_not_assumed(tmp_path):
    """The reported peak counts rows actually pulled, one per open run.

    With the cap at 1 every row becomes its own run, so a merge that held one
    row per source reports exactly the run count. This is the degenerate setting
    that would expose a merge secretly buffering more than heapq's one-per-source.
    """

    units = _make_units(tmp_path / "units", n_units=2, rows_per_unit=5)
    original_cap = ChunkedRun._TRACKING_SORT_RUN_ROWS
    ChunkedRun._TRACKING_SORT_RUN_ROWS = 1
    try:
        result = _merge(tmp_path, "streamed", units, False)
    finally:
        ChunkedRun._TRACKING_SORT_RUN_ROWS = original_cap

    assert result["tracking_merge_peak_resident_rows"] == result["tracking_rows"]


def test_comment_lines_are_still_dropped(tmp_path):
    """``iter_tsv`` skips ``#`` lines wherever they sit, as ``read_tsv`` does.

    LRAA writes a version comment, and optionally an XW warning, above the
    tracking header. The merged table has never carried them; a streaming reader
    that stopped skipping them would emit one as a data row.
    """

    units = _make_units(tmp_path / "units", n_units=2, rows_per_unit=6)
    rows = _unit_tracking_rows(0, 2, 6)
    with gzip.open(units[0]["quant_prefix"] + ".quant.tracking.gz", "wt") as fh:
        print("# LRAA version 0.0.0", file=fh)
        print("\t".join(TRACKING_HEADER), file=fh)
        for i, row in enumerate(rows):
            if i == 2:
                print("# a comment in the middle", file=fh)
            print("\t".join(row), file=fh)

    original_cap = ChunkedRun._TRACKING_SORT_RUN_ROWS
    ChunkedRun._TRACKING_SORT_RUN_ROWS = 3
    try:
        streamed = _merge(tmp_path, "streamed", units, False)
    finally:
        ChunkedRun._TRACKING_SORT_RUN_ROWS = original_cap
    reference = _merge(tmp_path, "reference", units, False, reference=True)

    assert _tracking_bytes(streamed) == _tracking_bytes(reference)
    assert b"#" not in _tracking_bytes(streamed)
