#!/usr/bin/env python3

"""Stage 6 must not translate coordinates that are already whole-contig.

``merge_and_translate`` and ``merge_discovery_gtf`` exist to undo the extractor's
rebasing: a chunk of a cut contig sits on a mini contig starting at 1, so its
coordinates need its chunk ``offset`` added back. A chunk that spans its WHOLE
contig has offset 0, and then the shift adds 0 to every coordinate and the splice
hash is recomputed from an unchanged string -- work whose entire output is its
input.

That case is the rule, not the exception: 171 of the 475 chunks of the shipped
HG002 partition already spanned their whole contig, and the one-chunk-per-contig
whole-genome mode makes it true of every chunk. So the skip has to be proven
SAFE rather than plausible, which is what the byte-identity test here does: the
same input through both paths, compared as bytes.

The three properties, in the order they matter:

1. mixed offsets still translate, to the exact values the arithmetic produces;
2. all-zero offsets produce the translating path's bytes exactly;
3. the skip is decided by the offsets -- forced on where an offset is non-zero,
   it corrupts the table, which is why it is a condition and not a rewrite.
"""

import contextlib
import gzip
import importlib.util
import os
import re
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "pylib"))

import ChunkedRun  # noqa: E402
import Util_funcs  # noqa: E402


def _load_lraa():
    """LRAA's own read-assignment summary writers.

    Loaded rather than reimplemented because stage 6 now REQUIRES a summary from
    every unit (``ChunkedRun.merge_read_assignment_summaries``), so these fixtures
    have to produce one -- and a hand-copied 28-column schema would be free to
    drift from the one LRAA writes, which is exactly the disagreement the merge's
    field-name check exists to catch.
    """

    loader = SourceFileLoader("lraa_driver_summary_fixture", str(REPO / "LRAA"))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


lraa = _load_lraa()


def _write_read_assignment_summary(prefix, reads_total):
    """The per-unit summary a real quant unit leaves beside its quant.expr.

    Written through LRAA's worker writer and then its own merger, because that is
    how a unit's file comes to hold a ``worker`` row AND its own ``TOTAL``: stage 6
    has to skip the latter or it counts every unit twice.
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

# chunk-local, as a chunk worker emits them
MULTI_EXONS = "chrT:(+)[[100, 200], [300, 400]]"
MULTI_INTRONS = "chrT:(+)[[201, 299]]"
MONO_EXONS = "chrT:(+)[[600, 700]]"
# LRAA prints splice_hash_code as the digest of the introns string it prints
# alongside (Quantify.report_quant_results), which is what makes the recompute an
# identity when the introns string does not move.
MULTI_HASH = Util_funcs.get_hash_code(MULTI_INTRONS)


def _write_unit(root, unit_id, offset):
    """One quant unit's chunk-local artifacts: expr, tracking, model gtf, and the
    read-assignment summary stage 6 requires of every unit."""

    prefix = str(root / unit_id)
    multi = "t.{}.multi".format(unit_id)
    mono = "t.{}.mono".format(unit_id)
    gene = "g.{}".format(unit_id)

    with open(prefix + ".quant.expr", "wt") as fh:
        # A real unit's quant.expr opens with these; the merge reads the version
        # off them to stamp the merged artifact and refuses a unit carrying none.
        print("# LRAA version test", file=fh)
        print("# LRAA CMD: LRAA --contig chrT", file=fh)
        print("\t".join(EXPR_HEADER), file=fh)
        print(
            "\t".join(
                [gene, multi, "7", "9.0", "0.750", "1.000", "750000.000",
                 MULTI_EXONS, MULTI_INTRONS, MULTI_HASH, "90.000"]
            ),
            file=fh,
        )
        # monoexonic: no introns, and splice_hash_code IS the transcript id
        print(
            "\t".join(
                [gene, mono, "3", "3.0", "0.250", "1.000", "250000.000",
                 MONO_EXONS, "", mono, "30.000"]
            ),
            file=fh,
        )

    with open(prefix + ".quant.tracking", "wt") as fh:
        print("\t".join(TRACKING_HEADER), file=fh)
        print(
            "\t".join(
                [gene, multi, MULTI_HASH, "2", "MP1",
                 "read.{}.a".format(unit_id), "1.000000", "1.000000"]
            ),
            file=fh,
        )
        print(
            "\t".join(
                [gene, mono, mono, "1", "MP2",
                 "read.{}.b".format(unit_id), "1.000000", "1.000000"]
            ),
            file=fh,
        )

    attrs = 'gene_id "{}"; transcript_id "{}";'.format(gene, multi)
    with open(prefix + ".gtf", "wt") as fh:
        print("# LRAA version comment", file=fh)
        for feature, start, end in (
            ("transcript", "100", "400"),
            ("exon", "100", "200"),
            ("exon", "300", "400"),
        ):
            print(
                "\t".join(
                    ("chrT", "LRAA", feature, start, end, ".", "+", ".", attrs)
                ),
                file=fh,
            )

    # Two tracking rows above, so two reads: the merged TOTAL is then predictable
    # rather than an opaque constant.
    _write_read_assignment_summary(prefix, 2)

    return {"unit_id": unit_id, "offset": offset, "quant_prefix": prefix}


def _units(root, offsets):
    os.makedirs(root, exist_ok=True)
    return [
        _write_unit(root, "c{}".format(i), offset)
        for i, offset in enumerate(offsets)
    ]


@contextlib.contextmanager
def _forced_guard(value):
    """Pin ``coords_already_whole_contig`` so a path can be exercised on demand.

    Both merge functions look the name up on the module, so this reaches the
    decision itself without duplicating either function. Forcing True on non-zero
    offsets is exactly the mutation "make the fast path unconditional".
    """

    original = ChunkedRun.coords_already_whole_contig
    ChunkedRun.coords_already_whole_contig = lambda units: value
    try:
        yield
    finally:
        ChunkedRun.coords_already_whole_contig = original


def _merge(outdir, units, discovery, forced=None):
    os.makedirs(outdir, exist_ok=True)
    if forced is None:
        return ChunkedRun.merge_and_translate(
            str(outdir), units, discovery=discovery
        )
    with _forced_guard(forced):
        return ChunkedRun.merge_and_translate(
            str(outdir), units, discovery=discovery
        )


def _expr_rows(result):
    with open(result["quant_expr"], "rt") as fh:
        rows = [
            line.rstrip("\n").split("\t")
            for line in fh
            if not line.startswith("#")  # provenance comments precede the column row
        ]
    return rows[0], rows[1:]


def _artifact_bytes(result, discovery):
    """Every merged artifact as bytes, keyed by name.

    The tracking file is compared DECOMPRESSED: gzip stores an mtime, so two runs
    a second apart differ in the container while agreeing on every byte of the
    table. The gtf's leading comment is held out because it is REQUIRED to differ
    -- it reports which path ran -- and is asserted on its own.
    """

    blobs = {
        "quant.expr": Path(result["quant_expr"]).read_bytes(),
        "tpm_audit": Path(result["tpm_chunk_local_audit"]).read_bytes(),
        "quant.tracking": gzip.decompress(
            Path(result["quant_tracking"]).read_bytes()
        ),
    }
    if discovery:
        with open(result["gtf"], "rt") as fh:
            body = [line for line in fh if not line.startswith("#")]
        blobs["gtf_body"] = "".join(body).encode()
    return blobs


def _gtf_comment(result):
    """Every leading comment of the merged gtf, joined.

    Not ``readline()``: the merged header now opens with ``# LRAA version`` so that
    line 1 answers "which build made this" identically for chunked and unchunked
    outputs, and the merge-provenance line follows it. Callers assert with ``in``,
    so joining keeps them about content rather than position.
    """

    with open(result["gtf"], "rt") as fh:
        return "\n".join(
            line.rstrip("\n") for line in fh if line.startswith("#")
        )


# ------------------------------------------------------- 1. mixed offsets shift


def test_mixed_offsets_translate_to_exactly_the_values_the_shift_produces(tmp_path):
    """A cut contig's pieces still get their offset back, in every coordinate field.

    Asserted against the literal results of ``+ 5000``, not against a second
    implementation of the shift, so a test that agrees with a broken merge cannot
    be written by copying the merge. Unit c0 is offset 0 and unit c1 is offset
    5000 in the same run, which is what a partially cut genome looks like.
    """

    units = _units(tmp_path / "units", (0, 5000))
    result = _merge(tmp_path / "run", units, discovery=True)

    header, rows = _expr_rows(result)
    col = {name: i for i, name in enumerate(header)}
    by_id = {row[col["transcript_id"]]: row for row in rows}

    # No unit prefix: each of these ids names one model in one unit, so nothing
    # collides and the merge carries the id through as it found it. The prefix
    # appears only where two units name different models the same.
    c0 = by_id["t.c0.multi"]
    c1 = by_id["t.c1.multi"]

    # offset 0: unmoved
    assert c0[col["exons"]] == "chrT:(+)[[100, 200], [300, 400]]"
    assert c0[col["introns"]] == "chrT:(+)[[201, 299]]"
    assert c0[col["splice_hash_code"]] == MULTI_HASH

    # offset 5000: every integer inside the brackets, and nothing outside them
    assert c1[col["exons"]] == "chrT:(+)[[5100, 5200], [5300, 5400]]"
    assert c1[col["introns"]] == "chrT:(+)[[5201, 5299]]"
    shifted_hash = Util_funcs.get_hash_code("chrT:(+)[[5201, 5299]]")
    assert c1[col["splice_hash_code"]] == shifted_hash
    assert shifted_hash != MULTI_HASH

    # the monoexonic row has no introns to shift and keeps its id as its hash
    mono = by_id["t.c1.mono"]
    assert mono[col["introns"]] == ""
    assert mono[col["splice_hash_code"]] == "t.c1.mono"
    assert mono[col["exons"]] == "chrT:(+)[[5600, 5700]]"

    # tracking follows the expr recompute, per unit
    with gzip.open(result["quant_tracking"], "rt") as fh:
        track = [line.rstrip("\n").split("\t") for line in fh]
    tcol = {name: i for i, name in enumerate(track[0])}
    hashes = {
        row[tcol["transcript_id"]]: row[tcol["transcript_splice_hash_code"]]
        for row in track[1:]
    }
    assert hashes["t.c0.multi"] == MULTI_HASH
    assert hashes["t.c1.multi"] == shifted_hash

    # gtf columns 4 and 5, both, every row
    with open(result["gtf"], "rt") as fh:
        gtf = [line.rstrip("\n").split("\t") for line in fh if line[0] != "#"]
    assert [(row[3], row[4]) for row in gtf] == [
        ("100", "400"),
        ("100", "200"),
        ("300", "400"),
        ("5100", "5400"),
        ("5100", "5200"),
        ("5300", "5400"),
    ]

    assert result["coordinates_translated"] is True
    assert "coordinates translated to the whole-contig frame" in _gtf_comment(result)


# --------------------------------------- 2. all-zero offsets: byte-for-byte same


@pytest.mark.parametrize("discovery", (False, True))
def test_all_zero_offsets_emit_the_translating_paths_bytes(tmp_path, discovery):
    """The equivalence the skip rests on, asserted as bytes rather than argued.

    One input, two paths: the skip (all offsets 0, which is what one chunk per
    contig produces for every chunk) and the shift it replaces, reached by pinning
    the decision. Every merged artifact must agree byte for byte -- the table, the
    TPM audit, the tracking rows, and the model gtf's records.
    """

    units = _units(tmp_path / "units", (0, 0))

    skipped = _merge(tmp_path / "skip", units, discovery=discovery)
    shifted = _merge(tmp_path / "shift", units, discovery=discovery, forced=False)

    assert skipped["coordinates_translated"] is False
    assert shifted["coordinates_translated"] is True

    skipped_blobs = _artifact_bytes(skipped, discovery)
    shifted_blobs = _artifact_bytes(shifted, discovery)
    assert sorted(skipped_blobs) == sorted(shifted_blobs)
    for name in skipped_blobs:
        assert skipped_blobs[name] == shifted_blobs[name], name

    # and the tables are not vacuously equal
    _, rows = _expr_rows(skipped)
    assert len(rows) == 4

    # the shift did run in the control, on rows that could have moved
    assert shifted["splice_hash_codes_recomputed"] == 2
    # nothing was recomputed on the skipping path, and the report says so
    assert skipped["splice_hash_codes_recomputed"] == 0

    if discovery:
        comment = _gtf_comment(skipped)
        assert "and none were translated" in comment
        assert "coordinates translated to the whole-contig frame" not in comment
        assert (
            "coordinates translated to the whole-contig frame"
            in _gtf_comment(shifted)
        )


# ------------------------------------------ 3. the skip is a condition, not a rule


def test_the_skip_is_decided_by_the_offsets_and_would_corrupt_if_unconditional(
    tmp_path,
):
    """Fails if the fast path stops being conditional, in either way it could.

    Delete the condition and the mixed-offset expectations above break; make the
    DECISION always true and this breaks, because the two runs below then produce
    the same bytes. Non-zero offsets skipped is not a subtle degradation: the
    piece keeps its chunk-local coordinates under an absolute contig name, so two
    pieces of one contig land on top of each other.
    """

    units = _units(tmp_path / "units", (0, 5000))

    correct = _merge(tmp_path / "correct", units, discovery=True)
    unconditional = _merge(tmp_path / "bad", units, discovery=True, forced=True)

    correct_blobs = _artifact_bytes(correct, discovery=True)
    bad_blobs = _artifact_bytes(unconditional, discovery=True)
    assert bad_blobs["quant.expr"] != correct_blobs["quant.expr"]
    assert bad_blobs["gtf_body"] != correct_blobs["gtf_body"]

    # what "corrupt" means here: both units' models on the same coordinates
    header, rows = _expr_rows(unconditional)
    col = {name: i for i, name in enumerate(header)}
    exons = [row[col["exons"]] for row in rows if "multi" in row[col["transcript_id"]]]
    assert exons == [MULTI_EXONS, MULTI_EXONS]

    # the decision itself, read off the units and nothing else
    assert ChunkedRun.coords_already_whole_contig(units) is False
    assert ChunkedRun.coords_already_whole_contig(
        [{"offset": 0}, {"offset": 0}]
    ) is True
    assert ChunkedRun.coords_already_whole_contig([{"offset": 1}]) is False
    # a single cut piece is enough to put the whole merge back on the shift
    assert ChunkedRun.coords_already_whole_contig(
        [{"offset": 0}, {"offset": 0}, {"offset": 12345}]
    ) is False


def test_nothing_but_the_units_can_decide_the_skip():
    """It is derived, so there must be nothing to set.

    A knob here would let a caller ask for coordinates to be left chunk-local on a
    cut contig -- the corruption test 3 demonstrates -- and would let a resumed run
    disagree with the run that wrote the units. Held on the two surfaces a knob
    could appear on: the merge's own parameters, and the global config the CLI
    writes into.
    """

    import inspect

    import LRAA_Globals

    assert list(
        inspect.signature(ChunkedRun.merge_and_translate).parameters
    ) == ["outdir", "units", "discovery"]
    assert list(
        inspect.signature(ChunkedRun.merge_discovery_gtf).parameters
    ) == ["merged_dir", "units"]
    assert list(
        inspect.signature(ChunkedRun.coords_already_whole_contig).parameters
    ) == ["units"]

    assert [
        key
        for key in LRAA_Globals.config
        if re.search(r"translat|whole_contig_frame", key)
    ] == []
