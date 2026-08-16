#!/usr/bin/env python3

"""Tests for the tabix-indexed GTF path in extract_contig_region_inputs.py.

Extraction used to read the whole GTF once per chunk, so its cost scaled with
(chunk count x GTF size) rather than with genome size. The indexed path replaces
that scan with a region fetch, and these tests hold it to the only contract that
makes the substitution legitimate: for every gene the region needs, the indexed
answer is IDENTICAL to the answer the full scan gives, span and lines alike.

A truncated gene span is the failure that matters. Cut admissibility is decided
by comparing a gene's span against the boundary, so a span short by one exon does
not raise -- it silently declares an inadmissible cut admissible, and the locus
then falls in the gap between two chunks and vanishes from the run.
"""

import importlib.util
import os
from importlib.machinery import SourceFileLoader
from pathlib import Path

import pytest


def _load_extractor():
    path = (
        Path(__file__).resolve().parents[1]
        / "util"
        / "misc"
        / "extract_contig_region_inputs.py"
    )
    loader = SourceFileLoader("extract_contig_region_inputs_tabix_under_test", str(path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


extractor = _load_extractor()

CONTIG = "chrT"


def _row(feature, lend, rend, strand, attrs, contig=CONTIG):
    return "\t".join(
        (contig, "test", feature, str(lend), str(rend), ".", strand, ".", attrs)
    )


def _transcript_rows(gene_id, transcript_id, strand, exons, with_gene_row=True):
    """Rows for one transcript. ``with_gene_row`` controls whether a locus-spanning
    ``gene`` feature is present, which is what decides whether a naive region fetch
    could recover the gene's full span on its own."""

    lend = min(e[0] for e in exons)
    rend = max(e[1] for e in exons)
    attrs = 'gene_id "{}";'.format(gene_id)
    tx_attrs = 'gene_id "{}"; transcript_id "{}";'.format(gene_id, transcript_id)
    rows = []
    if with_gene_row:
        rows.append(_row("gene", lend, rend, strand, attrs))
    rows.append(_row("transcript", lend, rend, strand, tx_attrs))
    for exon_lend, exon_rend in exons:
        rows.append(_row("exon", exon_lend, exon_rend, strand, tx_attrs))
    return rows


def _write_gtf(path, rows):
    with open(path, "wt") as ofh:
        for row in rows:
            print(row, file=ofh)
    return str(path)


def _snapshot(annotation):
    """The state a consumer can observe, in a form two loaders must agree on."""

    return {
        gene_id: (
            gene.strand,
            gene.lend,
            gene.rend,
            sorted(gene.transcript_ids),
            sorted(line.rstrip("\n") for line in gene.lines),
        )
        for gene_id, gene in annotation.genes.items()
    }


def _transcript_snapshot(annotation):
    return {
        tid: (t.gene_id, t.strand, t.lend, t.rend, sorted(t.exons))
        for tid, t in annotation.transcripts.items()
    }


@pytest.fixture
def plain_gtf(tmp_path):
    """Three well-separated loci on + and one on -."""

    rows = []
    rows += _transcript_rows("gA", "gA.1", "+", [(1000, 1200), (1500, 1800)])
    rows += _transcript_rows("gB", "gB.1", "+", [(20000, 20500)])
    rows += _transcript_rows("gC", "gC.1", "-", [(40000, 40300), (41000, 41500)])
    return _write_gtf(tmp_path / "annot.gtf", rows)


def test_indexed_region_matches_full_scan_for_every_gene_it_returns(plain_gtf):
    """The substitution is only legitimate if the answers coincide exactly."""

    full = extractor.load_gtf(plain_gtf, CONTIG, "+")
    full_snapshot = _snapshot(full)

    for lend, rend in ((1, 50000), (900, 2000), (19000, 21000), (1500, 20100)):
        region = extractor.load_gtf_for_region(plain_gtf, CONTIG, "+", lend, rend)
        region_snapshot = _snapshot(region)
        assert region_snapshot, "region {}-{} returned nothing".format(lend, rend)
        for gene_id, observed in region_snapshot.items():
            assert observed == full_snapshot[gene_id]

        # every gene the region actually overlaps must be present, or a locus
        # straddling a boundary could go unnoticed
        for gene_id, (_, g_lend, g_rend, _, _) in full_snapshot.items():
            if g_lend <= rend and g_rend >= lend:
                assert gene_id in region_snapshot


def test_strand_filter_matches_full_scan(plain_gtf):
    indexed = extractor.load_gtf_for_region(plain_gtf, CONTIG, "-", 1, 50000)
    full = extractor.load_gtf(plain_gtf, CONTIG, "-")
    assert _snapshot(indexed) == _snapshot(full)
    assert set(indexed.genes) == {"gC"}


def test_transcript_models_match_full_scan(plain_gtf):
    indexed = extractor.load_gtf_for_region(plain_gtf, CONTIG, "+", 1, 50000)
    full = extractor.load_gtf(plain_gtf, CONTIG, "+")
    assert _transcript_snapshot(indexed) == _transcript_snapshot(full)


def test_gene_span_is_complete_when_transcripts_do_not_overlap(tmp_path):
    """The case widening exists for, and the one a region fetch alone gets wrong.

    One gene, two transcripts, a 39 kb gap between them and NO locus-spanning
    ``gene`` row. A fetch of the region around the far transcript returns only
    that transcript, making the gene look like it starts at 40000 when it starts
    at 1000 -- so a cut at 30000 would read as clear of the locus when it splits
    it. Widening by the longest gene span is what prevents that, and removing the
    widening turns this assertion red.
    """

    rows = _transcript_rows(
        "wide", "wide.1", "+", [(1000, 1200)], with_gene_row=False
    ) + _transcript_rows("wide", "wide.2", "+", [(40000, 40200)], with_gene_row=False)
    gtf = _write_gtf(tmp_path / "disjoint.gtf", rows)

    full = extractor.load_gtf(gtf, CONTIG, "+")
    assert full.genes["wide"].lend == 1000
    assert full.genes["wide"].rend == 40200

    indexed = extractor.load_gtf_for_region(gtf, CONTIG, "+", 39000, 41000)
    assert indexed.genes["wide"].lend == 1000
    assert indexed.genes["wide"].rend == 40200
    assert _snapshot(indexed)["wide"] == _snapshot(full)["wide"]


def test_max_gene_span_is_the_longest_locus(tmp_path):
    rows = _transcript_rows("short", "short.1", "+", [(100, 200)]) + _transcript_rows(
        "long", "long.1", "+", [(1000, 5000)]
    )
    gtf = _write_gtf(tmp_path / "spans.gtf", rows)
    indexed = extractor.ensure_gtf_index(gtf)
    assert indexed is not None
    assert indexed.max_gene_span == 5000 - 1000 + 1


def test_index_is_reused_rather_than_rebuilt(plain_gtf):
    first = extractor.ensure_gtf_index(plain_gtf)
    assert first is not None
    stat_before = os.stat(first.path)

    second = extractor.ensure_gtf_index(plain_gtf)
    assert second == first
    stat_after = os.stat(second.path)
    assert stat_after.st_mtime_ns == stat_before.st_mtime_ns
    assert stat_after.st_ino == stat_before.st_ino


def test_index_is_rebuilt_when_the_gtf_changes(tmp_path):
    gtf = _write_gtf(
        tmp_path / "mutating.gtf", _transcript_rows("g1", "g1.1", "+", [(100, 200)])
    )
    first = extractor.ensure_gtf_index(gtf)
    assert first.max_gene_span == 101

    _write_gtf(
        tmp_path / "mutating.gtf",
        _transcript_rows("g1", "g1.1", "+", [(100, 200)])
        + _transcript_rows("g2", "g2.1", "+", [(1000, 9000)]),
    )
    second = extractor.ensure_gtf_index(gtf)
    assert second.max_gene_span == 8001

    reloaded = extractor.load_gtf_for_region(gtf, CONTIG, "+", 1, 10000)
    assert set(reloaded.genes) == {"g1", "g2"}


def test_the_source_gtf_is_never_modified(plain_gtf):
    with open(plain_gtf, "rb") as fh:
        before = fh.read()
    extractor.ensure_gtf_index(plain_gtf)
    with open(plain_gtf, "rb") as fh:
        assert fh.read() == before


def test_index_lands_in_the_cache_dir_when_the_gtf_dir_is_unwritable(
    tmp_path, monkeypatch
):
    """An unwritable reference directory must not disable indexing."""

    gtf = _write_gtf(
        tmp_path / "ro.gtf", _transcript_rows("g1", "g1.1", "+", [(100, 200)])
    )
    cache = tmp_path / "cache"
    real_build = extractor._build_gtf_index

    def refuse_beside_the_gtf(gtf_filename, gz_path):
        if not str(gz_path).startswith(str(cache)):
            raise OSError("read-only file system")
        return real_build(gtf_filename, gz_path)

    monkeypatch.setattr(extractor, "_build_gtf_index", refuse_beside_the_gtf)

    indexed = extractor.ensure_gtf_index(gtf, cache_dir=str(cache))
    assert indexed is not None
    assert indexed.path.startswith(str(cache))


def test_falls_back_to_a_full_scan_when_no_index_can_be_built(plain_gtf, monkeypatch):
    """A GTF that cannot be indexed must still produce the right annotation.

    Losing the index costs time. Refusing the run, or quietly returning a
    truncated annotation, would cost correctness.
    """

    monkeypatch.setattr(extractor, "ensure_gtf_index", lambda *a, **k: None)
    fallback = extractor.load_gtf_for_region(plain_gtf, CONTIG, "+", 900, 2000)
    assert _snapshot(fallback) == _snapshot(extractor.load_gtf(plain_gtf, CONTIG, "+"))


def test_whole_contig_load_matches_full_scan(plain_gtf):
    for strand in ("+", "-", ""):
        indexed = extractor.load_gtf_for_contig(plain_gtf, CONTIG, strand)
        full = extractor.load_gtf(plain_gtf, CONTIG, strand)
        assert _snapshot(indexed) == _snapshot(full)
        assert _transcript_snapshot(indexed) == _transcript_snapshot(full)


def test_unknown_contig_yields_an_empty_annotation(plain_gtf):
    empty = extractor.load_gtf_for_region(plain_gtf, "chrAbsent", "+", 1, 1000)
    assert empty.genes == {}
    assert empty.transcripts == {}


def test_malformed_line_is_rejected_by_both_loaders(tmp_path):
    """The indexed path must not become the lenient one."""

    rows = _transcript_rows("g1", "g1.1", "+", [(100, 200)])
    rows.append(
        _row("exon", 300, 400, "+", 'note "no ids here";')
    )
    gtf = _write_gtf(tmp_path / "bad.gtf", rows)

    with pytest.raises(extractor.ExtractionError) as scan_err:
        extractor.load_gtf(gtf, CONTIG, "+")
    assert "neither gene_id nor transcript_id" in str(scan_err.value)

    with pytest.raises(extractor.ExtractionError) as indexed_err:
        extractor.load_gtf_for_region(gtf, CONTIG, "+", 1, 1000)
    assert "neither gene_id nor transcript_id" in str(indexed_err.value)
