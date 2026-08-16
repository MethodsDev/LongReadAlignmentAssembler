#!/usr/bin/env python3

"""Tests for the tabix-indexed GTF path in extract_contig_region_inputs.py.

Extraction used to read the whole GTF once per chunk, so its cost scaled with
(chunk count x GTF size) rather than with genome size. The indexed path replaces
that scan with a region fetch, and these tests hold it to the contract that makes
the substitution legitimate: the annotation a chunk sees is the same one the full
scan would have given it, over every coordinate the chunk's consumers consult.

"Every coordinate the consumers consult" is wider than the chunk. Cut selection
guarantees that a gene overlapping a chunk lies wholly inside it, so emission
needs nothing beyond the interval -- but ``admissibility_offenders`` tests each
boundary against ``cut +/- margin``, and a gene lying entirely outside the chunk,
within that margin, still blocks the cut. A fetch of the interval alone does not
fail on such a gene; it stops being able to see it, and reports zero offenders
where the scan reports one.
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
CONTIG_LENGTH = 100000


def _row(feature, lend, rend, strand, attrs, contig=CONTIG):
    return "\t".join(
        (contig, "test", feature, str(lend), str(rend), ".", strand, ".", attrs)
    )


def _transcript_rows(gene_id, transcript_id, strand, exons):
    lend = min(e[0] for e in exons)
    rend = max(e[1] for e in exons)
    attrs = 'gene_id "{}";'.format(gene_id)
    tx_attrs = 'gene_id "{}"; transcript_id "{}";'.format(gene_id, transcript_id)
    rows = [
        _row("gene", lend, rend, strand, attrs),
        _row("transcript", lend, rend, strand, tx_attrs),
    ]
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


def test_a_gene_outside_the_chunk_but_within_the_margin_still_blocks_the_cut(tmp_path):
    """The reason the fetch is widened at all, and the whole of that reason.

    ``near`` sits 50 bp past the chunk's right edge, so it overlaps nothing the
    chunk contains -- and it blocks the boundary anyway, because a cut must clear
    every locus by ``margin``. Fetch only the interval and the annotation has no
    record of it, so the check returns clean on a region the scan refuses.

    Set the widening to 0 and this goes red; it is what pins the margin.
    """

    rows = _transcript_rows("inside", "inside.1", "+", [(10000, 12000)])
    rows += _transcript_rows("near", "near.1", "+", [(20050, 20400)])
    gtf = _write_gtf(tmp_path / "near.gtf", rows)

    region = extractor.Region(CONTIG, "+", 5000, 20000)
    margin = extractor.DEFAULT_MARGIN
    index_path = extractor.ensure_gtf_index(gtf)

    scanned = extractor.load_gtf(gtf, CONTIG, "+")
    widened = extractor.load_gtf_region(
        index_path, gtf, CONTIG, "+", region.lend, region.rend, margin
    )
    unwidened = extractor.load_gtf_region(
        index_path, gtf, CONTIG, "+", region.lend, region.rend, 0
    )

    def offenders(annotation):
        return extractor.admissibility_offenders(
            annotation, region, CONTIG_LENGTH, margin
        )

    assert len(offenders(scanned)) == 1
    assert len(offenders(widened)) == len(offenders(scanned))
    assert offenders(widened) == offenders(scanned)
    # and the unwidened fetch is the failure this guards against
    assert offenders(unwidened) == []


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

        # every gene the region overlaps must be present and complete, since a
        # gene overlapping a chunk lies wholly inside it
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


def test_index_is_reused_rather_than_rebuilt(plain_gtf):
    first = extractor.ensure_gtf_index(plain_gtf)
    assert first is not None
    stat_before = os.stat(first)

    second = extractor.ensure_gtf_index(plain_gtf)
    assert second == first
    stat_after = os.stat(second)
    assert stat_after.st_mtime_ns == stat_before.st_mtime_ns
    assert stat_after.st_ino == stat_before.st_ino


def test_index_is_rebuilt_when_the_gtf_changes(tmp_path):
    """Staleness is size and mtime; a content hash would cost more than the scan."""

    gtf = _write_gtf(
        tmp_path / "mutating.gtf", _transcript_rows("g1", "g1.1", "+", [(100, 200)])
    )
    extractor.ensure_gtf_index(gtf)
    assert set(extractor.load_gtf_for_region(gtf, CONTIG, "+", 1, 10000).genes) == {"g1"}

    _write_gtf(
        tmp_path / "mutating.gtf",
        _transcript_rows("g1", "g1.1", "+", [(100, 200)])
        + _transcript_rows("g2", "g2.1", "+", [(1000, 9000)]),
    )
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

    index_path = extractor.ensure_gtf_index(gtf, cache_dir=str(cache))
    assert index_path is not None
    assert index_path.startswith(str(cache))


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
    rows.append(_row("exon", 300, 400, "+", 'note "no ids here";'))
    gtf = _write_gtf(tmp_path / "bad.gtf", rows)

    with pytest.raises(extractor.ExtractionError) as scan_err:
        extractor.load_gtf(gtf, CONTIG, "+")
    assert "neither gene_id nor transcript_id" in str(scan_err.value)

    with pytest.raises(extractor.ExtractionError) as indexed_err:
        extractor.load_gtf_for_region(gtf, CONTIG, "+", 1, 1000)
    assert "neither gene_id nor transcript_id" in str(indexed_err.value)
