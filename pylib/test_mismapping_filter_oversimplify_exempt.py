#!/usr/bin/env python3
# encoding: utf-8

"""The alignment-mismapping filter must NOT remove models on oversimplify contigs.

Oversimplify contigs (e.g. chrM) carry reference models forward whatever the reads did,
so every run over such a contig must report the same set. A coverage-driven filter that
removes different models in a sparsely covered cluster-guided input than in a dense one
makes merge_LRAA_GTFs refuse the disagreeing per-input record sets -- observed as a
cluster-guided scattered SC run failing at the final merge because one low-coverage
cluster's chrM lost a 0-read mito tRNA the others kept.

We patch the two detectors so a model on chrM is flagged, and assert run_mismapping_filter
drops it by default but keeps it when chrM is passed as an exempt (oversimplify) contig.
"""

import os
import AlignmentMismappingFilter as AMF


def _write(path, text):
    with open(path, "wt") as fh:
        fh.write(text)


def _make_inputs(tmp_path):
    gtf = str(tmp_path / "in.gtf")
    quant = str(tmp_path / "in.quant.expr")
    # two single-exon models: one on chrM (the flagged artifact), one on chr1 (control)
    _write(gtf,
        'chrM\tLRAA\ttranscript\t100\t400\t.\t+\t.\tgene_id "gM"; transcript_id "MT.art";\n'
        'chrM\tLRAA\texon\t100\t400\t.\t+\t.\tgene_id "gM"; transcript_id "MT.art";\n'
        'chr1\tLRAA\ttranscript\t100\t400\t.\t+\t.\tgene_id "g1"; transcript_id "c1.keep";\n'
        'chr1\tLRAA\texon\t100\t400\t.\t+\t.\tgene_id "g1"; transcript_id "c1.keep";\n')
    _write(quant,
        "gene_id\ttranscript_id\tall_reads\tTPM\n"
        "gM\tMT.art\t0.0\t0.000\n"
        "g1\tc1.keep\t50.0\t500000.000\n")
    return gtf, quant


def _run(tmp_path, monkeypatch, exempt):
    gtf, quant = _make_inputs(tmp_path)
    # flag the chrM model; leave the sequence detector empty
    monkeypatch.setattr(AMF, "_detect_mirror", lambda *a, **k: {"MT.art": ("c1.keep", 0.95, 0.004)})
    monkeypatch.setattr(AMF, "_detect_sequence", lambda *a, **k: {})
    out_gtf = str(tmp_path / f"out.{'exempt' if exempt else 'plain'}.gtf")
    out_quant = str(tmp_path / f"out.{'exempt' if exempt else 'plain'}.quant.expr")
    drop = AMF.run_mismapping_filter(
        gtf_in=gtf, quant_in=quant, genome_fasta="/does/not/matter",
        gtf_out=out_gtf, quant_out=out_quant,
        log_out=str(tmp_path / "log.txt"), workdir=str(tmp_path),
        exempt_contigs=({"chrM"} if exempt else None),
    )
    tids = set()
    for line in open(out_gtf):
        if "\ttranscript\t" in line:
            tids.add(line.split('transcript_id "')[1].split('"')[0])
    return drop, tids


def test_flagged_model_dropped_without_exemption(tmp_path, monkeypatch):
    drop, tids = _run(tmp_path, monkeypatch, exempt=False)
    assert "MT.art" in drop
    assert "MT.art" not in tids          # removed
    assert "c1.keep" in tids             # control untouched


def test_flagged_model_on_oversimplify_contig_is_retained(tmp_path, monkeypatch):
    drop, tids = _run(tmp_path, monkeypatch, exempt=True)
    assert "MT.art" not in drop          # exemption pulled it from the drop-set
    assert "MT.art" in tids              # carried forward
    assert "c1.keep" in tids


def test_exemption_only_spares_the_named_contig(tmp_path, monkeypatch):
    # chr1 is NOT exempt, so a chr1 flag is still dropped even when chrM is exempt.
    gtf, quant = _make_inputs(tmp_path)
    monkeypatch.setattr(AMF, "_detect_mirror", lambda *a, **k: {"c1.keep": ("MT.art", 0.95, 0.004)})
    monkeypatch.setattr(AMF, "_detect_sequence", lambda *a, **k: {})
    out_gtf = str(tmp_path / "out.gtf")
    drop = AMF.run_mismapping_filter(
        gtf_in=gtf, quant_in=quant, genome_fasta="/x",
        gtf_out=out_gtf, quant_out=str(tmp_path / "out.qe"),
        log_out=str(tmp_path / "log.txt"), workdir=str(tmp_path),
        exempt_contigs={"chrM"},
    )
    assert "c1.keep" in drop
