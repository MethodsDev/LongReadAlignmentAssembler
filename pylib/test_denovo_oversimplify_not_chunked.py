#!/usr/bin/env python3

"""De novo --oversimplify refuses a contig that gets cut.

An oversimplified contig is represented by one aggregate model per strand, minted for
the contig rather than discovered. Cut the contig and every chunk mints its own: the
extractor names each mini contig after the real one, so all of them carry the same id,
the merge's collision planner prefixes each with its unit id, and the contig emerges as
several unit-named aggregates. Two consequences, both observed:

  * identity encodes the cut geometry, so sibling cluster runs agree only while their
    cut plans do -- and merge_LRAA_GTFs carries an oversimplified contig forward
    verbatim, refusing inputs that disagree about it;
  * ``chrM_00_plus@g:chrM:+:OVSIMP`` does not match the anchored percent.mt default
    ``^(MT-|mt-|g:(chrM|MT|M):)``, so mitochondrial content silently reads zero.

Supporting it would mean coalescing per-chunk aggregates back into one model per
contig-strand. A contig nobody wants assembled is the opposite of a contig big enough
to need cutting, so it is refused instead -- at planning, before any chunk runs, so the
error names the cause rather than surfacing as a merge that rejects its own inputs.

Ref-guided is deliberately exempt: its models come from the annotation, which splits
across chunks deterministically.
"""

import os
import sys
import types

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if os.path.join(REPO_ROOT, "pylib") not in sys.path:
    sys.path.insert(0, os.path.join(REPO_ROOT, "pylib"))

import ChunkedRun


def _chunks(spec):
    """``{contig: n}`` as the chunk records the guard reads."""
    out = []
    for chrom, count in spec.items():
        for i in range(count):
            out.append({"chrom": chrom, "chunk_id": "{}_{:02d}".format(chrom, i)})
    return out


def _args(oversimplify="chrM", discovery=True, gtf=None):
    return types.SimpleNamespace(
        oversimplify=oversimplify, discovery=discovery, gtf=gtf
    )


def test_a_cut_denovo_oversimplified_contig_is_refused():
    with pytest.raises(ChunkedRun.PipelineError) as excinfo:
        ChunkedRun.refuse_chunked_denovo_oversimplify(
            _args(), _chunks({"chrM": 3, "chr1": 12})
        )

    message = str(excinfo.value)
    assert "chrM" in message
    assert "approx_MB_per_cut" in message, "the error must say how to proceed"
    assert "chr1" not in message, "only the oversimplified contig is at issue"


def test_an_uncut_denovo_oversimplified_contig_is_allowed():
    """The supported shape, and the default one: chrM is 16.5 kb against a 10 Mb cut."""
    ChunkedRun.refuse_chunked_denovo_oversimplify(
        _args(), _chunks({"chrM": 1, "chr1": 12})
    )


def test_a_cut_contig_nobody_oversimplified_is_allowed():
    ChunkedRun.refuse_chunked_denovo_oversimplify(
        _args(oversimplify="chrM"), _chunks({"chr1": 12})
    )


def test_ref_guided_may_cut_an_oversimplified_contig(tmp_path):
    """Its models come from the annotation, so the chunks partition a fixed set."""
    gtf = tmp_path / "annot.gtf"
    gtf.write_text(
        'chrM\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
    )
    ChunkedRun.refuse_chunked_denovo_oversimplify(
        _args(discovery=True, gtf=str(gtf)), _chunks({"chrM": 3})
    )


def test_an_empty_annotation_file_does_not_count_as_ref_guided(tmp_path):
    """``# no gtf records`` is what the workflow hands a chunk with no annotation.

    Reading its presence as "annotated" is the mistake that once let chrM fall through
    to ordinary discovery and come back as 75 spurious models.
    """
    gtf = tmp_path / "empty.gtf"
    gtf.write_text("# no gtf records\n")

    with pytest.raises(ChunkedRun.PipelineError):
        ChunkedRun.refuse_chunked_denovo_oversimplify(
            _args(discovery=True, gtf=str(gtf)), _chunks({"chrM": 2})
        )


def test_an_annotation_silent_about_this_contig_does_not_count_either(tmp_path):
    """Annotated elsewhere is not annotated here: chrM still takes the de novo path."""
    gtf = tmp_path / "annot.gtf"
    gtf.write_text(
        'chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
    )

    with pytest.raises(ChunkedRun.PipelineError):
        ChunkedRun.refuse_chunked_denovo_oversimplify(
            _args(discovery=True, gtf=str(gtf)), _chunks({"chrM": 2})
        )


def test_quant_only_is_not_subject_to_this():
    """No discovery, so no aggregate is minted and nothing can disagree."""
    ChunkedRun.refuse_chunked_denovo_oversimplify(
        _args(discovery=False), _chunks({"chrM": 4})
    )


def test_several_oversimplified_contigs_are_all_named():
    with pytest.raises(ChunkedRun.PipelineError) as excinfo:
        ChunkedRun.refuse_chunked_denovo_oversimplify(
            _args(oversimplify="chrM,chr5"), _chunks({"chrM": 2, "chr5": 9})
        )

    message = str(excinfo.value)
    assert "chrM" in message and "chr5" in message
