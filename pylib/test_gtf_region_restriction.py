#!/usr/bin/env python3

"""Region restriction must not invent isoforms.

The GTF reader filters by containment, and containment is a property of a whole
transcript.  Applied per feature row it silently TRUNCATES: the rows inside the
window are parsed and the rest are skipped, so a transcript straddling a bound
comes back with fewer exons than the annotation holds, indistinguishable
downstream from a real short isoform.  An exon-only GTF does this with no row
ever crossing a boundary.

Reachable from any --region run: LRAA passes these bounds positionally at
:2255-2260 and again on the worker path at :2810-2815.  The chunked pipeline is
safe for a different reason -- same-strand gene spans block absolutely and only
the gaps between islands are admissible cuts -- so these tests hold the reader
itself, which has no such guarantee.
"""

import logging
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "pylib") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "pylib"))

from Transcript import GTF_contig_to_transcripts


WINDOW = (1000, 2000)


def _exon(transcript_id, lend, rend, strand="+"):
    return (
        f"chr1\ttest\texon\t{lend}\t{rend}\t.\t{strand}\t.\t"
        f'gene_id "g_{transcript_id}"; transcript_id "{transcript_id}";'
    )


def _write_gtf(tmp_path, rows, name="in.gtf"):
    path = tmp_path / name
    path.write_text("\n".join(rows) + "\n")
    return str(path)


def _parse(path, window=WINDOW):
    contig_to_transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(
        path, "chr1", "+", window[0], window[1]
    )
    return {
        t.get_transcript_id(): t.get_coords()
        for t in contig_to_transcripts.get("chr1", [])
    }


def test_a_straddling_transcript_is_dropped_not_truncated(tmp_path):
    """The defect: an exon-only GTF where no single row crosses the boundary.

    t_strad has one exon inside the window and one outside.  Filtering rows keeps
    the inside exon and discards the other, yielding a single-exon transcript at
    (1500, 1600) that does not exist in the annotation.
    """

    path = _write_gtf(
        tmp_path,
        [
            _exon("t_in", 1100, 1200),
            _exon("t_strad", 1500, 1600),
            _exon("t_strad", 2500, 2600),
        ],
    )

    kept = _parse(path)

    assert "t_strad" not in kept, (
        "a transcript only partly inside the window must be dropped whole; keeping "
        "its contained exons invents an isoform the annotation does not contain"
    )
    assert kept == {"t_in": (1100, 1200)}


def test_a_contained_transcript_keeps_every_exon(tmp_path):
    path = _write_gtf(
        tmp_path,
        [_exon("t_multi", 1100, 1200), _exon("t_multi", 1700, 1800)],
    )

    assert _parse(path) == {"t_multi": (1100, 1800)}


def test_a_transcript_entirely_outside_is_filtered_without_a_warning(tmp_path, caplog):
    """Lying outside the window is correct filtering, not a loss to report.

    Only a partly-contained transcript is a dropped annotation; warning about
    every out-of-window transcript would make the signal useless on a whole-genome
    GTF restricted to one locus.
    """

    path = _write_gtf(
        tmp_path,
        [_exon("t_in", 1100, 1200), _exon("t_far", 50000, 50100)],
    )

    with caplog.at_level(logging.WARNING, logger="Transcript.Transcript"):
        kept = _parse(path)

    assert kept == {"t_in": (1100, 1200)}
    assert not [r for r in caplog.records if "dropped" in r.getMessage()]


def test_dropping_a_straddling_transcript_is_reported(tmp_path, caplog):
    """Nothing downstream can tell the locus was ever in the GTF, so say so."""

    path = _write_gtf(
        tmp_path,
        [_exon("t_strad", 1500, 1600), _exon("t_strad", 2500, 2600)],
    )

    with caplog.at_level(logging.WARNING):
        assert _parse(path) == {}

    messages = [r.getMessage() for r in caplog.records]
    assert any("t_strad" in m and "dropped" in m for m in messages), messages


def test_no_region_restriction_keeps_everything(tmp_path):
    """The filter must stay inert when no bounds are given."""

    path = _write_gtf(
        tmp_path,
        [_exon("t_a", 1100, 1200), _exon("t_b", 50000, 50100)],
    )

    contig_to_transcripts = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(path)
    kept = {
        t.get_transcript_id(): t.get_coords()
        for t in contig_to_transcripts.get("chr1", [])
    }
    assert kept == {"t_a": (1100, 1200), "t_b": (50000, 50100)}
