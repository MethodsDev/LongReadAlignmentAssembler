#!/usr/bin/env python3

"""A GTF read back in and re-exported writes each 3'-end attribute exactly once.

InternalPriming, PAS and PAS_offset are imported into _meta and also re-emitted by
to_GTF_format through their accessors (which fall back to _meta). Emitting both
duplicated the keys on every re-exported transcript line, which R's unnest_wider
refuses outright. A rescan that disagrees with the imported value must win, not be
written alongside it.
"""

import re

from Transcript import GTF_contig_to_transcripts


GTF_LINES = [
    'chr1\tLRAA\ttranscript\t100\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; '
    'PolyA "False"; TSS "False"; InternalPriming "False"; PAS "AATAAA"; PAS_offset "-28";',
    'chr1\tLRAA\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
    'chr1\tLRAA\texon\t300\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1";',
]


def _reexport_transcript_line(tmp_path, mutate=None):
    gtf = tmp_path / "in.gtf"
    gtf.write_text("\n".join(GTF_LINES) + "\n")
    (transcript,) = GTF_contig_to_transcripts.parse_GTF_to_Transcripts(str(gtf))["chr1"]
    transcript.set_simple_path(["n1"])
    if mutate is not None:
        mutate(transcript)
    return transcript.to_GTF_format().splitlines()[0]


def _values(line, key):
    return re.findall(r'\b{} "([^"]*)"'.format(key), line)


def test_imported_attributes_written_once(tmp_path):
    line = _reexport_transcript_line(tmp_path)
    assert _values(line, "InternalPriming") == ["False"]
    assert _values(line, "PAS") == ["AATAAA"]
    assert _values(line, "PAS_offset") == ["-28"]


def test_rescan_replaces_imported_value(tmp_path):
    def rescan(transcript):
        transcript.set_likely_internal_primed(True)
        transcript.set_polyA_signal(None, None)

    line = _reexport_transcript_line(tmp_path, rescan)
    assert _values(line, "InternalPriming") == ["True"]
    assert _values(line, "PAS") == ["none"]
    assert _values(line, "PAS_offset") == []
