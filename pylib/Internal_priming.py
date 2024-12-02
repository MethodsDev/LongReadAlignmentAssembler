import sys, os, re
from Transcript import Transcript


import logging

logger = logging.getLogger(__name__)


def filter_internally_primed_transcripts(
    transcripts, contig_seq_str, restrict_filter_to_monoexonic=True
):

    retained_transcripts = list()

    for transcript in transcripts:

        if restrict_filter_to_monoexonic and not transcript.is_monoexonic():
            retained_transcripts.append(transcript)
            continue

        # evaluate whether transcript looks internally primed.
        transcript_lend, transcript_rend = transcript.get_coords()
        strand = transcript.get_orient()

        if _looks_internally_primed(
            transcript_lend, transcript_rend, strand, contig_seq_str
        ):
            # gbye
            pass
        else:
            # keep
            retained_transcripts.append(transcript)

    return retained_transcripts


def _looks_internally_primed(
    transcript_lend, transcript_rend, strand, contig_seq_str, check_length=20
):

    if strand not in {"+", "-"}:
        raise ValueError("Strand must be '+' or '-'")

    target_base = "A" if strand == "+" else "T"
    target_polyA_motif = target_base * 6

    contig_length = len(contig_seq_str)

    if strand == "+":
        start = transcript_rend + 1
    else:
        start = transcript_lend - check_length - 1

    end = start + check_length

    # ensure coordinates within contig bounds
    start = max(1, start)
    end = min(end, contig_length)

    extracted_sequence = contig_seq_str[start - 1 : end].upper()

    has_flanking_polyA = (
        extracted_sequence.count(target_base) >= 12
        or target_polyA_motif in extracted_sequence
    )

    return has_flanking_polyA
