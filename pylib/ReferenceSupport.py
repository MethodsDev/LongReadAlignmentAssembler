#!/usr/bin/env python3

"""Reference identity of models in ref-guided runs.

Three things the filters need to know about a model besides its reads, all written to the
output GTF so a run guided by that GTF (a cluster run guided by the pooled init GTF) and
the cross-cluster merge can carry them forward:

- which guide models stand for reference transcripts (validate_guide_transcripts): a
  reference annotation, or an LRAA-produced GTF whose rows carry their own call;
- how many reads of the FULL bam carry exactly each reference intron chain
  (assign_ref_FSM_reads): ref_splice_pattern_reprieve's evidence, taken before depth
  normalization so a real full-length read cannot be thinned away;
- whether a single-exon model matches a reference single-exon transcript
  (assign_ref_supported_monoexonic): ref_supported_monoexonic's evidence.
"""

import bisect
import gzip
import logging
import re
import time
from collections import defaultdict

import LRAA_Globals
from Bam_alignment_extractor import Bam_alignment_extractor
from Pretty_alignment import Pretty_alignment

logger = logging.getLogger(__name__)


def reference_rules_active():
    """Whether either reference-identity rule is on (both are by default)."""
    return bool(
        LRAA_Globals.config["ref_splice_pattern_reprieve"]
        or LRAA_Globals.config["ref_supported_monoexonic"]
    )


def validate_guide_transcripts(guide_transcripts, gtf_label=""):
    """Classify a guide's models and refuse an LRAA-produced guide without its flags.

    A guide row whose GTF column 2 is "LRAA" was built from reads by an earlier run, not
    annotated, so it stands for a reference transcript only when that run said so:
    ref_splice_pattern on multi-exon rows, ref_supported on single-exon rows. A row
    without its flag (an LRAA GTF from before v0.50.0) cannot say, and treating it as a
    reference would turn every model of that run into one -- novel chains included -- so
    it is refused unless allow_unflagged_LRAA_guide is set, which restores exactly that.

    Used by the merge on its parsed inputs; LRAA itself checks its --gtf file once at
    startup (validate_guide_gtf_file). No-op when both rules are off. Returns
    (num_reference_annotation, num_LRAA, num_LRAA_unflagged), or None when there was
    nothing to check.
    """
    if not guide_transcripts or not reference_rules_active():
        return None

    num_ref = num_LRAA = num_unflagged = 0
    for transcript in guide_transcripts:
        if not transcript.is_LRAA_guide_model():
            num_ref += 1
            continue
        num_LRAA += 1
        if _lacks_flag(
            transcript.has_introns(),
            transcript.get_ref_splice_pattern() is not None,
            transcript.get_ref_supported() is not None,
        ):
            num_unflagged += 1

    return _report_guide_class(gtf_label, num_ref, num_LRAA, num_unflagged)


def validate_guide_gtf_file(gtf_path, exempt_contigs=()):
    """validate_guide_transcripts over a GTF file, streamed rather than parsed.

    LRAA calls it once, before any work unit is launched, so a refused guide stops the
    run at startup instead of at its first discovery unit. Only the rows with column 2
    "LRAA" are tracked (exon count and flags per transcript), so a full GENCODE guide
    costs one pass over the text. Rows on exempt_contigs (oversimplified contigs, which
    are carried rather than filtered) are skipped.
    """
    if not gtf_path or not reference_rules_active():
        return None

    opener = gzip.open if gtf_path.endswith(".gz") else open
    reference_ids = set()
    lraa_num_exons = defaultdict(int)
    lraa_flags = defaultdict(set)
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[0] in exempt_contigs:
                continue
            m = _TRANSCRIPT_ID.search(fields[8])
            if m is None:
                continue
            key = (fields[0], m.group(1))
            if fields[1] != "LRAA":
                reference_ids.add(key)
                continue
            if fields[2] == "exon":
                lraa_num_exons[key] += 1
            for flag in _FLAG_ATTRS.findall(fields[8]):
                lraa_flags[key].add(flag)

    num_unflagged = sum(
        _lacks_flag(
            num_exons > 1,
            "ref_splice_pattern" in lraa_flags[key],
            "ref_supported" in lraa_flags[key],
        )
        for key, num_exons in lraa_num_exons.items()
    )
    return _report_guide_class(
        gtf_path, len(reference_ids), len(lraa_num_exons), num_unflagged
    )


_TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]*)"')
_FLAG_ATTRS = re.compile(r'\b(ref_splice_pattern|ref_supported) "')


def _lacks_flag(is_multi_exon, has_ref_splice_pattern, has_ref_supported):
    """Whether an LRAA guide row lacks the flag the active rule needs for its kind."""
    if is_multi_exon:
        return (
            LRAA_Globals.config["ref_splice_pattern_reprieve"]
            and not has_ref_splice_pattern
        )
    return LRAA_Globals.config["ref_supported_monoexonic"] and not has_ref_supported


def _report_guide_class(gtf_label, num_ref, num_LRAA, num_unflagged):
    logger.info(
        "guide GTF %s: %d reference-annotation models, %d LRAA-produced models "
        "(%d without their ref_splice_pattern / ref_supported attribute)",
        gtf_label,
        num_ref,
        num_LRAA,
        num_unflagged,
    )

    if num_unflagged and not LRAA_Globals.config["allow_unflagged_LRAA_guide"]:
        raise RuntimeError(
            "guide GTF {} was written by LRAA (column 2 'LRAA') but {} of its models lack "
            "the ref_splice_pattern (multi-exon) / ref_supported (single-exon) attribute, "
            "so it cannot say which of them stand for reference transcripts. Regenerate it "
            "with LRAA v0.50.0 or later, or set allow_unflagged_LRAA_guide (via "
            "--config_update) to treat all of them as references as earlier versions "
            "did.".format(gtf_label, num_unflagged)
        )

    return num_ref, num_LRAA, num_unflagged


def record_carried_model_identity(transcripts):
    """Record reference identity on models an oversimplified contig writes unchanged.

    No filter runs there, so nothing else sets it, and every row LRAA writes must carry
    it or the GTF is refused as a guide. A reference annotation row is a reference
    transcript (ref_splice_pattern / ref_supported "True"); an LRAA guide row keeps what
    it was parsed with; a model built here (the de novo aggregate) is not a reference.
    """
    if not reference_rules_active():
        return
    for transcript in transcripts:
        if transcript.is_reference_annotation_model():
            is_reference = True
        elif transcript.is_LRAA_guide_model():
            continue
        else:
            is_reference = False
        if transcript.has_introns():
            if (
                LRAA_Globals.config["ref_splice_pattern_reprieve"]
                and transcript.get_ref_splice_pattern() is None
            ):
                transcript.set_ref_splice_pattern(is_reference)
        elif (
            LRAA_Globals.config["ref_supported_monoexonic"]
            and transcript.get_ref_supported() is None
        ):
            transcript.set_ref_supported(is_reference)


def assign_ref_supported_monoexonic(transcripts, guide_transcripts):
    """Set ref_supported on every single-exon model; return how many are supported.

    A single-exon model is supported when, on its strand, it and a reference single-exon
    transcript (guide.is_reference_single_exon_source()) each have at least
    ref_supported_monoexonic_min_reciprocal_overlap of their length in the overlap.
    Lying inside an exon of a multi-exon reference, spanning a reference intron, the
    opposite strand, or a model much longer or shorter than the reference is not
    support. With no guide every single-exon model is unsupported.
    """
    min_overlap = LRAA_Globals.config["ref_supported_monoexonic_min_reciprocal_overlap"]

    strand_to_refs = defaultdict(list)
    for guide in guide_transcripts or []:
        if guide.is_reference_single_exon_source():
            strand_to_refs[guide.get_strand()].append(tuple(guide.get_coords()))

    # Sorted by left end, with the longest reference length per strand, so each model
    # scans only the references that can reach it.
    strand_to_lends = dict()
    strand_to_max_len = dict()
    for strand, refs in strand_to_refs.items():
        refs.sort()
        strand_to_lends[strand] = [lend for lend, _ in refs]
        strand_to_max_len[strand] = max(rend - lend + 1 for lend, rend in refs)

    num_supported = 0
    for transcript in transcripts:
        if transcript.has_introns():
            continue
        lend, rend = transcript.get_coords()
        model_len = rend - lend + 1
        strand = transcript.get_strand()
        supported = False
        refs = strand_to_refs.get(strand, [])
        if refs:
            first = bisect.bisect_left(
                strand_to_lends[strand], lend - strand_to_max_len[strand]
            )
            for ref_lend, ref_rend in refs[first:]:
                if ref_lend > rend:
                    break
                overlap = min(rend, ref_rend) - max(lend, ref_lend) + 1
                if overlap <= 0:
                    continue
                if (
                    overlap / model_len >= min_overlap
                    and overlap / (ref_rend - ref_lend + 1) >= min_overlap
                ):
                    supported = True
                    break
        transcript.set_ref_supported(supported)
        num_supported += supported

    return num_supported


def pooled_ref_FSM_reads_from_guides(guide_transcripts):
    """{splice_hash: ref_FSM_reads} carried by the guide's reference-chain rows.

    An LRAA guide written with these counts (the init GTF guiding a cluster run) has
    already counted its reference chains over the pooled reads, which is the question
    the reprieve asks: is this annotated chain real in the sample. Returns None when any
    reference-chain guide row lacks the count (a reference annotation never carries it),
    and the caller then scans its own bam.
    """
    counts = dict()
    for guide in guide_transcripts or []:
        if not guide.is_reference_splice_pattern_source():
            continue
        num_reads = guide.get_ref_FSM_reads()
        if num_reads is None:
            return None
        splice_hash = guide.get_splice_hash()
        counts[splice_hash] = max(counts.get(splice_hash, 0), num_reads)
    return counts


def count_reads_by_reference_chain(
    bam_file, contig_acc, contig_strand, region_lend, region_rend, ref_splice_hashes
):
    """{splice_hash: reads of bam_file whose intron chain is exactly that chain}.

    Read selection and intron convention are quant's: Bam_alignment_extractor's
    iter_read_alignments (Util_funcs.quant_discard_reason: secondary, supplementary,
    unmapped, other strand, mapping quality, percent identity, intron length) and
    Pretty_alignment.get_splice_hashcode (pretty segments, so short deletions are
    merged and long ones split like an intron). A read is fetched when it overlaps the
    region, so a read starting before region_lend still counts. In a chunked run each
    chunk's bam lacks the reads severed at its cuts; those are not counted anywhere.
    """
    counts = dict.fromkeys(ref_splice_hashes, 0)
    if not counts or bam_file is None:
        return counts

    extractor = Bam_alignment_extractor(bam_file)
    for read in extractor.iter_read_alignments(
        contig_acc, contig_strand, region_lend, region_rend, pretty=False
    ):
        # no N and no D: no gap to become an intron, so no chain to match
        if not any(op in (2, 3) for op, _ in read.cigartuples):
            continue
        pretty_alignment = Pretty_alignment.get_pretty_alignment(read)
        splice_hash = pretty_alignment.get_splice_hashcode(contig_acc)
        if splice_hash in counts:
            counts[splice_hash] += 1

    return counts


def assign_ref_FSM_reads(
    transcripts,
    guide_transcripts,
    ref_splice_hashes,
    bam_file,
    contig_acc,
    contig_strand,
    region_lend,
    region_rend,
):
    """Set ref_FSM_reads on each model whose chain is a reference chain.

    The counts come from the guide when it carries them
    (pooled_ref_FSM_reads_from_guides), else from one pass over this work unit's full
    bam (count_reads_by_reference_chain). Models with a novel chain are left at None.
    Returns {splice_hash: count}.
    """
    start_time = time.time()
    counts = pooled_ref_FSM_reads_from_guides(guide_transcripts)
    source = "guide GTF (pooled)"
    if counts is None:
        counts = count_reads_by_reference_chain(
            bam_file,
            contig_acc,
            contig_strand,
            region_lend,
            region_rend,
            ref_splice_hashes,
        )
        source = "full bam"

    for transcript in transcripts:
        if transcript.has_introns() and transcript.has_novel_splice_pattern() is False:
            transcript.set_ref_FSM_reads(counts.get(transcript.get_splice_hash(), 0))

    logger.info(
        "[%s%s] ref_FSM_reads from %s: %d reference chains, %d carried by >= 1 read, %.1f s",
        contig_acc,
        contig_strand,
        source,
        len(counts),
        sum(1 for n in counts.values() if n > 0),
        time.time() - start_time,
    )

    return counts
