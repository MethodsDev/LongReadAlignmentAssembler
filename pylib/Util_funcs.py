#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
import subprocess
import logging
import string
import time
import pysam
import LRAA_Globals
from hashlib import blake2s

from collections import defaultdict

logger = logging.getLogger(__name__)


def retrieve_contig_seq_from_fasta_file(contig_acc, fasta_filename):

    # samtools faidx <fasta> <region>
    attempts = 3
    contig_seq_str = None
    last_error = None

    for attempt in range(1, attempts + 1):
        try:
            contig_seq_str = subprocess.check_output(
                "samtools faidx {} {}".format(fasta_filename, contig_acc),
                shell=True,
                encoding="utf-8",
            )
            break
        except subprocess.CalledProcessError as exc:
            last_error = exc
            logger.warning(
                "samtools faidx failed for %s (attempt %d/%d): %s",
                contig_acc,
                attempt,
                attempts,
                exc,
            )
            if attempt < attempts:
                time.sleep(1 * attempt)

    if contig_seq_str is None:
        raise last_error

    contig_seq_str = contig_seq_str.upper()
    contig_seq_str = contig_seq_str.split("\n")
    contig_seq_str = contig_seq_str[1:]
    contig_seq_str = "".join(contig_seq_str)

    contig_seq_str = re.sub("\\s", "", contig_seq_str)  # just in case

    return contig_seq_str


def coordpairs_overlap(coordset_A, coordset_B):

    A_lend, A_rend = coordset_A

    B_lend, B_rend = coordset_B

    if A_lend <= B_rend and A_rend >= B_lend:
        return True
    else:
        return False


def get_num_overlapping_bases(coordset_A, coordset_B):

    if not coordpairs_overlap(coordset_A, coordset_B):
        return 0

    coords = sorted([coordset_A[0], coordset_A[1], coordset_B[0], coordset_B[1]])
    overlap_len = coords[2] - coords[1] + 1

    return overlap_len


def get_read_name_include_sc_encoding(pysam_read_alignment):

    cell_barcode_tag = LRAA_Globals.config["cell_barcode_tag"]
    read_umi_tag = LRAA_Globals.config["read_umi_tag"]

    read = pysam_read_alignment
    if read.has_tag(cell_barcode_tag) and read.has_tag(read_umi_tag):
        cell_barcode = read.get_tag(cell_barcode_tag)
        umi = read.get_tag(read_umi_tag)
        read_name = "^".join([cell_barcode, umi, read.query_name])
        return read_name
    else:
        return read.query_name


# CIGAR operation 'N' = skipped region from the reference, ie. an intron.
# Introns are identified from the CIGAR rather than from read.get_blocks(),
# because get_blocks() reports an identical gap for a deletion ('D') as for an
# intron of the same length, and additionally splits blocks at insertions.
# The 'N' operations are exactly the introns, terminal and internal alike.
INTRON_CIGAR_OP = 3


def get_longest_intron_length(read):
    """length of the longest intron (CIGAR 'N' operation) of the read, or zero if unspliced"""

    cigartuples = read.cigartuples
    if cigartuples is None:
        return 0

    longest_intron_length = 0
    for cigar_op, cigar_op_length in cigartuples:
        if cigar_op == INTRON_CIGAR_OP and cigar_op_length > longest_intron_length:
            longest_intron_length = cigar_op_length

    return longest_intron_length


def has_disqualifying_long_intron(read, max_intron_length):
    """True when the read carries an intron the pipeline declines to model.

    An intron of exactly max_intron_length is acceptable; one base longer is not,
    matching the thresholds enforced in Pretty_alignment and Splice_graph.
    max_intron_length <= 0 disables intron length filtering.

    This is the single implementation of the intron-length rule.  It governs both
    the splice-graph input (util/separate_bam_by_strand.py) and read assignment
    (Bam_alignment_extractor), so the two always see one record set.
    """

    if max_intron_length <= 0:
        return False

    return get_longest_intron_length(read) > max_intron_length


def quant_discard_reason(
    read,
    contig_strand=None,
    max_intron_length=None,
    min_mapping_quality=None,
    min_per_id=None,
):
    """Why quantification would discard this alignment, or None if it keeps it.

    The full retention policy Bam_alignment_extractor applies, in one callable
    place.  Reasons are the keys its discard counter uses, so a caller can report
    them in the same vocabulary.

    This exists because the policy was previously only expressed inline inside the
    extractor's fetch loop, which meant nothing else could ask the question and
    every other consumer reimplemented an approximation instead.  Cut selection's
    ``retained_for_extraction`` is a deliberate superset -- it omits mapping
    quality and percent identity, which is sound for BOUNDING the cost of a cut,
    since over-counting only biases toward safer positions.  It is not sound for
    naming the alignments quantification will actually use: a read this predicate
    rejects is one no downstream stage ever sees, so attributing anything to it is
    a false positive.

    Percent identity is cheap here -- CIGAR stats plus the NM/nM tag, no sequence
    -- so there is no cost argument for leaving it out.
    """

    if max_intron_length is None:
        max_intron_length = LRAA_Globals.config["max_intron_length"]
    if min_mapping_quality is None:
        min_mapping_quality = int(LRAA_Globals.config["min_mapping_quality"])
    if min_per_id is None:
        min_per_id = LRAA_Globals.config["min_per_id"]

    if read.is_unmapped:
        return "unmapped"
    if read.reference_id < 0:
        return "no_chromosome"
    if contig_strand is not None:
        if read.is_forward and contig_strand != "+":
            return "wrong_strand"
        if read.is_reverse and contig_strand != "-":
            return "wrong_strand"
    if read.mapping_quality < min_mapping_quality:
        return "min_mapping_quality"
    if read.is_paired and not read.is_proper_pair:
        return "improper_pair"
    if read.is_duplicate:
        return "duplicate"
    if read.is_qcfail:
        return "qcfail"
    if read.is_supplementary:
        return "supplementary"
    if read.is_secondary:
        return "secondary"
    if has_disqualifying_long_intron(read, max_intron_length):
        return "long_intron"

    # Absent NM and nM there is nothing to measure, and the extractor keeps the
    # read rather than guessing; identical behaviour here.
    cigar_stats = read.get_cigar_stats()
    aligned_base_count = cigar_stats[0][0]
    if aligned_base_count == 0:
        aligned_base_count = cigar_stats[0][7] + cigar_stats[0][8]
    mismatch_count = None
    if read.has_tag("NM"):
        mismatch_count = int(read.get_tag("NM"))
    elif read.has_tag("nM"):
        mismatch_count = int(read.get_tag("nM"))
    if mismatch_count is not None and aligned_base_count > 0:
        per_id = 100 - (mismatch_count / aligned_base_count) * 100
        if per_id < min_per_id:
            return "low_perID"

    return None


def retained_for_quantification(read, contig_strand=None, **kwargs):
    """True when quantification would use this alignment.  See quant_discard_reason."""

    return quant_discard_reason(read, contig_strand, **kwargs) is None


def frac_base_composition(nuc_seq, nuc_base):

    assert len(nuc_seq) > 0, "Error, nuc_seq is empty"

    counter = 0

    for char in nuc_seq:
        if char.upper() == nuc_base.upper():
            counter += 1

    frac_base = counter / len(nuc_seq)

    return frac_base


def looks_internally_primed(
    contig_seq_str, three_prime_pos, strand, check_length=20
):
    """True when the genome immediately 3' of `three_prime_pos` is itself A-rich in
    transcript sense: the oligo-dT primer had a template to bind without any
    post-transcriptional tail, so the 3' end is an internal-priming artifact rather
    than a cleavage site.

    `three_prime_pos` is the 1-based genomic coordinate of the 3'-most transcribed base
    -- a transcript's rend on '+' and its lend on '-', or equally a candidate PolyA
    site's coordinate. The test fires on either a diffuse A-run (>= 12 of the 20
    downstream bases) or a hexamer of A's within the first 8, unchanged from the rule
    that has always been applied at transcript filtering.

    Single definition, called from two stages: the PolyA candidate gate in
    Splice_graph._incorporate_PolyA_objects, which stops an A-rich read-derived
    candidate becoming a vertex, and
    TranscriptFiltering.filter_internally_primed_transcripts, which judges the emitted
    model's own terminus. Copying the rule into either caller is what would let the two
    drift apart.
    """

    if strand not in {"+", "-"}:
        raise ValueError("Strand must be '+' or '-'")

    target_base = "A" if strand == "+" else "T"
    target_polyA_motif = target_base * 6

    contig_length = len(contig_seq_str)

    if strand == "+":
        start = three_prime_pos + 1
        end = start + check_length - 1
    else:
        end = three_prime_pos - 1
        start = end - check_length + 1

    # ensure coordinates within contig bounds
    start = max(1, start)
    end = min(end, contig_length)

    extracted_long_sequence = contig_seq_str[start - 1 : end].upper()
    extracted_short_sequence = (
        extracted_long_sequence[-8:] if strand == "-" else extracted_long_sequence[0:8]
    )

    return (
        extracted_long_sequence.count(target_base) >= 12
        or target_polyA_motif in extracted_short_sequence
    )


def get_hash_code(input_string):
    hash_object = blake2s(digest_size=11)
    hash_object.update(input_string.encode("utf-8"))
    hex_digest = hash_object.hexdigest()
    hex_digest = str(hex_digest)

    return hex_digest


def file_identity_token(path):
    """Short digest of which file this is and whether it has changed since.

    The resolved path separates two inputs that happen to share a basename; size
    and modification time catch one replaced in place. Contents are deliberately
    not hashed: these bams run to gigabytes, and reading one on every startup
    would cost more than the step being cached, to cover a case the stat pair
    already covers.
    """
    resolved = os.path.realpath(path)
    try:
        stat_result = os.stat(resolved)
        identity = "{}:{}:{}".format(
            resolved, stat_result.st_size, stat_result.st_mtime_ns
        )
    except OSError:
        # A missing input is the caller's problem to report; naming it by path
        # alone keeps this from raising first and obscuring that.
        identity = resolved

    return get_hash_code(identity)[:12]


def splice_graph_norm_cache_stem(
    base_root, normalize_max_cov_level, source_bam, max_intron_length
):
    """Name for a normalized bam and its work directory.

    Everything that determines the contents belongs in the name, because a hit is
    taken on trust and nothing downstream can audit it.

    The method matters most: reads selected by an algorithm no longer in the tree
    are undetectable once written, since a bam from the read-start-binning era
    carries no weight tag and an untagged read correctly weighs 1.

    The source is identified by content fingerprint rather than by name. base_root
    is only a basename, so without this a second input called sample.bam from a
    different directory, or the same path regenerated, would land on the cache
    built for the first.

    The target depth is here for the same reason, and it also scopes the work
    directory, whose checkpoints would otherwise let a run at one depth reuse the
    strand split and merge performed for another.

    The intron cap decides which records the strand split emits at all, so two
    runs differing only in --max_intron_length produce different normalized bams
    and must not share a name.
    """
    return "{}.norm_{}.maxintron_{}.{}.{}".format(
        base_root,
        normalize_max_cov_level,
        max_intron_length,
        LRAA_Globals.SPLICE_GRAPH_NORMALIZATION_METHOD,
        file_identity_token(source_bam),
    )


def available_cpus():
    """Number of CPUs this process may actually run on.

    `os.cpu_count()` reports the machine's CPUs, not the ones the process is allowed to
    use. Under a cpuset -- every container runtime, Slurm with `--cpus-per-task`, taskset
    -- those differ, and sizing a thread pool from `cpu_count()` oversubscribes by the
    ratio between them: 4x for 4 cpus granted out of 16 present.

    `os.sched_getaffinity` is the affinity-aware count but exists only on Linux, so it is
    probed rather than assumed. `os.cpu_count()` can itself return None when the count is
    indeterminable, hence the floor of 1.
    """

    get_affinity = getattr(os, "sched_getaffinity", None)
    if get_affinity is not None:
        try:
            return max(1, len(get_affinity(0)))
        except OSError:
            pass

    return max(1, os.cpu_count() or 1)


# Canonical polyadenylation signal hexamers, in transcript sense.  These two account
# for the large majority of characterised human cleavage sites; the many
# single-substitution variants are deliberately excluded, so a positive call means the
# canonical signal rather than a permissive match.
POLYA_SIGNAL_HEXAMERS = ("AATAAA", "ATTAAA")

# Transcript-sense bounds of the EXTRACTED REGION, relative to the cleavage site.  A
# motif must fall entirely inside it, so a 6-mer in -40..-10 has first-base offsets of
# -40..-15.  Note that is NOT the same set as "every start 10 to 30 nt upstream".
#
# Deliberately broader than the textbook 10-30 nt spacing.  This is an ANNOTATION, and
# PAS_offset is emitted alongside the motif, so a scan that reports a motif at -38 loses
# nothing: a consumer wanting the canonical spacing filters on the offset.  Narrowing
# the scan instead would make an out-of-window motif indistinguishable from no motif at
# all, which is the one thing an annotation should not do.
#
# -40..-10 also matches the offline analysis these numbers are compared against,
# investigations/PIP_work/pip_profile.py, which extracts sense(pos, -40, -10) and
# substring-searches it.
#
# For scale, on 3,073 chr20 ref-free termini against 9,219 strand-matched random genomic
# positions, canonical first-base offsets are ~26-62x more frequent than background at
# -23..-17, ~5x by -15, and indistinguishable from background by -14 and beyond -41.
# That is an enrichment estimate against a random-position null, not a precision or
# recall measurement -- nothing here carries validated true-positive labels -- so it
# describes where signal concentrates, and does not by itself justify a narrower default.
POLYA_SIGNAL_WINDOW = (-40, -10)

_DNA_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _genomic_spelling_on_minus(motif):
    """How a transcript-sense motif is written in the reference on the minus strand.

    Scanning for this instead of reverse-complementing the extracted slice keeps the
    reference read-only and matches what looks_internally_primed does when it counts T
    rather than A on '-'.
    """
    return "".join(reversed([_DNA_COMPLEMENT[base] for base in motif]))


def find_polyA_signal(
    contig_seq_str,
    three_prime_pos,
    strand,
    window=POLYA_SIGNAL_WINDOW,
    hexamers=POLYA_SIGNAL_HEXAMERS,
):
    """Look for a canonical polyadenylation signal upstream of a 3' terminus.

    Returns (hexamer, offset) for the match closest to the site, or (None, None).
    `hexamer` is always reported in transcript sense, so a minus-strand hit reads
    AATAAA rather than its genomic spelling TTTATT.  `offset` is the signed distance
    from `three_prime_pos` to the hexamer's first base in transcript sense, so it is
    negative: -21 means the signal begins 21 nt upstream of the cleavage site.

    `window` bounds the extracted region and the hexamer must lie entirely within it,
    so a 6-mer in the default -40..-10 window has first-base offsets of -40..-15.

    This is evidence ABOUT a 3' end and is deliberately independent of
    `looks_internally_primed`, which interrogates the opposite side of the same site.
    A terminus can carry both -- a genuine signal upstream and an A-rich stretch
    downstream -- and the two belong in the output separately rather than collapsed
    into one verdict.

    `three_prime_pos` is the 1-based genomic coordinate of the 3'-most transcribed
    base: a transcript's rend on '+' and its lend on '-'.
    """

    if strand not in {"+", "-"}:
        raise ValueError("Strand must be '+' or '-'")

    lo, hi = window
    if lo > hi:
        raise ValueError("window must be (low, high) with low <= high")

    # The scan reads fixed-width k-mers and derives the containment bound from that one
    # width, so a mixed-width motif set would silently mis-place the shorter members.
    hexamers = tuple(hexamers)
    if not hexamers:
        raise ValueError("hexamers must be a non-empty sequence of motifs")
    widths = {len(motif) for motif in hexamers}
    if len(widths) != 1:
        raise ValueError(
            f"all motifs must share one length; got lengths {sorted(widths)}"
        )
    if any(set(motif) - set("ACGT") for motif in hexamers):
        raise ValueError("motifs must be uppercase and contain only A, C, G or T")

    if strand == "+":
        patterns = {motif: motif for motif in hexamers}
    else:
        patterns = {_genomic_spelling_on_minus(motif): motif for motif in hexamers}

    width = widths.pop()

    # Genomic span that can hold a hexamer whose transcript-sense first base falls in
    # [lo, hi].  On '+' that first base is the leftmost genomic base; on '-' it is the
    # rightmost, so the span shifts by the hexamer width.
    # Only fully contained hexamers count, so the first base can lie no later than
    # hi - width + 1 and the region needs no widening beyond the window itself.
    if strand == "+":
        start, end = three_prime_pos + lo, three_prime_pos + hi
    else:
        start, end = three_prime_pos - hi, three_prime_pos - lo

    start = max(1, start)
    end = min(len(contig_seq_str), end)
    if end - start + 1 < width:
        return (None, None)

    region = contig_seq_str[start - 1 : end].upper()

    best_hexamer, best_offset = None, None
    for i in range(0, len(region) - width + 1):
        hexamer = patterns.get(region[i : i + width])
        if hexamer is None:
            continue

        genomic_start = start + i
        if strand == "+":
            offset = genomic_start - three_prime_pos
        else:
            offset = three_prime_pos - (genomic_start + width - 1)

        # fully contained: both the first and last base inside the window
        if not (lo <= offset and offset + width - 1 <= hi):
            continue

        # closest to the site wins, i.e. the least negative offset
        if best_offset is None or offset > best_offset:
            best_hexamer, best_offset = hexamer, offset

    return (best_hexamer, best_offset)


def resolve_polyA_signal_settings(config=None):
    """Validate and normalise the PAS annotation settings from config.

    Returns (motifs, window) as a tuple of uppercase motifs and an (lo, hi) pair.
    Raises ValueError with a message naming the offending key, so a bad organism
    configuration fails at startup rather than part-way through a run.
    """

    if config is None:
        config = LRAA_Globals.config

    raw_motifs = config.get("polyA_signal_motifs", POLYA_SIGNAL_HEXAMERS)
    if isinstance(raw_motifs, str):
        raw_motifs = [token for token in re.split(r"[\s,]+", raw_motifs.strip()) if token]
    motifs = tuple(str(motif).strip().upper() for motif in raw_motifs)

    if not motifs:
        raise ValueError("polyA_signal_motifs must list at least one motif")
    widths = {len(motif) for motif in motifs}
    if len(widths) != 1:
        raise ValueError(
            "polyA_signal_motifs must all share one length, since the search derives "
            f"its containment bound from that length; got lengths {sorted(widths)} "
            f"for {list(motifs)}"
        )
    for motif in motifs:
        if set(motif) - set("ACGT"):
            raise ValueError(
                f"polyA_signal_motifs entry {motif!r} must contain only A, C, G or T"
            )

    raw_window = config.get("polyA_signal_window", POLYA_SIGNAL_WINDOW)
    if isinstance(raw_window, str):
        raw_window = [t for t in re.split(r"[\s,]+", raw_window.strip()) if t]
    try:
        window = tuple(int(bound) for bound in raw_window)
    except (TypeError, ValueError):
        raise ValueError(
            f"polyA_signal_window must be two integers; got {raw_window!r}"
        )
    if len(window) != 2:
        raise ValueError(
            f"polyA_signal_window must be exactly two integers (low, high); got {list(window)}"
        )
    lo, hi = window
    if lo > hi:
        raise ValueError(
            f"polyA_signal_window must be (low, high) with low <= high; got {list(window)}"
        )
    width = widths.pop()
    if hi - lo + 1 < width:
        raise ValueError(
            f"polyA_signal_window {list(window)} spans {hi - lo + 1} nt, too narrow to "
            f"contain a motif of length {width}"
        )

    return motifs, window
