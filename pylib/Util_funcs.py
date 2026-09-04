#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
import subprocess
import logging
import string
import time
import pysam
import LRAA_Globals
import RdnaMask
from hashlib import blake2s

from collections import defaultdict

logger = logging.getLogger(__name__)


def aligned_base_count(read):
    """Read bases sitting in aligned columns: M + = + X.

    A CIGAR conventionally uses either M or the =/X pair, but mixing them in one
    record is legal and minimap2 emits it -- measured at 0.32% of alignments on
    XP132160.ucsc.bam, e.g. `627=1M25=`. Taking M alone and falling back to =/X
    only when M is exactly zero therefore under-counts those records to the M
    portion, which is a silent error wherever the result is a denominator.

    Summing all three is correct for every convention at once rather than
    special-casing two of them, and needs no knowledge of whether the aligner ran
    with --eqx -- which cannot be relied on regardless: XP132160.ucsc.bam is =/X
    with a mixed tail while PBMCs_pbio.aligned.sorted.bam is 99.82% M.

    Excludes I, D and N: inserted bases occupy no reference column, deleted and
    skipped bases occupy no read base. Callers needing alignment columns rather
    than read bases add I and D themselves.
    """
    stats = read.get_cigar_stats()[0]
    return int(stats[0] + stats[7] + stats[8])


def substitution_count(read):
    """Substituted bases only, indels excluded, or None when unmeasurable.

    Callers that want "how many read bases did this alignment get wrong" need
    substitutions alone, and the two mismatch tags reach that differently:

      NM  edit distance, so inserted and deleted bases are included and must be
          subtracted out.
      nM  STAR's tag, already mismatches only, so nothing is subtracted.

    Subtracting indels from nM destroys real substitutions rather than correcting
    anything: `10M5I10M` with nM:i:2 has two substituted bases, and NM-style
    subtraction reports zero, scoring a wrong alignment as perfect. The tag read
    therefore has to determine the arithmetic, which is why this is one function
    rather than the same three lines at each call site.

    FLOOR AT ZERO, deliberately, and it should be unreachable. NM counts every
    inserted and deleted base, so NM >= I + D holds for any self-consistent
    record, with equality when nothing is substituted. Reaching the floor means
    the NM tag and the CIGAR disagree -- the record is internally inconsistent,
    not merely low quality.

    Resolving that to zero substitutions is the permissive reading, chosen to
    match how the rest of this module treats an alignment it cannot assess: a
    missing tag returns None and callers pass the read rather than fail it. A
    malformed record is no more informative than an unmeasurable one, so it gets
    the same benefit of the doubt instead of a fabricated substitution count.
    Deliberately NOT logged: this runs per alignment over tens of millions of
    records, and a warning there would either be dropped or would bury the log.
    A caller needing to know should count the condition itself, which it can,
    since NM, I and D are all still available to it.
    """
    stats = read.get_cigar_stats()[0]
    if read.has_tag("NM"):
        edits = int(read.get_tag("NM"))
        return max(0, edits - int(stats[1]) - int(stats[2]))
    if read.has_tag("nM"):
        # nM needs no subtraction, so the floor here only guards a negative tag.
        return max(0, int(read.get_tag("nM")))
    return None


def alignment_edit_count(read):
    """Edit distance -- substitutions plus inserted and deleted bases -- or None.

    The complement of substitution_count(): same tag disambiguation, opposite
    adjustment. NM already IS an edit distance and is returned as-is. nM counts
    mismatches only, so I and D are ADDED to reach the same quantity.

    Without that addition a score built on nM charges nothing for indels, so two
    candidate alignments of one read differing only in indel content tie -- which
    defeats the comparison the score exists to make. `10M5I10M` scored 18 under
    nM:i:2 and 13 under the equivalent NM:i:7, for the same alignment. Producer
    agreement does not save it: both candidates carry nM, but each has its own I
    and D, and neither is reflected.
    """
    stats = read.get_cigar_stats()[0]
    if read.has_tag("NM"):
        return max(0, int(read.get_tag("NM")))
    if read.has_tag("nM"):
        return max(0, int(read.get_tag("nM")) + int(stats[1]) + int(stats[2]))
    return None


def alignment_per_id(read):
    """Percent identity of an alignment, or None when it cannot be determined.

    The single definition shared by everything that filters on identity. Coverage
    normalization has to apply the same floor as alignment extraction, because depth
    measured over reads the consumer rejects yields acceptance probabilities -- and so
    XW weights -- for a coverage level nothing downstream sees. Two implementations
    that agree today drift tomorrow, and the failure is silent: the bam still looks
    normalized, its weights are just wrong.

    None means no mismatch tag was present, which callers treat as passing rather than
    failing -- absence of evidence is not evidence of a bad alignment.
    """
    # The denominator must be the same set of alignment columns that NM counts
    # edits over, or the ratio is not a fraction and the result is not bounded.
    # NM is an edit distance: mismatched bases PLUS inserted bases PLUS deleted
    # bases. So the denominator is M + = + X + I + D. N is excluded deliberately
    # -- a reference skip is an intron, not an edit, and NM does not count it.
    #
    # Two distinct defects were measured here, both silent. Figures below are from
    # the FINAL remeasurement with this definition in place, comparing all three
    # candidate denominators on the same samples -- earlier drafts of this comment
    # quoted intermediate numbers taken before defect 2 was found.
    #
    # 1. Aligned bases were taken as M alone, falling back to =/X only when M was
    #    exactly zero. A CIGAR conventionally uses either M or the =/X pair, but
    #    mixing them in one record is legal and minimap2 emits it. On
    #    XP132160.ucsc.bam (46,123,848 records, 1-in-50 sample = 922,477) 0.3209%
    #    of alignments carry a few M ops among many =/X ops, e.g. `627=1M25=`:
    #    653 aligned bases and 99.85% identity, scored as 1 aligned base and
    #    0.0%, so every consumer applying a floor discarded it. Values reached
    #    -13900%.
    #
    # 2. Excluding I and D while NM counts them left the result unbounded even for
    #    records with a single CIGAR convention. On PBMCs_pbio.aligned.sorted.bam,
    #    99.82% M-op records with ZERO mixed CIGARs, 0.0320% of alignments still
    #    came out negative from indel-heavy alignments whose NM exceeded their
    #    aligned-base count. Fixing defect 1 alone did not touch these.
    #
    # Effect of the two fixes together, apparent-vs-genuine below an 80% floor on
    # SBX and a 97% floor on PacBio:
    #
    #                   below floor    negatives    median      min
    #   SBX   old         0.5280%       0.2866%     99.5305   -13900.00
    #         M+=+X       0.2089%       0.0212%     99.5338     -308.38
    #         final       0.1767%       0.0000%     99.5345       18.75
    #   PB    old         1.7673%       0.0320%     99.8473     -514.42
    #         M+=+X       1.7673%       0.0320%     99.8473     -514.42
    #         final       1.7392%       0.0000%     99.8476       14.00
    #
    # So 0.3513% of all SBX alignments and 0.0281% of PacBio's were being wrongly
    # discarded, and two thirds of everything appearing to fail SBX's floor was
    # these defects rather than poor alignment. Medians move by ten-thousandths,
    # so the tail is repaired without shifting the central distribution.
    #
    # The two corpora use OPPOSITE conventions -- SBX is =/X with a mixed tail,
    # PacBio is M -- so neither can be assumed and the presence of --eqx on the
    # aligner's @PG line does not settle it.
    #
    # This is GAP-AWARE percent identity: every edit, gaps included, over the
    # alignment's full length, gaps included.
    #
    #   100 - (mismatches + inserted + deleted) / (M + = + X + I + D) * 100
    #
    # Both parts of that have to move together. The defect was a numerator that
    # counted gaps over a denominator that did not, which is neither the gap-aware
    # nor the gap-excluded definition and is why the result could leave [0, 100].
    #
    # NM and nM reach the numerator differently and alignment_edit_count() is what
    # reconciles them: NM already IS mismatches + I + D, while nM counts
    # MISMATCHES ONLY and has I + D added. Using nM raw over the aligned bases
    # alone would compute a gap-EXCLUDED identity -- a different metric under the
    # same name, silently selected by which aligner wrote the bam. Measured on
    # `10M5I10M`: 90.0 gap-excluded against 72.0 gap-aware, for one alignment.
    #
    # The result is therefore tag-independent, and agrees exactly with
    # IsoformReadRescue._gap_aware_identity(), which reaches the same quantity by
    # the algebraically equivalent route of matched bases over span. That equality
    # is asserted in test_alignment_per_id_mixed_cigar.py; it is the only thing
    # keeping two modules' definitions from drifting apart again.
    #
    # Bounded by arithmetic rather than by input: mismatches <= M+=+X, insertions
    # = I and deletions = D, so the numerator cannot exceed the denominator under
    # either tag. One caveat belongs to nM's own definition -- STAR documents it
    # as mismatches per PAIRED read, summed across both mates, while this CIGAR
    # describes one mate, so for a proper pair the numerator can exceed this
    # record's alignment length. Not guarded against, because one record cannot
    # resolve it; quant_discard_reason() already rejects improper pairs, and this
    # codebase's inputs are long single-end reads.
    stats = read.get_cigar_stats()[0]
    edit_count = alignment_edit_count(read)
    if edit_count is None:
        return None
    alignment_length = aligned_base_count(read) + int(stats[1]) + int(stats[2])
    if alignment_length == 0:
        return None
    return 100 - (edit_count / alignment_length) * 100


def alignment_passes_per_id(read, min_per_id):
    """Whether this alignment clears the identity floor. No floor, or no tag, passes."""
    if not min_per_id or min_per_id <= 0:
        return True
    per_id = alignment_per_id(read)
    return per_id is None or per_id >= min_per_id


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
    rdna_mask=None,
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

    ``rdna_mask`` follows the same "None reads the config default" convention as
    every other threshold here: None looks up ``config["rdna_mask_intervals"]``,
    the {contig: IntervalTree} RdnaMask.build_rdna_mask_bed populates at startup
    (or None, when masking found nothing or was disabled), rather than requiring
    every caller to thread it through explicitly. See RdnaMask.read_overlaps_mask.
    """

    if max_intron_length is None:
        max_intron_length = LRAA_Globals.config["max_intron_length"]
    if min_mapping_quality is None:
        min_mapping_quality = int(LRAA_Globals.config["min_mapping_quality"])
    if min_per_id is None:
        min_per_id = LRAA_Globals.config["min_per_id"]
    if rdna_mask is None:
        rdna_mask = LRAA_Globals.config.get("rdna_mask_intervals")

    if read.is_unmapped:
        return "unmapped"
    if read.reference_id < 0:
        return "no_chromosome"
    # Truthiness, not `is not None`.  The other two members of this predicate
    # family -- extract_contig_region_inputs.retained_for_extraction and its
    # load_gtf -- already read a falsy strand as "any orientation", which is what
    # an already strand-split bam wants.  This one alone demanded an exact match,
    # so a caller passing "" got wrong_strand for EVERY record while the looser
    # predicate beside it kept them all.  That is not a stricter subset, it is an
    # empty one, and it silently emptied the severed-alignment bam for every run
    # of the chunked pipeline, which omits --strand by design.
    if contig_strand:
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
    if rdna_mask and RdnaMask.read_overlaps_mask(read, rdna_mask):
        return "rdna_masked"
    if has_disqualifying_long_intron(read, max_intron_length):
        return "long_intron"

    # Absent NM and nM there is nothing to measure, and the extractor keeps the
    # read rather than guessing; identical behaviour here.
    #
    # Delegated to alignment_per_id rather than recomputed. This was an inline
    # copy of the formula, which is how the floor came to differ by code path:
    # the copy kept the M-only aligned-base count after the shared definition
    # was corrected, so the same read could be discarded here and kept by the
    # extractor. If it must ever diverge again, it needs a different NAME, not a
    # second implementation of this one.
    per_id = alignment_per_id(read)
    if per_id is not None and per_id < min_per_id:
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


def format_splice_compatible_id_set(id_sets, transcript_id):
    """Render one transcript's splice-compat-containment/contained-by entry.

    Shared by Quantify.report_quant_results (the default path) and
    StreamingQuant.write_expr (the --stream_reads path) so the two format the
    splice_compat_contained/splice_contained_by columns identically -- both are sets
    of transcript ids, whose iteration order is unspecified, so sorted() rather than
    set order makes the column reproducible between runs and byte-identical between
    the two quant paths. Absent from the map (no relationship recorded) prints "",
    not "{}", matching how neither path emits an entry for it in the first place.
    """
    if transcript_id not in id_sets:
        return ""
    return "{" + ",".join(sorted(id_sets[transcript_id])) + "}"


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
    base_root,
    normalize_max_cov_level,
    source_bam,
    max_intron_length,
    min_per_id,
    min_mapping_quality,
    depth_window,
    random_seed,
    window_origin,
    scope,
    rdna_mask_bed,
):
    """Name for a normalized bam and its work directory.

    Everything that determines the contents belongs in the name, because a hit is taken on
    trust and nothing downstream can audit it -- and because the driver returns on seeing
    this stem's checkpoint, it never runs the utility and so never consults the utility's
    own finer-grained token. Anything missing here is invisible.

    Every parameter is therefore required rather than defaulted: a caller that forgets one
    would key a bam by a name that does not describe it, and the omission would surface as
    a silently reused cache rather than an error.

    The method matters most: reads selected by an algorithm no longer in the tree are
    undetectable once written, since a bam from the read-start-binning era carries no
    weight tag and an untagged read correctly weighs 1.

    The source is identified by content fingerprint rather than by name. base_root is only
    a basename, so without this a second input called sample.bam from a different
    directory, or the same path regenerated, would land on the cache built for the first.

    The target depth is here for the same reason, and it also scopes the work directory,
    whose checkpoints would otherwise let a run at one depth reuse the strand split and
    merge performed for another.

    The intron cap decides which records the strand split emits at all, so two runs
    differing only in --max_intron_length produce different normalized bams and must not
    share a name.

    The identity floor and the mapping-quality floor decide which records count toward
    depth and which get written. The identity floor is why a HiFi run must not share a
    cache with a default one: 97 against 80 admits a materially different read population.
    The mapping-quality floor is whichever value the consumer will actually enforce, which
    under the quant stage is --min_mapping_quality_for_final_quant rather than
    --min_mapping_quality.

    The window, the seed and the window origin change the contents too -- the seed salts
    the per-read acceptance draw, while the window and its origin both move the
    depth-window boundaries and so decide which reads survive. The driver holds them fixed
    today, which is exactly why omitting them would be easy and would go unnoticed until
    someone changed a default.

    ``scope`` names WHICH PART of ``source_bam`` was read: ``None`` for the whole file,
    else an iterable of contig names the read was restricted to. Gated to render inert
    (``"none"``, the same convention ``window_origin`` uses for "unset") rather than simply
    omitted when unrestricted, so today's whole-source callers keep today's names. Without
    this field, normalizing a contig directly from a shared whole-genome bam would collide
    on the exact name a whole-bam normalization of the same file already uses -- two
    different outputs sharing one cache key -- which is silent in exactly the way this
    function's own docstring says every other field must not be. Each name is trusted
    pre-sanitized, the same contract ``base_root`` already has; sorted here so the caller's
    argument order is not another axis the key would need to track separately.

    ``rdna_mask_bed`` is the BED of excluded regions RdnaMask.build_rdna_mask_bed
    wrote (or None, when masking is disabled or found nothing for this genome):
    it decides which reads pass_2 counts as depth and which it writes, exactly as
    min_per_id and min_mapping_quality do, so a bam normalized under one mask
    state must not be reused under another. "none" rather than omitted, for the
    same reason window_origin renders "none": a caller that forgot to pass a
    real path must not collide with one that legitimately has no mask.
    """
    if scope is None:
        scope_token = "none"
    else:
        sorted_scope = sorted(scope)
        if not sorted_scope:
            raise ValueError(
                "scope must be None (unrestricted) or a non-empty iterable of contig "
                "names; an empty iterable names neither"
            )
        joined = "+".join(sorted_scope)
        # Bounded rather than always literal: --restrict_to_chromosomes can name dozens
        # of contigs, and a stem long enough to exceed a filesystem's path-component
        # limit would fail every caller rather than only the caller that hit it. The
        # common case -- one contig from --contig or --region -- stays fully readable.
        scope_token = (
            joined
            if len(joined) <= 80
            else "{}contigs.{}".format(
                len(sorted_scope), get_hash_code(joined)[:12]
            )
        )
    return "{}.norm_{}.maxintron_{}.{}.pid{}.mapq{}.w{}.s{}.o{}.scope{}.rdna{}.{}".format(
        base_root,
        normalize_max_cov_level,
        max_intron_length,
        LRAA_Globals.SPLICE_GRAPH_NORMALIZATION_METHOD,
        ("%g" % float(min_per_id or 0)),
        int(min_mapping_quality or 0),
        int(depth_window),
        int(random_seed),
        # "none" rather than 0: an unset origin anchors the grid per contig on the first
        # aligned base seen there, which is a different placement from the absolute grid
        # at 0, and collapsing the two here would let them share one cached bam.
        ("none" if window_origin is None else int(window_origin)),
        scope_token,
        ("none" if rdna_mask_bed is None else file_identity_token(rdna_mask_bed)),
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


def paths_name_one_file(path_a, path_b):
    """Whether two paths are ONE file on disk.

    Device plus inode, via ``os.path.samefile``. That catches everything path
    comparison catches AND a HARD LINK, which ``os.path.realpath`` equality does not:
    a hard link is the same file under a second name, and realpath reports the two
    names as distinct because there is no link to resolve.

    What this does NOT establish is that the two paths hold different DATA. A
    byte-identical copy has its own inode and is reported here as two files, which is
    the truth about the filesystem and not about the contents. Deciding contents needs
    the contents, and ``file_identity_token`` above is not a substitute on two counts:
    it hashes resolved path plus size plus mtime rather than bytes, and it includes the
    resolved path, so two distinct paths can never produce one token however identical
    their bytes. Callers guarding against "the same data supplied twice" must say that
    a copy defeats them rather than implying it does not.
    """

    try:
        return os.path.samefile(path_a, path_b)
    except OSError:
        # One path is gone or unstattable. Fall back to the weaker comparison rather
        # than raising from inside a guard: whatever needs the file reports its absence
        # with the context to act on, and a guard that raises first hides that.
        return os.path.realpath(path_a) == os.path.realpath(path_b)

def intron_chain_from_simple_path(simple_path):
    """The intron nodes of a splice-graph path, in order.

    Introns are nodes of the gene's own graph, so identity of node sequence is identity
    of intron chain -- which makes a full-splice-match test an exact tuple comparison
    without resolving either side back to coordinates. Shared by the reporting paths
    and TranscriptFiltering so "FSM" means one thing: the read's intron chain is
    EXACTLY the isoform's, not a prefix, not a subset, not merely compatible.

    Only valid within a single gene's graph. Pretty_alignment.get_splice_hashcode
    carries the coordinate form where stability across graphs matters.
    """
    return tuple(node for node in simple_path if node.startswith("I:"))


def intron_chain_from_segments(segments):
    """The intron chain implied by a pretty alignment's exon segments.

    Gaps BETWEEN consecutive segments, as (start, end) inclusive coordinates. Valid
    because these are LRAA's own post-processed segments -- the same ones the splice
    graph is built from, with alignment noise already merged -- so a gap here is an
    intron rather than an indel artifact. Raw CIGAR gaps would need the small-gap merge
    applied first and must not be passed here.

    Returns () for a single-segment (monoexonic) alignment.
    """
    ordered = sorted(segments, key=lambda s: s[0])
    return tuple(
        (ordered[i - 1][1] + 1, ordered[i][0] - 1) for i in range(1, len(ordered))
    )
