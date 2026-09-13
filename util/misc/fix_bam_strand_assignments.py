#!/usr/bin/env python3

"""Rewrite read alignment strands to the transcribed orientation.

The bam counterpart to fix_gtf_strand_assignments.py.  An aligner sets a read's
reverse bit from which genomic strand the read SEQUENCE matched, which for an
unstranded cDNA library says nothing about which strand was TRANSCRIBED.  The
introns do say it: a canonical donor/acceptor pair reads GT..AG on the
transcribed strand and CT..AC on the other one, so a spliced read carries its
own transcribed orientation regardless of how the aligner happened to place it.

A read's own splice motifs are the only direct evidence and are used first.
When they do not decide -- the read is unspliced, its introns are noncanonical,
or its canonical introns split evenly -- and a --gtf is supplied, the read is
given the orientation of the annotated transcript it shares the most exonic
bases with.  That is weaker evidence: it reports where the read landed rather
than what it is, so it is never allowed to overrule the motifs.

Every input record is written exactly once, in input order.  A record is left
alone unless the evidence names the orientation opposite its flag; then bit
0x10 is set to the transcribed orientation and the ORIGINAL aligned strand is
recorded in a tag (XD by default), so the change is both visible and
reversible.

What a flip does and does not mean:

  - SEQ, QUAL, CIGAR and POS are untouched, so the alignment stays valid
    against the genome and the file stays coordinate sorted.  Only the claim
    the flag makes about the source molecule changes.
  - After this tool, 0x10 means "the transcript is on the minus strand", not
    "SEQ is the reverse complement of the read as sequenced".  Anything that
    reverse complements SEQ by the flag to recover the original read must use
    the XD tag instead.
  - minimap2's ts:A is read-relative, so flipping the flag without flipping ts
    would change the genomic orientation ts encodes.  ts is flipped alongside,
    leaving that orientation where it was.
  - SA:Z strand fields describe the aligner's mapping and are left as found.
"""

import sys, os, re
import argparse
import gzip
import logging
from collections import defaultdict

import intervaltree as itree
import pysam

sys.path.insert(
    0,
    os.path.sep.join(
        [os.path.dirname(os.path.realpath(__file__)), "..", "..", "pylib"]
    ),
)

import Util_funcs
from GenomeFeature import Intron
from Pretty_alignment import Pretty_alignment


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


BAM_CREF_SKIP = 3  # the cigar N operation, the only source of an intron

DEFAULT_FLIP_TAG = "XD"

# Every record lands in exactly one of these, so they sum to num_records.
OUTCOMES = (
    "unchanged_unmapped",
    "unchanged_paired",
    "strand_agree",
    "strand_flipped",
    "strand_uncertain",
)

# Why a record's own splice motifs did not decide it.  These sum to the number
# of records the annotation, when supplied, was asked about.
SPLICE_UNDECIDED_REASONS = (
    "unspliced",
    "no_canonical_intron",
    "conflicting_introns",
)


def main():

    parser = argparse.ArgumentParser(
        description="fix bam read alignment strands to the transcribed orientation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input_bam", type=str, required=True, help="input bam filename"
    )

    parser.add_argument(
        "--output_bam",
        type=str,
        required=True,
        help="output bam file, strands set to the transcribed orientation",
    )

    parser.add_argument(
        "--genome", type=str, required=True, help="genome fasta file (faidx indexed)"
    )

    parser.add_argument(
        "--gtf",
        type=str,
        required=False,
        default=None,
        help="optional annotation; a read whose own splice motifs do not decide "
        "it takes the strand of the transcript it most overlaps",
    )

    parser.add_argument(
        "--flip_tag",
        type=str,
        default=DEFAULT_FLIP_TAG,
        help="tag set on flipped records, holding the original aligned strand",
    )

    args = parser.parse_args()

    if len(args.flip_tag) != 2:
        sys.exit("Error - --flip_tag must be a two character bam tag")

    contig_to_exon_itree = None
    if args.gtf is not None:
        contig_to_exon_itree = build_contig_exon_itrees(args.gtf)

    counters = fix_bam_strand_assignments(
        args.input_bam,
        args.output_bam,
        args.genome,
        contig_to_exon_itree,
        args.flip_tag,
    )

    report_counters(counters)

    index_if_coordinate_sorted(args.output_bam)

    logger.info("Done.")

    sys.exit(0)


def new_counters():

    counters = defaultdict(int)

    counters["num_records"] = 0
    counters["num_records_spliced"] = 0
    counters["num_inferred_by_splice_dinucs"] = 0
    counters["num_inferred_by_annot_overlap"] = 0
    counters["num_ts_tags_flipped"] = 0
    for outcome in OUTCOMES:
        counters["num_records_{}".format(outcome)] = 0
    for reason in SPLICE_UNDECIDED_REASONS:
        counters["num_splice_undecided_{}".format(reason)] = 0

    return counters


def fix_bam_strand_assignments(
    input_bam_filename,
    output_bam_filename,
    genome_fasta,
    contig_to_exon_itree=None,
    flip_tag=DEFAULT_FLIP_TAG,
):
    """Stream every record, flipping the strand bit where the evidence disagrees.

    Records are handled one at a time and written immediately, so memory does
    not grow with the bam.  Bam_alignment_extractor is deliberately not the
    intake here: it DISCARDS records (unmapped, secondary, supplementary, low
    identity), materializes a list, and lightens away the pysam record a
    rewrite has to write back out.  A bam repair tool that dropped reads would
    quietly shrink the library, so the pysam iterator is the intake and
    Pretty_alignment is built per record for its intron structure alone.
    """

    counters = new_counters()

    contig_seq = None
    contig_seq_acc = None
    contigs_loaded = set()
    warned_paired = False

    bamfile_reader = pysam.AlignmentFile(input_bam_filename, "rb")
    bamfile_writer = pysam.AlignmentFile(
        output_bam_filename, "wb", template=bamfile_reader
    )

    for read in bamfile_reader:

        counters["num_records"] += 1

        if read.is_unmapped or read.reference_id < 0:
            counters["num_records_unchanged_unmapped"] += 1
            bamfile_writer.write(read)
            continue

        if read.is_paired:
            # flipping 0x10 here would leave the mate's 0x20 on the old strand,
            # and this tool never holds both mates at once
            if not warned_paired:
                warned_paired = True
                logger.warning(
                    "paired records found (eg. {}) and left untouched: flipping one "
                    "mate's 0x10 would strand the other mate's 0x20".format(
                        read.query_name
                    )
                )
            counters["num_records_unchanged_paired"] += 1
            bamfile_writer.write(read)
            continue

        contig_acc = read.reference_name
        alignment_segments = None
        transcribed_orient = "?"

        if is_spliced_cigar(read.cigartuples):

            if contig_acc != contig_seq_acc:
                if contig_acc in contigs_loaded:
                    logger.warning(
                        "contig {} revisited: input is not grouped by contig, so "
                        "the genome sequence is being re-read".format(contig_acc)
                    )
                contig_seq = Util_funcs.retrieve_contig_seq_from_fasta_file(
                    contig_acc, genome_fasta
                )
                contig_seq_acc = contig_acc
                contigs_loaded.add(contig_acc)

            pretty_alignment = Pretty_alignment.get_pretty_alignment(read)
            alignment_segments = pretty_alignment.get_pretty_alignment_segments()
            introns = pretty_alignment.get_introns()

            if introns:
                counters["num_records_spliced"] += 1
                transcribed_orient, reason = infer_orient_via_splice_motifs(
                    introns, contig_seq
                )
                if transcribed_orient != "?":
                    counters["num_inferred_by_splice_dinucs"] += 1
            else:
                # every N was shorter than read_aln_gap_merge_int and merged back
                reason = "unspliced"
        else:
            reason = "unspliced"

        if transcribed_orient == "?":

            counters["num_splice_undecided_{}".format(reason)] += 1

            if contig_to_exon_itree is not None:
                if alignment_segments is None:
                    alignment_segments = unspliced_alignment_segments(read)
                transcribed_orient = infer_orient_via_annotation(
                    alignment_segments, contig_to_exon_itree.get(contig_acc)
                )
                if transcribed_orient != "?":
                    counters["num_inferred_by_annot_overlap"] += 1

        if transcribed_orient == "?":
            counters["num_records_strand_uncertain"] += 1
        else:
            aligned_orient = "-" if read.is_reverse else "+"
            if transcribed_orient == aligned_orient:
                counters["num_records_strand_agree"] += 1
            else:
                flip_record_strand(read, aligned_orient, flip_tag, counters)
                counters["num_records_strand_flipped"] += 1

        bamfile_writer.write(read)

    bamfile_writer.close()
    bamfile_reader.close()

    return counters


def is_spliced_cigar(cigartuples):
    """Whether this cigar can hold an intron at all.

    Only N skips reference without the read: a deletion keeps its block
    adjacent to its neighbors, so Pretty_alignment merges it back into the
    surrounding exon and no intron ever comes of it.  Checking here keeps the
    monoexonic majority from paying for a Pretty_alignment, which reads the
    query sequence to measure soft clipping this tool has no use for.
    """

    if not cigartuples:
        return False

    for opcode, _ in cigartuples:
        if opcode == BAM_CREF_SKIP:
            return True

    return False


def unspliced_alignment_segments(read):
    """The single exon block of a record with no N, without building a Pretty.

    Equivalent to what Pretty_alignment would return for such a cigar: with
    nothing to skip the reference, every block merges into one spanning the
    whole alignment, 1-based and inclusive.
    """

    return [(read.reference_start + 1, read.reference_end)]


def infer_orient_via_splice_motifs(introns, contig_seq):
    """Majority orientation among this read's canonical introns.

    Intron.check_canonical_splicing is the same classifier the splice graph
    admits junctions with, so a read is oriented by the motifs LRAA already
    believes, rather than by a fourth private copy of the dinucleotide sets.

    Returns the orientation, or "?" with the reason nothing was decided.
    """

    num_top_strand = 0
    num_bottom_strand = 0

    for intron_lend, intron_rend in introns:

        orient = Intron.check_canonical_splicing(intron_lend, intron_rend, contig_seq)

        if orient == "+":
            num_top_strand += 1
        elif orient == "-":
            num_bottom_strand += 1

    if num_top_strand > num_bottom_strand:
        return "+", None

    if num_bottom_strand > num_top_strand:
        return "-", None

    if num_top_strand > 0:
        return "?", "conflicting_introns"

    return "?", "no_canonical_intron"


def build_contig_exon_itrees(gtf_file):
    """contig -> interval tree of annotated exons, each holding (strand, id).

    Exons rather than transcript spans, so a read sitting in the intron of one
    gene and the exon of another is credited to the one it is transcribed from.
    Strandless ('.') annotation names no orientation and is not loaded.
    """

    logger.info("-building exon itrees from: " + gtf_file)

    contig_to_exon_itree = dict()
    num_exons = 0

    opener = gzip.open if gtf_file.endswith(".gz") else open

    with opener(gtf_file, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue

            vals = line.rstrip().split("\t")
            if len(vals) < 9:
                continue

            if vals[2] != "exon":
                continue

            strand = vals[6]
            if strand not in ("+", "-"):
                continue

            m = re.search('transcript_id "([^"]+)"', vals[8])
            if m is None:
                raise RuntimeError(
                    "Error, couldn't extract transcript_id from line: {}".format(line)
                )

            contig_acc = vals[0]
            lend = int(vals[3])
            rend = int(vals[4])

            if contig_acc not in contig_to_exon_itree:
                contig_to_exon_itree[contig_acc] = itree.IntervalTree()

            contig_to_exon_itree[contig_acc][lend : rend + 1] = (m.group(1), strand)
            num_exons += 1

    logger.info(
        "-loaded {} exons across {} contigs.".format(
            num_exons, len(contig_to_exon_itree)
        )
    )

    return contig_to_exon_itree


def infer_orient_via_annotation(alignment_segments, exon_itree):
    """The strand of the transcript this read shares the most exonic bases with.

    Bases, and per transcript, rather than a count of overlapping exon
    records: a read lying over one exon of a five exon transcript would
    otherwise be outvoted by the five exon transcript on the other strand it
    barely touches.  A read whose best transcript on each strand overlaps it
    equally is not decided.
    """

    if exon_itree is None:
        return "?"

    bases_per_transcript = defaultdict(int)
    strand_of_transcript = dict()

    for seg_lend, seg_rend in alignment_segments:
        for exon in exon_itree[seg_lend : seg_rend + 1]:
            transcript_id, strand = exon.data
            # itree intervals are half open, so exon.end is one past its last base
            overlap = min(seg_rend, exon.end - 1) - max(seg_lend, exon.begin) + 1
            bases_per_transcript[transcript_id] += overlap
            strand_of_transcript[transcript_id] = strand

    best_bases = {"+": 0, "-": 0}
    for transcript_id, num_bases in bases_per_transcript.items():
        strand = strand_of_transcript[transcript_id]
        if num_bases > best_bases[strand]:
            best_bases[strand] = num_bases

    if best_bases["+"] > best_bases["-"]:
        return "+"

    if best_bases["-"] > best_bases["+"]:
        return "-"

    return "?"


def flip_record_strand(read, orig_aligned_orient, flip_tag, counters):
    """Set 0x10 to the transcribed orientation, keeping the record readable back.

    The original strand goes into the flip tag rather than a bare flag, because
    the flip is only reversible if the strand it replaced survives, and a
    consumer asking "was this read moved" and "where was it" should not need
    two tags.
    """

    read.is_reverse = orig_aligned_orient == "+"

    read.set_tag(flip_tag, orig_aligned_orient, value_type="A")

    if read.has_tag("ts"):
        ts = read.get_tag("ts")
        if ts in ("+", "-"):
            # ts is relative to the read, so it moves with the read's strand
            read.set_tag("ts", "-" if ts == "+" else "+", value_type="A")
            counters["num_ts_tags_flipped"] += 1

    return


def report_counters(counters):

    num_records = counters["num_records"]

    def as_pct(count):
        if num_records < 1:
            return 0.0
        return 100.0 * count / num_records

    logger.info("-records examined: {}".format(num_records))
    for counter_name in sorted(counters.keys()):
        if counter_name == "num_records":
            continue
        count = counters[counter_name]
        logger.info("-{}: {} ({:.2f}%)".format(counter_name, count, as_pct(count)))

    return


def index_if_coordinate_sorted(bam_filename):
    """Index the output when the input ordering survived, which it always does.

    A strand flip changes no coordinate, so a coordinate sorted input yields a
    coordinate sorted output.  An input in any other order cannot be indexed
    and is left as written.
    """

    with pysam.AlignmentFile(bam_filename, "rb", check_sq=False) as reader:
        sort_order = reader.header.to_dict().get("HD", {}).get("SO")

    if sort_order != "coordinate":
        logger.info("-output not indexed: header sort order is {}".format(sort_order))
        return

    pysam.index(bam_filename)
    logger.info("-output indexed.")

    return


if __name__ == "__main__":
    main()
