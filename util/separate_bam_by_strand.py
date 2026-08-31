#!/usr/bin/env python3

import sys, os, re
import argparse
import multiprocessing
import pysam
import logging
import shutil
import subprocess
import tempfile
import time
import intervaltree as itree
import gzip
from collections import defaultdict

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../pylib"])
)
from Pretty_alignment import Pretty_alignment
import Util_funcs
import LRAA_Globals

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# The intron-length rule lives in pylib/Util_funcs.py, which is the single
# implementation shared with pylib/Bam_alignment_extractor.py so that the
# splice-graph input and read assignment always see one record set.  Re-bound here
# because it is part of this module's filtering vocabulary.
get_longest_intron_length = Util_funcs.get_longest_intron_length

# reasons for discarding an input record, in the order they are tested.
# Each discarded record is counted under exactly one reason, so that the
# discard counts plus the written records sum to the records read.
DISCARD_REASONS = (
    "unmapped",
    "no_chromosome",
    "secondary",
    "supplementary",
    "duplicate",
    "qcfail",
    "long_intron",
)


def discard_counter_name(discard_reason):
    return "num_records_discarded_" + discard_reason


def get_discard_reason(read, max_intron_length):
    """returns the reason for discarding the read, or None if the read is retained.

    max_intron_length <= 0 disables intron length filtering.
    """

    if read.is_unmapped:
        return "unmapped"

    if read.reference_id < 0:
        return "no_chromosome"

    if read.is_secondary:
        return "secondary"

    if read.is_supplementary:
        return "supplementary"

    if read.is_duplicate:
        return "duplicate"

    if read.is_qcfail:
        return "qcfail"

    if Util_funcs.has_disqualifying_long_intron(read, max_intron_length):
        return "long_intron"

    return None


# the records with no coordinate at all, fetched as their own scope by the
# fan-out below.  "*" is htslib's own spelling of it, accepted by
# pysam.AlignmentFile.fetch and by `samtools view`.
UNPLACED_SCOPE = "*"


def new_counters():
    """the zeroed record accounting a single pass fills in.

    One definition, so a pass, a merge of passes and the report cannot disagree
    about which counters exist.
    """

    counters = {
        # general stats
        "num_records": 0,
        "num_forward": 0,
        "num_reverse": 0,
        "num_neither": 0,
        # infer stats
        "num_inferred_by_splice_dinucs": 0,
        "num_inferred_by_annot_overlap": 0,
        "num_records_strand_flipped": 0,
        "num_records_strand_uncertain": 0,
    }

    # discard stats, one counter per discard reason
    for discard_reason in DISCARD_REASONS:
        counters[discard_counter_name(discard_reason)] = 0

    return counters


def merge_counters(counter_dicts):
    """the accounting of several disjoint passes, as one dict.

    Every counter here counts records and the passes partition the input, so the
    merge is the ELEMENTWISE SUM -- one sum per counter, including one per
    discard reason.  The per-reason breakdown is summed rather than collapsed
    into a total because the report attributes every drop to the filter that
    made it, and a total cannot say which filter that was.

    Adopting any single pass's dict instead -- the last one to finish, say --
    would leave ``num_records`` at that pass's record count, and every
    percentage in ``build_report`` is computed against it.
    """

    merged = new_counters()

    for counters in counter_dicts:
        unexpected = sorted(set(counters) - set(merged))
        missing = sorted(set(merged) - set(counters))
        if unexpected or missing:
            raise ValueError(
                "cannot merge record accounting: unexpected counter(s) {}, "
                "missing counter(s) {}".format(unexpected, missing)
            )
        for counter_name, count in counters.items():
            merged[counter_name] += count

    return merged


def available_cpu_count():
    """the CPUs this process is actually allowed to run on."""

    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return os.cpu_count() or 1


def resolve_num_workers(requested, num_tasks=None):
    """workers to run, given the request and the work there is to divide."""

    workers = requested if requested and requested > 0 else available_cpu_count()
    if num_tasks is not None:
        workers = min(workers, num_tasks)
    return max(1, workers)


def main():

    parser = argparse.ArgumentParser(
        description="separate bam into strand-specific bam files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--bam", type=str, required=True, help="input bam filename")

    parser.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="output prefix: files named ${output_prefix}.${strand}.bam",
    )

    parser.add_argument(
        "--infer_read_orient",
        action="store_true",
        default=False,
        help="infer read orientation based on splicing evidence or overlapping reference annotations (in that order)",
    )

    parser.add_argument(
        "--gtf",
        type=str,
        required=False,
        help="reference annotation used for inferring transcribed orientation of read (required if --infer_read_orient",
    )

    parser.add_argument(
        "--genome",
        type=str,
        required=False,
        help="genome fasta file, required if --infer_read_orient",
    )

    parser.add_argument(
        "--max_intron_length",
        type=int,
        default=LRAA_Globals.config["max_intron_length"],
        help="alignments containing any intron (CIGAR 'N' operation) longer than this are discarded; set to 0 or a negative value to disable intron length filtering",
    )

    parser.add_argument(
        "--contig",
        type=str,
        default=None,
        help="restrict the split to one contig, fetched by name from an indexed "
        "bam. Written for the chunked pipeline's whole-contig chunks, whose bam "
        "is the SOURCE rather than an extracted copy and therefore still holds "
        "every other contig. A NAME rather than a region on purpose: the only "
        "restriction anyone needs here is a whole contig, and a name has no "
        "1-based/0-based boundary to get wrong. Omitted, every record in the "
        "file is read, which is what stage 1 does and what this tool has always "
        "done.",
    )

    parser.add_argument(
        "--num_workers",
        type=int,
        default=0,
        help="concurrent per-contig passes over the whole bam. The whole-file "
        "sweep is one thread no matter the machine: measured at 915 s over "
        "48.1 M records with 27 of 28 cores idle, and it does not shrink with "
        "the contig list either, because there is no contig list to shrink. The "
        "fan-out reads each reference through the index instead, in parallel, "
        "and concatenates the per-contig output in header order. 0 (the "
        "default) means one worker per CPU this process may run on; 1 means the "
        "single whole-file sweep, unchanged, and is the way out if the input "
        "cannot be indexed. Ignored with --contig, which is already one pass "
        "over one reference.",
    )

    args = parser.parse_args()

    input_bam_filename = args.bam
    output_prefix = args.output_prefix

    infer_read_orient_flag = args.infer_read_orient
    genome_fasta = args.genome
    gtf_file = args.gtf
    max_intron_length = args.max_intron_length
    contig = args.contig
    num_workers = args.num_workers

    #########
    ### begin

    chrom_to_itree = None

    if infer_read_orient_flag:
        if genome_fasta is None:
            sys.exit("Error - with --infer_read_orient need --genome")

        if gtf_file:
            chrom_to_itree = build_chrom_itrees(gtf_file)

    top_strand_bam_filename = output_prefix + ".+.bam"
    bottom_strand_bam_filename = output_prefix + ".-.bam"

    output_bam_files = (top_strand_bam_filename, bottom_strand_bam_filename)

    if max_intron_length > 0:
        logger.info(
            "-discarding alignments having any intron longer than {}".format(
                max_intron_length
            )
        )
    else:
        logger.info("-intron length filtering disabled")

    if contig:
        logger.info("-restricting the split to contig {}".format(contig))
        counters = split_bam_by_strand(
            input_bam_filename,
            top_strand_bam_filename,
            bottom_strand_bam_filename,
            max_intron_length,
            infer_read_orient_flag=infer_read_orient_flag,
            genome_fasta=genome_fasta,
            chrom_to_itree=chrom_to_itree,
            contig=contig,
        )

    elif num_workers == 1:
        logger.info("-one worker requested: single whole-file pass")
        counters = split_bam_by_strand(
            input_bam_filename,
            top_strand_bam_filename,
            bottom_strand_bam_filename,
            max_intron_length,
            infer_read_orient_flag=infer_read_orient_flag,
            genome_fasta=genome_fasta,
            chrom_to_itree=chrom_to_itree,
        )

    else:
        if not ensure_bam_index(
            input_bam_filename, threads=resolve_num_workers(num_workers)
        ):
            sys.exit(
                "Error - {0} has no index and one could not be built, so the "
                "per-contig fan-out cannot address it: it reaches a reference "
                "by name through the index. Index it where it lives "
                "(samtools index {0}), stage it somewhere writable first, or "
                "pass --num_workers 1 to take the single-threaded whole-file "
                "pass deliberately.".format(input_bam_filename)
            )
        counters = split_bam_by_strand_parallel(
            input_bam_filename,
            top_strand_bam_filename,
            bottom_strand_bam_filename,
            max_intron_length,
            infer_read_orient_flag=infer_read_orient_flag,
            genome_fasta=genome_fasta,
            chrom_to_itree=chrom_to_itree,
            num_workers=num_workers,
        )

    if counters["num_records"] == 0:
        logger.warning(
            "No records detected for {}{}".format(
                input_bam_filename,
                " on contig {}".format(contig) if contig else "",
            )
        )
    else:
        logger.info(build_report(counters, infer_read_orient_flag))

    # index the bams
    index_threads = max(resolve_num_workers(num_workers) - 1, 0)
    for output_bam_file in output_bam_files:
        subprocess.check_call(
            ["samtools", "index", "-@", str(index_threads), output_bam_file]
        )

    sys.exit(0)


def split_bam_by_strand(
    input_bam_filename,
    top_strand_bam_filename,
    bottom_strand_bam_filename,
    max_intron_length,
    infer_read_orient_flag=False,
    genome_fasta=None,
    chrom_to_itree=None,
    contig=None,
):
    """writes the retained records of the input bam to the strand-specific bam files
    and returns the record accounting as a dict of counters.

    Only primary, non-supplementary alignments passing the intron length
    criterion are retained; every other input record is discarded here and so
    excluded from all downstream processing.

    ``contig`` restricts the pass to one reference, fetched from the index, or to
    the records with no coordinate at all when it is ``UNPLACED_SCOPE``. The
    counters then carry that scope's records as their denominator rather than the
    file's -- which is the honest reading, since nothing else was looked at.  A
    coordinate fetch does return the unmapped records PLACED on the reference,
    which are discarded under ``unmapped`` exactly as the whole-file sweep
    discards them; the coordinate-less ones live in ``UNPLACED_SCOPE``, so the
    two kinds of scope together account for every record in the file.  Whole-file
    iteration is unchanged and is what a single-worker run does.
    """

    bamfile_reader = pysam.AlignmentFile(input_bam_filename, "rb")

    top_strand_bamfile_writer = pysam.AlignmentFile(
        top_strand_bam_filename, "wb", template=bamfile_reader
    )
    bottom_strand_bamfile_writer = pysam.AlignmentFile(
        bottom_strand_bam_filename, "wb", template=bamfile_reader
    )

    chrom_seq = None
    prev_chrom = None

    counters = new_counters()

    reader = bamfile_reader.fetch(contig) if contig else bamfile_reader
    for read in reader:

        counters["num_records"] += 1

        discard_reason = get_discard_reason(read, max_intron_length)
        if discard_reason is not None:
            counters[discard_counter_name(discard_reason)] += 1
            continue

        chrom = bamfile_reader.get_reference_name(read.reference_id)

        strand = "+" if read.is_forward else "-"
        init_strand = strand

        if infer_read_orient_flag:

            if prev_chrom is None or prev_chrom != chrom:
                prev_chrom = chrom
                chrom_seq = Util_funcs.retrieve_contig_seq_from_fasta_file(
                    chrom, genome_fasta
                )

            pretty_alignment = Pretty_alignment.get_pretty_alignment(read)

            # first try by dinuc splice sites of spliced introns from read
            strand = infer_spliced_orient(pretty_alignment, chrom_seq)
            if strand != "?":
                counters["num_inferred_by_splice_dinucs"] += 1
            elif chrom_to_itree is not None:
                # try by annotation mapping
                strand = infer_transcribed_orient_via_annotation_mapping(
                    pretty_alignment, chrom_to_itree[chrom]
                )
                if strand != "?":
                    counters["num_inferred_by_annot_overlap"] += 1

            # fix strand setting in the read alignment record
            if strand != "?" and strand != init_strand:
                read.is_reverse = True if strand == "-" else False
                counters["num_records_strand_flipped"] += 1

        if strand == "?":
            counters["num_records_strand_uncertain"] += 1
            # set to aligned orientation
            strand = init_strand

        # write strand-specific records
        if strand == "+":
            top_strand_bamfile_writer.write(read)
            counters["num_forward"] += 1

        elif strand == "-":
            bottom_strand_bamfile_writer.write(read)
            counters["num_reverse"] += 1

        else:
            counters["num_neither"] += 1

    top_strand_bamfile_writer.close()
    bottom_strand_bamfile_writer.close()
    bamfile_reader.close()

    return counters


def ensure_bam_index(bam_filename, threads=1):
    """makes sure the input carries an index; returns whether one is available.

    The fan-out reaches a reference by NAME, which is an index lookup, so the
    index is load-bearing here rather than incidental.  Built rather than merely
    demanded because neither whole-bam caller guarantees one: stage 1 of the
    chunked pipeline hands over the raw input, and a bam localized by a WDL task
    carries no ``.bai`` at all unless the task declared one.  Building it is one
    loud pass; refusing outright would send both of them back to the
    single-threaded sweep this exists to remove.

    Returns False rather than raising when the index cannot be written -- an
    input in a read-only directory is the honest case -- and leaves the caller to
    say what that costs.
    """

    with pysam.AlignmentFile(bam_filename, "rb") as reader:
        if reader.has_index():
            return True

    logger.info(
        "-{} has no index; building one, which the per-contig fan-out needs to "
        "address a reference by name".format(bam_filename)
    )
    try:
        subprocess.check_call(
            ["samtools", "index", "-@", str(max(threads - 1, 0)), bam_filename]
        )
    except (subprocess.CalledProcessError, OSError) as err:
        logger.warning("-could not index {}: {}".format(bam_filename, err))
        return False

    return True


def require_coordinate_sorted(bam_filename):
    """refuses an input whose header says it is not in coordinate order.

    The fan-out's output is coordinate-sorted BY CONSTRUCTION and only by
    construction: each part is one reference's records in the order an index
    iterator returns them, which is ascending position for a coordinate-sorted
    file, and the parts are concatenated in header order.  Both halves of that
    rest on the input being coordinate-sorted -- from an unsorted input the parts
    themselves would be unsorted, and nothing later would say so.

    Sortedness matters to the consumers, not just to samtools: cut selection
    fetches these bams by region (util/misc/select_contig_cut_points.py), where a
    misplaced record is simply absent from the answer, and the depth
    normalization that reads them anchors each contig's window grid on the first
    aligned base it sees there (normalize_bam_by_strand.py:_reject_read_before_anchor).

    Only an explicit contradiction is refused.  A header carrying no ``SO`` at
    all is common and says nothing either way; there the index stands in for the
    claim, since samtools will not build one over an unsorted bam, and the
    ``samtools index`` of the OUTPUT that ends every run rejects a concatenation
    that came out in the wrong order.
    """

    with pysam.AlignmentFile(bam_filename, "rb") as reader:
        sort_order = reader.header.to_dict().get("HD", {}).get("SO")

    if sort_order is not None and sort_order != "coordinate":
        raise RuntimeError(
            "{} declares sort order {!r}, and the per-contig fan-out reads it "
            "through the index and concatenates the parts in header order, which "
            "is only coordinate order for a coordinate-sorted input. Sort it "
            "(samtools sort) or pass --num_workers 1 for the whole-file "
            "pass, which writes records in the order it reads "
            "them.".format(bam_filename, sort_order)
        )


def index_record_counts(bam_filename):
    """``[(scope, num_records)]`` in header order, plus the file's record total.

    The scopes PARTITION the file: each reference carries its own records --
    mapped, and the unmapped ones placed on it, which a coordinate fetch does
    return -- and ``UNPLACED_SCOPE`` carries the records with no coordinate at
    all.  Nothing else can exist in a bam, so summing the scopes reproduces the
    file's record count, which is what the returned total is checked against
    here before any work starts.

    Counts come from the index, so a reference holding nothing costs nothing to
    skip -- a whole-genome header carries hundreds of them -- and the caller can
    start the biggest scope first.  Order is HEADER order, which for a
    coordinate-sorted bam is also the order the output must be concatenated in.
    """

    with pysam.AlignmentFile(bam_filename, "rb") as reader:
        stats = reader.get_index_statistics()
        num_unplaced = reader.nocoordinate
        num_records = reader.mapped + reader.unmapped
        references = list(reader.references)

    per_reference = {entry.contig: entry.total for entry in stats}

    scopes = [(ref, per_reference[ref]) for ref in references if per_reference.get(ref)]
    if num_unplaced:
        # last, which is where a coordinate-sorted bam keeps them
        scopes.append((UNPLACED_SCOPE, num_unplaced))

    accounted = sum(count for _, count in scopes)
    if accounted != num_records:
        raise RuntimeError(
            "the index of {} accounts for {} record(s) across {} scope(s) but "
            "reports {} in the file. The fan-out would read a different record "
            "set than a whole-file pass, so it refuses to start; reindex the bam "
            "(samtools index) or pass --num_workers 1.".format(
                bam_filename, accounted, len(scopes), num_records
            )
        )

    return scopes, num_records


# set once per worker process, so a fork inherits the annotation itrees rather
# than pickling them per task
_split_worker_context = dict()


def _init_split_worker(context):
    _split_worker_context.update(context)


def _split_one_scope(task):
    return split_scope(_split_worker_context, task)


def split_scope(context, task):
    """one scope's pass, as the worker runs it."""

    started = time.time()
    counters = split_bam_by_strand(
        context["input_bam_filename"],
        task["+"],
        task["-"],
        context["max_intron_length"],
        infer_read_orient_flag=context["infer_read_orient_flag"],
        genome_fasta=context["genome_fasta"],
        chrom_to_itree=context["chrom_to_itree"],
        contig=task["scope"],
    )
    logger.info(
        "-split {}: {} record(s) in {:.1f}s".format(
            task["scope"], counters["num_records"], time.time() - started
        )
    )

    return task["scope"], counters


def verify_part_order(tasks, part_key, label=None):
    """refuses to concatenate parts that would not come out coordinate-sorted.

    Costs one record read per part rather than a pass over the output.  That is
    enough because of what a part IS: one reference's records in the order an
    index iterator returned them, so its first record's reference names the whole
    part and the records within it already ascend.  Parts whose references ascend
    therefore concatenate to a coordinate-sorted bam.

    This check has to exist because ``samtools index`` does NOT catch the failure
    it guards.  Measured: an index over a two-reference bam whose references
    appear in reverse header order is built with exit 0 and no warning, and the
    resulting file then looks fine to a region query while being, in file order,
    unsorted -- so a merge or any streaming consumer downstream reads it wrongly.
    Within a reference htslib does refuse ("Unsorted positions on sequence #N"),
    which is the half of the property this design cannot get wrong anyway.

    ``part_key`` names the entry in each task holding that part's filename, and
    ``label`` names the SERIES being checked in the refusals.  Both are arguments
    because the same guarantee is owed by two fan-outs over the same partition:
    this split, whose parts are keyed by strand, and the per-contig depth
    normalization in normalize_bam_by_strand.py, whose units carry one part each.
    One implementation rather than two, because two copies of an ordering
    assertion drift and the drifting copy is the one that stops catching
    anything.
    """

    if label is None:
        label = part_key

    previous_reference_id = None
    previous_scope = None

    for task in tasks:
        with pysam.AlignmentFile(task[part_key], "rb") as reader:
            first = next(reader.fetch(until_eof=True), None)

        if first is None:
            # nothing of this scope survived here: no records to place
            continue

        if first.reference_id != task["reference_id"]:
            raise RuntimeError(
                "the {} part for scope {} holds a record on reference id {} "
                "rather than {}: the parts do not describe the scopes they were "
                "split from, so concatenating them would not reproduce the "
                "input's order.".format(
                    label, task["scope"], first.reference_id, task["reference_id"]
                )
            )

        if previous_reference_id is not None and first.reference_id < previous_reference_id:
            raise RuntimeError(
                "the {} parts are not in header order: scope {} (reference id "
                "{}) follows scope {} (reference id {}). Concatenating them would "
                "produce a bam that is not coordinate-sorted, which samtools "
                "would index without complaint and every consumer fetching it by "
                "region would read wrongly.".format(
                    label,
                    task["scope"],
                    first.reference_id,
                    previous_scope,
                    previous_reference_id,
                )
            )

        previous_reference_id = first.reference_id
        previous_scope = task["scope"]


def concatenate_bam_parts(output_bam_filename, part_filenames):
    """``samtools cat`` of the parts, in the order given.

    ``cat`` rather than ``merge``: every part was written from one template, so
    they share the input's header byte for byte and this is a copy of BGZF blocks
    under the first part's header.  A merge would instead re-derive the header and
    uniquify each inherited @PG record once per input, which is the pathology
    documented at the merge in normalize_bam_by_strand.py -- here it would run
    over one part per reference.

    ``-b`` because a whole-genome header can hold more references, and so more
    parts, than an argv can carry.  ``--no-PG`` because a concatenation of one
    file's own halves is not a step the header should grow a record for.
    """

    listing = output_bam_filename + ".parts.fofn"
    with open(listing, "wt") as fh:
        for part_filename in part_filenames:
            fh.write(part_filename + "\n")
    try:
        subprocess.check_call(
            ["samtools", "cat", "--no-PG", "-b", listing, "-o", output_bam_filename]
        )
    finally:
        os.remove(listing)


def split_bam_by_strand_parallel(
    input_bam_filename,
    top_strand_bam_filename,
    bottom_strand_bam_filename,
    max_intron_length,
    infer_read_orient_flag=False,
    genome_fasta=None,
    chrom_to_itree=None,
    num_workers=0,
):
    """the whole-bam split as one concurrent pass per scope, concatenated.

    Equivalent to the whole-file sweep rather than merely similar to it:

      - the scopes PARTITION the file, checked against the index's own record
        total before any work starts and against the merged accounting
        afterwards.  So no record is read twice and none is missed, unmapped
        records included: the placed ones inside their reference's scope, the
        coordinate-less ones in ``UNPLACED_SCOPE``.  Both are counted and then
        discarded under ``unmapped``, in the same place the sweep does it, which
        is what keeps ``num_records`` -- the denominator of every percentage in
        the report -- the file's record count and not something smaller.
      - the pass is the same function over the same records.  Its only
        cross-record state is ``prev_chrom``, which gates loading a contig
        sequence and resets on every change of chrom; a scope holds one chrom, so
        a worker loads once.
      - the parts are concatenated in HEADER order, never in completion order:
        the scopes are numbered by header position when the tasks are built and
        the concatenation walks that numbering, while the workers are handed the
        biggest scope first for load balance.  Each part holds one reference's
        records in the order an index iterator returns them, which is ascending
        position for the coordinate-sorted input ``require_coordinate_sorted``
        insists on, and coordinate-less records go in the last part.  So the
        concatenation is coordinate-sorted in the input's own reference order,
        which is what the sweep produced and what a consumer fetching these bams
        by region needs.  ``verify_part_order`` asserts it against the parts
        themselves before they are joined, because nothing later would: measured,
        ``samtools index`` builds an index over references in reverse header
        order with exit 0 and no warning.

    Each worker opens its own reader and its own pair of writers: a pysam file
    handle cannot be shared across processes, and the per-scope output has to be
    a separate file anyway to be concatenated in an order the completion order
    does not decide.
    """

    require_coordinate_sorted(input_bam_filename)
    scopes, expected_num_records = index_record_counts(input_bam_filename)

    if not scopes:
        # an empty bam still owes the caller two bams carrying its header
        with pysam.AlignmentFile(input_bam_filename, "rb") as reader:
            for output_bam_filename in (
                top_strand_bam_filename,
                bottom_strand_bam_filename,
            ):
                pysam.AlignmentFile(
                    output_bam_filename, "wb", template=reader
                ).close()
        return new_counters()

    num_workers = resolve_num_workers(num_workers, len(scopes))
    logger.info(
        "-splitting {} record(s) across {} scope(s) on {} worker(s)".format(
            expected_num_records, len(scopes), num_workers
        )
    )

    # alongside the output, so the parts share its filesystem and its disk budget
    part_dir = tempfile.mkdtemp(
        prefix=".{}.parts.".format(os.path.basename(top_strand_bam_filename)),
        dir=os.path.dirname(os.path.abspath(top_strand_bam_filename)),
    )

    try:
        tasks = list()
        with pysam.AlignmentFile(input_bam_filename, "rb") as reader:
            for index, (scope, num_records) in enumerate(scopes):
                # numbered rather than named: a reference name may hold characters
                # a filename may not, and the number is the header position the
                # parts are concatenated in
                tasks.append(
                    {
                        "scope": scope,
                        "num_records": num_records,
                        # -1 for the unplaced scope, which is what a
                        # coordinate-less record's reference_id is
                        "reference_id": (
                            -1
                            if scope == UNPLACED_SCOPE
                            else reader.get_tid(scope)
                        ),
                        "+": os.path.join(part_dir, "part_{:05d}.+.bam".format(index)),
                        "-": os.path.join(part_dir, "part_{:05d}.-.bam".format(index)),
                    }
                )

        context = {
            "input_bam_filename": input_bam_filename,
            "max_intron_length": max_intron_length,
            "infer_read_orient_flag": infer_read_orient_flag,
            "genome_fasta": genome_fasta,
            "chrom_to_itree": chrom_to_itree,
        }

        # biggest scope first, so the longest pass is not the one that starts last.
        # This is the ONLY place completion order is allowed to differ from header
        # order; the concatenation below walks `tasks`, not this.
        scheduled = sorted(tasks, key=lambda task: task["num_records"], reverse=True)

        started = time.time()
        per_scope_counters = list()

        if num_workers == 1:
            for task in scheduled:
                per_scope_counters.append(split_scope(context, task)[1])
        else:
            # fork inherits the annotation itrees for free; the fallback pickles
            # the context once per worker, which is still once rather than per
            # task
            start_method = (
                "fork" if "fork" in multiprocessing.get_all_start_methods() else None
            )
            pool_context = multiprocessing.get_context(start_method)
            with pool_context.Pool(
                processes=num_workers,
                initializer=_init_split_worker,
                initargs=(context,),
            ) as pool:
                for _scope, counters in pool.imap_unordered(
                    _split_one_scope, scheduled
                ):
                    per_scope_counters.append(counters)

        counters = merge_counters(per_scope_counters)

        if counters["num_records"] != expected_num_records:
            raise RuntimeError(
                "the fan-out read {} record(s) of {} where the index says {} "
                "live in the file. Every percentage in the report divides by "
                "that count, so a short read would misdescribe the split rather "
                "than fail it; refusing instead.".format(
                    counters["num_records"], input_bam_filename, expected_num_records
                )
            )

        logger.info(
            "-{} scope(s) split in {:.1f}s; concatenating".format(
                len(scopes), time.time() - started
            )
        )

        for strand, output_bam_filename in (
            ("+", top_strand_bam_filename),
            ("-", bottom_strand_bam_filename),
        ):
            verify_part_order(tasks, strand)
            concatenate_bam_parts(
                output_bam_filename, [task[strand] for task in tasks]
            )

    finally:
        shutil.rmtree(part_dir, ignore_errors=True)

    return counters


def build_report(counters, infer_read_orient_flag=False):
    """record accounting summary, every count reported against the number of
    records read as denominator"""

    num_records = counters["num_records"]

    def as_pct_of_total(count):
        return "{} = {:.1f}% of {}".format(
            count, count / num_records * 100, num_records
        )

    report_vals = ["Num input bam records: {}".format(num_records)]

    num_records_discarded = 0
    for discard_reason in DISCARD_REASONS:
        count = counters[discard_counter_name(discard_reason)]
        num_records_discarded += count
        report_vals.append(
            "Num records discarded as {}: {}".format(
                discard_reason, as_pct_of_total(count)
            )
        )

    report_vals += [
        "Num records discarded (all reasons): {}".format(
            as_pct_of_total(num_records_discarded)
        ),
        "Num records retained: {}".format(
            as_pct_of_total(num_records - num_records_discarded)
        ),
        "Num top strand: {}".format(as_pct_of_total(counters["num_forward"])),
        "Num bottom strand: {}".format(as_pct_of_total(counters["num_reverse"])),
        "Num neither strand and ignored: {}".format(
            as_pct_of_total(counters["num_neither"])
        ),
    ]

    if infer_read_orient_flag:
        report_vals += [
            "Num read orientations inferred by dinuc splice sites: {}".format(
                as_pct_of_total(counters["num_inferred_by_splice_dinucs"])
            ),
            "Num read orientations inferred by annotation overlap: {}".format(
                as_pct_of_total(counters["num_inferred_by_annot_overlap"])
            ),
            "Num read orientations flipped: {}".format(
                as_pct_of_total(counters["num_records_strand_flipped"])
            ),
            "Num reads with uncertain orientation: {}".format(
                as_pct_of_total(counters["num_records_strand_uncertain"])
            ),
        ]

    return "\n".join(report_vals)


def infer_spliced_orient(pretty_alignment, contig_seq):

    introns_coordsets = pretty_alignment.get_introns()

    if len(introns_coordsets) < 1:
        return "?"

    return majority_vote_intron_orient(introns_coordsets, contig_seq)


def majority_vote_intron_orient(intron_coordsets, contig_seq):

    splice_dinucs_top_strand = {"GTAG", "GCAG", "ATAC"}
    splice_dinucs_bottom_strand = {
        "CTAC",
        "CTGC",
        "GTAT",
    }  # revcomp of top strand dinucs

    orient_counts = {"+": 0, "-": 0}

    for intron_coordset in intron_coordsets:
        intron_lend, intron_rend = intron_coordset
        dinuc_left = contig_seq[intron_lend - 1] + contig_seq[intron_lend - 1 + 1]
        dinuc_right = contig_seq[intron_rend - 1 - 1] + contig_seq[intron_rend - 1]
        dinuc_combo = dinuc_left + dinuc_right

        if dinuc_combo in splice_dinucs_top_strand:
            orient_counts["+"] += 1
        elif dinuc_combo in splice_dinucs_bottom_strand:
            orient_counts["-"] += 1

    # check tie or not match
    if orient_counts["+"] == orient_counts["-"]:
        return "?"
    elif orient_counts["+"] > orient_counts["-"]:
        return "+"
    else:
        return "-"


def build_chrom_itrees(gtf_file):

    logger.info("-building chrom itrees from: " + gtf_file)

    chrom_to_itree = defaultdict(lambda: itree.IntervalTree())

    if re.search(gtf_file, "\\.gz"):
        opener = gzip.open
    else:
        opener = open

    with opener(gtf_file, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue

            line = line.rstrip()
            vals = line.split("\t")
            if len(vals) < 8:
                continue

            chrom = vals[0]
            lend = int(vals[3])
            rend = int(vals[4])
            feature_type = vals[2]

            if feature_type != "exon":
                continue

            strand = vals[6]

            chrom_to_itree[chrom][lend : rend + 1] = strand

    return chrom_to_itree


def infer_transcribed_orient_via_annotation_mapping(pretty_alignment, chrom_itree):

    alignment_segments = pretty_alignment.get_pretty_alignment_segments()

    strand_counter = defaultdict(int)

    for segment in alignment_segments:
        lend, rend = segment

        for overlapping_exon in chrom_itree[lend : rend + 1]:
            strand = overlapping_exon.data
            strand_counter[strand] += 1

    if strand_counter["+"] == strand_counter["-"]:
        return "?"
    elif strand_counter["+"] > strand_counter["-"]:
        return "+"
    else:
        return "-"


if __name__ == "__main__":
    main()
