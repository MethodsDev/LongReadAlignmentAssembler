#!/usr/bin/env python3

import sys, os, re
import subprocess
import argparse
import hashlib
import logging
import pysam
from array import array
from collections import defaultdict
from statistics import median
sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../pylib"]))
from Pipeliner import Pipeliner, Command
import LRAA_Globals
import Util_funcs

# Named in every artifact this writes, so a cache from a different method cannot
# be mistaken for a current one. Defined in LRAA_Globals because the driver keys
# its own cache on it too.
METHOD = LRAA_Globals.SPLICE_GRAPH_NORMALIZATION_METHOD

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)




def main():

    parser = argparse.ArgumentParser(description="normalize bam in a strand-specific manner", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--input_bam", type=str, required=True, help="input bam filename")
    parser.add_argument("--output_bam", type=str, required=True, help="output for normalied bam file")
    parser.add_argument("--normalize_max_cov_level", type=int, default=1000, help="target read depth; coverage above this is thinned, coverage below it is left alone (default: 1000)")
    parser.add_argument("--depth_window", type=int, default=100, help="resolution in bases at which read depth is measured")
    parser.add_argument("--random_seed", type=int, default=42, help="random seed for reproducible sampling (default: 42)")
    parser.add_argument("--max_intron_length", type=int, default=LRAA_Globals.config["max_intron_length"], help="alignments containing any intron longer than this are discarded during the strand split; set to 0 or a negative value to disable intron length filtering")
    parser.add_argument("--input_is_single_strand", action="store_true", help="the input is already orientation-pure (one strand only), as a per-chunk input from the partitioned pipeline is: normalize it directly, skipping the strand split and the merge that would otherwise write one populated bam plus one empty one and glue them back together")
    parser.add_argument("--window_origin", type=int, default=None, help="absolute reference coordinate (0-based) that position 0 of this input maps to: 0 for a whole-contig bam, the rebase offset for a rebased chunk. Aligns the depth-window grid to the absolute coordinate grid, so the same absolute locus lands in the same window whether it is normalized whole or as part of a chunk. Default (unset): the grid is anchored per contig on the first aligned base seen there, making window boundaries a function of which records happen to be in the input")

    args = parser.parse_args()
    
    input_bam_filename = args.input_bam
    output_bam_filename = args.output_bam
    normalize_max_cov_level = args.normalize_max_cov_level
    depth_window = args.depth_window
    random_seed = args.random_seed
    input_is_single_strand = args.input_is_single_strand
    window_origin = args.window_origin
    if window_origin is not None and window_origin < 0:
        parser.error("--window_origin is an absolute reference coordinate and cannot be negative")
    # 0 and any negative value all mean "no intron filtering"; canonicalizing to
    # 0 keeps them from hashing to different checkpoint tokens for one behavior.
    max_intron_length = max(args.max_intron_length, 0)

    # sampling is keyed on a hash of the read name, so a read's fate depends on
    # neither its position in the file nor the order reads are visited in
    logger.info(f"Using random seed: {random_seed}")

    split_token, run_token = compute_tokens(
        input_bam_filename,
        output_bam_filename,
        max_intron_length,
        normalize_max_cov_level,
        depth_window,
        random_seed,
        window_origin=window_origin,
        input_is_single_strand=input_is_single_strand,
    )

    pipeliner = Pipeliner("__chckpts")

    # what to normalize, where it goes, and a checkpoint naming what produced it
    norm_jobs = list()

    if input_is_single_strand:
        # The input carries one orientation already, so the split would write a
        # copy of it plus an empty bam and the merge would glue those back
        # together: one pass over every record and two files, for nothing.
        logger.info(
            "input declared single-strand: normalizing {} directly, no strand split and no merge".format(
                input_bam_filename
            )
        )
        if max_intron_length > 0:
            # said out loud because the intron filter lives in the split, and the
            # split is precisely what this path skips
            logger.info(
                "--max_intron_length {} is enforced by the strand split only, which this path skips: long-intron alignments must already have been removed upstream".format(
                    max_intron_length
                )
            )
        norm_jobs.append(
            (
                input_bam_filename,
                output_bam_filename,
                # the normalized bam is the final output here, so it cannot carry
                # the run token in its own name; the checkpoint carries it instead
                f"{output_bam_filename}.norm_{normalize_max_cov_level}.{run_token}.ok",
            )
        )

    else:
        # first separate input bam into strand-specific files
        SS_output_prefix = f"{os.path.basename(input_bam_filename)}.{split_token}.SS"

        scriptdir = os.path.abspath(os.path.dirname(__file__))
        cmd = " ".join([os.path.join(scriptdir, "separate_bam_by_strand.py"),
                        "--bam {}".format(input_bam_filename),
                        "--output_prefix {}".format(SS_output_prefix),
                        "--max_intron_length {}".format(max_intron_length)])


        pipeliner.add_commands([Command(cmd, f"sep_by_strand.{split_token}.ok")])

        pipeliner.run()

        for SS_bam_file in [SS_output_prefix + x for x in (".+.bam", ".-.bam")]:
            norm_bam_filename = f"{SS_bam_file}.norm_{normalize_max_cov_level}.{run_token}.bam"
            norm_jobs.append((SS_bam_file, norm_bam_filename, norm_bam_filename + ".ok"))

    ## run normalizations
    norm_bam_files = list()

    for source_bam_file, norm_bam_filename, norm_bam_checkpoint in norm_jobs:

        if not os.path.exists(norm_bam_checkpoint):
            sift_bam(source_bam_file, norm_bam_filename, normalize_max_cov_level, depth_window, random_seed, window_origin)
            subprocess.check_call("touch {}".format(norm_bam_checkpoint), shell=True)

        norm_bam_files.append(norm_bam_filename)

    if input_is_single_strand:
        # sift_bam already wrote the final output; there is nothing to merge, and
        # a checkpoint named for a merge that never ran would misdescribe it
        index_checkpoint = f"index_single_strand_norm.{run_token}.ok"
    else:
        # merge the norm SS bam filenames into the final output file
        cmd = f"samtools merge -f {output_bam_filename} " + " ".join(norm_bam_files)
        pipeliner.add_commands([Command(cmd, f"SS_merge.{run_token}.ok")])
        index_checkpoint = f"index_merged.{run_token}.ok"

    cmd = f"samtools index {output_bam_filename}"
    pipeliner.add_commands([Command(cmd, index_checkpoint)])

    pipeliner.run()

    logger.info("Done.  See SS-normalized bam: {}".format(output_bam_filename))

    sys.exit(0)


def compute_tokens(
    input_bam_filename,
    output_bam_filename,
    max_intron_length,
    normalize_max_cov_level,
    depth_window,
    random_seed,
    window_origin=None,
    input_is_single_strand=False,
):
    """Checkpoint tokens for the two stages of this script.

    Every artifact here is checkpointed, and a checkpoint is trusted on sight,
    so its name has to carry whatever determines the contents. Two tokens,
    because the two stages depend on different things.
    """

    # The strand split depends on which bam came in -- and only on its basename
    # otherwise, so two inputs of the same name from different directories would
    # collide here.
    source_token = Util_funcs.file_identity_token(input_bam_filename)

    # The window origin moves the depth-window boundaries, so it decides which
    # reads survive sampling; --input_is_single_strand decides whether the split
    # stage runs at all and so what stands in for its output. Both therefore have
    # to name every artifact downstream of them. They are appended only when set
    # so that a default invocation hashes the same string it hashed before these
    # options existed and keeps hitting its existing caches. Default and non-
    # default cannot collide, because a non-default run appends text no default
    # run ever produces.
    variant_fields = list()
    if input_is_single_strand:
        variant_fields.append("single_strand")
    if window_origin is not None:
        variant_fields.append("window_origin={}".format(window_origin))

    # The split also drops the records failing the intron length criterion, so
    # the threshold determines which records it emits and has to name the split's
    # outputs and checkpoint alongside the input identity.
    split_token = Util_funcs.get_hash_code(
        "|".join([source_token, str(max_intron_length)] + variant_fields)
    )[:12]

    # Sampling, merge and index consume the split's output, so they depend on the
    # threshold too, plus every flag that changes which reads are kept or where
    # they land. The driver holds the window and seed fixed, but both are exposed
    # on this command line.
    run_token = Util_funcs.get_hash_code(
        "|".join(
            [
                source_token,
                str(max_intron_length),
                str(normalize_max_cov_level),
                str(depth_window),
                str(random_seed),
                METHOD,
                os.path.realpath(output_bam_filename),
            ]
            + variant_fields
        )
    )[:12]

    return split_token, run_token


def _read_variate(read_name, random_seed):
    """A uniform [0,1) draw determined by the read name, not by visit order."""
    digest = hashlib.blake2b(
        "{}:{}".format(random_seed, read_name).encode("utf-8"), digest_size=8
    ).digest()
    return int.from_bytes(digest, "big") / float(1 << 64)


def _window_span(lend, rend, depth_window, anchor):
    """Windows touched by the half-open interval [lend, rend).

    Offsets are measured from `anchor`, so window k covers [anchor + k*W,
    anchor + (k+1)*W). Where the anchor comes from decides what the partition is
    a function of, and there are only two sane choices:

    the input's first aligned base (the default, and all this used to do)
        Translating a locus translates the anchor with it, so a read falls in the
        same window as before and the depth estimate is unchanged. But the
        boundaries are then a function of which records are in the input:
        normalize a chunk of a contig and its first read is not the contig's
        first read, so the same locus bins differently than it did whole.

    the caller's coordinate grid (see `window_origin` on sift_bam)
        Boundaries sit at absolute multiples of W regardless of contents, so the
        same absolute locus lands in the same window whether it is normalized
        whole or inside any chunk. That equivalence is the reason normalization
        can move inside a chunk at all.
    """
    return range((lend - anchor) // depth_window, (rend - 1 - anchor) // depth_window + 1)


def _reject_read_before_anchor(read, contig, anchor, bam_file):
    """Refuse a record the default anchor cannot bin.

    With no `window_origin` the grid is anchored on the first aligned base seen
    on a contig, which is that contig's minimum only if the input is
    coordinate-sorted. Out of order, a read starting before the anchor yields a
    negative window -- which Python resolves against the array tail, corrupting
    an unrelated locus rather than failing. Skipping the read instead would
    undercount depth just as silently, so the precondition is enforced.

    Checked against the records rather than the header's `SO`, which costs one
    comparison per read and is not something a mislabelled bam can talk its way
    past. Given an explicit origin the anchor is at or below zero and no aligned
    position can precede it, so this never fires and the ordering of the input
    stops mattering.
    """
    if read.reference_start >= anchor:
        return
    raise ValueError(
        "{}: depth windows on {} are anchored at {}, the first aligned base seen "
        "there, but read {} starts at {}. This requires coordinate-sorted input. "
        "Sort the bam, or pass --window_origin to anchor on an absolute grid "
        "whose boundaries do not depend on visit order.".format(
            bam_file, contig, anchor, read.query_name, read.reference_start
        )
    )


def _grow_windows(bases, window):
    """Extend a per-contig depth array until `window` is addressable.

    The initial size comes from the header's reference length, but an alignment
    may run past it -- coordinate-remapped bams routinely carry introns longer
    than the contig they were rebased onto. Growing keeps every window a read
    occupies measurable, so pass 1 and pass 2 partition the same reads over the
    same grid. Doubling keeps a monotonically advancing scan amortized O(1).
    """
    target = max(window + 1, 2 * len(bases))
    bases.extend(array("q", bytes(8 * (target - len(bases)))))


def sift_bam(
    SS_bam_file,
    norm_bam_filename,
    normalize_max_cov_level,
    depth_window,
    random_seed,
    window_origin=None,
):
    """Thin coverage toward a target depth, recording each read's sampling weight.

    Two sequential passes over the strand-specific bam:

      pass 1  measure read depth per window and count junction support exactly
      pass 2  keep each read with probability p = min(1, target / local_depth)
              and record 1/p in the XW tag

    Acceptance follows local depth, so reads are retained untouched wherever
    coverage already sits below the target and thinned only where it does not.
    That is what makes this a normalization rather than a downsampling.

    A read carrying a junction supported by fewer than the target number of
    reads is always kept. Such a junction cannot be why depth is over target,
    and retaining it in full makes its support exact rather than estimated --
    which matters because the splice graph decides whether to keep a junction by
    comparing its support against the strongest junction sharing its exon
    island.

    Recording 1/p lets a consumer recover unbiased counts by summing weights
    instead of reads. Some such correction is required of any scheme whose
    acceptance rate varies along the genome. Sampling read starts per bin, as
    this did previously, varies just as much and recorded nothing: a junction at
    a crowded TSS could only be supported by reads starting in the one bin that
    was sampled hardest, while a junction further into the gene collected reads
    from many sparse bins that survived whole. Measured on a locus at 15,000x,
    that deflated the 5'-most junction's apparent frequency by 30% and dropped a
    transcript supported by 147 full-length reads.

    `window_origin` is the absolute reference coordinate that position 0 of this
    input maps to: 0 for a whole-contig bam, the rebase offset for a chunk
    carrying rebased coordinates. Supplying it puts the window boundaries at
    absolute multiples of `depth_window`, so a chunk measures depth over the same
    windows a whole-contig run measures and a read at a shared locus faces the
    same acceptance probability in both. Left unset, the grid is anchored per
    contig on the first aligned base seen there, exactly as this always did --
    which is correct for a whole input and wrong for a chunk of one, because the
    chunk's first read is not the contig's first read.
    """

    window_bases = dict()  # contig -> aligned bases per window
    contig_anchor = dict()  # contig -> window origin, in this input's coordinates
    junction_support = defaultdict(int)

    # Only the origin's position within a window matters: shifting the grid by a
    # whole number of windows renumbers the windows without moving any boundary.
    # Reducing it keeps window indices non-negative and the per-contig array the
    # size it always was, even for a chunk rebased a hundred megabases in.
    grid_anchor = None if window_origin is None else -(window_origin % depth_window)
    if grid_anchor is not None:
        logger.info(
            "depth windows anchored on the caller's grid: origin {} (offset {} within a {} bp window)".format(
                window_origin, -grid_anchor, depth_window
            )
        )

    with pysam.AlignmentFile(SS_bam_file, "rb") as reader:
        contig_length = dict(zip(reader.references, reader.lengths))

        overhanging = 0
        max_overhang = 0

        for read in reader.fetch(until_eof=True):
            if read.is_unmapped:
                continue

            contig = read.reference_name
            bases = window_bases.get(contig)
            if bases is None:
                if grid_anchor is None:
                    # the contig's minimum aligned position, given sorted input;
                    # _reject_read_before_anchor holds every later read to it
                    contig_anchor[contig] = read.reference_start
                else:
                    contig_anchor[contig] = grid_anchor
                bases = array("q", bytes(8 * (contig_length[contig] // depth_window + 2)))
                window_bases[contig] = bases
            anchor = contig_anchor[contig]
            _reject_read_before_anchor(read, contig, anchor, SS_bam_file)

            overhang = (read.reference_end or 0) - contig_length[contig]
            if overhang > 0:
                overhanging += 1
                max_overhang = max(max_overhang, overhang)

            for block_lend, block_rend in read.get_blocks():
                for window in _window_span(block_lend, block_rend, depth_window, anchor):
                    if window >= len(bases):
                        _grow_windows(bases, window)
                    window_lend = anchor + window * depth_window
                    covered_lend = max(block_lend, window_lend)
                    covered_rend = min(block_rend, window_lend + depth_window)
                    bases[window] += covered_rend - covered_lend

            for junction in _read_junctions(read):
                junction_support[(contig, junction)] += 1

    if overhanging:
        logger.warning(
            "{} alignment(s) extend past the reference length declared in the bam header, "
            "by up to {} bp. Their depth is measured over the windows they actually occupy; "
            "the records themselves are left untouched.".format(overhanging, max_overhang)
        )

    logger.info(
        "measured depth over {} contig(s), {} distinct junction(s)".format(
            len(window_bases), len(junction_support)
        )
    )

    kept = 0
    total = 0

    with pysam.AlignmentFile(SS_bam_file, "rb") as reader:
        with pysam.AlignmentFile(norm_bam_filename, "wb", template=reader) as writer:
            for read in reader.fetch(until_eof=True):
                if read.is_unmapped:
                    continue
                total += 1

                contig = read.reference_name
                probability = _acceptance_probability(
                    read,
                    contig,
                    window_bases[contig],
                    contig_anchor[contig],
                    junction_support,
                    normalize_max_cov_level,
                    depth_window,
                )

                if probability < 1.0 and _read_variate(read.query_name, random_seed) >= probability:
                    continue

                read.set_tag("XW", 1.0 / probability, value_type="f")
                writer.write(read)
                kept += 1

    logger.info(
        "{}: retained {} of {} reads targeting depth {}".format(
            os.path.basename(norm_bam_filename), kept, total, normalize_max_cov_level
        )
    )

    return


def _read_junctions(read):
    """Introns implied by the alignment, as half-open reference intervals."""
    junctions = []
    position = read.reference_start
    for code, length in read.cigartuples or ():
        if code in (0, 2, 7, 8):  # M D = X consume reference
            position += length
        elif code == 3:  # N
            junctions.append((position, position + length))
            position += length
    return junctions


def _acceptance_probability(
    read,
    contig,
    window_bases,
    anchor,
    junction_support,
    normalize_max_cov_level,
    depth_window,
):
    if normalize_max_cov_level <= 0:
        return 1.0

    for junction in _read_junctions(read):
        if junction_support[(contig, junction)] < normalize_max_cov_level:
            # scarce junction: keep every read carrying it so its support is exact
            return 1.0

    depths = []
    for block_lend, block_rend in read.get_blocks():
        for window in _window_span(block_lend, block_rend, depth_window, anchor):
            if not 0 <= window < len(window_bases):
                # pass 1 walked these same reads with this same anchor and grew
                # the array to cover every window they touch, so an unaddressable
                # window here means the two passes disagree. Skipping would drop
                # this read's depth from its own acceptance decision silently.
                raise AssertionError(
                    "{}: window {} outside the {} measured for it in pass 1 "
                    "(read {} block [{}, {}), anchor {})".format(
                        contig, window, len(window_bases), read.query_name,
                        block_lend, block_rend, anchor,
                    )
                )
            depths.append(window_bases[window] / depth_window)
    if not depths:
        return 1.0

    # the median keeps a read that merely clips a narrow peak from being judged
    # by that peak
    local_depth = median(depths)
    if local_depth <= normalize_max_cov_level:
        return 1.0

    return normalize_max_cov_level / local_depth
    
        
        

    
if __name__=='__main__':
    main()
