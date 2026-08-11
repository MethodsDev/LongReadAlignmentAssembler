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

    args = parser.parse_args()
    
    input_bam_filename = args.input_bam
    output_bam_filename = args.output_bam
    normalize_max_cov_level = args.normalize_max_cov_level
    depth_window = args.depth_window
    random_seed = args.random_seed

    # sampling is keyed on a hash of the read name, so a read's fate depends on
    # neither its position in the file nor the order reads are visited in
    logger.info(f"Using random seed: {random_seed}")

    # Every artifact below is checkpointed, and a checkpoint is trusted on sight,
    # so its name has to carry whatever determines the contents. Two tokens,
    # because the two stages depend on different things.
    #
    # The strand split depends only on which bam came in -- and only on its
    # basename otherwise, so two inputs of the same name from different
    # directories would collide here.
    source_token = Util_funcs.file_identity_token(input_bam_filename)

    # Sampling, merge and index additionally depend on every flag that changes
    # which reads are kept or where they land. The driver holds the window and
    # seed fixed, but both are exposed on this command line.
    run_token = Util_funcs.get_hash_code(
        "|".join(
            [
                source_token,
                str(normalize_max_cov_level),
                str(depth_window),
                str(random_seed),
                METHOD,
                os.path.realpath(output_bam_filename),
            ]
        )
    )[:12]

    pipeliner = Pipeliner("__chckpts")

    # first separate input bam into strand-specific files
    SS_output_prefix = f"{os.path.basename(input_bam_filename)}.{source_token}.SS"

    scriptdir = os.path.abspath(os.path.dirname(__file__))
    cmd = " ".join([os.path.join(scriptdir, "separate_bam_by_strand.py"),
                    "--bam {}".format(input_bam_filename),
                    "--output_prefix {}".format(SS_output_prefix)])


    pipeliner.add_commands([Command(cmd, f"sep_by_strand.{source_token}.ok")])

    pipeliner.run()
    
    ## run normalizations
    SS_bam_files = [SS_output_prefix + x for x in (".+.bam", ".-.bam") ]

    SS_norm_bam_files = list()
    
    for SS_bam_file in SS_bam_files:

        norm_bam_filename = f"{SS_bam_file}.norm_{normalize_max_cov_level}.{run_token}.bam"

        norm_bam_checkpoint = norm_bam_filename + ".ok"

        if not os.path.exists(norm_bam_checkpoint):
            sift_bam(SS_bam_file, norm_bam_filename, normalize_max_cov_level, depth_window, random_seed)
            subprocess.check_call("touch {}".format(norm_bam_checkpoint), shell=True)
        
        SS_norm_bam_files.append(norm_bam_filename)
        

    # merge the norm SS bam filenames into the final output file

    cmd = f"samtools merge -f {output_bam_filename} " + " ".join(SS_norm_bam_files)
    pipeliner.add_commands([Command(cmd, f"SS_merge.{run_token}.ok")])

    cmd = f"samtools index {output_bam_filename}"
    pipeliner.add_commands([Command(cmd, f"index_merged.{run_token}.ok")])

    pipeliner.run()

    logger.info("Done.  See SS-normalized bam: {}".format(output_bam_filename))

    sys.exit(0)


def _read_variate(read_name, random_seed):
    """A uniform [0,1) draw determined by the read name, not by visit order."""
    digest = hashlib.blake2b(
        "{}:{}".format(random_seed, read_name).encode("utf-8"), digest_size=8
    ).digest()
    return int.from_bytes(digest, "big") / float(1 << 64)


def _window_span(lend, rend, depth_window):
    """Windows touched by the half-open interval [lend, rend)."""
    return range(lend // depth_window, (rend - 1) // depth_window + 1)


def sift_bam(
    SS_bam_file, norm_bam_filename, normalize_max_cov_level, depth_window, random_seed
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
    """

    window_bases = dict()  # contig -> aligned bases per window
    junction_support = defaultdict(int)

    with pysam.AlignmentFile(SS_bam_file, "rb") as reader:
        contig_length = dict(zip(reader.references, reader.lengths))

        for read in reader.fetch(until_eof=True):
            if read.is_unmapped:
                continue

            contig = read.reference_name
            bases = window_bases.get(contig)
            if bases is None:
                bases = array("q", bytes(8 * (contig_length[contig] // depth_window + 2)))
                window_bases[contig] = bases

            for block_lend, block_rend in read.get_blocks():
                windows = _window_span(block_lend, block_rend, depth_window)
                for window in windows:
                    covered_lend = max(block_lend, window * depth_window)
                    covered_rend = min(block_rend, (window + 1) * depth_window)
                    bases[window] += covered_rend - covered_lend

            for junction in _read_junctions(read):
                junction_support[(contig, junction)] += 1

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
    read, contig, window_bases, junction_support, normalize_max_cov_level, depth_window
):
    if normalize_max_cov_level <= 0:
        return 1.0

    for junction in _read_junctions(read):
        if junction_support[(contig, junction)] < normalize_max_cov_level:
            # scarce junction: keep every read carrying it so its support is exact
            return 1.0

    depths = []
    for block_lend, block_rend in read.get_blocks():
        for window in _window_span(block_lend, block_rend, depth_window):
            if window < len(window_bases):
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
