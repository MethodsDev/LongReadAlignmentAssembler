#!/usr/bin/env python3

"""Plan how to cut a coordinate-sorted, indexed bam into shards of similar size.

Prints one line per shard, in the order the shards must be appended to reproduce the
bam's own order. A shard is one or more CHUNKS, joined with commas (a comma cannot be
in a contig name), each chunk being tab separated fields:

    contigs <TAB> name name name <TAB> N      whole contigs (small ones grouped, in header
                                              order)
    range   <TAB> name <TAB> S <TAB> E <TAB> N   reads of one contig whose 1-based start
                                              POS is in [S, E]
    contigs <TAB> * <TAB> N                  the reads with no coordinate, if any

N is the estimated number of records, last on every line, so a consumer can start the
largest chunks first.

A chunk is the unit one process classifies; a shard is what one VM runs, its chunks in
parallel on its cores, so it should hold several. Chunks follow each other in the bam's
order and are packed into shards in that order, so appending the chunks' outputs in
order, then the shards', gives the order of a single run.

A contig holding more than --max_reads_per_chunk records is cut into ranges of about
that many. A read belongs to the range holding its start position, so every read is in
exactly one shard (a read may extend past the end of its range, and a range's query
returns reads that began before it; the consumer keeps only POS >= S). Nothing needs to
be unspanned: classifying a read depends only on the read and its contig's annotation.

The cut positions come from the index, not from a pass over the reads. The bai's linear
index gives, for each 16 kb window, the file offset of the first record overlapping it.
That is too coarse on its own (long spliced reads make many windows share one offset), so
each distinct offset is opened and the record found there is read: (file offset,
position) pairs, exact, a few hundred per contig. Compressed bytes alone are a poor
stand-in for read counts, though: reads of a very highly expressed locus compress far
better than the rest, so an equal share of bytes can hold several times the reads (the
first Terra run: one chunk of 2.25 million reads against 575 thousand planned). So the
few blocks after each sample are decoded as well and the records counted, which gives the
records per compressed byte right there; the reads between two samples are estimated from
that density, and cuts are placed at equal shares of the estimated reads.
"""

import argparse
import math
import os
import struct
import sys

import pysam
from concurrent.futures import ProcessPoolExecutor

LINEAR_WINDOW = 16384
# uncompressed bytes per full BGZF block (htslib's BGZF_BLOCK_SIZE)
BGZF_BLOCK_BYTES = 65280
# blocks decoded after each sample to measure the local density of records
DENSITY_BLOCKS = 4


def read_linear_indexes(bai_path):
    """Per reference, the list of linear-index virtual offsets (BAI, not CSI)."""
    with open(bai_path, "rb") as fh:
        if fh.read(4) != b"BAI\1":
            return None
        (n_ref,) = struct.unpack("<i", fh.read(4))
        indexes = []
        for _ in range(n_ref):
            (n_bin,) = struct.unpack("<i", fh.read(4))
            for _ in range(n_bin):
                _bin, n_chunk = struct.unpack("<Ii", fh.read(8))
                fh.seek(16 * n_chunk, os.SEEK_CUR)
            (n_intv,) = struct.unpack("<i", fh.read(4))
            indexes.append(
                list(struct.unpack("<%dQ" % n_intv, fh.read(8 * n_intv)))
                if n_intv
                else []
            )
    return indexes


def sample_points(bam, linear_index, contig, contig_end_coffset):
    """(compressed offset, record start position, records per compressed byte) along one
    contig, in order."""
    points = []
    seen = set()
    for voffset in linear_index:
        if voffset == 0 or voffset in seen:
            continue
        seen.add(voffset)
        bam.seek(voffset)
        try:
            first = next(bam)
        except StopIteration:
            continue
        if first.reference_name != contig:
            continue

        # Count the records over the next few blocks. The first block is entered part
        # way (the offset points at a record), so the compressed bytes it contributes
        # are scaled by the fraction of it that is left.
        first_coffset = voffset >> 16
        entered = min((voffset & 0xFFFF) / BGZF_BLOCK_BYTES, 1.0)
        records = 1
        next_block = None  # offset of the block after the first
        distinct_blocks = 1
        last_block = first_coffset
        stop_coffset = None
        while True:
            coffset = bam.tell() >> 16
            if coffset != last_block:
                if next_block is None:
                    next_block = coffset
                last_block = coffset
                distinct_blocks += 1
                if distinct_blocks > DENSITY_BLOCKS:
                    stop_coffset = coffset
                    break
            try:
                read = next(bam)
            except StopIteration:
                break
            if read.reference_name != contig:
                break
            records += 1
        density = 0.0
        if next_block is not None:
            end_coffset = stop_coffset if stop_coffset is not None else last_block
            covered = (end_coffset - first_coffset) - entered * (next_block - first_coffset)
            if covered > 0:
                density = records / covered
        # file offset of the block, plus how far into it: both move forward together
        # with the position, so the pair is monotone
        points.append((first_coffset, first.reference_start + 1, density))
    points.sort()
    return points, contig_end_coffset


def sample_points_in_worker(bam_path, linear_index, contig, contig_end_coffset):
    # one handle per call: a handle cannot be shared between processes
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        return sample_points(bam, linear_index, contig, contig_end_coffset)


def cut_positions(points, end_coffset, contig_length, pieces):
    """Positions at which to start pieces 2..n, from equal shares of the estimated reads."""
    if len(points) < 2:
        step = math.ceil(contig_length / pieces)
        return [1 + step * i for i in range(1, pieces)]

    # the span after the last sample runs to the contig's end, at the last density
    xs = [c for c, _, _ in points] + [end_coffset]
    ps = [p for _, p, _ in points] + [contig_length + 1]
    ds = [d for _, _, d in points]
    ds.append(ds[-1])
    if not any(ds):
        # no density measured anywhere (every sample inside one block): bytes instead
        ds = [1.0] * len(xs)

    # estimated records from the first sample up to each sample: the compressed bytes
    # between two samples times the mean density at their ends
    cumulative = [0.0]
    for i in range(len(xs) - 1):
        cumulative.append(
            cumulative[-1] + max(xs[i + 1] - xs[i], 0) * (ds[i] + ds[i + 1]) / 2
        )
    total = cumulative[-1]
    if total <= 0:
        step = math.ceil(contig_length / pieces)
        return [1 + step * i for i in range(1, pieces)]

    cuts = []
    segment = 0
    for i in range(1, pieces):
        target = total * i / pieces
        while segment + 1 < len(xs) - 1 and cumulative[segment + 1] <= target:
            segment += 1
        r0, r1 = cumulative[segment], cumulative[segment + 1]
        p0, p1 = ps[segment], ps[segment + 1]
        fraction = (target - r0) / (r1 - r0) if r1 > r0 else 0.0
        cuts.append(int(p0 + fraction * (p1 - p0)))
    # strictly increasing and inside the contig
    cleaned = []
    for cut in cuts:
        if cut > (cleaned[-1] if cleaned else 1) and cut <= contig_length:
            cleaned.append(cut)
    return cleaned


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--bam", required=True)
    parser.add_argument(
        "--bai", default=None, help="default: the index pysam finds beside the bam"
    )
    parser.add_argument(
        "--max_reads_per_chunk",
        type=int,
        required=True,
        help="records one process classifies: contigs are grouped up to this, a larger "
        "contig is cut into ranges of about this",
    )
    parser.add_argument(
        "--max_reads_per_shard",
        type=int,
        required=True,
        help="records one VM classifies: consecutive chunks are packed up to this. "
        "Should be several times --max_reads_per_chunk, so a shard's cores stay busy",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="processes sampling the index, one contig each (default: the cores this "
        "process may use). The sampling is the slow part of planning, and the contigs "
        "are independent",
    )
    parser.add_argument(
        "--no_split",
        action="store_true",
        help="never cut a contig, only group them (also what happens without a .bai)",
    )
    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    stats = {s.contig: s.mapped + s.unmapped for s in bam.get_index_statistics()}
    unplaced = bam.nocoordinate
    contigs = [c for c in bam.references if stats.get(c, 0) > 0]
    lengths = dict(zip(bam.references, bam.lengths))

    linear = None
    if not args.no_split:
        linear = read_linear_indexes(args.bai or args.bam + ".bai")
        if linear is None:
            sys.stderr.write("not a .bai index; contigs are grouped, not cut\n")

    # where each contig's records end in the file: where the next contig with records
    # begins, the last one at the end of the file (less its 28-byte EOF block)
    def first_coffset(contig):
        ref_id = bam.get_tid(contig)
        offsets = [o for o in linear[ref_id] if o]
        return (min(offsets) >> 16) if offsets else None

    file_end = os.path.getsize(args.bam) - 28

    # The contigs to cut, each with where its records end in the file, sampled in
    # parallel before the planning below needs them.
    sampled = {}
    if linear is not None:
        to_cut = []
        for index, contig in enumerate(contigs):
            if stats[contig] > args.max_reads_per_chunk:
                end = file_end
                for later in contigs[index + 1 :]:
                    later_start = first_coffset(later)
                    if later_start is not None:
                        end = later_start
                        break
                to_cut.append((contig, end))
        threads = args.threads or len(os.sched_getaffinity(0))
        workers = max(1, min(threads, len(to_cut)))
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                contig: pool.submit(
                    sample_points_in_worker,
                    args.bam,
                    linear[bam.get_tid(contig)],
                    contig,
                    end,
                )
                for contig, end in to_cut
            }
            sampled = {contig: future.result() for contig, future in futures.items()}

    # (spec, estimated records), in the bam's order
    chunks = []
    group, held = [], 0

    def flush():
        nonlocal group, held
        if group:
            chunks.append(("contigs\t{}\t{}".format(" ".join(group), held), held))
        group, held = [], 0

    for index, contig in enumerate(contigs):
        n = stats[contig]
        if n > args.max_reads_per_chunk and linear is not None:
            flush()
            pieces = math.ceil(n / args.max_reads_per_chunk)
            points, end = sampled[contig]
            cuts = cut_positions(points, end, lengths[contig], pieces)
            bounds = [1] + cuts + [lengths[contig]]
            for start, stop in zip(bounds[:-1], bounds[1:]):
                last = stop == lengths[contig]
                estimate = n // (len(bounds) - 1)
                chunks.append(
                    (
                        "range\t{}\t{}\t{}\t{}".format(
                            contig, start, stop if last else stop - 1, estimate
                        ),
                        estimate,
                    )
                )
            continue
        if held > 0 and held + n > args.max_reads_per_chunk:
            flush()
        group.append(contig)
        held += n
    flush()
    if unplaced > 0:
        chunks.append(("contigs\t*\t{}".format(unplaced), unplaced))

    if not chunks:
        sys.exit("Error, {} has no reads".format(args.bam))

    # Pack consecutive chunks into shards, never reordering: the outputs are appended
    # in this order. As few shards as the cap allows, then evenly: a chunk goes to the
    # shard its midpoint falls in when the reads are laid end to end, so no shard is
    # left holding one small chunk while another holds a full load.
    total = sum(reads for _, reads in chunks)
    num_shards = max(1, math.ceil(total / args.max_reads_per_shard))
    per_shard = total / num_shards
    packed = [[] for _ in range(num_shards)]
    before = 0
    for spec, reads in chunks:
        packed[min(int((before + reads / 2) / per_shard), num_shards - 1)].append(spec)
        before += reads
    shards = [",".join(group) for group in packed if group]
    print("\n".join(shards))


if __name__ == "__main__":
    main()
