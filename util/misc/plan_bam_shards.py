#!/usr/bin/env python3

"""Plan how to cut a coordinate-sorted, indexed bam into shards of similar size.

Prints one line per shard, in the order the shards must be appended to reproduce the
bam's own order. A shard is one or more CHUNKS, joined with commas (a comma cannot be
in a contig name), each chunk being tab separated fields:

    contigs <TAB> name name name      whole contigs (small ones grouped, in header order)
    range   <TAB> name <TAB> S <TAB> E   reads of one contig whose 1-based start POS is
                                          in [S, E]
    contigs <TAB> *                   the reads with no coordinate, if any

A chunk is the unit one process classifies; a shard is what one VM runs, its chunks in
parallel on its cores, so it should hold several. Chunks follow each other in the bam's
order and are packed into shards in that order, so appending the chunks' outputs in
order, then the shards', gives the order of a single run.

A contig holding more than --max_reads_per_chunk records is cut into ranges of about
that many. A read belongs to the range holding its start position, so every read is in
exactly one shard (a read may extend past the end of its range, and a range's query
returns reads that began before it; the consumer keeps only POS >= S). Nothing needs to
be unspanned: classifying a read depends only on the read and its contig's annotation.

The cut positions come from the index, not from reading the bam. The bai's linear index
gives, for each 16 kb window, the file offset of the first record overlapping it. That
is too coarse on its own (long spliced reads make many windows share one offset), so
each distinct offset is opened and the position of the record found there is read:
(file offset, position) pairs, exact, a few hundred per contig. The compressed bytes
between them stand in for read counts, and the cuts are placed at equal shares of the
compressed bytes by interpolating between the pairs.
"""

import argparse
import math
import os
import struct
import sys

import pysam

LINEAR_WINDOW = 16384


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
    """(compressed offset, record start position) pairs along one contig, in order."""
    points = []
    seen = set()
    for voffset in linear_index:
        if voffset == 0 or voffset in seen:
            continue
        seen.add(voffset)
        bam.seek(voffset)
        try:
            read = next(bam)
        except StopIteration:
            continue
        if read.reference_name != contig:
            continue
        # file offset of the block, plus how far into it: both move forward together
        # with the position, so the pair is monotone
        points.append((voffset >> 16, read.reference_start + 1))
    points.sort()
    # a trailing sentinel so the last segment has an end
    return points, contig_end_coffset


def cut_positions(points, end_coffset, contig_length, pieces):
    """Positions at which to start pieces 2..n, from equal shares of compressed bytes."""
    if len(points) < 2:
        step = math.ceil(contig_length / pieces)
        return [1 + step * i for i in range(1, pieces)]

    start_coffset = points[0][0]
    total = max(end_coffset - start_coffset, 1)
    # the span after the last sample runs to the contig's end
    xs = [c for c, _ in points] + [end_coffset]
    ps = [p for _, p in points] + [contig_length + 1]

    cuts = []
    segment = 0
    for i in range(1, pieces):
        target = start_coffset + total * i / pieces
        while segment + 1 < len(xs) - 1 and xs[segment + 1] <= target:
            segment += 1
        x0, x1 = xs[segment], xs[segment + 1]
        p0, p1 = ps[segment], ps[segment + 1]
        fraction = (target - x0) / (x1 - x0) if x1 > x0 else 0.0
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
    # (spec, estimated records), in the bam's order
    chunks = []
    group, held = [], 0

    def flush():
        nonlocal group, held
        if group:
            chunks.append(("contigs\t" + " ".join(group), held))
        group, held = [], 0

    for index, contig in enumerate(contigs):
        n = stats[contig]
        if n > args.max_reads_per_chunk and linear is not None:
            flush()
            end = file_end
            for later in contigs[index + 1 :]:
                later_start = first_coffset(later)
                if later_start is not None:
                    end = later_start
                    break
            pieces = math.ceil(n / args.max_reads_per_chunk)
            ref_id = bam.get_tid(contig)
            points, end = sample_points(bam, linear[ref_id], contig, end)
            cuts = cut_positions(points, end, lengths[contig], pieces)
            bounds = [1] + cuts + [lengths[contig]]
            for start, stop in zip(bounds[:-1], bounds[1:]):
                last = stop == lengths[contig]
                chunks.append(
                    (
                        "range\t{}\t{}\t{}".format(
                            contig, start, stop if last else stop - 1
                        ),
                        n // (len(bounds) - 1),
                    )
                )
            continue
        if held > 0 and held + n > args.max_reads_per_chunk:
            flush()
        group.append(contig)
        held += n
    flush()
    if unplaced > 0:
        chunks.append(("contigs\t*", unplaced))

    if not chunks:
        sys.exit("Error, {} has no reads".format(args.bam))

    # Pack consecutive chunks into shards. Never reorder: the outputs are appended in
    # this order.
    shards = []
    current, current_reads = [], 0
    for spec, reads in chunks:
        if current and current_reads + reads > args.max_reads_per_shard:
            shards.append(",".join(current))
            current, current_reads = [], 0
        current.append(spec)
        current_reads += reads
    shards.append(",".join(current))
    print("\n".join(shards))


if __name__ == "__main__":
    main()
