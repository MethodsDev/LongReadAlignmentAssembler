#!/usr/bin/env python3

"""Collapse a BAM's accumulated @PG chain, preserving @HD/@SQ/@RG/@CO.

WHY THIS EXISTS

samtools appends one @PG record PER EXISTING CHAIN TIP on every write, and
`samtools merge` first concatenates all of its inputs' chains. On a header
carrying many parallel chains that is multiplicative, not additive.

Measured on XP132160.ucsc.bam, which arrives from upstream alignment with
34,976 @PG records across 5,824 parallel chains -- one per minimap2 alignment
shard, never collapsed. Through the cluster-guided path that became 2,727,296
records = 1,188,154,439 bytes of uncompressed SAM header text, in a merged BAM
of 3.68 GiB on disk. Consequences, all measured:

  * Resolving a region NAME to a tid forces htslib to parse the whole header,
    so each per-chromosome region query cost ~5 min of pure parsing regardless
    of how many alignments it returned (chr1 at 1.67M and chrY at 27K
    alignments measured within 2.2% of each other). ~27.5 h of one whole-genome
    run was header parsing.
  * `samtools view -b` copies the source header into every output, so each
    per-chromosome BAM inherited the whole chain: a chrY.bam holding 26,725
    alignments measured 87.9 MB on disk, and ~27 GB of 44 GB of per-chromosome
    outputs was duplicated header rather than alignment data.

Passing `--no-PG` at the writing sites caps further growth but cannot remove
inherited history: measured, it would leave 1,818,752 records / ~792 MB, still
20% of the file. Removing the accumulated chain requires rewriting the header,
which is what this does.

SAFETY

Dropping @PG definitions is only sound if no alignment record references one
via a PG:Z: tag. Verified on this corpus by full sequential scans -- 46,123,848
records of the input BAM and 17,030,856 of the merged BAM, 63,154,704 in total,
zero PG:Z: tags -- and LRAA itself never reads or writes PG tags. This script
re-checks rather than trusting that: it refuses to collapse if it finds a
PG:Z: tag, unless --force is given. A header with zero @PG records is valid
SAM, so no synthetic replacement record is required; --keep-provenance adds
one if you would rather retain a marker.

MEMORY

The header is read straight from the BGZF stream, never through `samtools
view -H` (which appends one record per chain tip, so it cannot report its own
input: a 34,976-record header reads as 40,800, and it has been observed to
hang at this size) and never through pysam's to_dict(). Peak RSS measured on
an 81.3 MB / 200,000-record header:

    streaming header count (this script)   16.3 MB   0.20x header
    samtools view -c -d PG (tag scan)      90.0 MB   1.11x
    samtools reheader -P -c (the collapse) 90.1 MB   1.11x
    pysam AlignmentFile open              104.2 MB   1.28x
    pysam header.to_dict()                469.0 MB   5.77x

KNOWN LIMIT, not a safety claim: ~1.11x the header text is irreducible here,
because htslib materializes the header text to rewrite it and `samtools
reheader` is what performs the collapse. At the worst case measured on this
corpus -- 1,188,154,439 bytes of header -- that implies roughly 1.3 GB of RAM
for the reheader step [INFERENCE: linear extrapolation from the table above,
not measured at that size]. Nothing this script does adds to that floor; the
point of the streaming reader is that the check for whether a collapse is even
needed does not multiply it.

ORDERING

Rewriting the header shifts every BGZF virtual offset, so any pre-existing
.bai is invalid afterwards -- a region query against the result fails with
"Invalid BGZF header at offset N". This script therefore reindexes by default;
pass --no-index only if the caller indexes afterwards itself.
"""

import argparse
import gzip
import logging
import os
import struct
import subprocess
import sys
import tempfile

logging.basicConfig(
    format="%(asctime)s : %(levelname)s : %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

PRESERVED = (b"@HD", b"@SQ", b"@RG", b"@CO")
CHUNK = 1 << 20


def _open_header(bam_filename):
    """Open a BAM and position at its header text, returning (stream, n_bytes).

    The header is read from the BGZF stream directly rather than through
    samtools or pysam, for two independent reasons. samtools appends one @PG
    record per existing chain tip to the header it emits, so it cannot report
    its own input faithfully (measured in the pinned image: 1 for an empty
    chain, 85 for an 84-record one). And both `samtools view -H` and pysam's
    to_dict() materialize the whole header at once, which on the headers this
    script exists to remove is a multi-GB allocation -- measured on an 85 MB
    header, peak RSS was 191 MB decoding and splitting the text and 481 MB via
    to_dict(), against 17 MB for the streaming path below.
    """
    fh = gzip.open(bam_filename, "rb")
    magic = fh.read(4)
    if magic != b"BAM\x01":
        fh.close()
        raise ValueError(f"{bam_filename} is not a BAM file (magic {magic!r})")
    return fh, struct.unpack("<i", fh.read(4))[0]


def _header_lines(fh, remaining):
    """Yield header lines as bytes, holding at most one chunk plus one line."""
    pending = b""
    while remaining > 0:
        chunk = fh.read(min(CHUNK, remaining))
        if not chunk:
            break
        remaining -= len(chunk)
        pending += chunk
        lines = pending.split(b"\n")
        pending = lines.pop()
        for line in lines:
            yield line
    # Header text may be NUL-padded and need not end in a newline.
    pending = pending.rstrip(b"\x00")
    if pending:
        yield pending


def count_pg(bam_filename):
    """@PG record count, exact and in bounded memory."""
    fh, remaining = _open_header(bam_filename)
    try:
        total = 0
        prev = b"\n"  # a newline notionally precedes the first line
        while remaining > 0:
            chunk = fh.read(min(CHUNK, remaining))
            if not chunk:
                break
            remaining -= len(chunk)
            window = prev + chunk
            total += window.count(b"\n@PG")
            prev = window[-3:]  # shortest overlap that can hide a b"\n@PG"
        return total
    finally:
        fh.close()


def count_pg_tagged_records(bam_filename, threads=1):
    """Number of alignments carrying a PG:Z: tag. Sequential, no region query.

    Delegated to samtools rather than pysam on measured grounds. On an 81.3 MB
    header, peak RSS was 90.0 MB for this call against 104.2 MB merely to open
    the file through pysam (and 469.0 MB via pysam's to_dict). More to the
    point, 90.0 MB is what the `samtools reheader` collapse immediately after
    it already costs, so this check adds no peak the operation was not paying
    anyway -- whereas a pure-Python record scan, though bounded, would spend
    minutes of interpreter time per 10M records to avoid a cost we still pay.
    """
    done = subprocess.run(
        ["samtools", "view", "-@", str(threads), "-c", "-d", "PG", bam_filename],
        capture_output=True,
        text=True,
        check=True,
    )
    return int(done.stdout.strip())


def write_collapsed_header(bam_filename, out_path, keep_provenance):
    """Write a header with every @PG line dropped, streaming, byte-preserving.

    Filtering the raw text rather than round-tripping pysam's to_dict()/
    from_dict() keeps @HD/@SQ/@RG/@CO exactly as they arrived -- field order and
    any tags pysam does not model included -- and never holds more than one
    chunk. Returns (kept, dropped).
    """
    fh, remaining = _open_header(bam_filename)
    kept = dropped = 0
    try:
        with open(out_path, "wb") as out:
            for line in _header_lines(fh, remaining):
                if line.startswith(b"@PG"):
                    dropped += 1
                    continue
                if line.startswith(b"@") or line:
                    out.write(line + b"\n")
                    kept += 1
            if keep_provenance:
                record = "\t".join(
                    [
                        "@PG",
                        "ID:lraa_collapse_pg",
                        "PN:collapse_bam_pg_header.py",
                        "CL:" + " ".join(sys.argv).replace("\t", " "),
                        "DS:prior @PG chain collapsed by LRAA "
                        "util/misc/collapse_bam_pg_header.py",
                    ]
                )
                out.write(record.encode("utf-8") + b"\n")
                kept += 1
    finally:
        fh.close()
    return kept, dropped


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input_bam", required=True, help="BAM to collapse")
    parser.add_argument(
        "--output_bam",
        help="output path; omit to rewrite --input_bam in place via a temp file",
    )
    parser.add_argument(
        "--keep-provenance",
        action="store_true",
        help="retain a single @PG record naming this collapse instead of none",
    )
    parser.add_argument(
        "--no-index",
        action="store_true",
        help="skip reindexing; only safe if the caller indexes the output itself",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="collapse even if a PG:Z: tag is found (leaves dangling references)",
    )
    parser.add_argument("--threads", type=int, default=1, help="samtools threads")
    args = parser.parse_args()

    in_bam = args.input_bam
    if not os.path.exists(in_bam):
        logger.error("no such file: %s", in_bam)
        sys.exit(1)

    before = count_pg(in_bam)
    logger.info("@PG records before: %d", before)
    if before == 0 and not args.keep_provenance:
        logger.info("nothing to collapse; leaving %s untouched", in_bam)
        sys.exit(0)

    n_tagged = count_pg_tagged_records(in_bam, threads=args.threads)
    if n_tagged:
        msg = (
            "%d alignment records carry a PG:Z: tag; collapsing @PG definitions "
            "would leave those references dangling" % n_tagged
        )
        if not args.force:
            logger.error("%s -- refusing (use --force to override)", msg)
            sys.exit(2)
        logger.warning("%s -- proceeding under --force", msg)
    else:
        logger.info("no PG:Z: tags on any record; collapse is reference-safe")

    in_place = args.output_bam is None
    out_bam = args.output_bam or in_bam + ".pgcollapsed.tmp.bam"

    hdr_fd, hdr_path = tempfile.mkstemp(suffix=".sam", dir=os.path.dirname(out_bam) or ".")
    os.close(hdr_fd)
    try:
        kept, dropped = write_collapsed_header(in_bam, hdr_path, args.keep_provenance)
        logger.info("header rewritten: %d lines kept, %d @PG dropped", kept, dropped)
        # -P so reheader does not append an @PG record of its own.
        with open(out_bam, "wb") as out_fh:
            subprocess.run(
                ["samtools", "reheader", "-P", hdr_path, in_bam],
                stdout=out_fh,
                check=True,
            )
    finally:
        os.path.exists(hdr_path) and os.remove(hdr_path)

    if in_place:
        os.replace(out_bam, in_bam)
        out_bam = in_bam
        # A .bai from before the rewrite is invalid: the header length changed, so
        # every BGZF virtual offset moved. Remove it so a stale index cannot be
        # picked up if --no-index was passed.
        for stale in (in_bam + ".bai", os.path.splitext(in_bam)[0] + ".bai"):
            os.path.exists(stale) and os.remove(stale)

    after = count_pg(out_bam)
    logger.info(
        "@PG records after: %d (removed %d) -> %s", after, before - after, out_bam
    )

    if not args.no_index:
        subprocess.run(
            ["samtools", "index", "-@", str(args.threads), out_bam], check=True
        )
        logger.info("reindexed %s", out_bam)
    else:
        logger.warning(
            "--no-index given; %s has no valid index until the caller builds one",
            out_bam,
        )

    logger.info("Done.")


if __name__ == "__main__":
    main()
