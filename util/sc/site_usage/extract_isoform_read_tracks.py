#!/usr/bin/env python3

"""Alignment blocks of a sample of reads assigned to chosen isoforms, per cell cluster,
for drawing isoform structures with their supporting reads (a read-track figure).

Reads come from an LRAA quant.tracking file (read_name CB^UMI^molecule/N, is_unique,
is_FSM per assignment); by default only reads uniquely assigned to the isoform as a
full splice match. Each read's alignment is fetched from the BAM (primary alignment,
QNAME = the read name's last ^-field) and written as its aligned blocks (CIGAR M/=/X/D
runs, split at N), so introns and the read's own 5' / 3' ends are drawn as aligned.

Up to --max_reads reads are sampled per (isoform, cluster), at random with --seed; with
--proportional, --max_reads reads per cluster are sampled from all chosen isoforms
together, so each isoform's share of a cluster's reads is kept (the switch between
clusters shows in the drawing).

Output tsv, one row per aligned block:
  transcript_id, cluster, read_name, read_start, read_end, strand, block_start, block_end

With --ends_output, also counts, per cluster, the 5' ends (--end 5, TSS) or 3' ends
(--end 3, PolyA) of ALL reads in the region on the isoforms' strand (primary alignments
of cells in the kept clusters, any isoform or none): cluster, pos, reads. That is the
read-end density the sampled reads are drawn from.
"""

import argparse
import collections
import csv
import gzip
import logging
import random
import sys

import pysam

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s", level=logging.INFO)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tracking", required=True, help="LRAA quant.tracking(.gz), may be pre-filtered to the gene")
    parser.add_argument("--transcripts", required=True,
                        help="comma-separated transcript ids as in the tracking file (e.g. t:chr11:+:comp-569:iso-10)")
    parser.add_argument("--bam", required=True, help="indexed BAM the reads were quantified from")
    parser.add_argument("--region", required=True, help="chrom:start-end covering the isoforms")
    parser.add_argument("--cell_clusters", required=True, help="cell_barcode <tab> cluster (header skipped)")
    parser.add_argument("--clusters", default=None, help="comma-separated cluster numbers to keep (default all)")
    parser.add_argument("--max_reads", type=int, default=20, help="reads sampled per isoform and cluster")
    parser.add_argument("--proportional", action="store_true",
                        help="sample --max_reads per cluster across the isoforms, keeping their read shares")
    parser.add_argument("--all_reads", action="store_true",
                        help="any read assigned to the isoform, not only unique full-splice matches")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--ends_output", default=None, help="optional tsv of per-cluster read-end counts")
    parser.add_argument("--end", choices=["5", "3"], default="5", help="read end to count for --ends_output")
    parser.add_argument("--strand", choices=["+", "-"], default=None,
                        help="transcript strand for --ends_output (default: the sampled reads' strand)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    want_tx = set(args.transcripts.split(","))
    keep_clusters = set(args.clusters.split(",")) if args.clusters else None
    cluster_of = {}
    for line in open(args.cell_clusters):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[1].strip().lstrip("-").isdigit():
            cluster_of[f[0]] = f[1].strip()

    groups = collections.defaultdict(list)  # (tx, cluster) -> [molecule name]
    opener = gzip.open if args.tracking.endswith(".gz") else open
    with opener(args.tracking, "rt") as fh:
        rows = csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")
        for r in rows:
            tx = r["transcript_id"].split("@")[-1]
            if tx not in want_tx:
                continue
            if not args.all_reads and not (r["is_unique"] == "1" and r["is_FSM"] == "1"):
                continue
            parts = r["read_name"].split("^")
            cl = cluster_of.get(parts[0])
            if cl is None or (keep_clusters and cl not in keep_clusters):
                continue
            groups[(tx, cl)].append(parts[-1])

    rng = random.Random(args.seed)
    chosen = {}
    if args.proportional:
        by_cluster = collections.defaultdict(list)
        for (tx, cl), names in groups.items():
            by_cluster[cl].extend((n, tx) for n in sorted(set(names)))
        for cl, pool in sorted(by_cluster.items()):
            pool.sort()
            pick = rng.sample(pool, min(args.max_reads, len(pool)))
            for n, tx in pick:
                chosen[n] = (tx, cl)
            tally = collections.Counter(tx for _, tx in pick)
            logger.info("cluster %s: %d reads, %d sampled: %s", cl, len(pool), len(pick), dict(tally))
    else:
        for key, names in sorted(groups.items()):
            names = sorted(set(names))
            pick = rng.sample(names, min(args.max_reads, len(names)))
            for n in pick:
                chosen[n] = key
            logger.info("%s cluster %s: %d reads, %d sampled", key[0], key[1], len(names), len(pick))
    for key, names in sorted(groups.items()):
        logger.info("  %s cluster %s: %d reads", key[0], key[1], len(set(names)))

    chrom, span = args.region.split(":")
    start, end = (int(x.replace(",", "")) for x in span.split("-"))
    found = set()
    with pysam.AlignmentFile(args.bam) as bam, open(args.output, "wt") as ofh:
        w = csv.writer(ofh, delimiter="\t", lineterminator="\n")
        w.writerow(["transcript_id", "cluster", "read_name", "read_start", "read_end", "strand",
                    "block_start", "block_end"])
        for read in bam.fetch(chrom, start - 1, end):
            if read.is_secondary or read.is_supplementary or read.query_name not in chosen:
                continue
            if read.query_name in found:
                continue
            found.add(read.query_name)
            tx, cl = chosen[read.query_name]
            strand = "-" if read.is_reverse else "+"
            if read.has_tag("ts") and read.get_tag("ts") == "-":
                strand = "+" if strand == "-" else "-"
            for b0, b1 in blocks(read):
                w.writerow([tx, cl, read.query_name, read.reference_start + 1, read.reference_end, strand, b0, b1])
    if args.ends_output:
        strands = collections.Counter()
        with open(args.output) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                strands[r["strand"]] += 1
        strand = args.strand or (strands.most_common(1)[0][0] if strands else "+")
        counts = collections.Counter()
        with pysam.AlignmentFile(args.bam) as bam:
            for read in bam.fetch(chrom, start - 1, end):
                if read.is_secondary or read.is_supplementary or not read.has_tag("CB"):
                    continue
                cl = cluster_of.get(read.get_tag("CB"))
                if cl is None or (keep_clusters and cl not in keep_clusters):
                    continue
                s = "-" if read.is_reverse else "+"
                if read.has_tag("ts") and read.get_tag("ts") == "-":
                    s = "+" if s == "-" else "-"
                if s != strand:
                    continue
                five = (read.reference_start + 1) if s == "+" else read.reference_end
                three = read.reference_end if s == "+" else (read.reference_start + 1)
                pos = five if args.end == "5" else three
                if start <= pos <= end:
                    counts[(cl, pos)] += 1
        with open(args.ends_output, "wt") as ofh:
            print("cluster\tpos\treads", file=ofh)
            for (cl, pos), n in sorted(counts.items(), key=lambda x: (int(x[0][0]), x[0][1])):
                print(f"{cl}\t{pos}\t{n}", file=ofh)
        logger.info("read-end counts (%s' ends, %s strand) written to %s", args.end, strand, args.ends_output)

    missing = set(chosen) - found
    if missing:
        logger.warning("%d sampled reads not found in the BAM region (e.g. %s)", len(missing), sorted(missing)[:3])
    logger.info("wrote %d reads", len(found))


def blocks(read):
    """aligned blocks, 1-based inclusive, split only at N (deletions stay inside a block)"""
    out, pos, cur = [], read.reference_start, None
    for op, n in read.cigartuples:
        if op in (0, 2, 7, 8):  # M D = X
            if cur is None:
                cur = pos
            pos += n
        elif op == 3:  # N
            if cur is not None:
                out.append((cur + 1, pos))
                cur = None
            pos += n
    if cur is not None:
        out.append((cur + 1, pos))
    return out


if __name__ == "__main__":
    main()
