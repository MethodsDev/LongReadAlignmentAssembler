#!/usr/bin/env python3
"""
Where reads start and end along a gene's main annotated transcript, per cell cluster: a
check for 5' truncation (reads that end at the polyA site but start short of the TSS).

For each gene (by gene_name in a reference gtf), reads overlapping the gene on its
transcribed strand that pass LRAA's read filters (as alt_termini_read_check.py applies them)
are kept if every intron they carry is an intron of one of the gene's annotated
transcripts (+/- 2 bp). The gene's main transcript is the one whose introns contain the
most full-splice-match reads. Each read's 5' and 3' ends are projected onto the main
transcript (distance in nt from its 5' end, along its exons; ends in an intron or beyond the
transcript are left as NA or clamped to 0 / length), so a read starting at the TSS has
five_prime_tx_pos ~ 0 and one ending at the polyA site has three_prime_tx_pos ~ tx_len.

Outputs:
  <prefix>.reads.tsv.gz    one row per read
  <prefix>.lengths.tsv.gz  aligned read lengths of a sample of all filtered reads, per cluster
                           (--length_sample reads, taken across the genome)
"""

import argparse
import bisect
import collections
import gzip
import logging
import os
import random
import sys

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "diff_iso_usage"))
from alt_termini_read_check import (make_read_filter, parse_cell_clusters, read_geometry,  # noqa: E402
                                    transcript_strand, register_contig)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

INTRON_TOL = 2


def parse_ref_gtf(gtf, genes):
    """gene -> {transcript_id: (contig, strand, sorted exons)}"""
    want = set(genes)
    tx = collections.defaultdict(dict)
    for line in open(gtf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        a = f[8]
        if 'gene_name "' not in a:
            continue
        g = a.split('gene_name "', 1)[1].split('"', 1)[0]
        if g not in want:
            continue
        tid = a.split('transcript_id "', 1)[1].split('"', 1)[0]
        rec = tx[g].setdefault(tid, (f[0], f[6], []))
        rec[2].append((int(f[3]), int(f[4])))
    for g in tx:
        for tid, (c, s, ex) in tx[g].items():
            ex.sort()
    return tx


def introns_of(exons):
    return [(exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1)]


def chain_in(read_introns, tx_introns):
    return all(any(abs(a - x) <= INTRON_TOL and abs(b - y) <= INTRON_TOL for x, y in tx_introns)
               for a, b in read_introns)


def tx_coord(pos, exons, strand):
    """distance (nt) of genomic pos from the transcript's 5' end along its exons; None if
    pos is in an intron; clamped to [0, len] if beyond the transcript"""
    tx_len = sum(e - s + 1 for s, e in exons)
    off = 0
    for s, e in exons:
        if s <= pos <= e:
            d = off + (pos - s)
            return d if strand == "+" else tx_len - 1 - d
        off += e - s + 1
    if pos < exons[0][0]:
        return 0 if strand == "+" else tx_len - 1
    if pos > exons[-1][1]:
        return tx_len - 1 if strand == "+" else 0
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True, help="reference annotation with gene_name attributes")
    ap.add_argument("--genes", required=True, help="comma-separated gene names")
    ap.add_argument("--cell_clusters", required=True, help="cell_barcode <tab> cluster number")
    ap.add_argument("--genome", required=True, help="genome fasta (strand corroboration)")
    ap.add_argument("--HiFi", action="store_true")
    ap.add_argument("--rdna_mask_bed", default=None)
    ap.add_argument("--length_sample", type=int, default=300000,
                    help="filtered reads sampled across the genome for the per-cluster length distribution")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--output_prefix", required=True)
    args = ap.parse_args()

    keep = make_read_filter(HiFi=args.HiFi, rdna_mask_bed=args.rdna_mask_bed)
    cell_to_cluster = parse_cell_clusters(args.cell_clusters)
    genome = pysam.FastaFile(args.genome)
    bam = pysam.AlignmentFile(args.bam)
    genes = args.genes.split(",")
    txs = parse_ref_gtf(args.gtf, genes)

    with gzip.open(f"{args.output_prefix}.reads.tsv.gz", "wt") as out:
        out.write("\t".join(["gene", "main_transcript", "tx_len", "cluster", "read_name", "aligned_len",
                             "n_introns", "five_prime_tx_pos", "three_prime_tx_pos", "first_exon_len"]) + "\n")
        for g in genes:
            if g not in txs:
                logger.warning("%s: not in the gtf", g)
                continue
            recs = txs[g]
            contig, strand = next(iter(recs.values()))[:2]
            lo = min(ex[0][0] for _, _, ex in recs.values())
            hi = max(ex[-1][1] for _, _, ex in recs.values())
            tx_introns = {t: introns_of(ex) for t, (_, _, ex) in recs.items()}
            register_contig(genome, contig)
            reads = []
            for r in bam.fetch(contig, lo - 1, hi):
                if not r.has_tag("CB") or r.get_tag("CB") not in cell_to_cluster:
                    continue
                if not keep(r) or transcript_strand(r, genome) != strand:
                    continue
                lend, rend, introns = read_geometry(r)
                if rend < lo or lend > hi:
                    continue
                fits = [t for t, ti in tx_introns.items() if chain_in(introns, ti)]
                if not fits:
                    continue
                reads.append((r.query_name, cell_to_cluster[r.get_tag("CB")], r.query_alignment_length,
                              lend, rend, introns, fits))
            # main transcript: the one whose full intron chain most reads reproduce exactly
            # (full-splice matches; a subset vote would favour long models that contain every
            # chain); ties and genes without such reads fall back to subset votes, then shorter
            tx_len = {t: sum(e - s + 1 for s, e in ex) for t, (_, _, ex) in recs.items()}
            fsm = collections.Counter(t for *_, introns, fits in reads if introns for t in fits
                                      if len(introns) == len(tx_introns[t]))
            subset = collections.Counter(t for *_, fits in reads for t in fits)
            main_t = max(subset, key=lambda t: (fsm[t], subset[t], -tx_len[t]))
            ex = recs[main_t][2]
            first_exon_len = (ex[0][1] - ex[0][0] + 1) if strand == "+" else (ex[-1][1] - ex[-1][0] + 1)
            n = 0
            for name, cl, alen, lend, rend, introns, fits in reads:
                if main_t not in fits:
                    continue
                five, three = (lend, rend) if strand == "+" else (rend, lend)
                f5, f3 = tx_coord(five, ex, strand), tx_coord(three, ex, strand)
                out.write("\t".join(map(str, [g, main_t, tx_len[main_t], cl, name, alen, len(introns),
                                              "" if f5 is None else f5, "" if f3 is None else f3,
                                              first_exon_len])) + "\n")
                n += 1
            logger.info("%s: main transcript %s (%d nt), %d reads", g, main_t, tx_len[main_t], n)

    # per-cluster aligned read lengths over a genome-wide sample: windows spread over the genome
    random.seed(args.seed)
    contigs = [(c, l) for c, l in zip(bam.references, bam.lengths) if c.startswith("chr") and "_" not in c
               and c not in ("chrM",)]
    total = sum(l for _, l in contigs)
    per_window = 2000
    n_windows = max(1, args.length_sample // per_window * 4)
    cum = []
    acc = 0
    for c, l in contigs:
        acc += l
        cum.append(acc)
    lengths = []
    for _ in range(n_windows):
        if len(lengths) >= args.length_sample:
            break
        x = random.randrange(total)
        i = bisect.bisect_right(cum, x)
        c = contigs[i][0]
        start = x - (cum[i - 1] if i else 0)
        k = 0
        for r in bam.fetch(c, start, start + 200000):
            if r.reference_start < start or not r.has_tag("CB") or r.get_tag("CB") not in cell_to_cluster:
                continue
            if not keep(r):
                continue
            lengths.append((cell_to_cluster[r.get_tag("CB")], r.query_alignment_length))
            k += 1
            if k >= per_window:
                break
    with gzip.open(f"{args.output_prefix}.lengths.tsv.gz", "wt") as out:
        out.write("cluster\taligned_len\n")
        for cl, l in lengths:
            out.write(f"{cl}\t{l}\n")
    logger.info("length sample: %d reads", len(lengths))


if __name__ == "__main__":
    main()
