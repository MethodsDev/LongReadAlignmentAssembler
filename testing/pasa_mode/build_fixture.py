#!/usr/bin/env python3
"""
Build testing/pasa_mode/fixture: synthetic staggered reads for --pasa_mode.

One contig per case, random sequence with GT..AG introns. Reads are transcript
fragments written on their transcript's strand, with no polyA tails:

  tiled_plus, tiled_minus  6 exons x 200 bp; 450 bp reads every 150 bp (2-3 exons
                           each, consecutive reads share an intron) -> one full model
  exon_skip                5 exons; staggered reads from exon 3 inclusion and skip
                           isoforms -> both full isoforms
  exon_overlap_only        two read pairs overlapping only inside exon 2, no shared
                           intron -> two models, not chained
  unspliced_bridge         spliced reads at both ends of a 2 kb exon, joined only by
                           an unspliced read -> not chained
  overhang                 read ends running 3-15 bp into introns -> full intron chain,
                           no exon ending inside an intron

Writes fixture/genome.fa(.fai), fixture/reads.bam(.bai) and fixture/expected.tsv
(case, strand, intron chain of each expected model).

usage: build_fixture.py   (needs pysam; deterministic)
"""

import os
import random

import pysam

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixture")
FLANK = 1000
rng = random.Random(7)


def rand_seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def locus(exon_lens, intron_len=500):
    """+ strand sequence and 1-based exon coords."""
    seq, exons = rand_seq(FLANK), []
    for i, el in enumerate(exon_lens):
        exons.append((len(seq) + 1, len(seq) + el))
        seq += rand_seq(el)
        if i < len(exon_lens) - 1:
            seq += "GT" + rand_seq(intron_len - 4) + "AG"
    return seq + rand_seq(FLANK), exons


def window(exons, a, b):
    """Genome segments of transcript bases [a, b) for + strand exons."""
    segs, pos = [], 0
    for lend, rend in exons:
        x, y = max(a, pos), min(b, pos + rend - lend + 1)
        if x < y:
            segs.append((lend + x - pos, lend + y - pos - 1))
        pos += rend - lend + 1
    return segs


def tiles(exons, read_len, step, copies=2, jitter=10):
    tlen = sum(r - l + 1 for l, r in exons)
    starts = sorted(set(list(range(0, tlen - read_len + 1, step)) + [tlen - read_len]))
    reads = []
    for s in starts:
        for _ in range(copies):
            a = max(0, s + rng.randint(-jitter, jitter))
            b = min(tlen, s + read_len + rng.randint(-jitter, jitter))
            reads.append(window(exons, a, b))
    return reads


def introns(exons):
    return [(exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1)]


def build():
    cases = []  # (name, strand, seq, reads, expected intron chains)

    for name, strand in (("tiled_plus", "+"), ("tiled_minus", "-")):
        seq, ex = locus([200] * 6)
        reads = tiles(ex, 450, 150)
        expected = [introns(ex)]
        if strand == "-":
            L = len(seq)
            flip = lambda segs: sorted((L - r + 1, L - l + 1) for l, r in segs)
            seq, reads = revcomp(seq), [flip(r) for r in reads]
            expected = [flip(introns(ex))]
        cases.append((name, strand, seq, reads, expected))

    seq, ex = locus([200] * 5)
    skip = [e for i, e in enumerate(ex) if i != 2]
    cases.append(("exon_skip", "+", seq, tiles(ex, 450, 150) + tiles(skip, 450, 150),
                  [introns(ex), introns(skip)]))

    seq, ex = locus([400, 400, 400])
    reads = [window(ex, 300, 700)] * 2 + [window(ex, 500, 900)] * 2
    cases.append(("exon_overlap_only", "+", seq, reads, [introns(ex[:2]), introns(ex[1:])]))

    seq, ex = locus([200, 2000, 200])
    reads = [window(ex, 0, 500)] * 2 + [window(ex, 1900, 2400)] * 2 + [window(ex, 300, 2100)] * 2
    cases.append(("unspliced_bridge", "+", seq, reads, [introns(ex[:2]), introns(ex[1:])]))

    seq, ex = locus([200] * 5)
    ends = {r for _, r in ex[:-1]}
    reads = []
    for segs in tiles(ex, 450, 150, jitter=0):
        if segs[-1][1] in ends:
            segs = segs[:-1] + [(segs[-1][0], segs[-1][1] + rng.randint(3, 15))]
        reads.append(segs)
    cases.append(("overhang", "+", seq, reads, [introns(ex)]))
    return cases


def main():
    os.makedirs(OUT, exist_ok=True)
    cases = build()
    fa = os.path.join(OUT, "genome.fa")
    with open(fa, "w") as fh:
        for name, _, seq, _, _ in cases:
            fh.write(f">{name}\n" + "".join(seq[i:i + 60] + "\n" for i in range(0, len(seq), 60)))
    pysam.faidx(fa)

    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": c[0], "LN": len(c[2])} for c in cases]}
    bam_path = os.path.join(OUT, "reads.bam")
    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for tid, (name, strand, seq, reads, _) in enumerate(cases):
            for i, segs in sorted(enumerate(reads), key=lambda x: x[1][0][0]):
                a = pysam.AlignedSegment()
                a.query_name = f"{name}.r{i}"
                a.reference_id = tid
                a.reference_start = segs[0][0] - 1
                cig = []
                for k, (l, r) in enumerate(segs):
                    if k:
                        cig.append((3, l - segs[k - 1][1] - 1))
                    cig.append((0, r - l + 1))
                a.cigartuples = cig
                a.query_sequence = "".join(seq[l - 1:r] for l, r in segs)
                a.flag = 16 if strand == "-" else 0
                a.mapping_quality = 60
                a.set_tag("NM", 0)
                bam.write(a)
    pysam.index(bam_path)

    with open(os.path.join(OUT, "expected.tsv"), "w") as fh:
        fh.write("#case\tstrand\tintron_chain\n")
        for name, strand, _, _, expected in cases:
            for chain in expected:
                fh.write(f"{name}\t{strand}\t{','.join(f'{l}-{r}' for l, r in chain)}\n")


if __name__ == "__main__":
    main()
