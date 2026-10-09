"""pytest: alt_termini_read_check.py counts the reads LRAA would use, on LRAA's strand.

One + strand gene: two isoforms sharing the intron 501-1000 (GT..AG), differing at the TSS
(101 vs 301). Cluster 1 cells start mostly at 101, cluster 2 cells mostly at 301. Extra
reads exercise the read filters and the ts-flip corroboration:
  - a supplementary alignment (never counted),
  - a 93%-identity read (counted only without --HiFi),
  - a ts:A:- read aligned reverse: the flip to + is corroborated by the GT..AG junction,
    so it counts (it is a + strand read sequenced antisense),
  - a ts:A:- read aligned forward: the flip to - is NOT corroborated (the junction is not
    canonical on -), so it stays + and counts; trusting ts alone dropped it.
"""

import csv
import os
import random
import subprocess
import sys

import pysam

HERE = os.path.dirname(os.path.realpath(__file__))
SCRIPT = os.path.join(HERE, "alt_termini_read_check.py")

CONTIG_LEN = 1500
INTRON = (501, 1000)
EXON2_END = 1300


def _contig():
    rng = random.Random(1)
    seq = [rng.choice("CG") for _ in range(CONTIG_LEN)]  # no A-runs: no polyA-terminal discards
    seq[INTRON[0] - 1:INTRON[0] + 1] = list("GT")
    seq[INTRON[1] - 2:INTRON[1]] = list("AG")
    return "".join(seq)


def _read(header, contig, name, start, cb, flag=0, nm=0, ts=None):
    """spliced read from `start` (1-based) through the intron to the end of exon 2"""
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.flag = flag
    a.reference_id = 0
    a.reference_start = start - 1
    a.mapping_quality = 60
    e1 = INTRON[0] - start
    e2 = EXON2_END - INTRON[1]
    a.cigartuples = [(0, e1), (3, INTRON[1] - INTRON[0] + 1), (0, e2)]
    a.query_sequence = contig[start - 1:INTRON[0] - 1] + contig[INTRON[1]:EXON2_END]
    a.query_qualities = pysam.qualitystring_to_array("I" * (e1 + e2))
    tags = [("CB", cb), ("NM", nm)]
    if ts:
        tags.append(("ts", ts, "A"))
    a.set_tags(tags)
    return a


def _inputs(tmp_path):
    contig = _contig()
    fa = tmp_path / "g.fa"
    fa.write_text(">chrT\n" + contig + "\n")
    pysam.faidx(str(fa))

    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                              "SQ": [{"SN": "chrT", "LN": CONTIG_LEN}]})
    reads = []
    n = 0

    def add(start, cb, **kw):
        nonlocal n
        n += 1
        reads.append(_read(header, contig, f"r{n}", start, cb, **kw))

    for i in range(8):
        add(101, "cA")
    for i in range(2):
        add(301, "cA")
    for i in range(2):
        add(101, "cB")
    for i in range(8):
        add(301, "cB")
    add(301, "cA", flag=2048)              # supplementary
    add(301, "cA", nm=50)                  # ~93% identity
    add(101, "cB", flag=16, ts="-")        # antisense-sequenced + read, flip corroborated
    add(101, "cA", ts="-")                 # ts says -, junction not canonical on -: stays +
    reads.sort(key=lambda r: r.reference_start)
    unsorted = tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for r in reads:
            out.write(r)
    pysam.index(str(unsorted))

    gtf = tmp_path / "t.gtf"
    with open(gtf, "w") as fh:
        for tid, tss in (("G^dom", 101), ("G^alt", 301)):
            for s, e in ((tss, INTRON[0] - 1), (INTRON[1] + 1, EXON2_END)):
                fh.write(f'chrT\tLRAA\texon\t{s}\t{e}\t.\t+\t.\tgene_id "G"; transcript_id "{tid}";\n')
    clusters = tmp_path / "clusters.tsv"
    clusters.write_text("cell_barcode\tcluster\ncA\t1\ncB\t2\n")
    cand = tmp_path / "cand.tsv"
    cand.write_text("gene_symbol\talt_terminus\tdominant_transcript_ids\talternate_transcript_ids\tcluster_A\tcluster_B\n"
                    "G\tTSS\tG^dom\tG^alt\tCluster_1\tCluster_2\n")
    return fa, unsorted, gtf, clusters, cand


def _run(tmp_path, *extra):
    fa, bam, gtf, clusters, cand = _inputs(tmp_path)
    out = tmp_path / ("out" + "_".join(extra) + ".tsv")
    subprocess.run([sys.executable, SCRIPT, "--candidates", str(cand), "--gtf", str(gtf), "--bam", str(bam),
                    "--cell_clusters", str(clusters), "--genome_fa", str(fa), "--output", str(out), *extra],
                   check=True, capture_output=True)
    return next(csv.DictReader(open(out), delimiter="\t"))


def test_hifi_filters_and_strand(tmp_path):
    r = _run(tmp_path, "--HiFi")
    # cluster 1: 8 + the uncorroborated-ts read at 101; 2 at 301 (supplementary and 93% read dropped)
    assert (int(r["reads_dom_A"]), int(r["reads_alt_A"])) == (9, 2)
    # cluster 2: 2 + the antisense-sequenced read at 101; 8 at 301
    assert (int(r["reads_dom_B"]), int(r["reads_alt_B"])) == (3, 8)


def test_default_identity_floor_keeps_93pct_read(tmp_path):
    r = _run(tmp_path)
    assert (int(r["reads_dom_A"]), int(r["reads_alt_A"])) == (9, 3)
    assert (int(r["reads_dom_B"]), int(r["reads_alt_B"])) == (3, 8)
