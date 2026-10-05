"""site_read_support_to_sparse_matrix.py: per-cell read support for LRAA's sites.

Pins the rules that make the counts LRAA's site support broken down by cell: an end
counts within int(max_dist / 2) = 25 nt of a site and not beyond; a soft clip at the
end disqualifies it unless LRAA strips it (a polyA tail at the 3' end, untemplated G's
at the 5' end); reads LRAA discards (secondary) and reads without a cell barcode are
not counted; minus-strand reads have their TSS at the alignment's right end.
"""
import gzip
import os
import subprocess
import sys

import pysam
import pytest
from scipy.io import mmread

SCRIPT = os.path.join(os.path.dirname(os.path.realpath(__file__)), "site_read_support_to_sparse_matrix.py")

CONTIG_LEN = 10000


def _write_bam(path, reads):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": CONTIG_LEN}]}
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for i, (start, cigar, seq, cb, flag) in enumerate(sorted(reads, key=lambda r: r[0])):
            a = pysam.AlignedSegment(out.header)
            a.query_name = f"r{i}"
            a.reference_id = 0
            a.reference_start = start  # 0-based
            a.cigarstring = cigar
            a.query_sequence = seq
            a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
            a.flag = flag
            a.mapping_quality = 60
            if cb:
                a.set_tag("CB", cb)
            out.write(a)
    pysam.index(path)


def _write_bed(path, rows, polyA=False):
    with open(path, "wt") as ofh:
        ofh.write("#chrom\tstart\tend\tname\tscore\tstrand\tsupport\tn\ttranscript_ids" +
                  ("\tpas\tpas_offset\tinternal_priming" if polyA else "") + "\n")
        for kind, pos, strand in rows:
            f = ["chr1", str(pos - 1), str(pos), f"{kind}:chr1:{pos}:{strand}", "10", strand, "10", "1", "t1"]
            if polyA:
                f += ["AATAAA", "-20", "False"]
            ofh.write("\t".join(f) + "\n")


def _run(tmp_path, reads):
    bam = str(tmp_path / "reads.bam")
    _write_bam(bam, reads)
    _write_bed(tmp_path / "tss.bed", [("TSS", 1001, "+"), ("TSS", 5000, "-")])
    _write_bed(tmp_path / "polya.bed", [("PolyA", 2000, "+")], polyA=True)
    clusters = tmp_path / "clusters.tsv"
    clusters.write_text("cell_barcode\tcluster\nA\t0\nB\t1\nC\t1\n")
    prefix = str(tmp_path / "out")
    subprocess.check_call([sys.executable, SCRIPT, "--bam", bam, "--TSS_bed", str(tmp_path / "tss.bed"),
                           "--PolyA_bed", str(tmp_path / "polya.bed"), "--cell_clusters", str(clusters),
                           "--CPU", "1", "--output_prefix", prefix])
    return prefix


def _matrix(prefix, kind):
    d = f"{prefix}.{kind}-sparseM"
    m = mmread(gzip.open(os.path.join(d, "matrix.mtx.gz"))).tocsr()
    sites = [l.strip() for l in gzip.open(os.path.join(d, "features.tsv.gz"), "rt")]
    cells = [l.strip() for l in gzip.open(os.path.join(d, "barcodes.tsv.gz"), "rt")]
    return {(sites[i], cells[j]): m[i, j] for i, j in zip(*m.nonzero())}


def test_site_support_rules(tmp_path):
    seq = lambda n, base="C": base * n
    reads = [
        # 1-based start 1001 = the TSS; end 2000 = the PolyA site
        (1000, "1000M", seq(1000), "A", 0),
        # TSS end 25 nt off: counted; 26 nt off: not
        (1025, "975M", seq(975), "B", 0),
        (1026, "974M", seq(974), "B", 0),
        # 5' soft clip of non-G bases: TSS end disqualified, PolyA end still counted
        (1000, "4S996M", "AAAA" + seq(996), "C", 0),
        # 5' clip of 3 untemplated G's is stripped: TSS counted
        (1000, "3S1000M", "GGG" + seq(1000), "C", 0),
        # 3' polyA tail clip is stripped: PolyA counted
        (1000, "1000M15S", seq(1000) + "A" * 15, "C", 0),
        # secondary alignment and a read without a cell barcode: never counted
        (1000, "1000M", seq(1000), "A", 256),
        (1000, "1000M", seq(1000), None, 0),
        # minus-strand read: TSS at its right end (5000)
        (4000, "1000M", seq(1000), "A", 16),
    ]
    prefix = _run(tmp_path, reads)

    tss = _matrix(prefix, "TSS")
    assert tss == {("TSS:chr1:1001:+", "A"): 1, ("TSS:chr1:1001:+", "B"): 1,
                   ("TSS:chr1:1001:+", "C"): 2, ("TSS:chr1:5000:-", "A"): 1}

    polya = _matrix(prefix, "PolyA")
    # read 1 (A), the 25/26-nt-offset reads end at 2000 too (B x2), the clipped-5' read and
    # the G-clipped read (C x2), the tail-stripped read (C)
    assert polya == {("PolyA:chr1:2000:+", "A"): 1, ("PolyA:chr1:2000:+", "B"): 2,
                     ("PolyA:chr1:2000:+", "C"): 3}

    cl = [l.rstrip("\n").split("\t") for l in open(prefix + ".TSS.cluster_counts.tsv")]
    assert cl[0] == ["site", "0", "1"]
    assert cl[1] == ["TSS:chr1:1001:+", "1", "3"]

    summary = dict(l.rstrip("\n").split("\t") for l in open(prefix + ".site_read_support.summary.tsv"))
    assert summary["reads_discarded:secondary"] == "1"
    assert summary["reads_without_cell_barcode"] == "1"
    assert summary["TSS_ends_soft_clipped"] == "1"
