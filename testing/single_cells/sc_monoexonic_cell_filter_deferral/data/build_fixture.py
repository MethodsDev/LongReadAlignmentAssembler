#!/usr/bin/env python3
"""Builds the chrTest genome/bam fixture committed alongside this script.

Regenerates chrTest.genome.fa, chrTest.full.bam and chrTest.sg.bam byte-for-byte
(mapping quality, CIGAR, tags and read order are all fixed, so re-running this is a
no-op against the committed files -- kept for documentation/audit, not because the
test needs to regenerate anything at runtime).

Three 300bp monoexonic loci on one synthetic contig, each far enough apart (2000bp)
to fall into its own connected component and never share a read with another:

  locus A (1000-1300, g:chrTest:+:comp-1): the case this fixture exists to pin down.
    chrTest.full.bam carries 10 reads from 8 distinct cells (CELL_A01..CELL_A08).
    chrTest.sg.bam -- standing in for a coverage-normalized splice-graph bam, wired in
    via --bam_for_sg/--no_norm rather than relying on the real normalizer's probabilistic
    window sampling -- carries only 4 of those reads, from 2 distinct cells (CELL_A01,
    CELL_A02): enough reads (4) to clear min_reads_novel_isoform's default of 2, but
    fewer distinct cells (2) than min_monoexonic_supporting_cells's default of 5.
    A filter that judged this locus's cell count from splice-graph-bam evidence would
    wrongly drop it; the full bam's 8 distinct cells clear the bar easily.

  locus B (3000-3300, g:chrTest:+:comp-2): control -- genuinely under-supported.
    Identical 3-read/2-distinct-cell (CELL_B01, CELL_B02) content in BOTH bams, so it
    must stay filtered under any cell-count evidence source, proving a fix does not
    just relax the filter wholesale.

  locus C (5000-5300, g:chrTest:+:comp-3): control -- always kept, byte-identical
    input in both bams (10 reads / 8 distinct cells, CELL_C01..CELL_C08), so its
    component's quant state must be untouched by anything scoped to locus A/B.

Genome sequence is "ACGT" repeated: exactly 25% A content and no run over 1 base, so
no locus's downstream context can ever look internally primed (Util_funcs.looks_
internally_primed fires at >=12 A's in the 20bp window or a 6-mer run -- this pattern
has neither, anywhere on the contig).
"""

import os
import pysam

HERE = os.path.dirname(os.path.abspath(__file__))

CONTIG = "chrTest"
CONTIG_LEN = 8000
READ_LEN = 300

LOCUS_A_START = 1000
LOCUS_B_START = 3000
LOCUS_C_START = 5000


def build_genome():
    pattern = "ACGT"
    seq = (pattern * (CONTIG_LEN // len(pattern) + 1))[:CONTIG_LEN]
    assert len(seq) == CONTIG_LEN
    path = os.path.join(HERE, "chrTest.genome.fa")
    with open(path, "w") as f:
        f.write(f">{CONTIG}\n")
        for i in range(0, len(seq), 70):
            f.write(seq[i : i + 70] + "\n")
    return seq


def make_read(name, start0, cb, umi, seq):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = seq
    a.flag = 0
    a.reference_id = 0
    a.reference_start = start0
    a.mapping_quality = 60
    a.cigar = [(0, READ_LEN)]
    a.next_reference_id = -1
    a.next_reference_start = -1
    a.template_length = 0
    a.query_qualities = pysam.qualitystring_to_array("I" * READ_LEN)
    a.tags = [("NM", 0), ("CB", cb), ("XM", umi)]
    return a


def write_bam(path, reads, header):
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for r in sorted(reads, key=lambda x: x.reference_start):
            out.write(r)
    pysam.index(path)


def main():
    genome_seq = build_genome()
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}],
    }

    def locus_seq(start0):
        return genome_seq[start0 : start0 + READ_LEN]

    full_reads = []
    sg_reads = []

    # locus A: full bam = 10 reads / 8 distinct cells; sg bam = 4 reads / 2 distinct cells
    a_seq = locus_seq(LOCUS_A_START)
    a_cells = [
        "CELL_A01", "CELL_A01", "CELL_A02", "CELL_A02", "CELL_A03",
        "CELL_A04", "CELL_A05", "CELL_A06", "CELL_A07", "CELL_A08",
    ]
    a_reads = [
        make_read(f"readA{i:02d}", LOCUS_A_START, cb, f"UMI_A{i:03d}", a_seq)
        for i, cb in enumerate(a_cells, 1)
    ]
    full_reads.extend(a_reads)
    sg_reads.extend(a_reads[0:4])  # CELL_A01 x2 + CELL_A02 x2 only

    # locus B: identical in both bams; 3 reads / 2 distinct cells (stays filtered)
    b_seq = locus_seq(LOCUS_B_START)
    b_cells = ["CELL_B01", "CELL_B01", "CELL_B02"]
    b_reads = [
        make_read(f"readB{i:02d}", LOCUS_B_START, cb, f"UMI_B{i:03d}", b_seq)
        for i, cb in enumerate(b_cells, 1)
    ]
    full_reads.extend(b_reads)
    sg_reads.extend(b_reads)

    # locus C: identical in both bams; 10 reads / 8 distinct cells (always kept, unaffected)
    c_seq = locus_seq(LOCUS_C_START)
    c_cells = [
        "CELL_C01", "CELL_C01", "CELL_C02", "CELL_C02", "CELL_C03",
        "CELL_C04", "CELL_C05", "CELL_C06", "CELL_C07", "CELL_C08",
    ]
    c_reads = [
        make_read(f"readC{i:02d}", LOCUS_C_START, cb, f"UMI_C{i:03d}", c_seq)
        for i, cb in enumerate(c_cells, 1)
    ]
    full_reads.extend(c_reads)
    sg_reads.extend(c_reads)

    write_bam(os.path.join(HERE, "chrTest.full.bam"), full_reads, header)
    write_bam(os.path.join(HERE, "chrTest.sg.bam"), sg_reads, header)

    print(
        "wrote {} full-bam reads, {} sg-bam reads".format(
            len(full_reads), len(sg_reads)
        )
    )


if __name__ == "__main__":
    main()
