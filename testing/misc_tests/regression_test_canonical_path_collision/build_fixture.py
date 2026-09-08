#!/usr/bin/env python3
"""Shrink the chr21 canonical-path-collision repro to a locus-sized fixture.

Window is coordinate-SHIFTED to offset 0 so the fixture carries a 40 kb contig
instead of a 6.7 Mb one. The collision is a property of one read-sharing
component (comp-1119, chr21:6,286,342-6,313,221), so translating every
coordinate by a constant preserves it; only the literal path string in the
error message changes, and the test asserts on the exception type, not on it.
"""
import sys, pysam

SRC = sys.argv[1]        # repro dir
DST = sys.argv[2]        # fixture dir
START, END = 6_280_000, 6_320_000        # 1-based inclusive
OFF = START - 1                           # subtract from every 1-based coord
NEWLEN = END - START + 1
CTG = "chr21"

# ---- fasta ----
fa = pysam.FastaFile(f"{SRC}/chunk.fa")
seq = fa.fetch(CTG, START - 1, END)       # 0-based half-open
assert len(seq) == NEWLEN, (len(seq), NEWLEN)
with open(f"{DST}/locus.fa", "w") as out:
    out.write(f">{CTG}\n")
    for i in range(0, len(seq), 60):
        out.write(seq[i:i + 60] + "\n")

# ---- bams ----
def shift_bam(src, dst):
    bam = pysam.AlignmentFile(src, "rb")
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"},
           "SQ": [{"SN": CTG, "LN": NEWLEN}]}
    kept = dropped = 0
    out = pysam.AlignmentFile(dst, "wb", header=hdr)
    for r in bam.fetch(CTG, START - 1, END):
        # fully inside only: a read clipped by the window would change the
        # splice graph the component is built from
        if r.reference_start < START - 1 or r.reference_end > END:
            dropped += 1
            continue
        a = pysam.AlignedSegment(out.header)
        a.query_name = r.query_name
        a.query_sequence = r.query_sequence
        a.flag = r.flag
        a.reference_id = 0
        a.reference_start = r.reference_start - OFF
        a.mapping_quality = r.mapping_quality
        a.cigar = r.cigar
        a.next_reference_id = -1
        a.next_reference_start = -1
        a.template_length = 0
        a.query_qualities = r.query_qualities
        a.tags = r.tags
        out.write(a)
        kept += 1
    out.close(); bam.close()
    pysam.index(dst)
    print(f"  {dst}: kept {kept}, dropped {dropped} (window-clipped)")

shift_bam(f"{SRC}/chunk.strand.+.bam", f"{DST}/locus.strand.+.bam")
shift_bam(f"{SRC}/chunk.plus.norm.bam", f"{DST}/locus.plus.norm.bam")

# ---- gtf ----
kept = 0
tx_in = set()
lines = []
for ln in open(f"{SRC}/chunk.strand.+.gtf"):
    if ln.startswith("#"):
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 9 or f[0] != CTG:
        continue
    s, e = int(f[3]), int(f[4])
    if s < START or e > END:
        continue
    f[3], f[4] = str(s - OFF), str(e - OFF)
    lines.append("\t".join(f) + "\n")
    kept += 1
    if f[2] == "transcript":
        tx_in.add(f[8])
with open(f"{DST}/locus.gtf", "w") as out:
    out.writelines(lines)
print(f"  locus.gtf: {kept} lines, {len(tx_in)} transcripts")
print(f"  locus.fa: contig {CTG} length {NEWLEN} (window {START}-{END}, offset -{OFF})")
