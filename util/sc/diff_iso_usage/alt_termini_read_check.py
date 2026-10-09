#!/usr/bin/env python3

"""Check alt-termini DTU calls against where the reads actually start or end.

Two isoforms that share a splice pattern and differ only at one terminus are
told apart by the EM mostly through reads compatible with both, so a significant
DTU call can rest on how those reads were apportioned rather than on reads that
reach either terminus. This script goes back to the alignments: for each
candidate pair it collects the reads carrying the pair's intron next to the
varying terminus, takes each read's 5' end (alt TSS) or 3' end (alt PolyA), and
counts, per cell cluster, how many land at the dominant vs the alternate
terminus.

A call is borne out by the reads when the per-cluster share of reads ending at
the dominant terminus moves between the pair's two clusters in the same
direction, and by a comparable amount, as the model's within-pair share
dominant_pi / (dominant_pi + alternate_pi).

Two kinds of read are counted: those carrying the intron next to the varying
terminus (+/- 2 bp), and, unless --spliced_only, unspliced reads lying wholly in
the terminal exon region beyond that intron (the last exon's 3' UTR for an alt
PolyA, the first exon for an alt TSS). Many reads ending at a distal polyA site
start inside a long 3' UTR and never reach the last intron, so counting only
intron-carrying reads understates distal-site use. Reads ending at an unmodelled
site (an internally primed A-run, say) are still counted: they come from an
expressed isoform, only their end cannot be trusted as a site.

Reads are kept only on the pair's transcript strand. minimap2's ts tag gives the
transcript strand relative to the read, so the genomic transcript strand is the
read's alignment orientation flipped when ts is "-"; reads without the tag are
taken as oriented to the transcript, as long reads from oriented libraries are.

read_frac_elsewhere, the share of reads ending away from both modelled termini,
is computed over the intron-carrying reads only: reads in a long terminal exon
also end at unmodelled sites (internally primed A-runs among them), which says
nothing about whether the modelled termini are where the pair's reads end.
"""

import argparse
import collections
import csv
import logging
import sys

import pysam

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s", level=logging.INFO)

MIN_READ_INTRON_LEN = 30  # shorter alignment gaps are deletions, not introns
INTRON_TOLERANCE = 2


def main():

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--candidates", required=True,
                        help="tsv: gene_symbol, alt_terminus (TSS|PolyA), dominant_transcript_ids, "
                             "alternate_transcript_ids, cluster_A, cluster_B, plus any other columns (passed through)")
    parser.add_argument("--gtf", required=True, help="LRAA gtf holding the candidate transcripts (segmented is fine)")
    parser.add_argument("--bam", required=True, help="aligned reads with CB cell-barcode tags")
    parser.add_argument("--cell_clusters", required=True,
                        help="cell_barcode <tab> cluster number, as given to the DTU test (a header line is skipped)")
    parser.add_argument("--genome_fa", required=True, help="genome fasta, for the A-content downstream of PolyA sites")
    parser.add_argument("--site_window", type=int, default=50,
                        help="a read end within this many bp of a terminus counts as ending there "
                             "(at most half the distance between the pair's two termini)")
    parser.add_argument("--spliced_only", action="store_true",
                        help="count only reads carrying the intron next to the varying terminus "
                             "(the behaviour before terminal-exon reads were counted)")
    parser.add_argument("--output", required=True, help="output tsv: the candidate columns plus read-level columns")
    args = parser.parse_args()

    candidates = list(csv.DictReader(open(args.candidates), delimiter="\t"))
    want = {c[k] for c in candidates for k in ("dominant_transcript_ids", "alternate_transcript_ids")}

    exons, strand_of, contig_of = parse_gtf_exons(args.gtf, want)
    missing = want - set(exons)
    if missing:
        sys.exit(f"transcripts not found in {args.gtf}: {sorted(missing)[:5]} ...")

    cell_to_cluster = parse_cell_clusters(args.cell_clusters)
    bam = pysam.AlignmentFile(args.bam)
    genome = pysam.FastaFile(args.genome_fa)

    out_rows = []
    for c in candidates:
        res = check_pair(c, exons, strand_of, contig_of, bam, genome, cell_to_cluster, args.site_window,
                         include_terminal_exon_reads=not args.spliced_only)
        out_rows.append({**c, **res})
        logger.info("%s %s: %d reads, dominant share %s -> %s",
                    c["gene_symbol"], c["alt_terminus"], res["n_reads"],
                    res["read_dom_share_A"], res["read_dom_share_B"])

    with open(args.output, "wt") as ofh:
        writer = csv.DictWriter(ofh, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)


def parse_gtf_exons(gtf, want):
    exons = collections.defaultdict(list)
    strand_of, contig_of = {}, {}
    for line in open(gtf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        tid = f[8].split('transcript_id "', 1)[1].split('"', 1)[0]
        if tid not in want:
            continue
        exons[tid].append((int(f[3]), int(f[4])))
        strand_of[tid], contig_of[tid] = f[6], f[0]

    # a segmented gtf splits exons at shared boundaries; rejoin abutting segments
    merged = {}
    for tid, segs in exons.items():
        segs.sort()
        m = [list(segs[0])]
        for s, e in segs[1:]:
            if s <= m[-1][1] + 1:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        merged[tid] = m
    return merged, strand_of, contig_of


def parse_cell_clusters(filename):
    cell_to_cluster = {}
    for line in open(filename):
        f = line.rstrip("\n").split("\t")
        if len(f) < 2 or not f[1].strip().lstrip("-").isdigit():
            continue  # header
        cell_to_cluster[f[0]] = "Cluster_" + f[1].strip()
    return cell_to_cluster


def check_pair(c, exons, strand_of, contig_of, bam, genome, cell_to_cluster, window,
               include_terminal_exon_reads=True):

    kind = c["alt_terminus"]
    if kind not in ("TSS", "PolyA"):
        raise ValueError(f"alt_terminus must be TSS or PolyA, got {kind!r}")

    dom, alt = c["dominant_transcript_ids"], c["alternate_transcript_ids"]
    strand, contig = strand_of[dom], contig_of[dom]
    plus = strand == "+"

    dom_exons = exons[dom]
    introns = [(a[1] + 1, b[0] - 1) for a, b in zip(dom_exons, dom_exons[1:])]
    if not introns:
        raise ValueError(f"{dom} is monoexonic; nothing ties reads to it")

    # the intron next to the varying terminus, and the terminus coordinate for each isoform
    five_prime_side = kind == "TSS"
    key_intron = introns[0] if five_prime_side == plus else introns[-1]

    def terminus(tid):
        lo, hi = exons[tid][0][0], exons[tid][-1][1]
        return (lo if plus else hi) if five_prime_side else (hi if plus else lo)

    dom_pos, alt_pos = terminus(dom), terminus(alt)

    fetch_lo = max(0, min(dom_pos, alt_pos, key_intron[0]) - 200)
    fetch_hi = max(dom_pos, alt_pos, key_intron[1]) + 200

    ends = []  # (position, cluster, carries the key intron)
    for read in bam.fetch(contig, fetch_lo, fetch_hi):
        if read.is_secondary or read.is_supplementary or not read.has_tag("CB"):
            continue
        if transcript_strand(read) != strand:
            continue
        blocks = read.get_blocks()
        read_introns = [(a[1] + 1, b[0]) for a, b in zip(blocks, blocks[1:])
                        if b[0] - a[1] > MIN_READ_INTRON_LEN]
        carries_intron = any(abs(i[0] - key_intron[0]) <= INTRON_TOLERANCE and abs(i[1] - key_intron[1]) <= INTRON_TOLERANCE
                             for i in read_introns)
        if not carries_intron:
            if not include_terminal_exon_reads or read_introns:
                continue
            # unspliced, and wholly on the terminal-exon side of the key intron
            if five_prime_side == plus:
                in_terminal_exon = read.reference_end < key_intron[0]
            else:
                in_terminal_exon = read.reference_start + 1 > key_intron[1]
            if not in_terminal_exon:
                continue

        if five_prime_side:
            pos = read.reference_start + 1 if plus else read.reference_end
            beyond = pos < key_intron[0] if plus else pos > key_intron[1]
        else:
            pos = read.reference_end if plus else read.reference_start + 1
            beyond = pos > key_intron[1] if plus else pos < key_intron[0]
        if not beyond:
            continue  # ends inside the intron's flank on the wrong side: not a terminal-exon end

        ends.append((pos, cell_to_cluster.get(read.get_tag("CB")), carries_intron))

    # for closely spaced termini the window shrinks to half their separation, so a read end
    # can count toward only one of them
    pair_window = min(window, abs(dom_pos - alt_pos) // 2)

    def site(pos):
        if abs(pos - dom_pos) <= pair_window:
            return "dom"
        if abs(pos - alt_pos) <= pair_window:
            return "alt"
        return "other"

    sites = [(site(p), cl) for p, cl, _ in ends]
    n = len(sites)
    count = collections.Counter(s for s, _ in sites)
    spliced = [site(p) for p, _, carries in ends if carries]
    n_spliced = len(spliced)
    count_spliced = collections.Counter(spliced)

    def share(cluster):
        nd = sum(1 for s, cl in sites if cl == cluster and s == "dom")
        na = sum(1 for s, cl in sites if cl == cluster and s == "alt")
        return nd, na, (round(nd / (nd + na), 3) if nd + na else "")

    a_dom, a_alt, a_share = share(c["cluster_A"])
    b_dom, b_alt, b_share = share(c["cluster_B"])

    peaks = collections.Counter((p // 10) * 10 for p, _, _ in ends).most_common(3)

    res = {
        "dom_terminus_pos": dom_pos,
        "alt_terminus_pos": alt_pos,
        "n_reads": n,
        "n_reads_with_key_intron": n_spliced,
        "read_frac_at_dom": round(count["dom"] / n, 3) if n else "",
        "read_frac_at_alt": round(count["alt"] / n, 3) if n else "",
        # over the intron-carrying reads; see the module docstring
        "read_frac_elsewhere": round(count_spliced["other"] / n_spliced, 3) if n_spliced else "",
        "read_frac_elsewhere_all_reads": round(count["other"] / n, 3) if n else "",
        "reads_dom_A": a_dom, "reads_alt_A": a_alt, "read_dom_share_A": a_share,
        "reads_dom_B": b_dom, "reads_alt_B": b_alt, "read_dom_share_B": b_share,
        "top_read_end_peaks": ",".join(f"{p}:{k}" for p, k in peaks),
    }

    if kind == "PolyA":
        res["dom_downstream_A_of_20"] = downstream_A(genome, contig, dom_pos, plus)
        res["alt_downstream_A_of_20"] = downstream_A(genome, contig, alt_pos, plus)
    else:
        res["dom_downstream_A_of_20"] = res["alt_downstream_A_of_20"] = ""

    return res


def transcript_strand(read):
    """Genomic strand of the transcript a read came from: its alignment orientation, flipped
    when minimap2's ts tag says the read is antisense to the transcript; reads without the tag
    are taken as oriented to the transcript."""
    strand = "-" if read.is_reverse else "+"
    if read.has_tag("ts") and read.get_tag("ts") == "-":
        strand = "+" if strand == "-" else "-"
    return strand


def downstream_A(genome, contig, pos, plus):
    """A's (T's on '-') in the 20 genomic bases 3' of a PolyA site: the oligo-dT priming template."""
    if plus:
        seq = genome.fetch(contig, pos, pos + 20)
        return seq.upper().count("A")
    seq = genome.fetch(contig, max(0, pos - 21), pos - 1)
    return seq.upper().count("T")


if __name__ == "__main__":
    main()
