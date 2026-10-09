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

Reads are those LRAA itself would use, so these counts match the site-usage
counting (util/sc/site_read_support_to_sparse_matrix.py):
  - LRAA's read filters (Util_funcs.quant_discard_reason): secondary,
    supplementary and duplicate alignments, mapping quality, percent identity
    (--HiFi for LRAA's HiFi floor), over-long introns, and the rDNA mask
    (--rdna_mask_bed);
  - LRAA's transcribed strand (Util_funcs.transcribed_strand): the alignment
    orientation, flipped when minimap2's ts tag is "-" only if the read's own
    junctions carry canonical splice motifs on the flipped strand (checked
    against --genome_fa); a flip the motifs don't corroborate is not applied;
  - LRAA's read geometry (Pretty_alignment): the aligned span's ends as the
    read's 5' and 3' ends, and its reference skips (CIGAR N) > 30 bp as introns,
    deletions staying in the exon. Read straight from the CIGAR rather than by
    building a Pretty_alignment per read, which is most of a read's cost.
Only reads on the pair's transcript strand are counted.

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

import os

import pysam

sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "../../../pylib"]))

import LRAA_Globals  # noqa: E402
import RdnaMask  # noqa: E402
import Util_funcs  # noqa: E402

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
    parser.add_argument("--HiFi", action="store_true",
                        help="apply LRAA's HiFi read-identity floor (min_per_id {}), as the LRAA run given --HiFi "
                             "did".format(LRAA_Globals.HIFI_MIN_PER_ID))
    parser.add_argument("--min_per_id", type=float, default=None,
                        help="override the percent-identity floor (default: LRAA config, or the HiFi one)")
    parser.add_argument("--min_mapping_quality", type=int, default=None,
                        help="override the mapping-quality floor (default: LRAA config min_mapping_quality)")
    parser.add_argument("--rdna_mask_bed", default=None,
                        help="rDNA mask bed LRAA built for this genome: reads overlapping it are discarded, as LRAA "
                             "discards them")
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

    read_filter = make_read_filter(HiFi=args.HiFi, min_per_id=args.min_per_id,
                                   min_mapping_quality=args.min_mapping_quality, rdna_mask_bed=args.rdna_mask_bed)

    # a contig at a time: the strand check holds one contig's sequence
    out_rows = [None] * len(candidates)
    order = sorted(range(len(candidates)), key=lambda i: contig_of[candidates[i]["dominant_transcript_ids"]])
    for i in order:
        c = candidates[i]
        res = check_pair(c, exons, strand_of, contig_of, bam, genome, cell_to_cluster, args.site_window,
                         include_terminal_exon_reads=not args.spliced_only, read_filter=read_filter)
        out_rows[i] = {**c, **res}
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


def read_geometry(read):
    """(lend, rend, introns > MIN_READ_INTRON_LEN) of a read, 1-based, as Pretty_alignment
    gives them: introns are the CIGAR's reference skips (N, consecutive ones joined);
    deletions (D) stay inside the exon, however long."""
    introns = []
    pos = read.reference_start  # 0-based
    skip_start = None
    for op, n in read.cigartuples:
        if op == 3:  # N
            if skip_start is None:
                skip_start = pos
            pos += n
            continue
        if skip_start is not None:
            if pos - skip_start > MIN_READ_INTRON_LEN:
                introns.append((skip_start + 1, pos))
            skip_start = None
        if op in (0, 2, 7, 8):  # M, D, =, X consume the reference
            pos += n
    return read.reference_start + 1, read.reference_end, introns


def make_read_filter(HiFi=False, min_per_id=None, min_mapping_quality=None, rdna_mask_bed=None):
    """keep(read) -> True for a read LRAA would use (Util_funcs.quant_discard_reason), with
    the identity / mapping-quality floors and rDNA mask the site-usage counting applies."""
    if min_per_id is None:
        min_per_id = LRAA_Globals.HIFI_MIN_PER_ID if HiFi else LRAA_Globals.config["min_per_id"]
    if min_mapping_quality is None:
        min_mapping_quality = int(LRAA_Globals.config["min_mapping_quality"])
    # {} rather than None: None makes quant_discard_reason read the LRAA run's config mask
    rdna_mask = (RdnaMask.load_mask_bed(rdna_mask_bed) if rdna_mask_bed else None) or {}

    def keep(read):
        return Util_funcs.quant_discard_reason(read, None, min_mapping_quality=min_mapping_quality,
                                               min_per_id=min_per_id, rdna_mask=rdna_mask) is None
    return keep


def register_contig(genome, contig):
    """Hold the contig's sequence for transcript_strand's splice-motif check (one contig at
    a time; a contig missing from the fasta keeps every read's aligned strand)."""
    if genome is None or Util_funcs.contig_seq_for_strand_check(contig) is not None:
        return
    if contig in genome.references:
        Util_funcs.register_contig_seq_for_strand_check(contig, genome.fetch(contig).upper())


def check_pair(c, exons, strand_of, contig_of, bam, genome, cell_to_cluster, window,
               include_terminal_exon_reads=True, read_filter=None):
    """read_filter: keep(read) from make_read_filter (default: LRAA's non-HiFi filters)."""

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

    if read_filter is None:
        read_filter = make_read_filter()
    register_contig(genome, contig)

    ends = []  # (position, cluster, carries the key intron)
    for read in bam.fetch(contig, fetch_lo, fetch_hi):
        if not read.has_tag("CB") or transcript_strand(read, genome) != strand:
            continue
        read_lend, read_rend, read_introns = read_geometry(read)
        carries_intron = any(abs(i[0] - key_intron[0]) <= INTRON_TOLERANCE and abs(i[1] - key_intron[1]) <= INTRON_TOLERANCE
                             for i in read_introns)
        if not carries_intron:
            if not include_terminal_exon_reads or read_introns:
                continue
            # unspliced, and wholly on the terminal-exon side of the key intron
            if five_prime_side == plus:
                in_terminal_exon = read_rend < key_intron[0]
            else:
                in_terminal_exon = read_lend > key_intron[1]
            if not in_terminal_exon:
                continue

        if five_prime_side:
            pos = read_lend if plus else read_rend
            beyond = pos < key_intron[0] if plus else pos > key_intron[1]
        else:
            pos = read_rend if plus else read_lend
            beyond = pos > key_intron[1] if plus else pos < key_intron[0]
        if not beyond:
            continue  # ends inside the intron's flank on the wrong side: not a terminal-exon end
        # LRAA's read filters last: the costliest test, needed only by reads that would count
        if not read_filter(read):
            continue

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


def transcript_strand(read, genome=None):
    """Genomic strand of the transcript a read came from, as LRAA assigns it
    (Util_funcs.transcribed_strand): the alignment orientation, flipped for ts:A:- only when
    the read's junctions carry canonical motifs on the flipped strand. With `genome` (a
    pysam.FastaFile) the read's contig is registered for that check; without it, only an
    already-registered contig can corroborate a flip, and otherwise the aligned strand
    stands."""
    if genome is not None:
        register_contig(genome, read.reference_name)
    return Util_funcs.transcribed_strand(read)


def downstream_A(genome, contig, pos, plus):
    """A's (T's on '-') in the 20 genomic bases 3' of a PolyA site: the oligo-dT priming template."""
    if plus:
        seq = genome.fetch(contig, pos, pos + 20)
        return seq.upper().count("A")
    seq = genome.fetch(contig, max(0, pos - 21), pos - 1)
    return seq.upper().count("T")


if __name__ == "__main__":
    main()
