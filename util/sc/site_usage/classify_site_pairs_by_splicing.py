#!/usr/bin/env python3

"""Split site-switching events into alternative terminal usage and switches
coupled with alternative splicing, from the intron chains of the reads at each
site.

An event (annotate_site_usage_events.py) moves read ends from a lost site to a
gained site. Whether that is a change of terminus alone or comes with different
splicing is read off the reads themselves: the reads starting (TSS) or ending
(PolyA) at each site, assigned exactly as count_site_read_ends.py assigns them
(nearest site of the kind on the read's transcript strand, within its window),
and their introns (CIGAR N operations). Reads are pooled over all clustered
cells: the question is which splicing each site goes with, not where it is used.

For each (gene, gained site, lost site) pair, one site is the inner one, nearer
the gene body (the proximal PolyA, the downstream TSS), the other the outer one.
Reads from the outer site pass the inner site's position on their way into the
gene, so they show whether the two sites share a terminal exon:

  inner adjacent intron  the inner site's most common intron next to the varying
                   terminus (the read's first intron for a TSS, its last for a
                   PolyA).
  outer share carrying it  among the outer site's reads spanning that intron
                   (reaching past both of its ends), the share carrying it.
  outer share spliced between sites  among the outer site's reads reaching
                   past the inner site or carrying an intron between the two, the
                   share with an intron between them (or over the inner site): the
                   inner site lies in an intron or an internal exon of the outer
                   site's transcripts, not on their terminal exon.
  splicing divergence  for sites sharing a terminal exon, over the introns seen
                   at either site, the largest difference between the two sites in
                   the share of reads carrying the intron among the reads spanning
                   it. Conditioning on spanning matters throughout: reads at a
                   distal PolyA site are longer and lose more of their 5' ends, so
                   raw intron frequencies would fall at the distal site and look
                   like a splicing difference. Introns spanned by fewer than
                   --min_spanning reads at either site are skipped.

  class            alt_splicing:terminal_exon: the outer reads splice between the
                     two sites, or don't carry its adjacent intron
                     (share < --min_adjacent_share): an alternative first or last
                     exon, an intronic PolyA, a retained intron.
                   alt_splicing:internal: same terminal exon, but an intron
                     further in differs by >= --min_divergence.
                   alt_terminal_usage: same terminal exon, same splicing: tandem
                     TSS / tandem 3' UTR.
                   unspliced_site: fewer than --min_spliced spliced reads at the
                     inner site (a monoexonic model, reads inside an intron): read
                     ends alone can't tell splicing from pre-mRNA there.
                   unresolved: too few reads, or too few outer reads reaching the
                     inner site's adjacent intron (long 3' UTRs whose reads start
                     inside the UTR).

Reads per site are capped at --max_reads_per_site (the first ones met), so highly
expressed genes don't dominate the run.
"""

import argparse
import collections
import csv
import logging
from multiprocessing import Pool

import numpy as np
import pandas as pd
import pysam

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s", level=logging.INFO)

G = {}


def main():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--events", required=True, nargs="+",
                        help="event tables from annotate_site_usage_events.py, one per kind, as KIND=path")
    parser.add_argument("--sites", required=True, help="site table from prep_site_table.py")
    parser.add_argument("--gene_spans", required=True, help="gene spans from prep_site_table.py")
    parser.add_argument("--bam", required=True, help="aligned reads with CB cell-barcode tags (indexed)")
    parser.add_argument("--cell_clusters", required=True, help="as given to count_site_read_ends.py")
    parser.add_argument("--max_reads_per_site", type=int, default=3000)
    parser.add_argument("--min_spliced", type=int, default=10)
    parser.add_argument("--min_spanning", type=int, default=10)
    parser.add_argument("--min_adjacent_share", type=float, default=0.5)
    parser.add_argument("--min_divergence", type=float, default=0.25)
    parser.add_argument("--CPU", type=int, default=8)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    pairs = []
    for spec in args.events:
        kind, path = spec.split("=", 1)
        e = pd.read_csv(path, sep="\t", usecols=["gene_key", "gene_symbol", "gained_site", "lost_site"])
        e = e.drop_duplicates(["gene_key", "gained_site", "lost_site"])
        e["kind"] = kind
        pairs.append(e)
    pairs = pd.concat(pairs, ignore_index=True)
    logger.info("%d site pairs in %d genes", len(pairs), pairs.gene_key.nunique())

    sites = list(csv.DictReader(open(args.sites), delimiter="\t"))
    spans = {s["gene_key"]: (int(s["start"]), int(s["end"]))
             for s in csv.DictReader(open(args.gene_spans), delimiter="\t")}

    G["bam"] = args.bam
    G["cells"] = set(parse_cells(args.cell_clusters))
    G["site_index"] = index_sites(sites)
    G["site_id"] = [s["site_id"] for s in sites]
    G["site_row"] = {s["site_id"]: i for i, s in enumerate(sites)}
    G["site_span"] = {s["site_id"]: (int(s["span_lo"]), int(s["span_hi"]), int(s["window"])) for s in sites}
    G["spans"] = spans
    G["args"] = args

    jobs = [(gk, g.kind.iloc[0] if g.kind.nunique() == 1 else None, g) for gk, g in pairs.groupby("gene_key")]
    # a gene can have events of both kinds; split them so each job reads one kind of end
    jobs = [(gk, kind, g[g.kind == kind]) for gk, _, g in jobs for kind in g.kind.unique()]
    # largest genes first, so the long ones don't finish last
    jobs.sort(key=lambda j: -(spans.get(j[0], (0, 0))[1] - spans.get(j[0], (0, 0))[0]))

    out = []
    with Pool(args.CPU) as pool:
        for i, rows in enumerate(pool.imap_unordered(classify_gene, jobs, chunksize=4)):
            out.extend(rows)
            if (i + 1) % 200 == 0:
                logger.info("%d / %d gene jobs", i + 1, len(jobs))

    out = pd.DataFrame(out)
    out.to_csv(args.output, sep="\t", index=False)
    for kind, g in out.groupby("kind"):
        logger.info("%s: %s", kind, g.splicing_class.value_counts().to_dict())


def parse_cells(filename):
    for line in open(filename):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[1].strip().lstrip("-").isdigit():
            yield f[0]


def index_sites(sites):
    """(kind, contig, strand) -> sorted span starts, span ends, windows, row numbers (as count_site_read_ends.py)"""
    idx = collections.defaultdict(list)
    for i, s in enumerate(sites):
        idx[(s["kind"], s["chrom"], s["strand"])].append((int(s["span_lo"]), int(s["span_hi"]), int(s["window"]), i))
    out = {}
    for k, v in idx.items():
        v.sort()
        out[k] = tuple(np.array(x) for x in zip(*v))
    return out


def nearest_site(index, pos):
    """identical to count_site_read_ends.nearest_site"""
    import bisect
    lo, hi, win, row = index
    i = bisect.bisect_right(lo, pos)
    best, best_d = None, None
    for j in (i - 1, i):
        if 0 <= j < len(lo):
            d = 0 if lo[j] <= pos <= hi[j] else min(abs(pos - lo[j]), abs(pos - hi[j]))
            if d <= win[j] and (best_d is None or d < best_d):
                best, best_d = row[j], d
    return best


def transcript_strand(read):
    strand = "-" if read.is_reverse else "+"
    if read.has_tag("ts") and read.get_tag("ts") == "-":
        strand = "+" if strand == "-" else "-"
    return strand


def read_introns(read):
    introns, pos = [], read.reference_start
    for op, length in read.cigartuples:
        if op == 3:  # N
            introns.append((pos + 1, pos + length))
            pos += length
        elif op in (0, 2, 7, 8):  # M D = X
            pos += length
    return tuple(introns)


def classify_gene(job):
    gene_key, kind, pairs = job
    args = G["args"]
    _, contig, strand = gene_key.split("|")
    plus = strand == "+"
    want = set(pairs.gained_site) | set(pairs.lost_site)
    want_rows = {G["site_row"][s]: s for s in want}

    lo, hi = G["spans"].get(gene_key, (None, None))
    for s in want:
        s_lo, s_hi, w = G["site_span"][s]
        lo = s_lo - w if lo is None else min(lo, s_lo - w)
        hi = s_hi + w if hi is None else max(hi, s_hi + w)

    index = G["site_index"][(kind, contig, strand)]
    reads = collections.defaultdict(list)  # site_id -> [(span_lo, span_hi, introns)]
    full = set()
    with pysam.AlignmentFile(G["bam"]) as bam:
        for read in bam.fetch(contig, max(0, lo - 1), hi):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or not read.has_tag("CB"):
                continue
            if transcript_strand(read) != strand or read.get_tag("CB") not in G["cells"]:
                continue
            five_prime = kind == "TSS"
            pos = (read.reference_start + 1 if plus else read.reference_end) if five_prime else \
                  (read.reference_end if plus else read.reference_start + 1)
            row = nearest_site(index, pos)
            sid = want_rows.get(row)
            if sid is None or sid in full:
                continue
            reads[sid].append((read.reference_start + 1, read.reference_end, read_introns(read)))
            if len(reads[sid]) >= args.max_reads_per_site:
                full.add(sid)
                if full == want:
                    break

    profiles = {}
    for s in want:
        s_lo, s_hi, _ = G["site_span"][s]
        profiles[s] = {**site_profile(reads[s], kind, plus), "pos": (s_lo + s_hi) // 2, "kind": kind, "plus": plus}
    rows = []
    for p in pairs.itertuples():
        rows.append({"kind": kind, "gene_key": gene_key, "gene_symbol": p.gene_symbol,
                     "gained_site": p.gained_site, "lost_site": p.lost_site,
                     **compare_sites(profiles[p.gained_site], profiles[p.lost_site], args)})
    return rows


def site_profile(reads, kind, plus):
    """per-site read counts, adjacent-intron counts, and the reads' spans and intron sets"""
    spliced = [r for r in reads if r[2]]
    adjacent = collections.Counter()
    for _, _, introns in spliced:
        # the intron next to the varying terminus: the transcript's first (TSS) or last (PolyA) intron
        first_in_transcript = introns[0] if plus else introns[-1]
        last_in_transcript = introns[-1] if plus else introns[0]
        adjacent[first_in_transcript if kind == "TSS" else last_in_transcript] += 1
    return {"n_reads": len(reads), "n_spliced": len(spliced), "adjacent": adjacent,
            "spans": np.array([(r[0], r[1]) for r in reads]) if reads else np.zeros((0, 2), dtype=int),
            "intron_sets": [set(r[2]) for r in reads]}


def intron_share(profile, intron):
    """share of the site's reads spanning the intron that carry it, and how many span it"""
    if not profile["n_reads"]:
        return None, 0
    spans = profile["spans"]
    spanning = (spans[:, 0] < intron[0]) & (spans[:, 1] > intron[1])
    n = int(spanning.sum())
    if not n:
        return None, 0
    carried = sum(1 for k, s in enumerate(profile["intron_sets"]) if spanning[k] and intron in s)
    return carried / n, n


def compare_sites(g, l, args):
    """g, l: the gained and lost sites' profiles"""
    # the inner site is the one nearer the gene body (proximal PolyA, downstream TSS); reads from the
    # outer site pass through the inner site's position on their way into the gene
    outer_is_gained = (g["pos"] > l["pos"]) == ((g["kind"] == "PolyA") == g["plus"])
    outer, inner = (g, l) if outer_is_gained else (l, g)

    res = {"gained_reads": g["n_reads"], "lost_reads": l["n_reads"],
           "gained_spliced_frac": frac(g["n_spliced"], g["n_reads"]),
           "lost_spliced_frac": frac(l["n_spliced"], l["n_reads"]),
           "inner_site": "gained" if outer is l else "lost",
           "inner_adjacent_intron": "", "outer_reads_spanning_it": 0, "outer_share_carrying_it": "",
           "outer_share_spliced_between_sites": "",
           "splicing_divergence": "", "divergent_intron": "",
           "divergent_intron_share_gained": "", "divergent_intron_share_lost": ""}

    if inner["n_spliced"] < args.min_spliced:
        # the inner site's reads are unspliced: a monoexonic model, or reads inside an intron or a
        # retained intron. Whether that is splicing or pre-mRNA can't be told from read ends.
        res["splicing_class"] = "unspliced_site" if inner["n_reads"] >= args.min_spliced else "unresolved"
        return res

    # outer reads with an intron between the two sites (or over the inner one): the two sites are not
    # on one terminal exon of the outer site's transcripts. The inner site sits in an intron of those
    # transcripts (an alternative first exon, an intronic or alternative last exon PolyA), or in one of
    # their internal exons (an alternative last exon, a downstream TSS of 5'-truncated reads)
    over, n_over = spliced_between_share(outer, inner["pos"])
    res["outer_share_spliced_between_sites"] = round(over, 3) if n_over else ""

    a_inner, _ = inner["adjacent"].most_common(1)[0]
    share, n_span = intron_share(outer, a_inner)
    res.update({"inner_adjacent_intron": fmt_intron(a_inner), "outer_reads_spanning_it": n_span,
                "outer_share_carrying_it": round(share, 3) if n_span else ""})

    if n_over >= args.min_spanning and over >= args.min_adjacent_share:
        res["splicing_class"] = "alt_splicing:terminal_exon"
        return res
    if n_span < args.min_spanning:
        res["splicing_class"] = "unresolved"  # too few outer reads reach the inner site's terminal intron
        return res
    if share < args.min_adjacent_share:
        res["splicing_class"] = "alt_splicing:terminal_exon"
        return res

    # same terminal exon: compare the splicing further in, over introns both sites' reads span
    seen = collections.Counter()
    for p in (g, l):
        c = collections.Counter(i for s in p["intron_sets"] for i in s)
        seen.update({i for i, n in c.items() if n >= 0.1 * p["n_spliced"]})
    best = (0.0, None, None, None)
    for intron in seen:
        sg, ng = intron_share(g, intron)
        sl, nl = intron_share(l, intron)
        if ng < args.min_spanning or nl < args.min_spanning:
            continue
        d = abs(sg - sl)
        if d > best[0]:
            best = (d, intron, sg, sl)
    res["splicing_divergence"] = round(best[0], 3)
    if best[1] is not None:
        res.update({"divergent_intron": fmt_intron(best[1]),
                    "divergent_intron_share_gained": round(best[2], 3), "divergent_intron_share_lost": round(best[3], 3)})
    res["splicing_class"] = "alt_splicing:internal" if best[0] >= args.min_divergence else "alt_terminal_usage"
    return res


def spliced_between_share(outer, inner_pos):
    """share of the outer site's informative reads with an intron between the two sites (or over
    the inner one). Informative: reaching past the inner site (the reads end at the outer site, so
    that spans both), or carrying an intron between the sites even if they stop short of the inner
    one, as reads of a distant alternative last / first exon do. Unspliced reads stopping short say
    nothing."""
    spans = outer["spans"]
    if not len(spans):
        return 0.0, 0
    lo, hi = sorted((outer["pos"], inner_pos))
    spanning = (spans[:, 0] < inner_pos) & (spans[:, 1] > inner_pos)
    n = over = 0
    for k, s in enumerate(outer["intron_sets"]):
        between = any(a <= hi and b >= lo for a, b in s)
        if spanning[k] or between:
            n += 1
            over += between
    return (over / n if n else 0.0), n


def frac(a, b):
    return round(a / b, 3) if b else ""


def fmt_intron(i):
    return f"{i[0]}-{i[1]}"


if __name__ == "__main__":
    main()
