#!/usr/bin/env python3

"""Build the TSS and PolyA site table for site-level usage testing.

Starts from LRAA's integrated TSS and PolyA site beds and, for each site:

  - assigns it to a gene: the gene symbols of the transcripts the bed names as
    carrying the site (cluster-guided sites), or else the symbols whose
    transcript span covers the site on its strand. Sites assigned to more than
    one symbol (read-through / cis-fusion models) or to none are kept in the
    table, so their reads are not counted toward a neighbour, but are marked
    competing=False and left out of the test.

  - optionally merges sites of the same gene, kind and strand lying within
    --merge_dist_TSS / --merge_dist_PolyA of each other (single linkage). Off by
    default (0: only sites at the same position merge): LRAA's own site
    definition already absorbs read ends within max_dist_between_alt_*_sites into
    one site, and integrate_TSS_PolyA_sites.py drops a basic site within that
    distance of a cluster-guided one, so the integrated beds hold no two sites of
    one gene closer than that. A merged site is reported at its best-supported
    member (cluster-guided sites ahead of basic ones, since basic support is a
    whole-sample value), keeps the span of its members, and carries the union of
    their transcripts.

  - gives each site a counting window: a read end within this many bp of the
    site's span is counted there, at most --window_TSS / --window_PolyA and at
    most half the gap to the neighbouring site of the same kind and strand, so
    that no read end can count toward two sites.

PolyA sites flagged as internally primed are kept for counting (their reads
should not spill onto a neighbouring site) but marked competing=False.
"""

import argparse
import bisect
import collections
import csv
import logging
import re
import sys

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s", level=logging.INFO)

TSS_COLS = ["chrom", "start0", "pos", "name", "score", "strand", "support", "n_tx", "transcript_ids", "source"]
POLYA_COLS = ["chrom", "start0", "pos", "name", "score", "strand", "support", "n_tx", "transcript_ids",
              "pas", "pas_offset", "internal_priming", "source"]

OUT_COLS = ["site_id", "kind", "chrom", "strand", "pos", "span_lo", "span_hi", "window", "n_merged",
            "merged_positions", "source", "support", "pas", "pas_offset", "internal_priming",
            "gene_symbol", "gene_key", "gene_assignment", "competing", "exclusion", "transcript_ids"]


def main():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--TSS_bed", required=True, help="LRAA integrated TSS site bed")
    parser.add_argument("--PolyA_bed", required=True, help="LRAA integrated PolyA site bed")
    parser.add_argument("--gtf", required=True,
                        help="LRAA gtf with gene symbols prefixed to transcript ids (SYMBOL^t:...), for gene spans")
    parser.add_argument("--gene_trans_map", required=True,
                        help="gene_transcript_splicehashcode.withGeneSymbols.tsv: transcript_id -> new_transcript_id / new_gene_id")
    parser.add_argument("--merge_dist_TSS", type=int, default=0)
    parser.add_argument("--merge_dist_PolyA", type=int, default=0)
    parser.add_argument("--window_TSS", type=int, default=50)
    parser.add_argument("--window_PolyA", type=int, default=25)
    parser.add_argument("--output", required=True, help="site table (tsv)")
    parser.add_argument("--gene_spans_output", required=True,
                        help="gene symbol spans (tsv), used to tally read ends at no site")
    args = parser.parse_args()

    merge_dist = {"TSS": args.merge_dist_TSS, "PolyA": args.merge_dist_PolyA}
    max_window = {"TSS": args.window_TSS, "PolyA": args.window_PolyA}

    tx_symbol = parse_transcript_symbols(args.gene_trans_map)
    spans = parse_symbol_spans(args.gtf)
    write_spans(spans, args.gene_spans_output)
    span_index = index_spans(spans)

    sites = []
    for kind, bed, cols in (("TSS", args.TSS_bed, TSS_COLS), ("PolyA", args.PolyA_bed, POLYA_COLS)):
        n = 0
        for row in read_bed(bed, cols):
            sites.append(assign_gene(kind, row, tx_symbol, span_index))
            n += 1
        logger.info("%s: %d sites read from %s", kind, n, bed)

    merged = merge_sites(sites, merge_dist)
    set_windows(merged, max_window)

    with open(args.output, "wt") as ofh:
        writer = csv.DictWriter(ofh, fieldnames=OUT_COLS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(merged)

    summ = collections.Counter((s["kind"], s["competing"], s["exclusion"]) for s in merged)
    for k, v in sorted(summ.items()):
        logger.info("%s competing=%s %s: %d sites", k[0], k[1], k[2] or "-", v)


def read_bed(filename, cols):
    for line in open(filename):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) != len(cols):
            sys.exit(f"{filename}: expected {len(cols)} columns, got {len(f)}: {line[:200]}")
        row = dict(zip(cols, f))
        row["pos"] = int(row["pos"])
        row["support"] = float(row["support"])
        yield row


def symbol_of(ident):
    return ident.split("^", 1)[0] if "^" in ident else None


def parse_transcript_symbols(filename):
    """transcript_id -> gene symbol: the transcript's own symbol, else its gene's."""
    tx_symbol = {}
    with open(filename) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sym = symbol_of(row["new_transcript_id"]) or symbol_of(row["new_gene_id"])
            if sym:
                tx_symbol[row["transcript_id"]] = sym
    return tx_symbol


def parse_symbol_spans(gtf):
    """(symbol, chrom, strand) -> [lo, hi] over the gtf's transcripts carrying that symbol."""
    spans = {}
    tid_re = re.compile(r'transcript_id "([^"]+)"')
    for line in open(gtf):
        f = line.split("\t", 9)
        if len(f) < 9 or f[2] != "transcript":
            continue
        sym = symbol_of(tid_re.search(f[8]).group(1))
        if sym is None:
            continue
        key = (sym, f[0], f[6])
        lo, hi = int(f[3]), int(f[4])
        if key in spans:
            spans[key][0] = min(spans[key][0], lo)
            spans[key][1] = max(spans[key][1], hi)
        else:
            spans[key] = [lo, hi]
    return spans


def write_spans(spans, filename):
    with open(filename, "wt") as ofh:
        print("\t".join(["gene_symbol", "chrom", "strand", "start", "end", "gene_key"]), file=ofh)
        for (sym, chrom, strand), (lo, hi) in sorted(spans.items(), key=lambda x: (x[0][1], x[1][0])):
            print("\t".join([sym, chrom, strand, str(lo), str(hi), gene_key(sym, chrom, strand)]), file=ofh)


def gene_key(sym, chrom, strand):
    # symbols recur at unrelated loci (paralog copies, PAR genes), so the key carries the locus
    return f"{sym}|{chrom}|{strand}"


def index_spans(spans):
    idx = collections.defaultdict(list)
    for (sym, chrom, strand), (lo, hi) in spans.items():
        idx[(chrom, strand)].append((lo, hi, sym))
    out = {}
    for k, v in idx.items():
        v.sort()
        out[k] = ([x[0] for x in v], v, max(x[1] - x[0] for x in v))
    return out


def symbols_covering(span_index, chrom, strand, pos):
    if (chrom, strand) not in span_index:
        return set()
    starts, v, longest = span_index[(chrom, strand)]
    i = bisect.bisect_right(starts, pos)
    lo_bound = pos - longest
    hits = set()
    j = i - 1
    while j >= 0 and v[j][0] >= lo_bound:
        if v[j][1] >= pos:
            hits.add(v[j][2])
        j -= 1
    return hits


def assign_gene(kind, row, tx_symbol, span_index):
    tids = [t for t in row["transcript_ids"].split(",") if t and t != "."]
    syms = set()
    if row["source"] == "cluster_guided":
        # basic sites name models of the initial catalog, whose ids are not in the map
        syms = {tx_symbol[t] for t in tids if t in tx_symbol}
    how = "transcripts"
    if not syms:
        syms = symbols_covering(span_index, row["chrom"], row["strand"], row["pos"])
        how = "span" if syms else "none"

    site = {
        "kind": kind, "chrom": row["chrom"], "strand": row["strand"], "pos": row["pos"],
        "source": row["source"], "support": row["support"],
        "pas": row.get("pas", ""), "pas_offset": row.get("pas_offset", ""),
        "internal_priming": row.get("internal_priming", ""),
        "transcript_ids": set(tids) if row["source"] == "cluster_guided" else set(),
        "gene_assignment": how,
    }
    if len(syms) == 1:
        sym = next(iter(syms))
        site["gene_symbol"], site["gene_key"], site["exclusion"] = sym, gene_key(sym, row["chrom"], row["strand"]), ""
    else:
        site["gene_symbol"] = ",".join(sorted(syms))
        site["gene_key"] = ""
        site["exclusion"] = "cross_gene" if syms else "no_gene"
    return site


def is_true(x):
    return str(x).strip().lower() in ("true", "1", "1.0")


def merge_sites(sites, merge_dist):
    """single-linkage merge of same-gene sites of one kind and strand within merge_dist"""
    groups = collections.defaultdict(list)
    for s in sites:
        # unassigned sites are never merged: each keeps its own reads
        key = (s["kind"], s["chrom"], s["strand"], s["gene_key"] or f"_unassigned_{id(s)}")
        groups[key].append(s)

    merged = []
    for (kind, chrom, strand, _), members in groups.items():
        members.sort(key=lambda s: s["pos"])
        clusters = [[members[0]]]
        for s in members[1:]:
            if s["pos"] - clusters[-1][-1]["pos"] <= merge_dist[kind]:
                clusters[-1].append(s)
            else:
                clusters.append([s])
        for c in clusters:
            rep = max(c, key=lambda s: (s["source"] == "cluster_guided", s["support"]))
            m = dict(rep)
            m["span_lo"], m["span_hi"] = c[0]["pos"], c[-1]["pos"]
            m["n_merged"] = len(c)
            m["merged_positions"] = ",".join(str(s["pos"]) for s in c)
            m["transcript_ids"] = ",".join(sorted(set().union(*(s["transcript_ids"] for s in c))))
            m["support"] = sum(s["support"] for s in c if s["source"] == rep["source"])
            ip = kind == "PolyA" and is_true(rep["internal_priming"])
            if ip and not m["exclusion"]:
                m["exclusion"] = "internal_priming"
            m["competing"] = not m["exclusion"]
            m["site_id"] = f"{kind}:{chrom}:{m['pos']}:{strand}"
            merged.append(m)

    merged.sort(key=lambda s: (s["kind"], s["chrom"], s["strand"], s["span_lo"]))
    return merged


def set_windows(merged, max_window):
    """window = max_window[kind], shrunk to half the gap to the neighbouring site's span (any gene)"""
    by_strand = collections.defaultdict(list)
    for s in merged:
        by_strand[(s["kind"], s["chrom"], s["strand"])].append(s)
    for (kind, _, _), ss in by_strand.items():
        ss.sort(key=lambda s: s["span_lo"])
        for i, s in enumerate(ss):
            w = max_window[kind]
            if i > 0:
                w = min(w, (s["span_lo"] - ss[i - 1]["span_hi"]) // 2)
            if i + 1 < len(ss):
                w = min(w, (ss[i + 1]["span_lo"] - s["span_hi"]) // 2)
            s["window"] = max(w, 0)


if __name__ == "__main__":
    main()
