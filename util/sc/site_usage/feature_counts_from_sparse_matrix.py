#!/usr/bin/env python3

"""Turn an LRAA single-cell feature x cell sparse matrix (EM-assigned isoform or
splice-pattern counts) into the count files and feature table site_usage_dexseq.R reads,
so the site-usage DEXSeq pseudo-replicate test can be run on other features.

Two modes:

  splice_pattern   (kind SplicePattern) features = the splice patterns (intron chains) of
                   the splice-pattern matrix (SYMBOL^hash ids), grouped by gene
                   (gene_key SYMBOL|chrom|strand). Unspliced / monoexonic entries (ids with
                   ':iso-') are left out, as the chi-square test was run with
                   --ignore_unspliced. A pattern whose isoforms carry more than one gene
                   symbol is kept in the table but marked competing False (not tested).

  isoform_termini  (kind IsoformTermini) features = (splice pattern, annotated TSS site,
                   annotated PolyA site) combinations: the EM counts of the multi-exon
                   isoforms sharing a splice pattern and the same TSS site and PolyA site
                   (sites from the site table's transcript_ids; 'none' when the isoform
                   carries no site at that end) are summed. The group tested (gene_key) is
                   the splice pattern, SYMBOL^hash|chrom|strand. Combinations with no
                   annotated site at either end are left out.

Outputs, under --output_prefix (the layout site_counts_from_lraa_support.py writes):
  .site_counts.mtx.gz      features x cells MatrixMarket (real: EM counts are fractional)
  .sites.tsv.gz            feature ids, in row order
  .barcodes.tsv.gz         clustered cell barcodes, in column order, with their cluster
  .features.tsv            the feature table (--sites for site_usage_dexseq.R): site_id, kind,
                           competing, gene_key, gene_symbol, chrom, strand, member isoforms, ...
  .cluster_library.tsv     all of the matrix's counts per cluster (every feature, before any
                           filter): the library for the annotator's CPM-based switch class
"""

import argparse
import collections
import gzip
import logging
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from site_counts_from_lraa_support import parse_cell_clusters  # noqa: E402

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s", level=logging.INFO)

KIND = {"splice_pattern": "SplicePattern", "isoform_termini": "IsoformTermini"}


def main():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mode", required=True, choices=sorted(KIND))
    parser.add_argument("--sparseM_dir", required=True,
                        help="10x-style dir (features.tsv.gz, barcodes.tsv.gz, matrix.mtx.gz): the splice-pattern "
                             "matrix (splice_pattern) or the isoform matrix (isoform_termini)")
    parser.add_argument("--gene_transcript_splicehash", required=True,
                        help="gene_transcript_splicehashcode.withGeneSymbols.tsv: isoform -> splice pattern, exon count")
    parser.add_argument("--cell_clusters", required=True, help="cell_barcode <tab> cluster; a header line is skipped")
    parser.add_argument("--sp_gtf", default=None,
                        help="splice-pattern collapsed gtf (splice_pattern mode): each pattern's main TSS / PolyA site")
    parser.add_argument("--sites", default=None,
                        help="site table from prep_site_table.py (isoform_termini mode): the TSS / PolyA site each "
                             "isoform carries (transcript_ids)")
    parser.add_argument("--output_prefix", required=True)
    args = parser.parse_args()

    kind = KIND[args.mode]
    barcodes, cluster_of = parse_cell_clusters(args.cell_clusters)
    m, features = read_sparse_matrix(args.sparseM_dir, barcodes)
    clusters = [cluster_of[b] for b in barcodes]
    write_library(m, clusters, args.output_prefix + ".cluster_library.tsv")

    tx = pd.read_csv(args.gene_transcript_splicehash, sep="\t", dtype=str)
    tx["symbol"] = [own_symbol(t) for t in tx.new_transcript_id]
    tx["chrom"] = tx.gene_id.str.split(":").str[1]
    tx["strand"] = tx.gene_id.str.split(":").str[2]
    tx["num_exons"] = tx.num_exons.astype(int)

    if args.mode == "splice_pattern":
        if not args.sp_gtf:
            sys.exit("--sp_gtf is required for --mode splice_pattern")
        table, rows, agg = splice_pattern_features(features, tx, args.sp_gtf)
    else:
        if not args.sites:
            sys.exit("--sites is required for --mode isoform_termini")
        table, rows, agg = isoform_termini_features(features, tx, args.sites)

    table.insert(1, "kind", kind)
    # rows of the output: the selected input rows, summed into features by the indicator agg
    out = (agg.T @ m[rows, :]).tocoo() if agg is not None else m[rows, :].tocoo()
    logger.info("%s: %d features (%d competing) in %d groups; %s of %s counts in clustered cells kept",
                kind, len(table), int(table.competing.sum()), table.loc[table.competing, "gene_key"].nunique(),
                f"{out.sum():,.0f}", f"{m.sum():,.0f}")

    write_mtx_real(args.output_prefix + ".site_counts.mtx.gz", out)
    with gzip.open(args.output_prefix + ".barcodes.tsv.gz", "wt") as ofh:
        for cb in barcodes:
            print(f"{cb}\t{cluster_of[cb]}", file=ofh)
    with gzip.open(args.output_prefix + ".sites.tsv.gz", "wt") as ofh:
        for sid in table.site_id:
            print(sid, file=ofh)
    table.to_csv(args.output_prefix + ".features.tsv", sep="\t", index=False)


def own_symbol(feature_id):
    """the gene symbol before '^' in an LRAA id; '' when there is none"""
    return feature_id.split("^", 1)[0] if "^" in feature_id else ""


def read_sparse_matrix(d, barcodes):
    """features x clustered-cells csr matrix (columns in the order of barcodes), and the feature ids"""
    features = [l.rstrip("\n").split("\t")[0] for l in gzip.open(os.path.join(d, "features.tsv.gz"), "rt")]
    cells = [l.rstrip("\n").split("\t")[0] for l in gzip.open(os.path.join(d, "barcodes.tsv.gz"), "rt")]
    with gzip.open(os.path.join(d, "matrix.mtx.gz"), "rt") as fh:
        header = fh.readline()
        if "coordinate" not in header:
            sys.exit(f"{d}/matrix.mtx.gz: not a coordinate MatrixMarket file")
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        nr, nc, nnz = (int(x) for x in line.split())
        trip = pd.read_csv(fh, sep=" ", header=None, names=["r", "c", "v"],
                           dtype={"r": np.int32, "c": np.int32, "v": np.float64}, engine="c")
    if nr != len(features) or nc != len(cells) or len(trip) != nnz:
        sys.exit(f"{d}: matrix dimensions do not match its features / barcodes")
    col_of = {cb: j for j, cb in enumerate(barcodes)}
    fcol = np.array([col_of.get(c, -1) for c in cells])
    c = fcol[trip.c.values - 1]
    keep = c >= 0
    m = sparse.csr_matrix((trip.v.values[keep], (trip.r.values[keep] - 1, c[keep])), shape=(nr, len(barcodes)))
    m.sum_duplicates()
    logger.info("%s: %d features x %d cells (%d clustered cells found of %d)", d, nr, nc,
                int((fcol >= 0).sum()), len(barcodes))
    return m, features


def write_library(m, clusters, filename):
    tot = collections.defaultdict(float)
    col_sums = np.asarray(m.sum(axis=0)).ravel()
    for cl, v in zip(clusters, col_sums):
        tot[cl] += v
    order = sorted(tot, key=lambda c: int(c.replace("Cluster_", "")))
    # one row, columns named by the bare cluster number: the layout of
    # site_read_support_to_sparse_matrix.py's cluster_counts.tsv that annotate_site_usage_events.py reads
    pd.DataFrame([[tot[c] for c in order]], index=["all_features"],
                 columns=[c.replace("Cluster_", "") for c in order]).to_csv(filename, sep="\t")


def write_mtx_real(filename, m):
    order = np.lexsort((m.row, m.col))
    with gzip.open(filename, "wt", compresslevel=4) as ofh:
        print("%%MatrixMarket matrix coordinate real general", file=ofh)
        print(f"{m.shape[0]} {m.shape[1]} {m.nnz}", file=ofh)
        for r, c, v in zip(m.row[order] + 1, m.col[order] + 1, m.data[order]):
            ofh.write(f"{r} {c} {v:.6g}\n")


## ---- splice patterns


def parse_sp_gtf(gtf):
    """transcript_id -> (chrom, strand, main TSS pos, main PolyA pos) from the splice-pattern gtf: the
    TSS_sites / PolyA_sites with most support, else the transcript terminus when TSS / PolyA "True", else None"""
    att_re = re.compile(r'(\S+) "([^"]*)"')
    info = {}
    for line in open(gtf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "transcript":
            continue
        a = dict(att_re.findall(f[8]))
        plus = f[6] == "+"
        lo, hi = int(f[3]), int(f[4])

        def main_site(end, five_prime):
            sites = a.get(f"{end}_sites")
            if sites:
                pos = [int(x) for x in sites.split(",")]
                sup = [float(x) for x in a.get(f"{end}_site_support", "").split(",")] if a.get(f"{end}_site_support") else [0] * len(pos)
                return pos[int(np.argmax(sup))]
            if a.get(end) == "True":
                return (lo if plus else hi) if five_prime else (hi if plus else lo)
            return None
        # a pattern of one isoform keeps that isoform's transcript_id, its pattern id in splice_pattern
        info[a.get("splice_pattern", a["transcript_id"])] = (f[0], f[6], main_site("TSS", True), main_site("PolyA", False))
    return info


def splice_pattern_features(features, tx, sp_gtf):
    sp = parse_sp_gtf(sp_gtf)
    members = tx.groupby("transcript_splice_hash_code").agg(
        isoforms=("transcript_id", lambda s: ",".join(sorted(s))),
        named_symbols=("symbol", lambda s: ",".join(sorted({x for x in s if x}))),
        num_exons=("num_exons", "max"), chrom=("chrom", "first"), strand=("strand", "first"))
    rows, recs = [], []
    n_unspliced = n_unmapped = 0
    for i, fid in enumerate(features):
        if ":iso-" in fid:
            n_unspliced += 1
            continue
        h = fid.split("^", 1)[-1]
        if h not in members.index:
            n_unmapped += 1
            continue
        mem = members.loc[h]
        sym = own_symbol(fid) or mem.named_symbols.split(",")[0]
        chrom, strand, tss, polya = sp.get(fid, (mem.chrom, mem.strand, None, None))
        ambiguous = "," in mem.named_symbols or (bool(mem.named_symbols) and sym != mem.named_symbols)
        rows.append(i)
        recs.append({"site_id": fid, "competing": bool(sym) and not ambiguous and mem.num_exons > 1,
                     "gene_key": f"{sym}|{chrom}|{strand}" if sym else "", "gene_symbol": sym,
                     "chrom": chrom, "strand": strand, "num_exons": int(mem.num_exons),
                     "main_TSS": "" if tss is None else tss, "main_PolyA": "" if polya is None else polya,
                     "in_sp_gtf": fid in sp, "isoform_symbols": mem.named_symbols, "transcript_ids": mem.isoforms})
    table = pd.DataFrame(recs)
    logger.info("splice patterns: %d unspliced features left out, %d not in the transcript map; "
                "%d of %d kept patterns not in the sp gtf; %d marked non-competing (isoforms of several genes)",
                n_unspliced, n_unmapped, int((~table.in_sp_gtf).sum()), len(table), int((~table.competing).sum()))
    return table, np.array(rows), None


## ---- terminal sites within splice patterns


def isoform_termini_features(features, tx, sites_file):
    sites = pd.read_csv(sites_file, sep="\t", keep_default_na=False, low_memory=False,
                        usecols=["site_id", "kind", "transcript_ids", "pas"])
    site_of = {"TSS": {}, "PolyA": {}}
    for sid, k, tids in zip(sites.site_id, sites.kind, sites.transcript_ids):
        for t in tids.split(","):
            if t:
                if t in site_of[k] and site_of[k][t] != sid:
                    logger.warning("%s carries two %s sites (%s, %s); keeping the first", t, k, site_of[k][t], sid)
                    continue
                site_of[k][t] = sid
    txi = tx.set_index("transcript_id")
    # a pattern's label: the gene symbol most of its isoforms carry
    sym_of_hash = tx[tx.symbol != ""].groupby("transcript_splice_hash_code").symbol.agg(
        lambda s: s.value_counts().index[0])
    n_isoforms_of_hash = tx.groupby("transcript_splice_hash_code").size()

    combo_of_row = {}
    n_mono = n_unmapped = n_nosite = 0
    for i, fid in enumerate(features):
        t = fid.split("^", 1)[-1]
        if t not in txi.index:
            n_unmapped += 1
            continue
        r = txi.loc[t]
        if r.num_exons < 2:
            n_mono += 1
            continue
        tss, pa = site_of["TSS"].get(t, "none"), site_of["PolyA"].get(t, "none")
        if tss == "none" and pa == "none":
            n_nosite += 1
            continue
        h = r.transcript_splice_hash_code
        combo_of_row[i] = (h, r.chrom, r.strand, tss, pa, t)
    logger.info("isoforms: %d in combinations; left out %d monoexonic, %d with no annotated site at either end, "
                "%d not in the transcript map", len(combo_of_row), n_mono, n_nosite, n_unmapped)

    # combinations, in order of chrom, pattern, sites; the indicator maps each selected input row to its combination
    sel = sorted(combo_of_row, key=lambda i: (combo_of_row[i][1], combo_of_row[i][0], combo_of_row[i][3], combo_of_row[i][4]))
    combos = {}
    jj = []
    for i in sel:
        h, chrom, strand, tss, pa, t = combo_of_row[i]
        sym = sym_of_hash.get(h, "")
        sp_id = f"{sym}^{h}" if sym else h
        fid = f"{sp_id}|{tss}|{pa}"
        if fid not in combos:
            combos[fid] = {"site_id": fid, "gene_key": f"{sp_id}|{chrom}|{strand}", "gene_symbol": sym,
                           "splice_pattern": sp_id, "chrom": chrom, "strand": strand, "TSS_site": tss,
                           "PolyA_site": pa, "transcript_ids": [], "n_pattern_isoforms": int(n_isoforms_of_hash[h]),
                           "col": len(combos)}
        combos[fid]["transcript_ids"].append(t)
        jj.append(combos[fid]["col"])
    agg = sparse.csr_matrix((np.ones(len(sel)), (np.arange(len(sel)), jj)), shape=(len(sel), len(combos)))

    table = pd.DataFrame(list(combos.values()))
    table["n_isoforms"] = table.transcript_ids.str.len()
    table["transcript_ids"] = table.transcript_ids.str.join(",")
    n_per_group = table.groupby("gene_key").site_id.transform("size")
    # a combination is testable when its pattern has another; the rest stay in the table, not tested
    table["competing"] = n_per_group >= 2
    table = table[["site_id", "competing", "gene_key", "gene_symbol", "splice_pattern", "chrom", "strand",
                   "TSS_site", "PolyA_site", "n_isoforms", "n_pattern_isoforms", "transcript_ids"]]
    logger.info("combinations: %d in %d splice patterns; %d patterns with >= 2 combinations",
                len(table), table.gene_key.nunique(), table.loc[table.competing, "gene_key"].nunique())
    return table, np.array(sel), agg


if __name__ == "__main__":
    main()
