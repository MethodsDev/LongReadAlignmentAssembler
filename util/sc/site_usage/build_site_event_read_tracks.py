#!/usr/bin/env python3

"""Read-track data for a set of site-switch events (annotate_site_usage_events.py rows),
for drawing each event's two isoforms over a sample of their reads plus each cluster's
read-end density (site_usage_funcs.R plot_isoform_read_tracks).

For each event:
  - the isoform drawn for each site is the one, among the isoforms carrying the site
    (site table transcript_ids), with the most unique full-splice-match reads summed
    over the cluster quantifications;
  - reads: up to --max_reads per cluster from the event's two clusters, sampled from the
    unique full-splice-match reads of the two isoforms together, so each isoform's share
    of a cluster's reads is kept (as extract_isoform_read_tracks.py --proportional);
  - read-end density: per cluster, the 5' (TSS events) or 3' (PolyA events) ends of all
    reads of the two clusters on the gene's strand within the two isoforms' span;
  - totals: unique full-splice-match reads per isoform and cluster.

The tracking file is read once for all events.

Outputs in --outdir: <tag>.reads.tsv, <tag>.ends.tsv, <tag>.totals.tsv per event, and
manifest.tsv (one row per event: tag, gene_symbol, kind, clusters, sites, the two
isoforms' tracking and gtf ids, region, files).
"""

import argparse
import collections
import csv
import gzip
import io
import logging
import os
import random
import re
import tarfile

import pysam

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s", level=logging.INFO)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--events", required=True,
                        help="tsv: tag, gene_symbol, kind (TSS|PolyA), gained_site, lost_site, cluster_A, cluster_B")
    parser.add_argument("--sites", required=True, help="site table from prep_site_table.py")
    parser.add_argument("--gtf", required=True, help="LRAA gtf with SYMBOL^ transcript ids")
    parser.add_argument("--cluster_quant_tar", required=True, help="tar.gz of per-cluster quant.expr (uniq_FSM_reads)")
    parser.add_argument("--tracking", required=True, help="merged per-cluster quant.tracking(.gz)")
    parser.add_argument("--bam", required=True)
    parser.add_argument("--cell_clusters", required=True, help="cell_barcode <tab> cluster (header skipped)")
    parser.add_argument("--max_reads", type=int, default=30)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    events = list(csv.DictReader(open(args.events), delimiter="\t"))
    sites = {r["site_id"]: r for r in csv.DictReader(open(args.sites), delimiter="\t")}
    fsm = parse_FSM(args.cluster_quant_tar)
    cluster_of = {}
    for line in open(args.cell_clusters):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[1].strip().lstrip("-").isdigit():
            cluster_of[f[0]] = f[1].strip()

    # choose the isoform per site
    for e in events:
        for side in ("gained", "lost"):
            tids = [t for t in sites[e[f"{side}_site"]]["transcript_ids"].split(",") if t]
            if not tids:
                raise SystemExit(f"{e['tag']}: site {e[side + '_site']} carries no isoform")
            e[f"{side}_tx"] = max(tids, key=lambda t: (fsm.get(t, 0), t))
        if e["gained_tx"] == e["lost_tx"]:
            logger.warning("%s: the same isoform carries both sites; skipped", e["tag"])
    events = [e for e in events if e["gained_tx"] != e["lost_tx"]]
    want = {e[k] for e in events for k in ("gained_tx", "lost_tx")}

    exons, gtf_id = parse_gtf(args.gtf, want)

    # one pass over the tracking file
    reads_by_tx = collections.defaultdict(list)  # tx -> [(cluster, molecule)] unique FSM
    opener = gzip.open if args.tracking.endswith(".gz") else open
    with opener(args.tracking, "rt") as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("gene_id"):
                continue
            f = line.rstrip("\n").split("\t")
            tx = f[1].split("@")[-1]
            if tx not in want or f[8] != "1" or f[9] != "1":
                continue
            parts = f[5].split("^")
            cl = cluster_of.get(parts[0])
            if cl is not None:
                reads_by_tx[tx].append((cl, parts[-1]))
    logger.info("tracking read: %d isoforms with unique FSM reads", len(reads_by_tx))

    rng = random.Random(args.seed)
    bam = pysam.AlignmentFile(args.bam)
    manifest = []
    for e in events:
        tag = e["tag"]
        clusters = [e["cluster_A"].replace("Cluster_", ""), e["cluster_B"].replace("Cluster_", "")]
        txs = [e["gained_tx"], e["lost_tx"]]
        chrom, strand = exons[txs[0]]["chrom"], exons[txs[0]]["strand"]
        lo = min(min(s for s, _ in exons[t]["exons"]) for t in txs) - 50
        hi = max(max(x for _, x in exons[t]["exons"]) for t in txs) + 50

        totals = collections.Counter()
        chosen = {}
        for cl in clusters:
            pool = sorted({(m, t) for t in txs for c, m in reads_by_tx[t] if c == cl})
            for t in txs:
                totals[(cl, t)] = len({m for c, m in reads_by_tx[t] if c == cl})
            for m, t in rng.sample(pool, min(args.max_reads, len(pool))):
                chosen[m] = (t, cl)

        ends = collections.Counter()
        found = set()
        with open(os.path.join(args.outdir, f"{tag}.reads.tsv"), "wt") as ofh:
            w = csv.writer(ofh, delimiter="\t", lineterminator="\n")
            w.writerow(["transcript_id", "cluster", "read_name", "read_start", "read_end", "strand",
                        "block_start", "block_end"])
            for read in bam.fetch(chrom, max(0, lo - 1), hi):
                if read.is_secondary or read.is_supplementary:
                    continue
                s = "-" if read.is_reverse else "+"
                if read.has_tag("ts") and read.get_tag("ts") == "-":
                    s = "+" if s == "-" else "-"
                if read.has_tag("CB") and s == strand:
                    cl = cluster_of.get(read.get_tag("CB"))
                    if cl in clusters:
                        five = (read.reference_start + 1) if s == "+" else read.reference_end
                        three = read.reference_end if s == "+" else (read.reference_start + 1)
                        pos = five if e["kind"] == "TSS" else three
                        if lo <= pos <= hi:
                            ends[(cl, pos)] += 1
                if read.query_name in chosen and read.query_name not in found:
                    found.add(read.query_name)
                    t, cl = chosen[read.query_name]
                    for b0, b1 in blocks(read):
                        w.writerow([t, cl, read.query_name, read.reference_start + 1, read.reference_end, s, b0, b1])
        if len(found) < len(chosen):
            logger.warning("%s: %d sampled reads not found in the region", tag, len(chosen) - len(found))
        with open(os.path.join(args.outdir, f"{tag}.ends.tsv"), "wt") as ofh:
            print("cluster\tpos\treads", file=ofh)
            for (cl, pos), n in sorted(ends.items()):
                print(f"{cl}\t{pos}\t{n}", file=ofh)
        with open(os.path.join(args.outdir, f"{tag}.totals.tsv"), "wt") as ofh:
            print("cluster\ttranscript_id\tn", file=ofh)
            for (cl, t), n in sorted(totals.items()):
                print(f"{cl}\t{t}\t{n}", file=ofh)

        manifest.append({"tag": tag, "gene_symbol": e["gene_symbol"], "kind": e["kind"],
                         "cluster_A": clusters[0], "cluster_B": clusters[1],
                         "gained_site": e["gained_site"], "lost_site": e["lost_site"],
                         "gained_pos": e["gained_site"].split(":")[2], "lost_pos": e["lost_site"].split(":")[2],
                         "gained_tx": txs[0], "lost_tx": txs[1],
                         "gained_gtf_id": gtf_id.get(txs[0], txs[0]), "lost_gtf_id": gtf_id.get(txs[1], txs[1]),
                         "gained_uniq_FSM": fsm.get(txs[0], 0), "lost_uniq_FSM": fsm.get(txs[1], 0),
                         "chrom": chrom, "strand": strand, "region_start": lo, "region_end": hi,
                         "n_reads_drawn": len(found)})
        logger.info("%s: %s / %s, %d reads drawn", tag, txs[0], txs[1], len(found))

    with open(os.path.join(args.outdir, "manifest.tsv"), "wt") as ofh:
        w = csv.DictWriter(ofh, fieldnames=list(manifest[0].keys()) if manifest else ["tag"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(manifest)


def parse_FSM(tar):
    fsm = collections.Counter()
    with tarfile.open(tar) as tf:
        for mem in tf.getmembers():
            if mem.name.endswith("quant.expr"):
                rows = csv.DictReader((l for l in io.TextIOWrapper(tf.extractfile(mem)) if not l.startswith("#")),
                                      delimiter="\t")
                for r in rows:
                    fsm[r["transcript_id"]] += float(r["uniq_FSM_reads"])
    return fsm


def parse_gtf(gtf, want):
    """exons and gtf transcript id (SYMBOL^id) of the wanted (bare) transcript ids"""
    exons, gtf_id = {}, {}
    tid_re = re.compile(r'transcript_id "([^"]+)"')
    for line in open(gtf):
        f = line.split("\t", 9)
        if len(f) < 9 or f[2] != "exon":
            continue
        full = tid_re.search(f[8]).group(1)
        bare = full.split("^", 1)[-1]
        if bare not in want:
            continue
        gtf_id[bare] = full
        d = exons.setdefault(bare, {"chrom": f[0], "strand": f[6], "exons": []})
        d["exons"].append((int(f[3]), int(f[4])))
    return exons, gtf_id


def blocks(read):
    """aligned blocks, 1-based inclusive, split only at N"""
    out, pos, cur = [], read.reference_start, None
    for op, n in read.cigartuples:
        if op in (0, 2, 7, 8):
            if cur is None:
                cur = pos
            pos += n
        elif op == 3:
            if cur is not None:
                out.append((cur + 1, pos))
                cur = None
            pos += n
    if cur is not None:
        out.append((cur + 1, pos))
    return out


if __name__ == "__main__":
    main()
