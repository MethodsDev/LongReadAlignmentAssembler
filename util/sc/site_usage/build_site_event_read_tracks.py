#!/usr/bin/env python3

"""Read-track data for a set of site-switch events (annotate_site_usage_events.py rows),
for drawing each event's two isoforms over a sample of their reads plus each cluster's
read-end density (site_usage_funcs.R plot_isoform_read_tracks).

For each event:
  - the isoform drawn for each site is the one, among the isoforms carrying the site
    (site table transcript_ids) with >= --min_uniq_FSM unique full-splice-match (FSM)
    reads (all of them if none has that many), with the most reads assigned in the
    cluster favouring that site: cluster_B for the gained site, cluster_A for the lost
    one -- so the pair drawn is the pair carrying the switch (ties: reads assigned over
    all clusters). Ranking by unique FSM reads alone favours short fragment models:
    full-length reads are shared among near-identical full-length models and so are
    rarely unique to any one of them. The manifest gives each isoform's reads in the two
    clusters and the gained isoform's share of the pair in each, to show whether the
    pair itself switches;
  - each isoform's reads in a cluster: its unique FSM reads, plus "compatible" reads --
    reads (any assignment) whose 5' (TSS events) or 3' (PolyA events) end lies within
    --site_tolerance of the isoform's site and whose alignment fits the model: introns
    a consecutive run of the model's introns, no block reaching into a model intron or
    past the model's far end;
  - drawn: up to --max_reads per cluster from the event's two clusters, split between the
    two isoforms in proportion to their reads there (FSM + compatible), each isoform's
    share filled from its unique FSM reads first, then from its compatible reads;
  - read-end density: per cluster, the 5' (TSS events) or 3' (PolyA events) ends of all
    reads of the two clusters on the gene's strand within the two isoforms' span;
  - totals: reads per isoform and cluster (FSM + compatible, and unique FSM alone).

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
    parser.add_argument("--min_uniq_FSM", type=float, default=5,
                        help="isoforms with fewer unique FSM reads are drawn only if no isoform at the site has this many")
    parser.add_argument("--site_tolerance", type=int, default=25,
                        help="a compatible read's terminus lies within this many bp of the site (LRAA: half the 50 bp site window)")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    events = list(csv.DictReader(open(args.events), delimiter="\t"))
    sites = {r["site_id"]: r for r in csv.DictReader(open(args.sites), delimiter="\t")}
    fsm, assigned, by_cluster = parse_FSM(args.cluster_quant_tar)
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
            ok = [t for t in tids if fsm.get(t, 0) >= args.min_uniq_FSM] or tids
            favoured = (e["cluster_B"] if side == "gained" else e["cluster_A"]).replace("Cluster_", "")
            in_cl = by_cluster.get(favoured, {})
            e[f"{side}_tx"] = max(ok, key=lambda t: (in_cl.get(t, 0), assigned.get(t, 0), fsm.get(t, 0), t))
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

        site_pos = {txs[0]: int(e["gained_site"].split(":")[2]), txs[1]: int(e["lost_site"].split(":")[2])}
        fsm_names = {t: {(c, m) for c, m in reads_by_tx[t] if c in clusters} for t in txs}
        fsm_set = {m for t in txs for _, m in fsm_names[t]}

        # one pass over the region: read-end density, the FSM reads' alignments, and the
        # compatible reads of each isoform
        ends = collections.Counter()
        aln = {}
        compatible = {t: set() for t in txs}
        for read in bam.fetch(chrom, max(0, lo - 1), hi):
            if read.is_secondary or read.is_supplementary:
                continue
            s = "-" if read.is_reverse else "+"
            if read.has_tag("ts") and read.get_tag("ts") == "-":
                s = "+" if s == "-" else "-"
            if s != strand or not read.has_tag("CB"):
                continue
            cl = cluster_of.get(read.get_tag("CB"))
            if cl not in clusters:
                continue
            five = (read.reference_start + 1) if s == "+" else read.reference_end
            three = read.reference_end if s == "+" else (read.reference_start + 1)
            pos = five if e["kind"] == "TSS" else three
            if lo <= pos <= hi:
                ends[(cl, pos)] += 1
            name = read.query_name
            if name in aln:
                continue
            blks = None
            if name in fsm_set:
                blks = blocks(read)
            else:
                for t in txs:
                    if abs(pos - site_pos[t]) <= args.site_tolerance:
                        b = blocks(read)
                        if fits_model(b, exons[t]["exons"], args.site_tolerance):
                            compatible[t].add((cl, name))
                            blks = b
                        break
            if blks is not None:
                aln[name] = (read.reference_start + 1, read.reference_end, s, blks)

        totals = collections.Counter()
        chosen = {}
        for cl in clusters:
            fsm_cl = {t: sorted(m for c, m in fsm_names[t] if c == cl and m in aln) for t in txs}
            comp_cl = {t: sorted(m for c, m in compatible[t] if c == cl) for t in txs}
            n = {t: len(fsm_cl[t]) + len(comp_cl[t]) for t in txs}
            for t in txs:
                totals[(cl, t)] = (n[t], len({m for c, m in fsm_names[t] if c == cl}))
            k = min(args.max_reads, sum(n.values()))
            k0 = round(k * n[txs[0]] / sum(n.values())) if k else 0
            for t, kt in ((txs[0], k0), (txs[1], k - k0)):
                take = rng.sample(fsm_cl[t], min(kt, len(fsm_cl[t])))
                take += rng.sample(comp_cl[t], min(kt - len(take), len(comp_cl[t])))
                for m in take:
                    chosen[m] = (t, cl, "uniq_FSM" if m in fsm_set else "compatible")

        found = set()
        with open(os.path.join(args.outdir, f"{tag}.reads.tsv"), "wt") as ofh:
            w = csv.writer(ofh, delimiter="\t", lineterminator="\n")
            w.writerow(["transcript_id", "cluster", "read_name", "read_start", "read_end", "strand",
                        "block_start", "block_end", "read_class"])
            for m, (t, cl, rc) in sorted(chosen.items(), key=lambda x: x[0]):
                if m not in aln:
                    continue
                found.add(m)
                rs, re_, s, blks = aln[m]
                for b0, b1 in blks:
                    w.writerow([t, cl, m, rs, re_, s, b0, b1, rc])
        if len(found) < len(chosen):
            logger.warning("%s: %d sampled reads not found in the region", tag, len(chosen) - len(found))
        with open(os.path.join(args.outdir, f"{tag}.ends.tsv"), "wt") as ofh:
            print("cluster\tpos\treads", file=ofh)
            for (cl, pos), n in sorted(ends.items()):
                print(f"{cl}\t{pos}\t{n}", file=ofh)
        with open(os.path.join(args.outdir, f"{tag}.totals.tsv"), "wt") as ofh:
            print("cluster\ttranscript_id\tn\tn_uniq_FSM", file=ofh)
            for (cl, t), (n, n_fsm) in sorted(totals.items()):
                print(f"{cl}\t{t}\t{n}\t{n_fsm}", file=ofh)

        manifest.append({"tag": tag, "gene_symbol": e["gene_symbol"], "kind": e["kind"],
                         "cluster_A": clusters[0], "cluster_B": clusters[1],
                         "gained_site": e["gained_site"], "lost_site": e["lost_site"],
                         "gained_pos": e["gained_site"].split(":")[2], "lost_pos": e["lost_site"].split(":")[2],
                         "gained_tx": txs[0], "lost_tx": txs[1],
                         "gained_gtf_id": gtf_id.get(txs[0], txs[0]), "lost_gtf_id": gtf_id.get(txs[1], txs[1]),
                         "gained_uniq_FSM": fsm.get(txs[0], 0), "lost_uniq_FSM": fsm.get(txs[1], 0),
                         "gained_reads_assigned": round(assigned.get(txs[0], 0), 1),
                         "lost_reads_assigned": round(assigned.get(txs[1], 0), 1),
                         **pair_usage(by_cluster, clusters, txs),
                         "chrom": chrom, "strand": strand, "region_start": lo, "region_end": hi,
                         "n_reads_drawn": len(found),
                         "n_compatible_drawn": sum(1 for m in found if chosen[m][2] == "compatible")})
        logger.info("%s: %s / %s, %d reads drawn", tag, txs[0], txs[1], len(found))

    with open(os.path.join(args.outdir, "manifest.tsv"), "wt") as ofh:
        w = csv.DictWriter(ofh, fieldnames=list(manifest[0].keys()) if manifest else ["tag"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(manifest)


def pair_usage(by_cluster, clusters, txs):
    """reads assigned to each isoform of the pair in each cluster, and the gained isoform's
    share of the pair there"""
    out = {}
    for side, t in zip(("gained", "lost"), txs):
        for ab, cl in zip("AB", clusters):
            out[f"{side}_reads_{ab}"] = round(by_cluster.get(cl, {}).get(t, 0), 1)
    for ab in "AB":
        tot = out[f"gained_reads_{ab}"] + out[f"lost_reads_{ab}"]
        out[f"gained_pair_frac_{ab}"] = round(out[f"gained_reads_{ab}"] / tot, 3) if tot else ""
    return out


def parse_FSM(tar):
    """unique FSM reads and all reads assigned per transcript, summed over the clusters,
    and all reads assigned per cluster (cluster from the file name: <prefix>.<cluster>.LRAA...)"""
    fsm, assigned = collections.Counter(), collections.Counter()
    by_cluster = collections.defaultdict(collections.Counter)
    with tarfile.open(tar) as tf:
        for mem in tf.getmembers():
            if mem.name.endswith("quant.expr"):
                m = re.search(r"\.([^./]+)\.LRAA[^/]*$", mem.name)
                cl = m.group(1) if m else None
                rows = csv.DictReader((l for l in io.TextIOWrapper(tf.extractfile(mem)) if not l.startswith("#")),
                                      delimiter="\t")
                for r in rows:
                    fsm[r["transcript_id"]] += float(r["uniq_FSM_reads"])
                    assigned[r["transcript_id"]] += float(r["all_reads"])
                    if cl is not None:
                        by_cluster[cl][r["transcript_id"]] += float(r["all_reads"])
    return fsm, assigned, by_cluster


def fits_model(read_blocks, model_exons, end_slack, intron_slack=3, edge_slack=10):
    """read structure compatible with the model: the read's introns are a consecutive run
    of the model's introns (+- intron_slack); its outer blocks stay inside the model exons
    they map to (+- edge_slack at internal exon edges, +- end_slack at the model's ends)"""
    ex = sorted(model_exons)
    m_introns = [(ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1)]
    r_introns = [(read_blocks[i][1] + 1, read_blocks[i + 1][0] - 1) for i in range(len(read_blocks) - 1)]
    if r_introns:
        first = [k for k, mi in enumerate(m_introns)
                 if abs(mi[0] - r_introns[0][0]) <= intron_slack and abs(mi[1] - r_introns[0][1]) <= intron_slack]
        if not first:
            return False
        k = first[0]
        if k + len(r_introns) > len(m_introns):
            return False
        for j, ri in enumerate(r_introns):
            mi = m_introns[k + j]
            if abs(mi[0] - ri[0]) > intron_slack or abs(mi[1] - ri[1]) > intron_slack:
                return False
        e_first, e_last = k, k + len(r_introns)
    else:
        hits = [k for k, (a, b) in enumerate(ex) if read_blocks[0][0] <= b and read_blocks[0][1] >= a]
        if len(hits) != 1:
            return False
        e_first = e_last = hits[0]
    lo_slack = end_slack if e_first == 0 else edge_slack
    hi_slack = end_slack if e_last == len(ex) - 1 else edge_slack
    return read_blocks[0][0] >= ex[e_first][0] - lo_slack and read_blocks[-1][1] <= ex[e_last][1] + hi_slack


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
