#!/usr/bin/env python3

"""Turn differential-usage pairs (two models, two clusters) into read-track events for
build_site_event_read_tracks.py, so a DTU comparison can be drawn over its reads the way
a site-usage switch is.

Input --pairs (tsv): tag, gene_symbol, gained_id, lost_id, cluster_A, cluster_B, where the
gained model's share rises from cluster_A to cluster_B. An id is either an LRAA isoform
model (SYMBOL^t:...:iso-N, as in --gtf) or a splice pattern (SYMBOL^<hash>, the
new_transcript_splice_hash_code of --splice_pattern_map). A splice pattern is drawn as
one of its member isoform models: among the members with >= --min_uniq_FSM unique
full-splice-match reads over all clusters (all members if none has that many), the one
with the most reads assigned in the cluster favouring it (cluster_B for the gained
pattern, cluster_A for the lost one); ties: reads over all clusters.

Each event gets a site kind and a site per side, from the site table (--sites, as used
by the site-usage analysis): the site whose counting window holds the drawn model's
5' end (TSS) or 3' end (PolyA). The kind is the end at which the two models sit at
different sites -- TSS first, then PolyA -- or, when they share both, the TSS (a splicing
switch at a shared terminus; build_site_event_read_tracks.py then requires a compatible
read to fit exactly one of the two models). A pair whose models reach no site at either
end is skipped and listed in --skipped.

Output --output: tag, gene_symbol, kind, gained_site, lost_site, cluster_A, cluster_B,
gained_tx, lost_tx (tracking ids, the part after SYMBOL^), plus gained_model / lost_model
(the drawn model ids) and the pair's ids as given.
"""

import argparse
import collections
import csv
import io
import logging
import re
import tarfile

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s", level=logging.INFO)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--gtf", required=True, help="LRAA gtf with SYMBOL^ transcript ids")
    ap.add_argument("--sites", required=True, help="site table (prep_site_table.py)")
    ap.add_argument("--cluster_quant_tar", required=True, help="tar.gz of per-cluster quant.expr")
    ap.add_argument("--splice_pattern_map", default=None,
                    help="gene_transcript_splicehashcode tsv (needed when ids are splice patterns)")
    ap.add_argument("--min_uniq_FSM", type=float, default=3)
    ap.add_argument("--output", required=True)
    ap.add_argument("--skipped", default=None, help="tsv of pairs that could not be placed (default: <output>.skipped.tsv)")
    args = ap.parse_args()

    pairs = list(csv.DictReader(open(args.pairs), delimiter="\t"))

    # pattern -> member tracking ids
    members = collections.defaultdict(list)
    if args.splice_pattern_map:
        for r in csv.DictReader(open(args.splice_pattern_map), delimiter="\t"):
            members[r["new_transcript_splice_hash_code"]].append(r["transcript_id"])

    # model exons, keyed by tracking id (the part after SYMBOL^)
    want = set()
    for p in pairs:
        for x in (p["gained_id"], p["lost_id"]):
            want.update(members.get(x, []) or [x.split("^", 1)[-1]])
    exons, strand_of, chrom_of = collections.defaultdict(list), {}, {}
    full_id = {}
    for line in open(args.gtf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m:
            continue
        tid = m.group(1)
        trk = tid.split("^", 1)[-1]
        if trk not in want:
            continue
        exons[trk].append((int(f[3]), int(f[4])))
        strand_of[trk], chrom_of[trk] = f[6], f[0]
        full_id[trk] = tid

    # reads per model and cluster, unique FSM reads per model
    reads = collections.defaultdict(float)
    fsm = collections.defaultdict(float)
    with tarfile.open(args.cluster_quant_tar) as tar:
        for mem in tar:
            if not mem.name.endswith("quant.expr"):
                continue
            cm = re.search(r"\.(\d+)\.LRAA\.", mem.name.split("/")[-1])
            if not cm:
                continue
            cl = "Cluster_" + cm.group(1)
            fh = io.TextIOWrapper(tar.extractfile(mem))
            hdr = None
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if hdr is None:
                    hdr = {k: i for i, k in enumerate(f)}
                    continue
                t = f[hdr["transcript_id"]]
                if t not in want:
                    continue
                reads[(t, cl)] += float(f[hdr["all_reads"]])
                reads[(t, "all")] += float(f[hdr["all_reads"]])
                fsm[t] += float(f[hdr["uniq_FSM_reads"]])

    # site windows
    by_strand = collections.defaultdict(list)
    for s in csv.DictReader(open(args.sites), delimiter="\t"):
        lo, hi = int(s["span_lo"]) - int(s["window"]), int(s["span_hi"]) + int(s["window"])
        by_strand[(s["kind"], s["chrom"], s["strand"])].append((int(s["pos"]), lo, hi, s["site_id"]))

    def site_at(kind, chrom, strand, pos):
        hits = [x for x in by_strand.get((kind, chrom, strand), []) if x[1] <= pos <= x[2]]
        return min(hits, key=lambda x: abs(x[0] - pos))[3] if hits else None

    def model_end(t, kind):
        xs = exons[t]
        plus = strand_of[t] == "+"
        first = min(s for s, _ in xs) if plus else max(e for _, e in xs)
        last = max(e for _, e in xs) if plus else min(s for s, _ in xs)
        return first if kind == "TSS" else last

    def pick(x, favouring):
        cands = [t for t in (members.get(x, []) or [x.split("^", 1)[-1]]) if t in exons]
        if not cands:
            return None
        ok = [t for t in cands if fsm[t] >= args.min_uniq_FSM] or cands
        return max(ok, key=lambda t: (reads[(t, favouring)], reads[(t, "all")]))

    out, skipped = [], []
    for p in pairs:
        cl = lambda x: "Cluster_" + x.replace("Cluster_", "")
        g = pick(p["gained_id"], cl(p["cluster_B"]))
        l = pick(p["lost_id"], cl(p["cluster_A"]))
        if g is None or l is None:
            skipped.append({**p, "reason": "model_not_in_gtf"})
            continue
        chrom, strand = chrom_of[g], strand_of[g]
        s = {k: (site_at(k, chrom, strand, model_end(g, k)), site_at(k, chrom, strand, model_end(l, k)))
             for k in ("TSS", "PolyA")}
        kind = next((k for k in ("TSS", "PolyA") if all(s[k]) and s[k][0] != s[k][1]), None) \
            or next((k for k in ("TSS", "PolyA") if all(s[k])), None)
        if kind is None:
            skipped.append({**p, "reason": "models_end_at_no_site"})
            continue
        out.append({"tag": p["tag"], "gene_symbol": p["gene_symbol"], "kind": kind,
                    "gained_site": s[kind][0], "lost_site": s[kind][1],
                    "cluster_A": p["cluster_A"], "cluster_B": p["cluster_B"],
                    "gained_tx": g, "lost_tx": l, "gained_model": full_id[g], "lost_model": full_id[l],
                    "gained_id": p["gained_id"], "lost_id": p["lost_id"]})

    cols = ["tag", "gene_symbol", "kind", "gained_site", "lost_site", "cluster_A", "cluster_B",
            "gained_tx", "lost_tx", "gained_model", "lost_model", "gained_id", "lost_id"]
    with open(args.output, "w") as ofh:
        w = csv.DictWriter(ofh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(out)
    sk = args.skipped or args.output + ".skipped.tsv"
    with open(sk, "w") as ofh:
        if skipped:
            w = csv.DictWriter(ofh, fieldnames=list(skipped[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(skipped)
    logger.info("%d events, %d skipped (%s)", len(out), len(skipped), sk)


if __name__ == "__main__":
    main()
