#!/usr/bin/env python3

"""Check site-level alternative-TSS events against FANTOM5 CAGE from sorted cell types.

Each event (site_usage_dexseq.R + annotate_site_usage_events.py: the gained TSS rises in
usage from cluster_A to cluster_B, the lost TSS falls) is matched to FANTOM5 CAGE peaks:
each of its two sites to the nearest peak within --tolerance nt on the same strand (peak
interval, not summit). The two clusters are mapped to FANTOM5 cell groups
(--cluster_groups), and each group's CAGE expression of a peak is the mean TPM over the
group's samples (--sample_groups: regexes on the TPM table's decoded sample names).

The CAGE share of the gained site, gained_TPM / (gained_TPM + lost_TPM), is computed in
cluster_A's and cluster_B's cell groups; the event agrees when that share rises from A to
B, as ours does. An event is not assessable when a site has no CAGE peak, a cluster has no
cell group or both clusters map to the same group, or a group has no CAGE signal at either
site.

Outputs:
  <output_prefix>.events.tsv    the events with cage, cage_share_A, cage_share_B,
                                cage_share_shift, gained_cage_peak, lost_cage_peak
  <output_prefix>.summary.tsv   per event set (high-confidence / all), unique
                                gene x site pair x cluster pair: counts per outcome, the
                                agreeing share of the assessable ones, and agree /
                                disagree counts at |CAGE shift| >= --min_shift
"""

import argparse
import bisect
import collections
import gzip
import logging
import re
import urllib.parse

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s", level=logging.INFO)

KEY = ["gene_key", "gained_site", "lost_site", "cluster_A", "cluster_B"]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--events", required=True, help="<prefix>.dexseq.TSS.events.tsv (annotate_site_usage_events.py)")
    p.add_argument("--cage_peaks", required=True, help="FANTOM5 hg38_fair+new_CAGE_peaks_phase1and2.bed.gz")
    p.add_argument("--cage_tpm", required=True, help="FANTOM5 hg38_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz")
    p.add_argument("--cluster_groups", required=True,
                   help="tsv: cluster <tab> group (cluster number or Cluster_<n>; '#' lines and a header skipped)")
    p.add_argument("--sample_groups", required=True,
                   help="tsv: group <tab> sample_regex (several lines per group allowed; '#' lines and a header "
                        "skipped); samples whose name mentions 'expanded' are excluded")
    p.add_argument("--tolerance", type=int, default=25, help="max distance (nt) from a site to a CAGE peak interval")
    p.add_argument("--min_shift", type=float, default=0.2, help="CAGE share shift reported separately")
    p.add_argument("--output_prefix", required=True)
    args = p.parse_args()

    ev = pd.read_csv(args.events, sep="\t", low_memory=False)
    cl_group = read_cluster_groups(args.cluster_groups)
    peaks = load_peaks(args.cage_peaks)
    sites = set(ev.gained_site) | set(ev.lost_site)
    site_peak = {}
    for s in sites:
        pk = nearest_peak(peaks, s, args.tolerance)
        if pk is not None:
            site_peak[s] = pk
    logger.info("%d of %d event sites within %d nt of a CAGE peak", len(site_peak), len(sites), args.tolerance)

    tpm = load_group_tpm(args.cage_tpm, read_sample_groups(args.sample_groups), set(site_peak.values()))

    def assess(r):
        gA, gB = cl_group.get(norm_cluster(r.cluster_A)), cl_group.get(norm_cluster(r.cluster_B))
        if r.gained_site not in site_peak or r.lost_site not in site_peak:
            return "no CAGE peak at a site", np.nan, np.nan
        if gA is None or gB is None or gA == gB:
            return "cell types not separable in FANTOM5", np.nan, np.nan
        g, lo = tpm[site_peak[r.gained_site]], tpm[site_peak[r.lost_site]]
        shA = g[gA] / (g[gA] + lo[gA]) if g[gA] + lo[gA] > 0 else np.nan
        shB = g[gB] / (g[gB] + lo[gB]) if g[gB] + lo[gB] > 0 else np.nan
        if np.isnan(shA) or np.isnan(shB):
            return "no CAGE signal in a cell type", shA, shB
        return ("agrees" if shB > shA else "disagrees"), shA, shB

    res = [assess(r) for r in ev.itertuples()]
    ev = ev.assign(cage=[x[0] for x in res], cage_share_A=[x[1] for x in res], cage_share_B=[x[2] for x in res],
                   gained_cage_peak=ev.gained_site.map(site_peak), lost_cage_peak=ev.lost_site.map(site_peak))
    ev["cage_share_shift"] = ev.cage_share_B - ev.cage_share_A
    ev.to_csv(f"{args.output_prefix}.events.tsv", sep="\t", index=False)

    rows = []
    sets = [("all", ev)]
    if "high_confidence" in ev.columns:
        sets.insert(0, ("high_confidence", ev[ev.high_confidence.astype(str) == "True"]))
    for name, d in sets:
        d = d.drop_duplicates(KEY)
        t = d[d.cage.isin(["agrees", "disagrees"])]
        big = t[t.cage_share_shift.abs() >= args.min_shift]
        row = {"event_set": name, "events": len(d)}
        row.update({f"n_{k.replace(' ', '_')}": int(v) for k, v in d.cage.value_counts().items()})
        row.update({"assessable": len(t), "frac_agree": round((t.cage == "agrees").mean(), 3) if len(t) else np.nan,
                    f"agree_shift_ge_{args.min_shift}": int((big.cage == "agrees").sum()),
                    f"disagree_shift_ge_{args.min_shift}": int((big.cage == "disagrees").sum())})
        rows.append(row)
        logger.info("%s: %d events, %d assessable, %.0f%% agree; |shift| >= %g: %d agree, %d disagree", name, len(d),
                    len(t), 100 * row["frac_agree"] if len(t) else float("nan"), args.min_shift,
                    row[f"agree_shift_ge_{args.min_shift}"], row[f"disagree_shift_ge_{args.min_shift}"])
    pd.DataFrame(rows).to_csv(f"{args.output_prefix}.summary.tsv", sep="\t", index=False)


def norm_cluster(c):
    c = str(c).strip()
    return c if c.startswith("Cluster_") else f"Cluster_{c}"


def read_table(path):
    out = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            out.append((f[0].strip(), f[1].strip()))
    return out[1:] if out and out[0][0] in ("cluster", "group") else out


def read_cluster_groups(path):
    return {norm_cluster(c): g for c, g in read_table(path) if g}


def read_sample_groups(path):
    groups = collections.defaultdict(list)
    for g, rx in read_table(path):
        groups[g].append(re.compile(rx))
    return groups


def load_peaks(path):
    """{(chrom, strand): ([starts], [(start, end, peak_id)])} from the FANTOM5 peak bed
    (1-based inclusive intervals; peak id = the name's last ';' field)"""
    idx = collections.defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            idx[(f[0], f[5])].append((int(f[1]) + 1, int(f[2]), f[3].split(";")[-1]))
    out = {}
    for k, v in idx.items():
        v.sort()
        out[k] = ([x[0] for x in v], v)
    return out


def nearest_peak(peaks, site_id, tol, max_peak_len=5000):
    """site_id 'TSS:chrom:pos:strand' -> the nearest peak id within tol nt, or None"""
    _, chrom, pos, strand = site_id.split(":")
    pos = int(pos)
    if (chrom, strand) not in peaks:
        return None
    starts, recs = peaks[(chrom, strand)]
    best = None
    for s, e, pid in recs[bisect.bisect_left(starts, pos - tol - max_peak_len):]:
        if s > pos + tol:
            break
        d = 0 if s <= pos <= e else min(abs(pos - s), abs(pos - e))
        if d <= tol and (best is None or d < best[0]):
            best = (d, pid)
    return best[1] if best else None


def load_group_tpm(path, sample_groups, want):
    """{peak_id: {group: mean TPM over the group's samples}} for the wanted peaks"""
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("##"):
                header = line.rstrip("\n").split("\t")
                break
        names = [urllib.parse.unquote(h)[4:] for h in header]   # 'tpm.' prefix
        gcols = {g: [i for i, n in enumerate(names) if any(rx.search(n) for rx in rxs) and "expanded" not in n]
                 for g, rxs in sample_groups.items()}
        for g, cs in gcols.items():
            logger.info("CAGE group %s: %d samples", g, len(cs))
            if not cs:
                raise SystemExit(f"no FANTOM5 sample matches group {g}")
        tpm = {}
        for line in fh:
            i = line.find("\t")
            pid = line[:i].split(";")[-1]
            if pid in want:
                f = line.rstrip("\n").split("\t")
                tpm[pid] = {g: float(np.mean([float(f[c]) for c in cs])) for g, cs in gcols.items()}
    return tpm


if __name__ == "__main__":
    main()
