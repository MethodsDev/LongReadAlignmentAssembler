#!/usr/bin/env python3

"""Classify read PolyA (3') ends by the evidence at their position: what share of reads end
at a known polyadenylation site, at a novel LRAA site, look internally primed, or none.

Input is the PolyA end histogram of site_read_support_to_sparse_matrix.py
(--polyA_end_histogram: chrom, strand, pos, residual_soft_clip 0/1/2/3 (= > 2), reads),
so the reads are exactly those the site-usage counts are drawn from (LRAA's read filters,
reads with a cell barcode), and each end is the 3' end LRAA takes, after it strips the
polyA tail from the soft clip.

Each end position is annotated with, within --tolerance nt (LRAA's end-to-site tolerance,
int(max_dist_between_alt_polyA_sites / 2)) on the same strand:
  ref_end      a reference transcript 3' end (--ref_gtf)
  polyasite    a PolyASite 2.0 cluster (any), and polyasite_frac_ok: one seen in
               >= --min_polyasite_frac of the atlas samples
  lraa_site    an LRAA PolyA site (--lraa_polyA_bed)
and with the genomic evidence at the end itself:
  pas          a canonical PAS hexamer (AATAAA / ATTAAA) 10-40 nt upstream
               (Util_funcs.find_polyA_signal)
  downstream_A A's (T's on '-') in the 20 genomic bases past the end

Categories, first match wins:
  known              ref_end or polyasite
  LRAA_site          lraa_site, not known, and not IP-like
  LRAA_site_IP_like  lraa_site, not known, IP-like (a site our evidence filter flags)
  IP_like            no site of any kind, IP-like
  other              none of the above
where IP-like = no PAS and > --max_downstream_A A's downstream: the complement of the
site-usage PolyA evidence filter (annotate_site_usage_events.py), less its PolyASite
term, which `known` already takes.

Outputs:
  <prefix>.evidence.tsv    reads by residual_soft_clip x every evidence combination
                           (ref_end, polyasite, polyasite_frac_ok, lraa_site, pas,
                           downstream_A): any other rule can be applied downstream
  <prefix>.categories.tsv  reads by residual_soft_clip x category
"""

import argparse
import collections
import gzip
import logging
import os
import sys

import numpy as np
import pandas as pd
import pysam

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../../pylib"))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

import LRAA_Globals  # noqa: E402
import Util_funcs  # noqa: E402
from annotate_site_usage_events import load_polyasite  # noqa: E402
from site_read_support_to_sparse_matrix import read_site_bed  # noqa: E402

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(asctime)s %(levelname)s %(message)s", level=logging.INFO)

CATEGORIES = ("known", "LRAA_site", "LRAA_site_IP_like", "IP_like", "other")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--histogram", required=True, help="site_read_support_to_sparse_matrix.py --polyA_end_histogram")
    parser.add_argument("--ref_gtf", required=True, help="reference annotation: transcript 3' ends are known sites")
    parser.add_argument("--polyasite_atlas", required=True, help="PolyASite 2.0 clusters bed(.gz)")
    parser.add_argument("--lraa_polyA_bed", required=True, help="LRAA PolyA site bed (as counted)")
    parser.add_argument("--genome", required=True, help="genome fasta (indexed)")
    parser.add_argument("--tolerance", type=int,
                        default=int(LRAA_Globals.config["max_dist_between_alt_polyA_sites"] / 2))
    parser.add_argument("--min_polyasite_frac", type=float, default=0.10)
    parser.add_argument("--max_downstream_A", type=int, default=7)
    parser.add_argument("--output_prefix", required=True)
    args = parser.parse_args()

    ends = pd.read_csv(args.histogram, sep="\t", dtype={"chrom": str, "strand": str})
    ends = ends.groupby(["chrom", "strand", "pos", "residual_soft_clip"], as_index=False)["reads"].sum()
    logger.info("%s reads at %s end positions", f"{ends.reads.sum():,}", f"{len(ends):,}")

    ref = point_index(gtf_transcript_ends(args.ref_gtf))
    lraa = point_index([(s["chrom"], s["strand"], s["pos"]) for s in read_site_bed(args.lraa_polyA_bed, "PolyA")])
    atlas = load_polyasite(args.polyasite_atlas)
    pas_all = interval_index(atlas, 0.0)
    pas_ok = interval_index(atlas, args.min_polyasite_frac)

    ev = annotate(ends, ref, lraa, pas_all, pas_ok, args.genome, args.tolerance)
    ev["category"] = categorize(ev, args.max_downstream_A)

    keys = ["residual_soft_clip", "ref_end", "polyasite", "polyasite_frac_ok", "lraa_site", "pas", "downstream_A"]
    ev.groupby(keys + ["category"], as_index=False)["reads"].sum().to_csv(
        f"{args.output_prefix}.evidence.tsv", sep="\t", index=False)
    cat = ev.groupby(["residual_soft_clip", "category"], as_index=False)["reads"].sum()
    cat.to_csv(f"{args.output_prefix}.categories.tsv", sep="\t", index=False)
    for clip, d in cat.groupby("residual_soft_clip"):
        logger.info("clip %s: %s", clip, ", ".join(f"{r.category} {r.reads:,}" for r in d.itertuples()))


def gtf_transcript_ends(gtf):
    out = []
    with (gzip.open(gtf, "rt") if gtf.endswith(".gz") else open(gtf)) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] != "transcript":
                continue
            out.append((f[0], f[6], int(f[4]) if f[6] == "+" else int(f[3])))
    return out


def point_index(points):
    """(chrom, strand) -> sorted positions"""
    idx = collections.defaultdict(list)
    for c, s, p in points:
        idx[(c, s)].append(p)
    return {k: np.unique(np.array(v, dtype=np.int64)) for k, v in idx.items()}


def interval_index(atlas, min_frac):
    """(chrom, strand) -> (sorted starts, running max of ends) for clusters seen in >= min_frac
    of samples: a position p lies within tol of one iff, over clusters starting at or before
    p + tol, the furthest end reaches p - tol"""
    out = {}
    for k, (_, v) in atlas.items():
        iv = [(s, e) for s, e, frac in v if frac >= min_frac]
        if iv:
            st = np.array([s for s, _ in iv], dtype=np.int64)
            en = np.maximum.accumulate(np.array([e for _, e in iv], dtype=np.int64))
            out[k] = (st, en)
    return out


def near_point(index, chrom, strand, pos, tol):
    a = index.get((chrom, strand))
    if a is None or not len(a):
        return np.zeros(len(pos), dtype=bool)
    i = np.searchsorted(a, pos)
    lo = np.abs(a[np.clip(i - 1, 0, len(a) - 1)] - pos)
    hi = np.abs(a[np.clip(i, 0, len(a) - 1)] - pos)
    return np.minimum(lo, hi) <= tol


def near_interval(index, chrom, strand, pos, tol):
    if (chrom, strand) not in index:
        return np.zeros(len(pos), dtype=bool)
    st, en = index[(chrom, strand)]
    i = np.searchsorted(st, pos + tol, side="right") - 1
    ok = i >= 0
    out = np.zeros(len(pos), dtype=bool)
    out[ok] = en[i[ok]] >= pos[ok] - tol
    return out


def annotate(ends, ref, lraa, pas_all, pas_ok, genome_fa, tol):
    parts = []
    with pysam.FastaFile(genome_fa) as fa:
        for (chrom, strand), d in ends.groupby(["chrom", "strand"], sort=False):
            d = d.copy()
            pos = d.pos.to_numpy(dtype=np.int64)
            d["ref_end"] = near_point(ref, chrom, strand, pos, tol)
            d["lraa_site"] = near_point(lraa, chrom, strand, pos, tol)
            d["polyasite"] = near_interval(pas_all, chrom, strand, pos, tol)
            d["polyasite_frac_ok"] = near_interval(pas_ok, chrom, strand, pos, tol)
            if chrom in fa.references:
                seq = fa.fetch(chrom).upper()
                d["pas"] = [Util_funcs.find_polyA_signal(seq, int(p), strand)[0] is not None for p in pos]
                d["downstream_A"] = [downstream_A(seq, int(p), strand) for p in pos]
            else:
                d["pas"] = False
                d["downstream_A"] = -1
            parts.append(d)
            logger.info("%s %s: %s positions", chrom, strand, f"{len(d):,}")
    return pd.concat(parts, ignore_index=True)


def downstream_A(seq, pos, strand):
    """A's in the 20 genomic bases past a 1-based 3' end, in transcript sense (T's on '-'),
    as annotate_site_usage_events.downstream_A"""
    if strand == "+":
        return seq[pos:pos + 20].count("A")
    return seq[max(0, pos - 21):pos - 1].count("T")


def categorize(ev, max_downstream_A):
    known = ev.ref_end | ev.polyasite
    ip_like = (~ev.pas) & (ev.downstream_A > max_downstream_A)
    return np.select([known, ev.lraa_site & ~ip_like, ev.lraa_site & ip_like, ip_like],
                     list(CATEGORIES[:4]), default="other")


if __name__ == "__main__":
    main()
