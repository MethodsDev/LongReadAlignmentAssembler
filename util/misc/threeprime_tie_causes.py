#!/usr/bin/env python3
"""Why do so many ambiguous reads get IDENTICAL 3'-agreement weights?

Measurement instrumentation, companion to threeprime_weight_distribution.py. That script
establishes THAT the median multi-compatible read gets a uniform weight vector at N=2 and
N=3. This one asks why, because the answer decides whether the feature is fixable:

  SHARED 3' END   the candidate transcripts terminate at the same 3' coordinate, so there
                  is no 3'-end difference to measure. The transcripts differ somewhere
                  else -- internal splicing, or the 5' end. 3'-end agreement is blind to
                  that BY CONSTRUCTION, and no reweighting or pseudocount recovers a
                  signal that was never present.

  DISTINCT 3' END the candidates end at different coordinates, so the signal exists, but
                  the read still scored equal distances to both. LRAA measures distance
                  as transcript sequence lying 3'-ward of the read's first shared node
                  (Quantify._get_simple_path_dist_to_termini), so this is the truncated
                  read case: the read stops short of the divergence. That IS a weighting
                  problem, and a differently shaped feature could do better.

The split is reported for tied reads and, as a control, for non-tied reads. The contrast
is the point: if tied reads are overwhelmingly shared-3'-end while non-tied reads are
overwhelmingly distinct-3'-end, the tie is structural rather than incidental.

3' coordinates come from the same GTF that was quantified, so they are the coordinates
the transcripts actually have. Distance itself is computed on splice-graph nodes, so
equal 3' coordinates are a proxy for equal distance -- a tight one, since a shared 3'
terminus is what makes the walk lengths equal, but a proxy.
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict

TXID = re.compile(r'transcript_id "([^"]+)"')
SPREAD_EPS = 1e-5


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def load_transcript_ends(gtf):
    """transcript_id -> (three_prime_coord, intron_chain_key, strand)."""
    exons = defaultdict(list)
    strand = {}
    with opener(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = TXID.search(f[8])
            if not m:
                continue
            tid = m.group(1)
            exons[tid].append((int(f[3]), int(f[4])))
            strand[tid] = f[6]

    out = {}
    for tid, ex in exons.items():
        ex.sort()
        s = strand[tid]
        three_prime = ex[-1][1] if s == "+" else ex[0][0]
        # Intron chain: the internal structure that 3'-end agreement cannot see.
        introns = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
        out[tid] = (three_prime, introns, s)
    return out


def load_read_groups(tracking):
    """read_name -> (list of (transcript_id, weight), assigned mass)."""
    members = defaultdict(list)
    mass = defaultdict(float)
    with opener(tracking) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("gene_id\t"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            members[f[5]].append((f[1], float(f[7])))
            mass[f[5]] += float(f[6])
    return members, mass


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tracking")
    ap.add_argument("gtf", help="the GTF that was quantified")
    ap.add_argument("--label", default="")
    ap.add_argument("--polyA_window", type=int, default=50,
                    help="max_dist_between_alt_polyA_sites; 3' ends closer than this are "
                         "the same site as far as the rest of LRAA is concerned")
    args = ap.parse_args()

    ends = load_transcript_ends(args.gtf)
    groups, read_mass = load_read_groups(args.tracking)
    if not groups:
        sys.exit(f"Error, no tracking rows parsed from {args.tracking}")

    # bucket -> [tied_count, nontied_count]
    cats = {"shared_3p_exact": [0, 0],
            f"within_{args.polyA_window}nt": [0, 0],
            "distinct_3p": [0, 0]}
    cat_mass = {k: [0.0, 0.0] for k in cats}
    # Of tied candidate sets that agree on the 3' end, does the ambiguity come from
    # internal splicing or only from the 5' end? Only the first is irreducible at the
    # read level: a read compatible with two isoforms differing by internal splicing
    # cannot span the junction that separates them, or it would match only one of them.
    tie_kind = {"internal_splicing": [0, 0.0], "five_prime_only": [0, 0.0]}
    unresolved = 0
    n_tied = n_nontied = 0
    tied_mass = nontied_mass = 0.0
    # Every read in the file, N=1 included -- the denominator Main asked for is the
    # whole library, not just the ambiguous part.
    total_library_mass = sum(read_mass.values())
    spread_by_cat = defaultdict(list)

    for readname, members in groups.items():
        if len(members) < 2:
            continue
        coords, chains = [], []
        for tid, _ in members:
            info = ends.get(tid)
            if info is None:
                coords = None
                break
            coords.append(info[0])
            chains.append(info[1])
        if coords is None:
            unresolved += 1
            continue

        ws = [w for _, w in members]
        m = read_mass[readname]
        lo, hi = min(ws), max(ws)
        tied = hi <= lo * (1.0 + SPREAD_EPS)
        idx = 0 if tied else 1
        if tied:
            n_tied += 1
            tied_mass += m
        else:
            n_nontied += 1
            nontied_mass += m

        spread = max(coords) - min(coords)
        if spread == 0:
            cat = "shared_3p_exact"
        elif spread <= args.polyA_window:
            cat = f"within_{args.polyA_window}nt"
        else:
            cat = "distinct_3p"
        cats[cat][idx] += 1
        cat_mass[cat][idx] += m
        if tied and cat != "distinct_3p":
            kind = "five_prime_only" if len(set(chains)) == 1 else "internal_splicing"
            tie_kind[kind][0] += 1
            tie_kind[kind][1] += m
        if lo > 0:
            spread_by_cat[cat].append(hi / lo)

    total_amb = n_tied + n_nontied
    if total_amb == 0:
        sys.exit("Error, no multi-compatible reads found")

    def pctlib(x):
        return 100.0 * x / total_library_mass if total_library_mass else float("nan")

    print(f"# {args.label}: total library mass {total_library_mass:.1f} over "
          f"{len(read_mass)} reads; {total_amb} multi-compatible "
          f"({n_tied} tied, {n_nontied} non-tied); {unresolved} skipped "
          f"(transcript absent from GTF)")
    print(f"# TIED MASS = {tied_mass:.1f} = {pctlib(tied_mass):.3f}% of total library mass"
          f"   (non-tied ambiguous = {pctlib(nontied_mass):.3f}%)")
    print("label\tcategory\tn_tied\tpct_of_tied\ttied_pct_of_library"
          "\tn_nontied\tpct_of_nontied\tmed_ratio_nontied")
    for cat in ("shared_3p_exact", f"within_{args.polyA_window}nt", "distinct_3p"):
        t, nt = cats[cat]
        r = sorted(spread_by_cat[cat])
        medr = r[len(r) // 2] if r else float("nan")
        print(f"{args.label}\t{cat}\t{t}\t"
              f"{100.0 * t / n_tied if n_tied else float('nan'):.2f}\t"
              f"{pctlib(cat_mass[cat][0]):.3f}\t{nt}\t"
              f"{100.0 * nt / n_nontied if n_nontied else float('nan'):.2f}\t{medr:.4f}")
    print("label\ttie_kind\tn_reads\tpct_of_tied\tpct_of_library_mass")
    for kind in ("internal_splicing", "five_prime_only"):
        n, mm = tie_kind[kind]
        print(f"{args.label}\t{kind}\t{n}\t"
              f"{100.0 * n / n_tied if n_tied else float('nan'):.2f}\t{pctlib(mm):.3f}")


if __name__ == "__main__":
    main()
