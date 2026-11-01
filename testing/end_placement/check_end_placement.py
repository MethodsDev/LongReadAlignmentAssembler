#!/usr/bin/env python3
"""Checks for testing/end_placement: where isoform ends land (v0.48.0 end rules, merge end
logic, SE weak-junction guard, internal TSS of an ISM isoform).

Coordinates are those of the extracted minigenomes (1-based from each region's start); the
regions and their origin are in README.md. Every check names the locus and what moved.

usage: check_end_placement.py <locus> <run_dir>
  locus: TIGD7 | TIGD7_merge | SPACER_revert | SEPHS2 | COMT
"""

import glob
import re
import sys


def read_models(gtf):
    """{transcript_id: dict(strand, exons, attrs)} from an LRAA GTF."""
    models = {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            m = models.setdefault(tid, {"strand": f[6], "exons": [], "attrs": {}})
            if f[2] == "exon":
                m["exons"].append((int(f[3]), int(f[4])))
            elif f[2] == "transcript":
                m["attrs"] = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
    for m in models.values():
        m["exons"].sort()
        ex = m["exons"]
        m["introns"] = tuple((ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1))
        m["five"] = ex[0][0] if m["strand"] == "+" else ex[-1][1]
        m["three"] = ex[-1][1] if m["strand"] == "+" else ex[0][0]
    return models


def read_quant(path):
    """{transcript_id: (uniq_reads, all_reads)}"""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("gene_id"):
                continue
            f = line.split("\t")
            out[f[1]] = (int(f[2]), float(f[3]))
    return out


FAILS = []


def check(ok, msg):
    print(("ok:   " if ok else "FAIL: ") + msg)
    if not ok:
        FAILS.append(msg)


def one(pattern):
    hits = glob.glob(pattern)
    if len(hits) != 1:
        sys.exit(f"FAIL: expected one file matching {pattern}, found {hits}")
    return hits[0]


def by_chain(models, strand, introns):
    return [m for m in models.values() if m["strand"] == strand and m["introns"] == tuple(introns)]


def near(a, b, tol=2):
    return a is not None and abs(a - b) <= tol


# TIGD7 (BT474 chr16, minus strand). Two first-exon acceptors (43213 and 42995) and three
# 5' ends: 43431 (unflagged, the reads' consensus inside node 43213-43728), 43728 (TSS site)
# and 43975. The defect fixed for v0.48.0: an extended end must stay at its read position
# inside the node it reaches, not run on to that node's far edge.
TIGD7_FORMS = {
    # introns (lend, rend) -> (5' end, TSS flag)
    ((40011, 43212),): (43431, "False"),
    ((38629, 43212),): (43728, "True"),
    ((40011, 42994),): (43975, "False"),
}


def check_TIGD7(run_dir):
    models = read_models(one(f"{run_dir}/*.LRAA.ref-free.gtf"))
    quant = read_quant(one(f"{run_dir}/*.LRAA.ref-free.quant.expr"))
    for introns, (five, tss) in TIGD7_FORMS.items():
        hits = by_chain(models, "-", introns)
        check(len(hits) == 1, f"TIGD7: one model with chain {introns} (found {len(hits)})")
        for m in hits:
            check(near(m["five"], five), f"TIGD7 {introns}: 5' end {m['five']}, expected {five}")
            check(m["attrs"].get("TSS") == tss, f"TIGD7 {introns}: TSS {m['attrs'].get('TSS')}, expected {tss}")
    dom = by_chain(models, "-", ((40011, 43212),))
    if dom:
        tid = [t for t, m in models.items() if m is dom[0]][0]
        uniq = quant.get(tid, (0, 0))[0]
        check(uniq >= 10, f"TIGD7 dominant form keeps its reads in quant (unique {uniq}, v0.47.0: 15)")


# Merge of the two halves of the TIGD7 reads, run separately at v0.47.0 (inputs in data/).
# Each merged model's 5' end is its input models' end: 43409/43410 (two inputs), 43431, 43431,
# 43975. v0.47.0's merge moved the first three to node edges (43975, 43975, 43212).
TIGD7_MERGE_FORMS = {
    ((40011, 43212),): 43409,
    ((38629, 43212),): 43431,
    ((38629, 42994),): 43431,
    ((40011, 42994),): 43975,
}


def check_TIGD7_merge(run_dir):
    models = read_models(one(f"{run_dir}/merged.gtf"))
    minus = [m for m in models.values() if m["strand"] == "-"]
    check(len(minus) == len(TIGD7_MERGE_FORMS), f"TIGD7 merge: {len(minus)} minus-strand models, expected {len(TIGD7_MERGE_FORMS)}")
    for introns, five in TIGD7_MERGE_FORMS.items():
        hits = by_chain(models, "-", introns)
        check(len(hits) == 1, f"TIGD7 merge: one model with chain {introns} (found {len(hits)})")
        for m in hits:
            check(near(m["five"], five), f"TIGD7 merge {introns}: 5' end {m['five']}, the input models' {five}")
            check(m["attrs"].get("TSS") == "False", f"TIGD7 merge {introns}: no TSS flag (no input has one)")


# SPACER revert (BT474 chr19:14,071,536-14,207,187, minus strand). The model's reads include a
# contained sibling's reads ending on exon segment 46154-46250, which its 5' exon (45967-46079)
# has no edge to. Extending the 5' end there put a SPACER in the model's path and the remap
# dropped it; the guard reverts the end into the exon and the model survives.
SPACER_CHAIN = (
    (21703, 22092), (22258, 25825), (25949, 26043), (26140, 26228), (26356, 29290),
    (29374, 31280), (31380, 35224), (35354, 35812), (35875, 45966),
)


def check_SPACER_revert(run_dir):
    models = read_models(one(f"{run_dir}/*.LRAA.ref-free.gtf"))
    hits = by_chain(models, "-", SPACER_CHAIN)
    check(len(hits) == 1, f"SPACER_revert: the 10-exon minus-strand model is reported (found {len(hits)})")
    for m in hits:
        check(45967 <= m["five"] <= 46079, f"SPACER_revert: 5' end {m['five']} inside its exon 45967-46079")
    logs = glob.glob(f"{run_dir}/*.chunked_work/logs/*") + glob.glob(f"{run_dir}/run.log")
    text = ""
    for p in logs:
        with open(p, errors="replace") as fh:
            text += fh.read()
    check("EXTENSION_REVERTED" in text, "SPACER_revert: the guard logged EXTENSION_REVERTED")
    check("Found SPACER in input transcript" not in text, "SPACER_revert: no model skipped for a SPACER")


# SEPHS2 (BT474 chr16, minus strand): single-exon TUSCO gene, ~1,400 unspliced full-length
# reads against a 2-read GT-AG junction. v0.47.0 reported only a spliced model (213 reads).
def check_SEPHS2(run_dir):
    models = read_models(one(f"{run_dir}/*.LRAA.ref-free.gtf"))
    quant = read_quant(one(f"{run_dir}/*.LRAA.ref-free.quant.expr"))
    se = [(t, m) for t, m in models.items() if m["strand"] == "-" and len(m["exons"]) == 1
          and m["exons"][0][0] <= 5100 and m["exons"][0][1] >= 7100]
    check(len(se) == 1, f"SEPHS2: one single-exon model over the gene (found {len(se)})")
    for t, m in se:
        uniq = quant.get(t, (0, 0))[0]
        check(uniq >= 1000, f"SEPHS2: single-exon model holds the unspliced reads (unique {uniq})")


# COMT (BT474 chr22, plus strand): S-COMT starts inside exon 2 of MB-COMT at an internal
# promoter (genomic chr22:19,962,553). Both TSSs must be called, the S form being an ISM of
# the long form. This is the capability the paper's internal-TSS study rests on.
def check_COMT(run_dir):
    models = read_models(one(f"{run_dir}/*.LRAA.ref-free.gtf"))
    plus = [m for m in models.values() if m["strand"] == "+" and m["attrs"].get("TSS") == "True"]
    long_forms = [m for m in plus if near(m["five"], 67355, 3)]
    short_forms = [m for m in plus if near(m["five"], 88037, 3)]
    check(bool(long_forms), "COMT: MB-COMT 5' end (local 67355) with a TSS")
    check(bool(short_forms), "COMT: S-COMT internal 5' end (local 88037) with a TSS")
    ism = [s for s in short_forms for lf in long_forms
           if s["introns"] and lf["introns"][-len(s["introns"]):] == s["introns"]
           and any(a <= 88037 <= b for a, b in lf["exons"])]
    check(bool(ism), "COMT: an S-COMT model is an ISM of an MB-COMT model, starting inside its exon")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    locus, run_dir = sys.argv[1], sys.argv[2]
    {
        "TIGD7": check_TIGD7,
        "TIGD7_merge": check_TIGD7_merge,
        "SPACER_revert": check_SPACER_revert,
        "SEPHS2": check_SEPHS2,
        "COMT": check_COMT,
    }[locus](run_dir)
    if FAILS:
        sys.exit(f"{locus}: {len(FAILS)} check(s) failed")
