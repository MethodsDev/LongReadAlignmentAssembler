#!/usr/bin/env python3
"""Checks for testing/ref_splice_pattern_reprieve (v0.50.0 reference identity).

Run after `make` has produced the run_* directories (see Makefile for what each run is).

usage: check_ref_identity.py <reads.bam>
"""

import collections
import os
import re
import sys

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../pylib"))
from Transcript import Transcript  # noqa: E402

GUIDE = "../sirvs/data/SIRVs1-7.annot.half.gtf"
UNFLAGGED_MSG = "allow_unflagged_LRAA_guide"

failures = []


def check(ok, msg):
    print(("ok:   " if ok else "FAIL: ") + msg)
    if not ok:
        failures.append(msg)


def read_models(gtf):
    """{transcript_id: dict(contig, strand, introns, attrs)} from a GTF."""
    models = {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            m = models.setdefault(tid, {"contig": f[0], "strand": f[6], "exons": [], "attrs": {}})
            if f[2] == "exon":
                m["exons"].append((int(f[3]), int(f[4])))
            elif f[2] == "transcript":
                m["attrs"] = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
    for m in models.values():
        ex = sorted(m["exons"])
        m["introns"] = tuple((ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1))
        m["chain"] = (m["contig"], m["strand"], m["introns"]) if m["introns"] else None
    return models


def chain_to(models, key):
    """{chain: attribute value} over the multi-exon models carrying the attribute."""
    return {m["chain"]: m["attrs"][key] for m in models.values() if m["chain"] and key in m["attrs"]}


def bam_chain_counts(bam):
    """Reads per intron chain, strand-agnostic: (contig, introns) -> count. Primary
    alignments, introns from N operations. The SIRV reads carry no long deletions and no
    antisense read repeats a SIRV chain, so this equals the strand-aware count LRAA takes."""
    counts = collections.Counter()
    for read in pysam.AlignmentFile(bam):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        pos, introns = read.reference_start, []
        for op, n in read.cigartuples:
            if op == 3:
                introns.append((pos + 1, pos + n))
            if op in (0, 2, 3, 7, 8):
                pos += n
        if introns:
            counts[(read.reference_name, tuple(introns))] += 1
    return counts


def log_text(path):
    with open(path) as fh:
        return fh.read()


def main():
    bam = sys.argv[1]

    guide_chains = {m["chain"] for m in read_models(GUIDE).values() if m["chain"]}
    rg = read_models("run_rg/rg.LRAA.ref-guided.gtf")
    me = [m for m in rg.values() if m["chain"]]
    se = [m for m in rg.values() if not m["chain"]]

    # 1. attributes of a run guided by a reference annotation
    check(
        all("splice_hash" in m["attrs"] and "ref_splice_pattern" in m["attrs"] for m in me),
        f"rg: all {len(me)} multi-exon models carry splice_hash and ref_splice_pattern",
    )
    check(
        all(
            m["attrs"]["splice_hash"]
            == Transcript(m["contig"], [list(e) for e in m["exons"]], m["strand"]).get_splice_hash()
            for m in me
        ),
        "rg: splice_hash is the hash of the model's intron chain",
    )
    check(all("ref_supported" in m["attrs"] for m in se), f"rg: all {len(se)} single-exon models carry ref_supported")
    check(
        all((m["attrs"]["ref_splice_pattern"] == "True") == (m["chain"] in guide_chains) for m in me),
        "rg: ref_splice_pattern is True exactly for the guide's chains",
    )
    num_true = sum(m["attrs"]["ref_splice_pattern"] == "True" for m in me)
    check(0 < num_true < len(me), f"rg: both reference ({num_true}) and novel ({len(me) - num_true}) chains present")
    check(
        all(("ref_FSM_reads" in m["attrs"]) == (m["attrs"]["ref_splice_pattern"] == "True") for m in me),
        "rg: ref_FSM_reads on reference-chain models only",
    )
    counts = bam_chain_counts(bam)
    wrong = [
        (m["chain"], m["attrs"]["ref_FSM_reads"], counts[(m["contig"], m["introns"])])
        for m in me
        if "ref_FSM_reads" in m["attrs"] and int(m["attrs"]["ref_FSM_reads"]) != counts[(m["contig"], m["introns"])]
    ]
    check(not wrong, f"rg: ref_FSM_reads equals the reads of the full bam carrying the chain ({len(wrong)} differ: {wrong[:3]})")
    collapsed = read_models("run_rg/rg.LRAA.ref-guided.splice_pattern_collapsed.gtf")
    check(
        all("ref_splice_pattern" in m["attrs"] for m in collapsed.values() if m["chain"]),
        "rg: the splice-pattern collapsed GTF keeps ref_splice_pattern",
    )

    # 1b. chunked pipeline: same counts for the same chains
    chunked = chain_to(read_models("run_rg_chunked/rg.LRAA.ref-guided.gtf"), "ref_FSM_reads")
    plain = chain_to(rg, "ref_FSM_reads")
    shared = set(chunked) & set(plain)
    check(
        shared and all(chunked[c] == plain[c] for c in shared),
        f"rg_chunked: ref_FSM_reads equal the unchunked run's for all {len(shared)} shared reference chains",
    )

    # 2. guided by run 1's LRAA GTF: only its True chains are references, pooled counts carried
    rg_true = {c for c, v in chain_to(rg, "ref_splice_pattern").items() if v == "True"}
    guided = read_models("run_guided/guided.LRAA.ref-guided.gtf")
    g_flag = chain_to(guided, "ref_splice_pattern")
    check(
        all((v == "True") == (c in rg_true) for c, v in g_flag.items()),
        "guided: ref_splice_pattern True exactly for the guide's True chains (its novel chains stay novel)",
    )
    g_reads = chain_to(guided, "ref_FSM_reads")
    check(
        g_reads and all(g_reads[c] == plain[c] for c in g_reads),
        f"guided: ref_FSM_reads are the guide's (pooled) counts for all {len(g_reads)} reference chains",
    )

    # 3-5. an LRAA guide without the attributes
    check(UNFLAGGED_MSG in log_text("run_unflagged/run.log"), "unflagged: discovery refuses the unflagged LRAA guide")
    # the Makefile also checks the run left nothing but its log
    check(
        UNFLAGGED_MSG in log_text("run_unflagged_chunked/run.log"),
        "unflagged: the chunked pipeline refuses it at startup, before any chunk work",
    )
    allowed = read_models("run_unflagged_allowed/allowed.LRAA.ref-guided.gtf")
    a_flag = chain_to(allowed, "ref_splice_pattern")
    rg_chains = {m["chain"] for m in me}
    check(
        all(v == "True" for c, v in a_flag.items() if c in rg_chains),
        "unflagged + allow_unflagged_LRAA_guide: every guide chain counts as a reference (pre-v0.50.0 behaviour)",
    )
    check(
        os.path.exists("run_unflagged_quant/quant.LRAA.quant-only.quant.expr")
        and UNFLAGGED_MSG not in log_text("run_unflagged_quant/run.log"),
        "unflagged: quant-only accepts it",
    )

    # 6-7. merge
    merged = read_models("run_merge/merged.gtf")
    m_flag = chain_to(merged, "ref_splice_pattern")
    merged_me = [m for m in merged.values() if m["chain"]]
    check(len(m_flag) == len(merged_me), f"merge: all {len(merged_me)} multi-exon models carry ref_splice_pattern")
    either_true = rg_true | {c for c, v in g_flag.items() if v == "True"}
    check(
        all((v == "True") == (c in either_true) for c, v in m_flag.items()),
        "merge: ref_splice_pattern True exactly for chains an input flagged True",
    )
    max_reads = {}
    for source in (plain, g_reads):
        for c, v in source.items():
            max_reads[c] = max(max_reads.get(c, 0), int(v))
    m_reads = chain_to(merged, "ref_FSM_reads")
    check(
        m_reads and all(int(v) == max_reads[c] for c, v in m_reads.items()),
        f"merge: ref_FSM_reads is the inputs' maximum for all {len(m_reads)} reference chains",
    )
    check(UNFLAGGED_MSG in log_text("run_merge_unflagged/merge.log"), "merge: refuses an unflagged LRAA input")

    if failures:
        print(f"{len(failures)} check(s) FAILED")
        sys.exit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
