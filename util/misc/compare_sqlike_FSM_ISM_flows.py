#!/usr/bin/env python3

"""Where FSM and ISM assignments move between run modes, at the splice-chain level.

The FSM and ISM bars differ substantially between basic and scg, and between ref-guided and
de novo.  A count alone cannot say whether a mode found new reference transcripts or merely
reconstructed known ones to full length: an FSM chain reproduces a reference transcript's
intron chain exactly, an ISM chain is contained within one, so extending a chain converts
ISM to FSM without any new reference transcript being involved.

Two keys, answering two different questions:

  chain   the splice hashcode, i.e. the thing the bars count.  Cross-tabulating a chain's
          category in mode A against its category in mode B decomposes each bar into the
          part inherited from the other mode and the part that is new to this one.

  ref     the reference transcript named in matching_isoforms.  Its status in a mode is FSM
          if any chain reproduces it exactly, else ISM if any chain is contained within it,
          else absent.  This is the user-facing question -- how many annotated transcripts
          does the mode recover in full, rather than how many models it emitted.

Category strings other than FSM/ISM collapse to "other" (a chain can also be NIC, NNIC,
antisense ...), and a key missing from a mode entirely is "absent".
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

# the four splice-chain modes, by the type labels of the regimes table
MODE_LABELS = ["denovo-basic-spC", "denovo-scg-spC", "refGuided-basic-spC", "refGuided-scg-spC"]


def read_modes(path):
    """(label, iso_cats, gene_prefixed) for the splice_chain rows of a regimes table, in
    MODE_LABELS order (collect_sqlike_read_support.py documents the table)."""
    rows = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["level"] == "splice_chain":
                rows[row["type"]] = (
                    row["type"],
                    row["iso_cats"],
                    row["gene_prefixed"].strip() in ("1", "True", "true"),
                )
    missing = [label for label in MODE_LABELS if label not in rows]
    if missing:
        raise SystemExit(f"{path}: no splice_chain rows for {missing}")
    return [rows[label] for label in MODE_LABELS]


# (comparison label, mode A, mode B) -- A is the baseline the bar is being read against
COMPARISONS = [
    ("denovo: basic -> scg", "denovo-basic-spC", "denovo-scg-spC"),
    ("refGuided: basic -> scg", "refGuided-basic-spC", "refGuided-scg-spC"),
    ("basic: denovo -> refGuided", "denovo-basic-spC", "refGuided-basic-spC"),
    ("scg: denovo -> refGuided", "denovo-scg-spC", "refGuided-scg-spC"),
]

ABSENT = "absent"
OTHER = "other"
STATUSES = ["FSM", "ISM", OTHER, ABSENT]


def read_mode(path, gene_prefixed):
    """chain hashcode -> category, and reference transcript -> best status (FSM beats ISM)."""
    chain_cat = {}
    fsm_refs = set()
    ism_refs = set()
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            feature = row["feature_name"]
            # scg feature names are GENE^hashcode, except at loci that got no gene symbol
            key = feature.rsplit("^", 1)[-1] if gene_prefixed else feature
            category = row["sqanti_cat"]
            chain_cat[key] = category
            if category not in ("FSM", "ISM"):
                continue
            matches = row["matching_isoforms"]
            if not matches:
                continue
            target = fsm_refs if category == "FSM" else ism_refs
            target.update(matches.split(","))
    # a reference transcript reproduced exactly by one chain and contained by another is
    # recovered in full, so FSM takes precedence
    ref_status = {ref: "FSM" for ref in fsm_refs}
    for ref in ism_refs:
        ref_status.setdefault(ref, "ISM")
    return chain_cat, ref_status


def collapse(category):
    return category if category in ("FSM", "ISM") else OTHER


def transitions(a_map, b_map, collapse_values):
    """Cross-tab of status in A against status in B over the union of keys."""
    counts = defaultdict(int)
    for key in a_map.keys() | b_map.keys():
        a = a_map.get(key)
        b = b_map.get(key)
        a = ABSENT if a is None else (collapse(a) if collapse_values else a)
        b = ABSENT if b is None else (collapse(b) if collapse_values else b)
        if a == ABSENT and b == ABSENT:
            continue
        counts[(a, b)] += 1
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--regimes",
        required=True,
        help="regimes table, as for collect_sqlike_read_support.py (its splice_chain rows)",
    )
    parser.add_argument(
        "--eval_dir",
        default="..",
        help="directory the regimes table's paths are relative to (the eval dir)",
    )
    parser.add_argument(
        "--output_prefix",
        required=True,
        help="written as <prefix>.transitions.tsv and <prefix>.ref_status.tsv",
    )
    args = parser.parse_args()
    modes = read_modes(args.regimes)

    chains = {}
    refs = {}
    for label, iso_cats_rel, gene_prefixed in modes:
        path = os.path.join(args.eval_dir, iso_cats_rel)
        chains[label], refs[label] = read_mode(path, gene_prefixed)
        n_fsm = sum(1 for status in refs[label].values() if status == "FSM")
        print(
            f"{label}: {len(chains[label])} chains, {len(refs[label])} reference transcripts "
            f"matched ({n_fsm} in full)",
            file=sys.stderr,
        )

    transitions_path = f"{args.output_prefix}.transitions.tsv"
    with open(transitions_path, "w", newline="") as fh:
        out = csv.writer(fh, delimiter="\t", lineterminator="\n")
        out.writerow(["comparison", "key", "from_mode", "to_mode", "from_status", "to_status", "n"])
        for label, mode_a, mode_b in COMPARISONS:
            for key_kind, maps in (("chain", chains), ("ref_transcript", refs)):
                counts = transitions(
                    maps[mode_a], maps[mode_b], collapse_values=(key_kind == "chain")
                )
                for from_status in STATUSES:
                    for to_status in STATUSES:
                        n = counts.get((from_status, to_status), 0)
                        if n:
                            out.writerow(
                                [label, key_kind, mode_a, mode_b, from_status, to_status, n]
                            )

    status_path = f"{args.output_prefix}.ref_status.tsv"
    all_refs = set()
    for ref_status in refs.values():
        all_refs.update(ref_status)
    with open(status_path, "w", newline="") as fh:
        out = csv.writer(fh, delimiter="\t", lineterminator="\n")
        mode_labels = [label for label, _, _ in modes]
        out.writerow(["ref_transcript"] + mode_labels)
        for ref in sorted(all_refs):
            out.writerow([ref] + [refs[label].get(ref, ABSENT) for label in mode_labels])

    print(f"wrote {transitions_path} and {status_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
