#!/usr/bin/env python3

"""Read support per SQANTI-like feature, for both model levels and both support metrics.

Motivation: in ref-guided mode the annotation is seeded with reference transcripts, so a
model can be categorized FSM -- its intron chain matches a reference transcript exactly --
while the library never supports it on its own.  Counting those models alongside observed
ones overstates what the data show, and the overstatement is specific to the ref-guided
regimes, which is exactly the comparison the plots make.

Three columns of LRAA's quant.expr answer three versions of "supported"
(Quantify.py:2152-2257):

  uniq_reads      reads compatible with exactly one isoform -- assigned to this model
                  without EM having to apportion them across competitors.  Defined for
                  monoexonic models too, so every SQANTI-like category is comparable.

  uniq_FSM_reads  the subset of those whose intron chain is exactly this model's.  Zero by
                  construction for a monoexonic model, which has no chain to reproduce, so
                  those features are reported n/a rather than unsupported.

  has_FSM_read    whether ANY assigned read reproduces the chain, exclusivity dropped.
                  New in v0.36.1.

`uniq_FSM_reads == 0` conflates two different
facts about a model, and the CHANGELOG entry for the column says so: "no read traverses
this chain" (the model was assembled out of partial reads) versus "reads do, but each also
fits another model" (the chain was observed whole, but the read is shared with a competitor
-- typically a model with the same chain and different termini).  So this script reports
the FSM metric as a three-way split:

  unique FSM read        uniq_FSM_reads >= 1
  shared FSM read only   uniq_FSM_reads == 0 and has_FSM_read == 1
  no FSM read            has_FSM_read == 0

which is still a partition of exactly the features the counts view plots, so the Rmd's
`check_support_partitions` invariant holds unchanged against the same totals.

`all_reads` (EM mass, which on its own supports nothing) is carried through to the
per-feature output as the looser fourth reading.

Two model levels are emitted, matching the two plot families of the SQLIKE gtf analysis
notebooks (util/sc/notebook_templates/sqlike_support_funcs.R draws them):

  model         one row per transcript model, joined to quant.expr by transcript_id.
  splice_chain  one row per splice-pattern-collapsed chain, joined by splice hashcode --
                the collapsed GTFs name each chain by the hashcode quant.expr carries in
                splice_hash_code -- aggregating support over the models sharing the chain.
                A monoexonic model has no chain to hash, so quant.expr puts its bare
                transcript id in that column; those features fall back to the transcript_id
                join, which is exact (stripping the chunk@ prefix instead would collide --
                comp numbering restarts per chunk).

Read counts aggregate by SUM across clusters and across the models sharing a chain;
`has_FSM_read` aggregates by OR, because it is a per-model flag and not a count.  Summing
it would be meaningless and max() is the same operation as "any cluster saw a whole read".

The scg SQLIKE feature names carry a `GENE^` prefix that the init names do not, except at
loci that got no gene symbol.  For the cluster-guided regimes support is aggregated over
the 14 per-cluster quant.expr files in the release EXPRs tarball, the cluster-guided
analogue of the single pooled init quant.

The inputs come from a regimes table (--regimes, tab-separated, with a header), one row per
(regime, model level), paths relative to --eval_dir:

  type           label used by the notebook, e.g. denovo-basic, refGuided-scg-spC
  level          model | splice_chain
  iso_cats       the regime's SQLIKE *.iso_cats.tsv
  quant          its quant.expr, or the scg cluster_quant.EXPRs.tar.gz; several comma-separated
  gene_prefixed  1 when the feature names carry a GENE^ prefix (the scg regimes), else 0

compare_sqlike_FSM_ISM_flows.py reads the splice_chain rows of the same table.
"""

import argparse
import csv
import gzip
import os
import sys
import tarfile
from collections import defaultdict



def read_regimes(path):
    """(type, level, iso_cats, quant paths, gene_prefixed) rows from a regimes table."""
    regimes = []
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["level"] not in ("model", "splice_chain"):
                raise SystemExit(f"{path}: level must be model or splice_chain: {row}")
            regimes.append(
                (
                    row["type"],
                    row["level"],
                    row["iso_cats"],
                    [q for q in row["quant"].split(",") if q],
                    row["gene_prefixed"].strip() in ("1", "True", "true"),
                )
            )
    if not regimes:
        raise SystemExit(f"no regimes in {path}")
    return regimes


UNIQ_SUPPORTED = "uniquely supported"
UNIQ_UNSUPPORTED = "no unique reads"
FSM_UNIQUE = "unique FSM read"
FSM_SHARED = "shared FSM read only"
FSM_NONE = "no FSM read"
FSM_UNMEASURED = "FSM unmeasured (NA)"
MONOEXONIC = "monoexonic (n/a)"
NOQUANT = "no quant record"

COUNT_COLUMNS = ("uniq_reads", "uniq_FSM_reads")
SUM_COLUMNS = ("all_reads",)
# per-model flag, so OR across clusters and across the models sharing a chain
FLAG_COLUMNS = ("has_FSM_read",)


def new_totals():
    totals = {column: defaultdict(int) for column in COUNT_COLUMNS}
    totals.update({column: defaultdict(float) for column in SUM_COLUMNS})
    totals.update({column: defaultdict(int) for column in FLAG_COLUMNS})
    return totals


def parse_measure(value):
    """quant.expr writes NA where a path computes no per-read chain comparison.

    The --oversimplify aggregate writer emits NA for uniq_FSM_reads and has_FSM_read on
    multi-exon models (LRAA:5632-5635): unmeasured, which is neither 0 nor a number. Return
    None for it rather than crashing on int(), so a run on that path reports "FSM unmeasured
    (NA)" instead of dying here.
    """
    if value in ("NA", "", "."):
        return None
    return int(value)


def accumulate_expr(handle, by_hash, by_tid, n_models, unmeasured):
    """Aggregate the read-support columns per splice hashcode and per transcript id.

    Columns are addressed by name: the trailing columns of quant.expr vary with run mode
    (splice-compatible containment columns are emitted only when that analysis ran, and
    has_FSM_read exists only from v0.36.1 on).
    """
    header = None
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if header is None:
            header = {name: i for i, name in enumerate(fields)}
            required = (
                ("splice_hash_code", "transcript_id")
                + COUNT_COLUMNS
                + SUM_COLUMNS
                + FLAG_COLUMNS
            )
            for column in required:
                if column not in header:
                    raise SystemExit(
                        f"quant.expr lacks required column {column}: {fields}"
                    )
            continue
        hashcode = fields[header["splice_hash_code"]]
        transcript_id = fields[header["transcript_id"]]
        for column in COUNT_COLUMNS:
            value = parse_measure(fields[header[column]])
            if value is None:
                unmeasured[column] += 1
                continue
            by_hash[column][hashcode] += value
            by_tid[column][transcript_id] += value
        for column in SUM_COLUMNS:
            value = float(fields[header[column]])
            by_hash[column][hashcode] += value
            by_tid[column][transcript_id] += value
        for column in FLAG_COLUMNS:
            value = parse_measure(fields[header[column]])
            if value is None:
                unmeasured[column] += 1
                continue
            # OR, not sum: the flag says "some read reproduced this chain", and that stays
            # true of a chain as soon as it is true of one of its models in one cluster
            by_hash[column][hashcode] = max(by_hash[column][hashcode], value)
            by_tid[column][transcript_id] = max(by_tid[column][transcript_id], value)
        n_models["hash"][hashcode] += 1
        n_models["tid"][transcript_id] += 1


def read_quant_sources(paths):
    by_hash, by_tid = new_totals(), new_totals()
    n_models = {"hash": defaultdict(int), "tid": defaultdict(int)}
    unmeasured = defaultdict(int)
    for path in paths:
        if path.endswith(".tar.gz"):
            with tarfile.open(path, "r:gz") as tar:
                members = [m for m in tar if m.isfile() and m.name.endswith(".quant.expr")]
                if not members:
                    raise SystemExit(f"no *.quant.expr members in {path}")
                for member in members:
                    stream = tar.extractfile(member)
                    accumulate_expr(
                        (line.decode() for line in stream),
                        by_hash,
                        by_tid,
                        n_models,
                        unmeasured,
                    )
                print(f"  {len(members)} cluster quant.expr from {path}", file=sys.stderr)
        else:
            with open(path) as handle:
                accumulate_expr(handle, by_hash, by_tid, n_models, unmeasured)
            print(f"  {path}", file=sys.stderr)
    for column, n in sorted(unmeasured.items()):
        print(f"  WARNING: {n} rows carry NA for {column}", file=sys.stderr)
    return by_hash, by_tid, n_models, unmeasured


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--regimes", required=True, help="regimes table (see the module docstring)"
    )
    parser.add_argument(
        "--eval_dir",
        default="..",
        help="directory the regimes table's paths are relative to (the eval dir)",
    )
    parser.add_argument(
        "--output_prefix",
        required=True,
        help="written as <prefix>.summary.tsv and <prefix>.features.tsv.gz",
    )
    args = parser.parse_args()
    regimes = read_regimes(args.regimes)

    summary_rows = []
    features_path = f"{args.output_prefix}.features.tsv.gz"
    quant_cache = {}

    with gzip.open(features_path, "wt", newline="") as features_fh:
        features = csv.writer(features_fh, delimiter="\t", lineterminator="\n")
        features.writerow(
            [
                "type",
                "level",
                "feature_name",
                "Category",
                "num_exon_segments",
                "uniq_read_support",
                "FSM_read_support",
                "uniq_reads",
                "uniq_FSM_reads",
                "has_FSM_read",
                "all_reads",
                "num_models",
            ]
        )

        for type_label, level, iso_cats_rel, quant_rels, gene_prefixed in regimes:
            print(f"{type_label} ({level}):", file=sys.stderr)
            iso_cats_path = os.path.join(args.eval_dir, iso_cats_rel)
            quant_paths = tuple(os.path.join(args.eval_dir, rel) for rel in quant_rels)

            # the two levels of a regime read the same quant, so parse it once
            if quant_paths not in quant_cache:
                quant_cache[quant_paths] = read_quant_sources(quant_paths)
            else:
                print(f"  (quant already parsed)", file=sys.stderr)
            by_hash, by_tid, n_models, unmeasured = quant_cache[quant_paths]

            counts = defaultdict(int)
            with open(iso_cats_path) as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    feature = row["feature_name"]
                    # scg feature names are GENE^key, except at loci that got no gene symbol
                    key = feature.rsplit("^", 1)[-1] if gene_prefixed else feature
                    category = row["sqanti_cat"]
                    monoexonic = int(row["num_exon_segments"]) == 1

                    if level == "splice_chain" and key in n_models["hash"]:
                        totals, counted = by_hash, n_models["hash"]
                    else:
                        totals, counted = by_tid, n_models["tid"]

                    if key not in counted:
                        uniq_support = NOQUANT
                        fsm_support = NOQUANT
                        uniq_reads = fsm_reads = has_fsm = all_reads = None
                        num_models = 0
                    else:
                        uniq_reads = totals["uniq_reads"][key]
                        fsm_reads = totals["uniq_FSM_reads"][key]
                        has_fsm = totals["has_FSM_read"][key]
                        all_reads = totals["all_reads"][key]
                        num_models = counted[key]
                        uniq_support = (
                            UNIQ_SUPPORTED if uniq_reads >= 1 else UNIQ_UNSUPPORTED
                        )
                        if monoexonic:
                            # no intron chain to reproduce, so unmeasurable rather than
                            # unsupported -- and both columns are 0 by construction
                            fsm_support = MONOEXONIC
                        elif unmeasured["has_FSM_read"] and not has_fsm and not fsm_reads:
                            # only reachable on an --oversimplify aggregate quant, where the
                            # flag is NA; without it "no FSM read" is not a claim we can make
                            fsm_support = FSM_UNMEASURED
                        elif fsm_reads >= 1:
                            fsm_support = FSM_UNIQUE
                        elif has_fsm:
                            fsm_support = FSM_SHARED
                        else:
                            fsm_support = FSM_NONE

                    counts[(category, "uniq_reads", uniq_support)] += 1
                    counts[(category, "FSM_reads", fsm_support)] += 1
                    features.writerow(
                        [
                            type_label,
                            level,
                            feature,
                            category,
                            row["num_exon_segments"],
                            uniq_support,
                            fsm_support,
                            "" if uniq_reads is None else uniq_reads,
                            "" if fsm_reads is None else fsm_reads,
                            "" if has_fsm is None else has_fsm,
                            "" if all_reads is None else f"{all_reads:.1f}",
                            num_models,
                        ]
                    )

            for (category, metric, support), count in sorted(counts.items()):
                summary_rows.append([type_label, level, category, metric, support, count])

    summary_path = f"{args.output_prefix}.summary.tsv"
    with open(summary_path, "w", newline="") as summary_fh:
        summary = csv.writer(summary_fh, delimiter="\t", lineterminator="\n")
        summary.writerow(["type", "level", "Category", "metric", "support", "Count"])
        summary.writerows(summary_rows)

    print(f"wrote {summary_path} and {features_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
