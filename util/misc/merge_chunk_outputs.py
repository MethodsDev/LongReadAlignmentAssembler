#!/usr/bin/env python3

"""Stage 6, standalone: merge per-chunk quant output back into one table.

This is the sixth of the six tools ``pylib/ChunkedRun.py`` drives. The other five
were already invocable on their own; this one only existed inside the driver, so
a workflow engine could scatter stages 1-5 across tasks and then had to bring
every chunk back to ONE machine to merge. That made the gather a serial tail on
an otherwise embarrassingly parallel run.

    stage 1  util/separate_bam_by_strand.py             orientation split
    stage 2  util/misc/select_contig_cut_points.py      cut selection
    stage 3  util/misc/extract_contig_region_inputs.py  chunk extraction
    stage 4  util/normalize_bam_by_strand.py            per-chunk normalization
    stage 5  LRAA, quant-only or discovery              per-chunk work
    stage 6  THIS SCRIPT                                merge + coordinate translation

ONE IMPLEMENTATION, NOT TWO. The merge itself stays in
``ChunkedRun.merge_and_translate`` and is called from here. That is deliberate
rather than lazy: ``pylib/test_chunked_merge_whole_contig_frame.py`` asserts the
signatures of ``merge_and_translate`` and ``merge_discovery_gtf`` by
``inspect.signature``, and ``pylib/test_chunked_tracking_streaming_merge.py``
pins the merged bytes against a COMMITTED artifact. Copying the logic here would
have created a second implementation that those tests do not cover, which is the
one outcome worth avoiding. So this file supplies the two things the driver had
and a task does not -- a file-based input contract and a CLI -- and nothing else.

THE MANIFEST is that contract, because a workflow task cannot be handed Python
dicts. It names only the three fields the merge actually reads off a unit
(verified against every ``unit[...]`` access in the merge path), plus an optional
group label this script consumes itself and never passes down:

    {
      "units": [
        {
          "unit_id":      "chr1_00_plus",   # namespaces ids under --discovery
          "quant_prefix": "/path/chunk_quant",  # + .quant.expr/.quant.tracking.gz/.gtf
          "offset":       0,                # rebase applied by extraction; 0 = uncut
          "group":        "chr1"            # OPTIONAL, for --group; ignored by the merge
        }
      ]
    }

UNIT ORDER IN THE MANIFEST IS LOAD-BEARING and is not re-derived here. The driver
sorts units into strand-first order before merging so that strandless and
strand-first runs produce merged tables that line up row for row
(``ChunkedRun.py:2060``). A caller assembling a manifest by hand, or a workflow
collecting scattered outputs in completion order, must preserve that order; this
script merges the units in the order given.

GROUPED MERGES, and the exact identity they give
------------------------------------------------
``--group`` merges only the units carrying that label, so N labels can run as N
parallel tasks. What you get back is NOT byte-identical to one global merge, and
the reason is worth stating because it decides how a workflow assembles them.

``quant.tracking`` is sorted on ``(read_name, transcript_id, gene_id)``. Read
names are not correlated with contig, so a global sort INTERLEAVES rows from every
chromosome, while grouped-then-concatenated output is contig-major with the sort
holding only inside each group. Same rows, same values, different order:
CONTENT-identical, not byte-identical.

Byte-identity is still reachable, cheaply, and this is the useful part: each
group's output is sorted on the GLOBAL key, restricted to its own rows. So the
final assembly of N grouped files is a k-way merge of N sorted streams --
``heapq.merge`` over N open handles, O(N) rows resident, no re-sort. That is the
same shape ``ChunkedRun._merge_tracking_streaming`` already uses to combine its
spilled runs. A workflow therefore wants: N parallel ``--group`` tasks, then one
cheap merge task, and the result is byte-identical to a single global merge.

``quant.expr`` needs a SEPARATE combining step and is NOT safe to concatenate:
its ``TPM`` is recomputed over the merged scope from each row's ``all_reads``
(see ``merge_and_translate``), so a per-group merge normalises TPM over that
GROUP rather than the run. Grouped ``quant.expr`` output is therefore correct
only within its group.

That combining step needs no new arithmetic, because ``merge_and_translate``
already reports ``merged_scope_total_all_reads`` for exactly this reason -- "so
nothing is hidden" -- and already preserves each row's ``all_reads`` untouched.
Summing partitions of a sum is associative, so the N groups' totals sum to the
run's true total, and ``ChunkedRun.rebase_tpm`` applied once more over the union
of their rows with that sum reproduces a single global merge's TPM exactly.
``--combine-expr`` below does this; VERIFIED byte-identical on the
strandless-parity gate corpus, 555 of 555 rows, by
``test_combine_expr_matches_global_merge``.

The grouped mode's real target is ``quant.tracking``, which is the table that
does not fit in memory; ``quant.expr`` is one row per model and would merge
anywhere on its own, EXCEPT for the TPM rebase, which is why the extra step
exists at all.
"""

import argparse
import json
import logging
import os
import sys

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "..", "..", "pylib"])
)

import ChunkedRun

logging.basicConfig(
    format="%(asctime)s : %(levelname)s : %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

# The three fields the merge path reads off a unit. Enforced rather than trusted:
# a manifest missing one would fail deep inside the merge with a KeyError naming a
# dict, which tells a workflow author nothing about which unit or which field.
REQUIRED_UNIT_FIELDS = ("unit_id", "quant_prefix", "offset")


def load_manifest(path):
    """Read and VALIDATE a unit manifest, returning the units in file order.

    Order is preserved deliberately; see the module docstring on why it is
    load-bearing. Validation is eager and names the offending index and field,
    because the alternative is a KeyError raised after the merge has already
    started writing.
    """

    with open(path, "rt") as fh:
        doc = json.load(fh)

    if not isinstance(doc, dict) or "units" not in doc:
        raise ValueError(
            "manifest {} must be an object carrying a 'units' array".format(path)
        )

    units = doc["units"]
    if not isinstance(units, list) or not units:
        raise ValueError("manifest {} carries no units".format(path))

    for i, unit in enumerate(units):
        if not isinstance(unit, dict):
            raise ValueError(
                "manifest {} unit {} is not an object".format(path, i)
            )
        missing = [f for f in REQUIRED_UNIT_FIELDS if f not in unit]
        if missing:
            raise ValueError(
                "manifest {} unit {} ({}) is missing: {}".format(
                    path, i, unit.get("unit_id", "unnamed"), ", ".join(missing)
                )
            )
        # An offset that arrived as a string would shift coordinates by string
        # concatenation or raise inside the merge; both are worse than failing here.
        if not isinstance(unit["offset"], int) or isinstance(unit["offset"], bool):
            raise ValueError(
                "manifest {} unit {} ({}) has a non-integer offset {!r}".format(
                    path, i, unit["unit_id"], unit["offset"]
                )
            )

    return units


def select_group(units, group):
    """Keep the units labelled ``group``, in manifest order.

    A group naming no unit is an error rather than an empty merge: in a scattered
    workflow it means the shard list and the manifest disagree, and an empty
    merged table is a silently wrong answer that looks like a successful task.
    """

    if group is None:
        return units

    kept = [u for u in units if u.get("group") == group]
    if not kept:
        seen = sorted({str(u.get("group")) for u in units})
        raise ValueError(
            "no unit in the manifest carries group {!r}; groups present: {}".format(
                group, ", ".join(seen)
            )
        )
    return kept


def write_manifest(path, units, group_of=None):
    """Emit a manifest for ``units``, so the driver and this script cannot drift.

    ``group_of`` maps a unit to its group label; the default groups by the contig
    embedded in ``unit_id``, which is what a per-chromosome scatter wants.
    """

    def default_group(unit):
        # unit_id is "<chunk_id>_<orientation>" and chunk_id is "<contig>_<NN>",
        # so the contig is everything before the last two underscore fields. Done
        # by rsplit rather than a regex because a contig name may itself carry
        # underscores (chrUn_KI270302v1, chr14_GL000009v2_random).
        parts = str(unit["unit_id"]).rsplit("_", 2)
        return parts[0] if len(parts) == 3 else str(unit["unit_id"])

    picker = group_of or default_group
    doc = {
        "units": [
            {
                "unit_id": u["unit_id"],
                "quant_prefix": u["quant_prefix"],
                "offset": u["offset"],
                "group": picker(u),
            }
            for u in units
        ]
    }
    with open(path, "wt") as ofh:
        json.dump(doc, ofh, indent=2, sort_keys=False)
        ofh.write("\n")
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--manifest",
        default=None,
        help="JSON manifest of units to merge; schema in this script's docstring. "
        "Required unless --combine-expr.",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="run directory; merged tables are written under <outdir>/merged/. "
        "Required unless --combine-expr.",
    )
    parser.add_argument(
        "--discovery",
        action="store_true",
        help="de novo run: namespace gene/transcript ids per unit and merge the "
        "per-chunk model GTFs. Quantifying, ids come from the supplied GTF and "
        "are already global, so this must be OFF or ids will be rewritten.",
    )
    parser.add_argument(
        "--group",
        default=None,
        help="merge only the units labelled this, so N labels can run as N "
        "parallel tasks. quant.tracking from this is BYTE-identical to a global "
        "merge once gathered with heapq.merge (each group is pre-sorted on the "
        "global key); quant.expr from this is correct only WITHIN the group -- "
        "combine it with --combine-expr, not by concatenating. See the module "
        "docstring.",
    )
    parser.add_argument(
        "--result",
        default=None,
        help="write the merge's own report here as JSON (row counts, whether "
        "coordinates were translated, peak resident tracking rows). For a "
        "--group run this is what --combine-expr later reads.",
    )
    parser.add_argument(
        "--combine-expr",
        action="store_true",
        help="COMBINE MODE: recombine N grouped merges' quant.expr into one with "
        "correct run-wide TPM. Takes --group-result (repeated) and --out; ignores "
        "--manifest/--outdir/--discovery/--group. See combine_grouped_expr's "
        "docstring for why this needs no new arithmetic.",
    )
    parser.add_argument(
        "--group-result",
        action="append",
        default=None,
        dest="group_results",
        help="(--combine-expr) a grouped run's --result JSON; repeat once per "
        "group. Each must carry quant_expr and merged_scope_total_all_reads, "
        "which every --group run's --result already does.",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="(--combine-expr) path for the combined quant.expr",
    )
    args = parser.parse_args(argv)

    if args.combine_expr:
        if not args.group_results:
            parser.error("--combine-expr requires at least one --group-result")
        if not args.out:
            parser.error("--combine-expr requires --out")
        logger.info(
            "stage 6 combine: recombining %d grouped quant.expr result(s) into %s",
            len(args.group_results),
            args.out,
        )
        combined = ChunkedRun.combine_grouped_expr(args.group_results, args.out)
        logger.info(
            "combine done: %s expr row(s) from %s group(s), "
            "merged_scope_total_all_reads=%s",
            combined["expr_rows"],
            combined["combined_from_groups"],
            combined["merged_scope_total_all_reads"],
        )
        if args.result:
            with open(args.result, "wt") as ofh:
                json.dump(combined, ofh, indent=2, sort_keys=True)
                ofh.write("\n")
        return 0

    if not args.manifest or not args.outdir:
        parser.error("--manifest and --outdir are required unless --combine-expr")

    units = select_group(load_manifest(args.manifest), args.group)
    logger.info(
        "stage 6: merging %d unit(s)%s into %s",
        len(units),
        "" if args.group is None else " from group {}".format(args.group),
        os.path.join(args.outdir, "merged"),
    )

    merged = ChunkedRun.merge_and_translate(
        args.outdir, units, discovery=args.discovery
    )

    logger.info(
        "stage 6 done: %s expr row(s), %s tracking row(s), coordinates %s, "
        "peak %s tracking row(s) resident",
        merged.get("expr_rows"),
        merged.get("tracking_rows"),
        "translated" if merged.get("coordinates_translated") else "already whole-contig",
        merged.get("tracking_merge_peak_resident_rows"),
    )

    if args.result:
        with open(args.result, "wt") as ofh:
            json.dump(merged, ofh, indent=2, sort_keys=True)
            ofh.write("\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
