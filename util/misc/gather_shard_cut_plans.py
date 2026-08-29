#!/usr/bin/env python3
# encoding: utf-8

"""Gather a scattered run's per-shard cut geometry into ONE shared chunk plan.

WHY THIS EXISTS. A chunked LRAA run selects its own cut positions inside its own
make-chunks phase, where selection OVERLAPS extraction -- ``run_prep_concurrently``
submits a contig's extractions the moment that contig's selection finishes. The
other way to obtain a shareable partition, ``ChunkedRun.py --emit_cut_plan``, pays
a second full selection in a task the run then waits on, which buys a file at the
cost of that overlap. So for a run that has ALREADY chosen its cuts the right move
is not to re-select them, it is to GATHER the ones it chose.

Every real chunked run therefore writes ``<work>/cuts/shard_cut_plan.json``
(``ChunkedRun.SHARD_CUT_PLAN_NAME``), which is a plan in the exact format
``--chunk_plan`` consumes, restricted to the contigs THAT shard processed. A
scattered workflow has one per shard; ``subwdls/LRAA_runner.wdl`` surfaces it as
``LRAA_shard_cut_plan`` and ``LRAA.wdl`` gathers them here.

WHAT THIS MUST GET RIGHT, and both halves are refusals rather than resolutions.

THE ENVELOPE HAS EXACTLY ONE SOURCE, keyed with the empty string.
``ChunkedRun.validate_cut_plan_geometry`` refuses anything else outright: a shared
plan is strandless geometry, and strand-first geometry is selected over per-run
stage-1 split bams and is not shareable at all. So the gather does not concatenate
``by_source`` lists -- N shards would produce N sources and a file its own
consumers reject. Every shard's contigs are merged UNDER one source, which is
sound because a strandless shard has exactly one source to begin with and because
contigs are disjoint across shards.

THE GEOMETRY PARAMETERS MUST BE IDENTICAL ACROSS SHARDS. The cuts themselves
cannot conflict -- disjoint contigs -- but the ``params`` block is a claim about
what CHOSE them, and it is the claim
``ChunkedRun.validate_cut_plan_geometry`` checks a consuming run against. Taking
the first shard's silently would produce a plan whose recorded geometry does not
match part of what produced it: applied by a consumer, the check would pass on the
recorded values and the run would then extract at bounds selected under different
ones. That is precisely the failure the check exists to catch, arriving with the
check reporting success. So a disagreement is named -- parameter, both values,
both shard paths -- and refused.

WHAT THE GATHERED PLAN DELIBERATELY DOES NOT CLAIM.

``num_total_reads`` is null and ``chunks_extracted`` false, exactly as
``--emit_cut_plan`` records them. The TPM denominator belongs to whoever
quantifies, and per-shard denominators are per-shard by construction. The chunk
directories the shards built are their own working files at their own paths, which
no consumer can address, so ``--only_chunk`` is refused by name (see
``ChunkedRun.rebuild_chunk_record``) rather than left to fail on an absent
manifest.

``by_source[0].bam`` and ``bam_identity`` are null. In a per-shard run each is
that shard's own per-contig bam, and there is no single file the gathered geometry
was selected from. Both are provenance that no consumer reproduces
(``write_chunk_plan``'s own comment says so), and stating one shard's bam for a
genome-wide plan would be a false claim. The real provenance is under
``gathered_from``, one entry per contributing shard.

``chunks`` IS rebuilt, through ``ChunkedRun.planned_chunks`` -- the same generator
stage 3 iterates -- over the merged selections in sorted-contig order. That is the
order ``ChunkedRun.enumerate_prep_contigs`` enumerates in (``sorted(lengths)``), so
the ids, regions and global ``order`` counters are the ones a single genome-wide
run would have assigned rather than per-shard counters restarting at zero.
"""

import argparse
import json
import os
import sys

sys.path.insert(
    0,
    os.path.sep.join(
        [os.path.dirname(os.path.realpath(__file__)), "..", "..", "pylib"]
    ),
)

import ChunkedRun  # noqa: E402  (path insert must precede the import)


class GatherError(Exception):
    """A refusal. Every message names the shard(s) and the values involved."""


def load_shard_plan(path):
    """One shard's plan, checked for the properties the merge depends on.

    Version and mode first, then the single-strandless-source envelope: a shard
    that is not strandless has nothing shareable to contribute, and saying so here
    names the file rather than letting the merged output be refused later by a
    consumer that cannot say which shard spoiled it.
    """

    plan = ChunkedRun.load_chunk_plan(path)
    geometry = plan.get("geometry")
    if not geometry:
        raise GatherError(
            "shard plan {} carries no geometry block, so it states no cut "
            "positions. Only a file written by a chunked run's make-chunks phase, "
            "--emit_cut_plan or --stop_after_make_chunks can be gathered.".format(
                path
            )
        )
    mode = geometry.get("mode")
    if mode != ChunkedRun.STRANDLESS_MODE:
        raise GatherError(
            "shard plan {} was selected in {} mode, and only {} geometry is "
            "shareable: strand-first selects over that run's own stage-1 "
            "orientation split, so its cut positions are not applicable to any "
            "other run.".format(path, mode, ChunkedRun.STRANDLESS_MODE)
        )
    by_source = geometry.get("by_source") or []
    if len(by_source) != 1 or by_source[0].get("key") != "":
        raise GatherError(
            "shard plan {} carries {} selection source(s) {}, and a strandless "
            "plan has exactly one keyed \"\". Merging this would give the gathered "
            "envelope more than one source, which ChunkedRun."
            "validate_cut_plan_geometry refuses.".format(
                path,
                len(by_source),
                [entry.get("key") for entry in by_source],
            )
        )
    if not geometry.get("params"):
        raise GatherError(
            "shard plan {} records no cut parameters, so there is nothing to "
            "check the other shards' against. A plan that does not state every "
            "parameter deciding a cut cannot be applied.".format(path)
        )
    if not by_source[0].get("selections"):
        raise GatherError(
            "shard plan {} names no contig selections, so it contributes no "
            "geometry. A shard that selected nothing did not run make-chunks; "
            "gathering it would silently drop its contigs from the "
            "plan.".format(path)
        )
    return plan


def _agree(name, values):
    """One value from every shard's claim about ``name``, or a named refusal.

    ``values`` is ``[(path, value), ...]``. The message carries every distinct
    value with the shards that stated it, because with more than two shards
    "these disagree" does not say which one to look at.
    """

    distinct = []
    for path, value in values:
        for seen_value, paths in distinct:
            if seen_value == value:
                paths.append(path)
                break
        else:
            distinct.append((value, [path]))
    if len(distinct) == 1:
        return distinct[0][0]
    raise GatherError(
        "the shards disagree about {}: {}. A gathered plan states ONE geometry, "
        "and a plan whose recorded parameters do not match what produced part of "
        "it would pass ChunkedRun.validate_cut_plan_geometry in a consuming run "
        "and then extract at bounds selected under different ones.".format(
            name,
            "; ".join(
                "{!r} from {}".format(value, ", ".join(paths))
                for value, paths in distinct
            ),
        )
    )


def check_geometry_agreement(loaded):
    """The ``params`` block every shard must state identically.

    Checked KEY BY KEY rather than as a dict comparison, so a disagreement names
    the parameter and both values instead of dumping two dicts. A key missing from
    some shards is itself a disagreement and is reported as one: the set of keys
    is what ``cut_geometry_params`` decided a cut on, so a shard that omits one
    was not selecting on the same rules.
    """

    names = set()
    for _path, plan in loaded:
        names.update(plan["geometry"]["params"])
    params = {}
    for name in sorted(names):
        params[name] = _agree(
            "the cut parameter {}".format(name),
            [
                (path, plan["geometry"]["params"].get(name, "<absent>"))
                for path, plan in loaded
            ],
        )
    return params


def merge_selections(loaded):
    """Every shard's per-contig selections under one source, in contig order.

    Sorted by contig NAME, which is the order a single genome-wide run enumerates
    in (``ChunkedRun.enumerate_prep_contigs`` returns ``sorted(lengths)``), so the
    rebuilt chunk ids and ``order`` counters below are the ones such a run would
    have assigned.

    A contig claimed by two shards is refused rather than deduplicated. Shards
    partition the genome, so the same contig in two of them means the shards were
    not the partition they are being gathered as -- and one contig has one
    partition, which ``selections_from_chunk_plan`` also refuses on the way in.
    """

    owner = {}
    selections = {}
    for path, plan in loaded:
        for selection in plan["geometry"]["by_source"][0]["selections"]:
            chrom = selection["chrom"]
            if chrom in owner:
                raise GatherError(
                    "contig {} is partitioned by two shards, {} and {}. One "
                    "contig has one partition; these shards are not a partition "
                    "of the genome.".format(chrom, owner[chrom], path)
                )
            owner[chrom] = path
            selections[chrom] = selection
    return [selections[chrom] for chrom in sorted(selections)], owner


def rebuilt_chunks(tag, selections, has_gtf):
    """The chunk entries a genome-wide run over these selections would plan.

    Built through ``ChunkedRun.planned_chunks`` on purpose: it is the generator
    stage 3 iterates and ``write_chunk_plan`` builds its own entries from, so the
    gathered plan cannot describe a partition LRAA would not build. The source
    tuple carries no bam and no parent token because neither is read for the
    fields taken here -- ``planned_chunks_for_selection`` uses ``key`` and ``tag``
    for the ids and regions, and nothing else.
    """

    source = ("", tag, None, None)
    return [
        {
            "chunk_id": planned["chunk_id"],
            "chrom": planned["chrom"],
            "strand": planned["key"],
            "strandless": not planned["key"],
            "region": planned["region"],
            "index": planned["index"],
            "order": planned["order"],
            "has_gtf": bool(has_gtf),
        }
        for planned in ChunkedRun.planned_chunks([source], {"": selections})
    ]


def shard_has_gtf(path, plan):
    """Whether the shard's chunks were extracted with an annotation.

    Read off the shard's own chunk entries, where ``write_chunk_plan`` records
    ``bool(args.gtf)``, and required to be uniform within the shard: it is one
    property of one run, so two values in one file means the file was assembled
    from two.
    """

    values = [(path, bool(entry.get("has_gtf"))) for entry in plan.get("chunks") or []]
    if not values:
        raise GatherError(
            "shard plan {} names no chunks, so it does not say whether its cuts "
            "were placed against an annotation.".format(path)
        )
    return _agree("has_gtf within {}".format(path), values)


def gather(paths, output_path):
    loaded = [(path, load_shard_plan(path)) for path in paths]
    params = check_geometry_agreement(loaded)
    tag = _agree(
        "the selection source tag",
        [(path, plan["geometry"]["by_source"][0].get("tag")) for path, plan in loaded],
    )
    discovery = _agree(
        "discovery mode",
        [(path, bool(plan.get("discovery"))) for path, plan in loaded],
    )
    lraa_suffix = _agree(
        "the LRAA output suffix",
        [(path, plan.get("lraa_suffix")) for path, plan in loaded],
    )
    has_gtf = _agree(
        "has_gtf", [(path, shard_has_gtf(path, plan)) for path, plan in loaded]
    )
    selections, owner = merge_selections(loaded)

    plan = {
        "version": ChunkedRun.CHUNK_PLAN_VERSION,
        "chunks_extracted": False,
        "num_total_reads": None,
        "discovery": discovery,
        "lraa_suffix": lraa_suffix,
        "geometry": {
            "mode": ChunkedRun.STRANDLESS_MODE,
            "params": params,
            "by_source": [
                {
                    "key": "",
                    "tag": tag,
                    "bam": None,
                    "bam_identity": None,
                    "selections": selections,
                }
            ],
        },
        "chunks": rebuilt_chunks(tag, selections, has_gtf),
        # Not part of the schema a consumer reads -- load_chunk_plan checks the
        # version and validate_cut_plan_geometry the geometry, and neither looks
        # here. It is what replaces the single ``bam`` a one-run plan carries: the
        # shards this geometry actually came from, and which contigs each brought.
        "gathered_from": [
            {
                "plan": path,
                "bam": plan_["geometry"]["by_source"][0].get("bam"),
                "bam_identity": plan_["geometry"]["by_source"][0].get("bam_identity"),
                "contigs": sorted(
                    chrom for chrom, source in owner.items() if source == path
                ),
            }
            for path, plan_ in loaded
        ],
    }
    with open(output_path, "wt") as fh:
        json.dump(plan, fh, indent=2, sort_keys=True)
        print("", file=fh)
    return plan


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--shard_plans",
        nargs="+",
        required=True,
        metavar="PATH",
        help="one shard's <work>/cuts/shard_cut_plan.json per scattered LRAA "
        "task, in any order (the merged contigs are sorted by name)",
    )
    parser.add_argument(
        "--output",
        required=True,
        metavar="PATH",
        help="where to write the gathered plan, ready to hand to "
        "LRAA --chunk_plan / ChunkedRun.py --chunk_plan",
    )
    args = parser.parse_args(argv)

    try:
        plan = gather(args.shard_plans, args.output)
    except (GatherError, ChunkedRun.PipelineError) as err:
        print("Error: {}".format(err), file=sys.stderr)
        return 1

    selections = plan["geometry"]["by_source"][0]["selections"]
    print(
        "gathered {} shard(s) into {}: {} contig(s), {} cut(s), {} chunk(s) of "
        "geometry. Apply it with --chunk_plan {}".format(
            len(args.shard_plans),
            args.output,
            len(selections),
            sum(len(selection["cuts"]) for selection in selections),
            len(plan["chunks"]),
            args.output,
        ),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
