#!/usr/bin/env python3

"""Collate N per-invocation read-assignment summaries into ONE cross-cluster table.

Each LRAA invocation writes one ``<prefix>.read_assignment.summary.tsv``: a
``worker`` row per work unit and one recomputed ``TOTAL``. On a cluster-guided
single-cell run there are N such files, one per cell cluster, and nothing brought
them together -- so the question a run is actually asked ("where did the reads
go, across the whole experiment, and which chunk is the outlier") could only be
answered by opening N files by hand.

THE TABLE HAS THREE LEVELS, and they are nested rather than parallel:

    chunk         one row per contributing quant unit of one cluster, carrying the
                  ``chunk_id`` that names the cut its reads came from. Its
                  fractions are chunk-local and are copied through unchanged --
                  they are already correct against that chunk's own denominator,
                  and they are what makes a bad chunk visible.
    cluster       one row per input file: that file's own ``TOTAL``, re-emitted.
    all_clusters  ONE row, aggregated over every cluster row.

ROWS OF DIFFERENT ``level`` MUST NEVER BE SUMMED. A chunk row's reads are already
counted in its cluster row, and every cluster row's reads are already counted in
the ``all_clusters`` row, so a sum across levels multiply-counts the same reads.
The output says so on its first line, because the file outlives this docstring.

WHAT THE ``all_clusters`` ROW IS, AND IS NOT
--------------------------------------------
Its ``reads_total`` is the SUM OF PER-CLUSTER PASS-1 POPULATIONS, each thinned
INDEPENDENTLY by coverage normalization. It is not a library size and it is not
comparable to a pooled-normalized or pseudobulk total, because thinning is
depth-dependent: normalizing one deep pooled bam removes far more than
normalizing N shallow per-cluster bams. MEASURED on the single-cell fixture:
~95,064 reads pooled-normalized against 198,415 summed over 32 clusters. Both
figures are correct and they are ~2x apart.

Left unlabelled, that number gets read as the library size, or diffed against the
pseudobulk file and reported as a defect. So it is labelled twice over: in the
output's header comment, and in a ``population`` column carrying the caller's
``--population`` label on every chunk and cluster row and ``sum_of_<label>`` on
the aggregate. The column is what makes the non-comparability machine-visible --
two collated files whose ``population`` labels differ must not have their totals
compared, and a consumer can now see that without parsing prose.

``--population`` has no default worth guessing at, so it defaults to
``unspecified``: a wrong label would be worse than an absent one, and the caller
that knows which normalization produced its inputs is the only one that can say.
Pass ``cluster_thinned`` for per-cluster phases and ``pooled_thinned`` for an
init/pseudobulk phase.

FRACTIONS ARE RECOMPUTED, NEVER AVERAGED. The ``all_clusters`` row's fractions
are summed numerator over summed denominator. Averaging the per-cluster fractions
would weight a 141-read cluster the same as a 16,167-read one, which is a
different and wrong number that looks entirely plausible next to the right one.

PSEUDOBULK IS NOT A CLUSTER. Do not pass an init/pseudobulk summary alongside the
per-cluster ones, for the reason given above: collated into one table, two totals
that legitimately disagree by 2x invite a reconciliation that does not exist.
Report pseudobulk as its own file: this script with that one input,
``--no_all_clusters``, and ``--population pooled_thinned``.

A MISSING OR UNREADABLE INPUT IS REFUSED, naming the file and its cluster, and no
output is written. This matches ``ChunkedRun.merge_read_assignment_summaries``,
which refuses a missing per-unit summary for the same reason: a collated total
that silently omits a cluster is an unstated subset of the experiment reported as
though it were the whole of it.
"""

import argparse
import collections
import csv
import os
import sys

sys.path.insert(
    0,
    os.path.sep.join(
        [os.path.dirname(os.path.realpath(__file__)), "..", "..", "pylib"]
    ),
)

import ChunkedRun  # noqa: E402  (path insert must precede the import)

LEVEL_CHUNK = "chunk"
LEVEL_CLUSTER = "cluster"
LEVEL_ALL_CLUSTERS = "all_clusters"

# Ahead of the inherited schema. ``row_type`` is REPLACED by ``level`` rather than
# kept beside it: two columns saying what a row is, disagreeing on vocabulary
# ("worker"/"TOTAL" against "chunk"/"cluster"/"all_clusters"), is one column too
# many. ``chunk_id`` moves up here because it is now an identifier of the row and
# not one of its counters.
LEADING_FIELDS = ("level", "cluster_id", "chunk_id")

# AFTER the inherited schema, deliberately: the columns from ``contig_acc`` onward
# stay contiguous and in their original order, so a reader who knows the
# per-invocation summary can read this table's tail unchanged.
POPULATION_FIELD = "population"

DEFAULT_POPULATION = "unspecified"

# The counters that sum. Imported rather than re-listed so this file and the
# per-invocation merger cannot drift into aggregating different column sets; a
# column outside it is left empty on an aggregate row exactly as that merger
# leaves it empty on its TOTAL.
SUMMABLE_FIELDS = ChunkedRun._SUMMARY_TOTAL_KEYS

# Read as ``key=count`` and summed BY KEY. Concatenating instead would repeat the
# same reason once per cluster and make the field unreadable.
REJECTIONS_FIELD = "rescue_alignment_rejections"

# Required of every input. Not the whole schema -- columns are handled by NAME
# throughout, so an input carrying extra counters collates fine -- but without
# these three there is no way to tell a worker row from a TOTAL or to compute a
# denominator.
REQUIRED_FIELDS = ("row_type", "contig_acc", "reads_total")

HEADER_COMMENT_LINES = (
    "# ROWS OF DIFFERENT level MUST NOT BE SUMMED. level is one of chunk | "
    "cluster | all_clusters, and they are nested: a chunk row's reads are already "
    "counted in its cluster row, and every cluster row's reads are already counted "
    "in the all_clusters row.",
    "# THE all_clusters ROW IS A SUM OF PER-CLUSTER PASS-1 POPULATIONS, EACH "
    "THINNED INDEPENDENTLY BY COVERAGE NORMALIZATION. Its reads_total is NOT a "
    "library size and MUST NOT be compared to a pooled-normalized or "
    "pseudobulk/init total: thinning is depth-dependent, so normalizing one deep "
    "pooled bam removes far more than normalizing N shallow per-cluster bams. "
    "Measured: ~95,064 reads pooled-normalized against 198,415 summed over 32 "
    "clusters -- both correct, ~2x apart.",
    "# The population column names which normalization produced these reads, so "
    "that non-comparability is machine-visible: two collated files whose "
    "population labels differ must not have their totals compared. The "
    "all_clusters row carries sum_of_<label>.",
    "# Every fraction on the all_clusters row is recomputed as summed numerator "
    "over summed denominator, never averaged across clusters. chunk rows carry "
    "chunk-local fractions.",
)


class CollationError(Exception):
    """An input this script refuses to collate around."""


def _read_summary(cluster_id, path):
    """One per-invocation summary as ``(fieldnames, worker_rows, total_row)``.

    Comment lines are dropped before parsing so a collated file, or any summary
    that grows a provenance header, can be read back by this same function.
    """

    if not os.path.exists(path):
        raise CollationError(
            "cluster {!r}: no read-assignment summary at {}. Every named cluster "
            "must contribute, so a collated total cannot be published without "
            "it".format(cluster_id, path)
        )

    try:
        with open(path, "rt", newline="") as fh:
            lines = [line for line in fh if not line.startswith("#")]
    except OSError as err:
        raise CollationError(
            "cluster {!r}: cannot read {}: {}".format(cluster_id, path, err)
        )

    reader = csv.DictReader(lines, delimiter="\t")
    fieldnames = list(reader.fieldnames or [])
    if not fieldnames:
        raise CollationError(
            "cluster {!r}: {} has no header row".format(cluster_id, path)
        )
    missing = [f for f in REQUIRED_FIELDS if f not in fieldnames]
    if missing:
        raise CollationError(
            "cluster {!r}: {} is missing column(s) {}, so its rows cannot be "
            "collated".format(cluster_id, path, ", ".join(missing))
        )

    workers = []
    totals = []
    for row in reader:
        if (row.get("row_type") or "").strip() == "TOTAL":
            totals.append(row)
        else:
            workers.append(row)

    # Exactly one, because the cluster row IS this row. Zero means the file was
    # truncated or written by something that is not LRAA's merger, and computing a
    # replacement from the worker rows would quietly paper over that.
    if len(totals) != 1:
        raise CollationError(
            "cluster {!r}: {} carries {} TOTAL row(s); exactly one is expected, "
            "and it is what the cluster-level row reports".format(
                cluster_id, path, len(totals)
            )
        )

    return fieldnames, workers, totals[0]


def _count(row, field, where):
    """One counter cell as an int, refusing a cell that is not one."""

    raw = (row.get(field) or "").strip()
    if not raw:
        return 0
    try:
        return int(raw)
    except ValueError:
        raise CollationError(
            "{}: column {} holds {!r}, which is not a count".format(where, field, raw)
        )


def _add_rejections(into, row):
    for item in (row.get(REJECTIONS_FIELD) or "").split(","):
        item = item.strip()
        if not item or "=" not in item:
            continue
        reason, _, count = item.partition("=")
        try:
            into[reason] = into.get(reason, 0) + int(count)
        except ValueError:
            continue


def _fill_fractions(out_row, fieldnames, totals):
    """Every ``frac_`` column, from the SUMMED counters.

    The three shapes are the same three the per-invocation merger recognises: a
    fraction of ``reads_total``, and the two ``_of_requested`` fractions whose
    denominator is ``reads_rescue_requested``. A ``frac_`` column matching none of
    them is left EMPTY rather than filled with a plausible zero.
    """

    denominator = totals.get("reads_total", 0)
    requested = totals.get("reads_rescue_requested", 0)

    def ratio(numerator, denom):
        if denom <= 0:
            return "0.000000"
        return "{:.6f}".format(float(numerator) / float(denom))

    for field in fieldnames:
        if not field.startswith("frac_"):
            continue
        counter = "reads_" + field[len("frac_") :]
        if counter in totals:
            out_row[field] = ratio(totals[counter], denominator)
        elif field == "frac_rescue_rescued_of_requested":
            out_row[field] = ratio(totals.get("reads_rescue_rescued", 0), requested)
        elif field == "frac_rescue_unrescued_of_requested":
            out_row[field] = ratio(totals.get("reads_rescue_unrescued", 0), requested)
        else:
            out_row[field] = ""


def collate(
    pairs, output_path, emit_all_clusters=True, population=DEFAULT_POPULATION
):
    """Write the collated table for ``pairs`` of ``(cluster_id, summary_path)``.

    Nothing is written until every input has been read and accepted, so a refusal
    leaves no partial file for a workflow to mistake for a result.
    """

    if not pairs:
        raise CollationError("no per-cluster summaries were given")

    population = (population or "").strip()
    if not population or "\t" in population:
        raise CollationError(
            "--population {!r} is not a usable label; it names the normalization "
            "these reads came from and is written on every row".format(population)
        )

    seen_clusters = collections.OrderedDict()
    for cluster_id, path in pairs:
        if not cluster_id.strip():
            raise CollationError(
                "a cluster id is empty; {} names no cluster".format(path)
            )
        if cluster_id in seen_clusters:
            raise CollationError(
                "cluster {!r} is named twice ({} and {}); two files claiming one "
                "cluster would be reported as one cluster's rows".format(
                    cluster_id, seen_clusters[cluster_id], path
                )
            )
        seen_clusters[cluster_id] = path

    in_fieldnames = None
    out_rows = []
    totals = {key: 0 for key in SUMMABLE_FIELDS}
    rejections = collections.OrderedDict()

    for cluster_id, path in pairs:
        fieldnames, workers, total_row = _read_summary(cluster_id, path)
        if in_fieldnames is None:
            in_fieldnames = fieldnames
        elif fieldnames != in_fieldnames:
            # Every column is handled by name, so a differing schema would place a
            # counter under the wrong heading rather than fail. Refuse instead.
            raise CollationError(
                "cluster {!r}: {} has {} column(s) against {} in the first summary "
                "read, so the two cannot share a table".format(
                    cluster_id, path, len(fieldnames), len(in_fieldnames)
                )
            )

        for row in workers:
            out = {f: row.get(f, "") for f in fieldnames}
            out["level"] = LEVEL_CHUNK
            out["cluster_id"] = cluster_id
            # Empty when the input predates the column. Not guessed at from
            # unit_id, which appends the orientation and so names a strand.
            out["chunk_id"] = (row.get("chunk_id") or "").strip()
            out[POPULATION_FIELD] = population
            out_rows.append(out)

        cluster_out = {f: total_row.get(f, "") for f in fieldnames}
        cluster_out["level"] = LEVEL_CLUSTER
        cluster_out["cluster_id"] = cluster_id
        # A cluster spans every chunk it was cut into, so no single chunk id is
        # true of it.
        cluster_out["chunk_id"] = ""
        cluster_out[POPULATION_FIELD] = population
        out_rows.append(cluster_out)

        # THE CLUSTER TOTAL, not the sum of its worker rows. It is the number that
        # cluster published, and on a non-chunked cluster it is the only row that
        # spans its whole invocation.
        where = "cluster {!r} TOTAL row of {}".format(cluster_id, path)
        for key in SUMMABLE_FIELDS:
            if key in fieldnames:
                totals[key] += _count(total_row, key, where)
        _add_rejections(rejections, total_row)

    out_fieldnames = (
        list(LEADING_FIELDS)
        + [f for f in in_fieldnames if f not in LEADING_FIELDS and f != "row_type"]
        + [POPULATION_FIELD]
    )

    if emit_all_clusters:
        all_row = {f: "" for f in out_fieldnames}
        all_row["level"] = LEVEL_ALL_CLUSTERS
        # Both empty: this row is every cluster and every chunk.
        all_row["cluster_id"] = ""
        all_row["chunk_id"] = ""
        # DERIVED from the row-level label rather than taken separately, so the
        # aggregate cannot claim a population its own rows do not.
        all_row[POPULATION_FIELD] = "sum_of_{}".format(population)
        if "contig_acc" in all_row:
            all_row["contig_acc"] = "TOTAL"
        if "contig_strand" in all_row:
            all_row["contig_strand"] = "."
        for key in SUMMABLE_FIELDS:
            if key in all_row:
                all_row[key] = str(totals[key])
        _fill_fractions(all_row, out_fieldnames, totals)
        if REJECTIONS_FIELD in all_row:
            all_row[REJECTIONS_FIELD] = ",".join(
                "{}={}".format(k, v) for k, v in sorted(rejections.items())
            )
        out_rows.append(all_row)

    with open(output_path, "wt", newline="") as ofh:
        for line in HEADER_COMMENT_LINES:
            ofh.write(line + "\n")
        # LF, not csv's default CRLF. The per-invocation summaries this reads are
        # CRLF by accident of ``csv.DictWriter``'s default, and a file whose
        # comment lines end LF while its rows end CRLF is worse than either --
        # ``cut -f$NF`` on it returns a trailing carriage return.
        writer = csv.DictWriter(
            ofh, fieldnames=out_fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in out_rows:
            writer.writerow({f: row.get(f, "") for f in out_fieldnames})

    counts = collections.Counter(row["level"] for row in out_rows)
    return {
        "path": output_path,
        "clusters": len(pairs),
        "population": population,
        "chunk_rows": counts[LEVEL_CHUNK],
        "cluster_rows": counts[LEVEL_CLUSTER],
        "all_clusters_row": counts[LEVEL_ALL_CLUSTERS],
        "reads_total": totals["reads_total"],
    }


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--summary",
        action="append",
        default=[],
        metavar="TSV",
        help="a per-invocation <prefix>.read_assignment.summary.tsv. Repeatable, "
        "and PAIRED BY POSITION with --cluster_id",
    )
    parser.add_argument(
        "--cluster_id",
        action="append",
        default=[],
        metavar="ID",
        help="the cluster the correspondingly-positioned --summary belongs to. "
        "Repeatable; there must be exactly as many as --summary",
    )
    parser.add_argument(
        "--output",
        required=True,
        metavar="TSV",
        help="the collated table. Written only once every input is accepted",
    )
    parser.add_argument(
        "--no_all_clusters",
        action="store_true",
        help="omit the aggregate all_clusters row. For a single-input table -- an "
        "init/pseudobulk phase reported on its own -- where that row would only "
        "restate the one cluster row",
    )
    parser.add_argument(
        "--population",
        default=DEFAULT_POPULATION,
        metavar="LABEL",
        help="which normalization produced these reads, written to the population "
        "column of every row ('sum_of_LABEL' on the all_clusters row). Pass "
        "cluster_thinned for per-cluster phases and pooled_thinned for an "
        "init/pseudobulk phase. Two collated files whose labels differ must not "
        "have their totals compared: coverage thinning is depth-dependent, so a "
        "pooled-normalized total and a sum of per-cluster totals legitimately "
        "disagree (measured ~95,064 against 198,415 over 32 clusters). Defaults "
        "to {!r} because a wrong label is worse than an absent one".format(
            DEFAULT_POPULATION
        ),
    )
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)

    # Positional pairing is checked here rather than trusted, because a workflow
    # building two arrays can drop one member of one of them and the result would
    # otherwise be every summary attributed to the wrong cluster.
    if len(args.summary) != len(args.cluster_id):
        print(
            "--summary was given {} time(s) and --cluster_id {}; they are paired by "
            "position and must be given equally often".format(
                len(args.summary), len(args.cluster_id)
            ),
            file=sys.stderr,
        )
        return 2

    pairs = list(zip(args.cluster_id, args.summary))

    try:
        result = collate(
            pairs,
            args.output,
            emit_all_clusters=not args.no_all_clusters,
            population=args.population,
        )
    except CollationError as err:
        print("refusing to collate: {}".format(err), file=sys.stderr)
        return 2

    print(
        "collated {} cluster(s) [{}]: {} chunk row(s), {} cluster row(s), {} "
        "all_clusters row, {} reads -> {}".format(
            result["clusters"],
            result["population"],
            result["chunk_rows"],
            result["cluster_rows"],
            result["all_clusters_row"],
            result["reads_total"],
            result["path"],
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
