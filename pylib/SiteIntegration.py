#!/usr/bin/env python3
# encoding: utf-8
"""Integrate two sets of TSS/PolyA site beds into one, preferring the first.

The single-cell pipeline reports boundary sites twice: from the initial (basic) catalog
and from the cluster-guided (scg) catalog. The cluster-guided calls are the deliverable,
but a basic call that no scg call reproduces is still evidence of a site, so the
integrated view keeps every primary (scg) site and supplements it with each basic site
that lies farther than the site-aggregation distance from every primary site of the
same type on the same contig and strand.

The distance is max_dist_between_alt_{TSS,polyA}_sites itself, inclusive: the window
over which LRAA's site definition (Splice_graph.aggregate_sites_within_window) absorbs
read ends into one site, so two sites called by one run are never closer than that. A
basic site within that distance of a scg site is the same site called twice, and keeping
it would put two sites closer together than either run would. (The half window,
max_dist / 2, is a different tolerance: how far a single read end may lie from a site
and still count as ending there. Using it here let basic sites 26-50 nt from a scg site
through as separate sites.)

Surviving supplement sites are NOT re-aggregated against each other; they are the basic
run's own distinct calls and already went through that run's site clustering.

Both inputs are the beds SplicePatternCollapse.write_site_bed produces, and the output
keeps those columns unchanged and appends `source` (cluster_guided | basic). Support and
transcript_ids are copied from whichever run called the site, so they are not comparable
across sources: basic support is a whole-sample value and its transcript_ids name models
of the initial catalog, not the cluster-guided one.
"""
import bisect
from collections import defaultdict

import LRAA_Globals
from SplicePatternCollapse import TSS_BED_COLUMNS, POLYA_BED_COLUMNS

PRIMARY_SOURCE = "cluster_guided"
SUPPLEMENT_SOURCE = "basic"

_WINDOW_CONFIG_KEY = {
    "TSS": "max_dist_between_alt_TSS_sites",
    "PolyA": "max_dist_between_alt_polyA_sites",
}

_COLUMNS = {"TSS": TSS_BED_COLUMNS, "PolyA": POLYA_BED_COLUMNS}


class SiteBedFormatError(Exception):
    """A bed does not carry the columns write_site_bed writes for its site type."""


def default_window(site_type):
    return LRAA_Globals.config[_WINDOW_CONFIG_KEY[site_type]]


def read_site_bed(path, site_type):
    """Rows of a site bed as lists of fields, validated against the expected header.

    The column header is the last '#'-prefixed line before the data; it must match
    exactly, because TSS and PolyA beds differ in width and a swapped pair of inputs
    would otherwise integrate silently into nonsense.
    """
    expected = _COLUMNS[site_type]
    header = None
    rows = []
    with open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                header = line[1:].split("\t")
                continue
            fields = line.split("\t")
            if len(fields) != len(expected):
                raise SiteBedFormatError(
                    "{}: row has {} fields, expected {} for a {} bed: {}".format(
                        path, len(fields), len(expected), site_type, line
                    )
                )
            rows.append(fields)
    if header != expected:
        raise SiteBedFormatError(
            "{}: column header {} does not match the {} bed columns {}".format(
                path, header, site_type, expected
            )
        )
    return rows


def _position(row):
    # write_site_bed: BED end == the 1-based site position.
    return int(row[2])


def integrate_sites(primary_rows, supplement_rows, window):
    """Return (integrated_rows, n_supplement_dropped).

    integrated_rows carry the source label appended and are sorted by
    (contig, position, strand), matching write_site_bed's order.
    """
    primary_positions = defaultdict(list)
    for row in primary_rows:
        primary_positions[(row[0], row[5])].append(_position(row))
    for positions in primary_positions.values():
        positions.sort()

    integrated = [row + [PRIMARY_SOURCE] for row in primary_rows]
    dropped = 0
    for row in supplement_rows:
        positions = primary_positions.get((row[0], row[5]), [])
        pos = _position(row)
        i = bisect.bisect_left(positions, pos - window)
        if i < len(positions) and positions[i] <= pos + window:
            dropped += 1
            continue
        integrated.append(row + [SUPPLEMENT_SOURCE])

    integrated.sort(key=lambda r: (r[0], _position(r), r[5]))
    return integrated, dropped


def write_integrated_bed(
    path, site_type, integrated, counts, window, primary_path, supplement_path
):
    with open(path, "wt") as ofh:
        ofh.write(
            "# integrated {} sites: all {} sites from {}, plus each {} site from {} "
            "farther than {} nt (config {}) from every {} site of the same "
            "contig and strand\n".format(
                site_type,
                PRIMARY_SOURCE,
                primary_path,
                SUPPLEMENT_SOURCE,
                supplement_path,
                window,
                _WINDOW_CONFIG_KEY[site_type],
                PRIMARY_SOURCE,
            )
        )
        ofh.write(
            "# counts: {primary} {p} kept, {supplement} {sk} kept, {supplement} {sd} "
            "dropped as within {hw} nt of a {primary} site\n".format(
                primary=PRIMARY_SOURCE,
                supplement=SUPPLEMENT_SOURCE,
                p=counts["primary"],
                sk=counts["supplement_kept"],
                sd=counts["supplement_dropped"],
                hw=window,
            )
        )
        ofh.write(
            "# reported_boundary_support and transcript_ids come from the run that "
            "called the site and are not comparable across sources: {} support is a "
            "whole-sample value and its transcript_ids name models of the initial "
            "catalog, not the cluster-guided one.\n".format(SUPPLEMENT_SOURCE)
        )
        ofh.write("#" + "\t".join(_COLUMNS[site_type] + ["source"]) + "\n")
        for row in integrated:
            ofh.write("\t".join(row) + "\n")


def integrate_site_beds(
    primary_bed, supplement_bed, site_type, output_bed, window=None
):
    """Integrate one site type; returns the counts written to the header."""
    if window is None:
        window = default_window(site_type)
    primary_rows = read_site_bed(primary_bed, site_type)
    supplement_rows = read_site_bed(supplement_bed, site_type)
    integrated, dropped = integrate_sites(primary_rows, supplement_rows, window)
    counts = {
        "site_type": site_type,
        "window": window,
        "primary": len(primary_rows),
        "supplement_total": len(supplement_rows),
        "supplement_kept": len(supplement_rows) - dropped,
        "supplement_dropped": dropped,
        "integrated": len(integrated),
    }
    write_integrated_bed(
        output_bed, site_type, integrated, counts, window, primary_bed, supplement_bed
    )
    return counts


# (counts key, summary column name): the file names the sources rather than the roles.
SUMMARY_COLUMNS = [
    ("site_type", "site_type"),
    ("window", "window"),
    ("primary", PRIMARY_SOURCE + "_sites"),
    ("supplement_total", SUPPLEMENT_SOURCE + "_sites"),
    ("supplement_kept", SUPPLEMENT_SOURCE + "_supplement_kept"),
    ("supplement_dropped", SUPPLEMENT_SOURCE + "_within_window_dropped"),
    ("integrated", "integrated_sites"),
]


def write_summary(path, all_counts):
    with open(path, "wt") as ofh:
        ofh.write("\t".join(name for _, name in SUMMARY_COLUMNS) + "\n")
        for counts in all_counts:
            ofh.write("\t".join(str(counts[key]) for key, _ in SUMMARY_COLUMNS) + "\n")
