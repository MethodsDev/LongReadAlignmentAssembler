#!/usr/bin/env python3
"""Asserts the chrTest fixture's expected outcome against one LRAA run's outputs.

Locked-in behavior (see data/build_fixture.py for the fixture's shape): locus A
(g:chrTest:+:comp-1, 8 distinct cells on the full bam, only 2 on the splice-graph
bam) must be RETAINED with its full 10-read support; locus B (g:chrTest:+:comp-2,
genuinely 2 distinct cells on both bams) must be FILTERED; locus C
(g:chrTest:+:comp-3, 8 distinct cells on both bams) must be RETAINED with its full
10-read support, byte-for-byte unaffected by whatever happened at A or B.

This is a regression test for the min_monoexonic_supporting_cells filter's cell-count
EVIDENCE SOURCE (LRAA:5726-5735's mp_all_counter, built from bam_file_for_quant, the
full un-normalized bam, unconditionally in both --stream_reads and default mode) --
not a test of the filter's own bar/threshold/reference-exemption logic, which is
already covered by pylib unit tests.
"""

import argparse
import gzip
import sys


def parse_quant_expr(path):
    rows = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("gene_id"):
                continue
            fields = line.split("\t")
            gene_id, transcript_id = fields[0], fields[1]
            rows[(gene_id, transcript_id)] = fields
    return rows


def parse_tracking_read_counts(path):
    """transcript_id -> count of distinct read_name rows."""
    opener = gzip.open if path.endswith(".gz") else open
    counts = {}
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("gene_id"):
                continue
            fields = line.split("\t")
            transcript_id = fields[1]
            counts[transcript_id] = counts.get(transcript_id, 0) + 1
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quant_expr", required=True)
    ap.add_argument("--tracking", required=True)
    ap.add_argument("--label", required=True, help="run label for error messages")
    args = ap.parse_args()

    rows = parse_quant_expr(args.quant_expr)

    locus_a = ("g:chrTest:+:comp-1", "t:chrTest:+:comp-1:iso-1")
    locus_b = ("g:chrTest:+:comp-2", "t:chrTest:+:comp-2:iso-1")
    locus_c = ("g:chrTest:+:comp-3", "t:chrTest:+:comp-3:iso-1")

    assert locus_a in rows, (
        "[{}] locus A (8 distinct cells on the full bam, 2 on the splice-graph bam) "
        "was filtered -- the min-cells filter is judging cell counts from thinned "
        "splice-graph-bam evidence instead of the full bam. Rows present: {}".format(
            args.label, sorted(rows)
        )
    )
    assert locus_b not in rows, (
        "[{}] locus B (genuinely 2 distinct cells on BOTH bams) survived -- the "
        "min-cells filter is not applying its threshold at all. Rows present: {}".format(
            args.label, sorted(rows)
        )
    )
    assert locus_c in rows, (
        "[{}] locus C (always-kept control, byte-identical in both bams) was "
        "filtered -- an unrelated component was disturbed. Rows present: {}".format(
            args.label, sorted(rows)
        )
    )

    all_reads_a = float(rows[locus_a][3])
    all_reads_c = float(rows[locus_c][3])
    assert all_reads_a == 10.0, (
        "[{}] locus A quant.expr all_reads = {}, expected 10.0 (the full bam's read "
        "count, not the splice-graph bam's 4)".format(args.label, all_reads_a)
    )
    assert all_reads_c == 10.0, (
        "[{}] locus C quant.expr all_reads = {}, expected 10.0".format(
            args.label, all_reads_c
        )
    )

    track_counts = parse_tracking_read_counts(args.tracking)
    assert track_counts.get("t:chrTest:+:comp-1:iso-1") == 10, (
        "[{}] locus A tracking file carries {} rows, expected 10 (one per full-bam "
        "read) -- tracking and quant.expr disagree about how many reads this "
        "isoform saw".format(
            args.label, track_counts.get("t:chrTest:+:comp-1:iso-1")
        )
    )
    assert track_counts.get("t:chrTest:+:comp-3:iso-1") == 10, (
        "[{}] locus C tracking file carries {} rows, expected 10".format(
            args.label, track_counts.get("t:chrTest:+:comp-3:iso-1")
        )
    )
    assert "t:chrTest:+:comp-2:iso-1" not in track_counts, (
        "[{}] locus B has {} tracking rows despite being filtered from quant.expr -- "
        "tracking and quant.expr disagree about which isoforms survived".format(
            args.label, track_counts.get("t:chrTest:+:comp-2:iso-1")
        )
    )

    print(
        "[{}] OK: locus A retained (10/10 reads, 8 distinct cells via the full bam), "
        "locus B filtered (2 distinct cells on both bams), locus C unaffected "
        "(10/10 reads)".format(args.label)
    )


if __name__ == "__main__":
    try:
        main()
    except AssertionError as e:
        print("FAIL: {}".format(e), file=sys.stderr)
        sys.exit(1)
