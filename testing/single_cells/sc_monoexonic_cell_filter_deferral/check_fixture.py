#!/usr/bin/env python3
"""Asserts the chrTest fixture's expected outcome against one LRAA run's outputs.

Locus A (g:chrTest:+:comp-1) has 8 distinct cells on the full bam and only 2 on the
splice-graph bam, so its fate depends on WHICH BAM the draft quantification read -- and
that is now the same bam the arm's final quant estimates from: the splice-graph bam under
--stream_reads, the full bam without it. Its expected fate is therefore an ARGUMENT here
(--locus_a) rather than a constant, and the Makefile passes the value for each arm. Locus
B (comp-2, genuinely 2 distinct cells on both bams) must be FILTERED either way, and locus
C (comp-3, 8 distinct cells on both) must be RETAINED with its full 10-read support either
way, unaffected by whatever happened at A or B.

This remains a regression test for the min_monoexonic_supporting_cells filter's cell-count
EVIDENCE SOURCE, not for the filter's own bar/threshold/reference-exemption logic, which
pylib unit tests already cover. What changed is that the evidence source is no longer
unconditionally the full bam: see the Makefile header for why counting cells in thinned
evidence is accepted, and for the narrow case where it still costs a real isoform.

Reported counts are full-bam in both arms regardless -- streaming sums pass 2 over it,
non-streaming quantifies it directly -- so a RETAINED locus A must show all 10 reads, not
the splice-graph bam's 4.
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
    ap.add_argument(
        "--locus_a",
        required=True,
        choices=("retained", "filtered"),
        help="expected fate of locus A, which differs by arm: its cells are counted "
        "from whichever bam that arm's draft quant reads",
    )
    args = ap.parse_args()

    rows = parse_quant_expr(args.quant_expr)

    locus_a = ("g:chrTest:+:comp-1", "t:chrTest:+:comp-1:iso-1")
    locus_b = ("g:chrTest:+:comp-2", "t:chrTest:+:comp-2:iso-1")
    locus_c = ("g:chrTest:+:comp-3", "t:chrTest:+:comp-3:iso-1")

    if args.locus_a == "retained":
        assert locus_a in rows, (
            "[{}] locus A (8 distinct cells on the full bam, 2 on the splice-graph bam) "
            "was filtered, but this arm's draft quant reads the FULL bam, so the gate "
            "should have seen 8 cells. Rows present: {}".format(
                args.label, sorted(rows)
            )
        )
    else:
        assert locus_a not in rows, (
            "[{}] locus A survived, but this arm's draft quant reads the SPLICE-GRAPH "
            "bam, where only 2 of its 8 cells are present -- below the default bar of "
            "5, so it should have been filtered. Either the gate's evidence source "
            "moved back to the full bam or the fixture's cell counts changed. Rows "
            "present: {}".format(args.label, sorted(rows))
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

    all_reads_c = float(rows[locus_c][3])
    assert all_reads_c == 10.0, (
        "[{}] locus C quant.expr all_reads = {}, expected 10.0".format(
            args.label, all_reads_c
        )
    )

    track_counts = parse_tracking_read_counts(args.tracking)

    if args.locus_a == "retained":
        # Reported counts come from the FULL bam in both arms -- streaming sums pass 2
        # over it, non-streaming quantifies it directly -- so a surviving locus A must
        # show all 10 of its reads, not the splice-graph bam's 4. That is the assertion
        # separating "the gate read the right bam" from "the gate read the right bam and
        # then reported the wrong one".
        all_reads_a = float(rows[locus_a][3])
        assert all_reads_a == 10.0, (
            "[{}] locus A quant.expr all_reads = {}, expected 10.0 (the full bam's read "
            "count, not the splice-graph bam's 4)".format(args.label, all_reads_a)
        )
        assert track_counts.get("t:chrTest:+:comp-1:iso-1") == 10, (
            "[{}] locus A tracking file carries {} rows, expected 10 (one per full-bam "
            "read) -- tracking and quant.expr disagree about how many reads this "
            "isoform saw".format(
                args.label, track_counts.get("t:chrTest:+:comp-1:iso-1")
            )
        )
    else:
        assert "t:chrTest:+:comp-1:iso-1" not in track_counts, (
            "[{}] locus A has {} tracking rows despite being filtered from quant.expr "
            "-- tracking and quant.expr disagree about which isoforms survived".format(
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
        "[{}] OK: locus A {} (as expected for this arm's draft-quant bam), locus B "
        "filtered (2 distinct cells on both bams), locus C unaffected "
        "(10/10 reads)".format(args.label, args.locus_a)
    )


if __name__ == "__main__":
    try:
        main()
    except AssertionError as e:
        print("FAIL: {}".format(e), file=sys.stderr)
        sys.exit(1)
