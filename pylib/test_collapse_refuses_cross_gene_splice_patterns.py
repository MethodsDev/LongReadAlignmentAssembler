"""`collapse_LRAA_GTF_by_splice_pattern.py` must not emit duplicate transcript_ids.

The collapsed transcript_id is minted from the splice pattern (plus any gene_symbol^
prefix), while isoforms are grouped by (gene_id, splice_pattern). So an intron pattern
carried by two gene_ids produced TWO groups that minted the SAME transcript_id, and both
were written: gtf has no uniqueness constraint, so nothing downstream objected.
MEASURED on chr8 of a 14-input single-cell merge before the fix -- 4 conflicting
patterns, 4 duplicate ids -- while the conflicts report named them and the run exited 0.

The conflict is refused rather than reconciled. Choosing a gene_id here would be a gene
REASSIGNMENT, and gene identity is decided by clustering: the expression matrices key
their features on the gene_id in the quant files (build_LRAA_expr_matrices.py:132-135),
so a gene_id invented by this script would appear in the collapsed gtf and nowhere else.
"""
import os
import re
import subprocess
import sys

import pytest

SCRIPT = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "..",
    "util",
    "collapse_LRAA_GTF_by_splice_pattern.py",
)


def _gtf_records(gene_id, transcript_id, exons):
    """A transcript row plus its exon rows, which is all the parser needs."""
    lines = []
    attrs = 'gene_id "{}"; transcript_id "{}";'.format(gene_id, transcript_id)
    lines.append(
        "\t".join(
            [
                "chr1",
                "LRAA",
                "transcript",
                str(exons[0][0]),
                str(exons[-1][1]),
                ".",
                "+",
                ".",
                attrs,
            ]
        )
    )
    for lend, rend in exons:
        lines.append(
            "\t".join(
                ["chr1", "LRAA", "exon", str(lend), str(rend), ".", "+", ".", attrs]
            )
        )
    return lines


# Same intron (1101,1999) either way; only the outer boundaries differ, which is
# exactly the case the collapse exists to merge.
SHARED_CHAIN_A = [(1000, 1100), (2000, 2100)]
SHARED_CHAIN_B = [(1050, 1100), (2000, 2200)]


def _write_gtf(path, records):
    path.write_text("\n".join(records) + "\n")
    return str(path)


def _run(gtf_path, out_path):
    return subprocess.run(
        [sys.executable, SCRIPT, "--gtf", gtf_path, "--output_gtf", str(out_path)],
        capture_output=True,
        text=True,
    )


def _transcript_ids(gtf_path):
    ids = []
    with open(gtf_path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            match = re.search(r'transcript_id\s+"([^"]*)"', cols[8])
            if match:
                ids.append(match.group(1))
    return ids


def test_same_splice_pattern_under_two_gene_ids_is_refused(tmp_path):
    gtf = _write_gtf(
        tmp_path / "conflict.gtf",
        _gtf_records("geneA", "geneA.1", SHARED_CHAIN_A)
        + _gtf_records("geneB", "geneB.1", SHARED_CHAIN_B),
    )
    out = tmp_path / "collapsed.gtf"

    result = _run(gtf, out)

    assert result.returncode != 0, "a duplicate-id-producing input must not collapse"
    assert "more than one gene_id" in (result.stdout + result.stderr)


def test_refusal_happens_before_the_output_gtf_is_written(tmp_path):
    """A partial gtf is worse than none: it looks like a collapse that worked."""
    gtf = _write_gtf(
        tmp_path / "conflict.gtf",
        _gtf_records("geneA", "geneA.1", SHARED_CHAIN_A)
        + _gtf_records("geneB", "geneB.1", SHARED_CHAIN_B),
    )
    out = tmp_path / "collapsed.gtf"

    _run(gtf, out)

    assert not out.exists(), "no output file may exist after a refusal"


def test_the_conflict_is_named_in_the_report(tmp_path):
    """The report is how a caller finds what to fix, so it must exist on refusal."""
    gtf = _write_gtf(
        tmp_path / "conflict.gtf",
        _gtf_records("geneA", "geneA.1", SHARED_CHAIN_A)
        + _gtf_records("geneB", "geneB.1", SHARED_CHAIN_B),
    )
    out = tmp_path / "collapsed.gtf"

    _run(gtf, out)

    report = tmp_path / "collapsed.gtf.gene_conflicts.tsv"
    assert report.exists()
    body = report.read_text().splitlines()[1:]
    assert len(body) == 2, "both gene_ids carrying the pattern must be listed"
    assert {row.split("\t")[3] for row in body} == {"geneA", "geneB"}


def test_the_same_pattern_within_one_gene_collapses(tmp_path):
    """The control: identical chains under ONE gene_id are what merging is for."""
    gtf = _write_gtf(
        tmp_path / "clean.gtf",
        _gtf_records("geneA", "geneA.1", SHARED_CHAIN_A)
        + _gtf_records("geneA", "geneA.2", SHARED_CHAIN_B),
    )
    out = tmp_path / "collapsed.gtf"

    result = _run(gtf, out)

    assert result.returncode == 0, result.stdout + result.stderr
    ids = _transcript_ids(str(out))
    assert len(ids) == 1, "the two isoforms must collapse to one"
    assert len(set(ids)) == len(ids)

    report = tmp_path / "collapsed.gtf.gene_conflicts.tsv"
    assert report.exists(), "the report is a declared output even with no conflicts"
    assert report.read_text().splitlines()[1:] == []


def test_distinct_patterns_keep_distinct_ids(tmp_path):
    """Guards the uniqueness self-check against firing on well-formed input."""
    gtf = _write_gtf(
        tmp_path / "clean.gtf",
        _gtf_records("geneA", "geneA.1", SHARED_CHAIN_A)
        + _gtf_records("geneB", "geneB.1", [(5000, 5100), (7000, 7100)])
        + _gtf_records("geneC", "geneC.1", [(9000, 9100)]),  # single-exon, passed through
    )
    out = tmp_path / "collapsed.gtf"

    result = _run(gtf, out)

    assert result.returncode == 0, result.stdout + result.stderr
    ids = _transcript_ids(str(out))
    assert len(ids) == 3
    assert len(set(ids)) == 3
