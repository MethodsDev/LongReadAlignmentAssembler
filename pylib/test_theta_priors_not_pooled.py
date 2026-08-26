#!/usr/bin/env python3

"""Pass 1 iterates the PRIORS bam, never the shared splice-graph bam.

WHAT THIS FILE ADDS OVER THE EXISTING COVERAGE
----------------------------------------------
`pylib/test_chunk_sg_bam.py` proves only PLUMBING -- that a distinct priors file
is sliced per chunk and arrives at stage 5 spelled `--bam_for_priors`. Those
assertions all hold if someone reorders `_first_pass_assignment_bam`
(LRAA:5658-5698) so the sg bam wins again: the flag would still be forwarded,
the slice would still exist on disk, and theta would quietly be pooled once more.

`pylib/test_chunk_priors_bam.py::test_pass_one_reads_the_priors_bam_when_there_is_one`
does hold the precedence, and MEASURED against a scratch tree whose
`_first_pass_assignment_bam` was patched to drop the priors branch, it fails
(`assert 4 == 8`). This file's first three tests restate that contract on a
fixture built for the purpose -- sg a STRICT SUPERSET of the priors bam and four
times its size, rather than the two overlapping subsets of that fixture -- so the
two counts cannot be confused and the separation is itself asserted.

The FOURTH test is coverage that exists nowhere else. MEASURED against a second
scratch tree, patched at LRAA:6283 so pass 2 apportions under a uniform prior
instead of the estimated theta -- the flag accepted, `reads_total` still correct,
the priors' CONTENTS ignored -- `test_pass_one_reads_the_priors_bam_when_there_is_one`
PASSES, this file's first three tests PASS, and only the fourth fails. That is
the regression class no other test in this repo can see.

THE REGRESSION BEING GUARDED
----------------------------
`_first_pass_assignment_bam` resolves three roles to three files:

    --bam             the library to quantify        (pass 2)
    --bam_for_sg      splice-graph evidence          (graph only, never pass 1)
    --bam_for_priors  pass-1 theta estimation

Before the priors input existed, cluster-guided single-cell runs
(`WDL/LRAA_quant_by_cluster.wdl`) handed every cluster ONE shared sg bam, and
under `--stream_reads` -- the default -- pass 1 read that shared bam. So theta
came from POOLED evidence and each cluster's ambiguous-read apportionment was a
function of every other cluster's expression, reported as cluster-local quant.

MEASURED on one cluster of the chr19 2 Mb fixture, holding `--bam` and
`--bam_for_sg` fixed and varying ONLY the priors: own-priors vs pooled differs on
152 of 654 transcripts in read assignment and 215 of 654 in TPM (max delta 31.4
reads); own-priors vs a DIFFERENT cluster's priors differs on 163 rows. Not
cosmetic.

WHAT IS ASSERTED, AND WHY reads_total IS THE READOUT
----------------------------------------------------
The read-assignment summary's `reads_total` is a count over the pass-1
population (LRAA:5687-5691, LRAA:5727-5732 -- the same
`Util_funcs.quant_discard_reason` retention policy `get_read_alignments`
applies). It is therefore an EXACT INTEGER identifying which bam that pass
opened, with no dependence on EM arithmetic or floating-point isoform
quantities. The fixture makes the sg bam a strict superset of the priors bam and
four times its size, so neither count can be mistaken for the other.

The fourth test is the one that catches a change which accepts the flag and then
ignores its CONTENTS: two priors bams of the SAME SIZE and different reads must
move the quantification, because that is what a prior IS. It is the local
analogue of the 163-row measurement above.
"""

import subprocess
import sys
from pathlib import Path

import pysam
import pytest

REPO = Path(__file__).resolve().parents[1]
LRAA = REPO / "LRAA"
CONTIG = "chr1"
CONTIG_LEN = 10000
READ_LEN = 100


# --------------------------------------------------------------------- fixtures
#
# One contig, one gene, two overlapping single-exon isoforms. The overlap is what
# makes a read AMBIGUOUS, which is what makes theta observable at all: an
# unambiguous read is assigned the same way under any prior, so a fixture without
# a shared zone would make the fourth test below vacuous no matter how the
# priors were resolved.
#
#     t_short  1001..............1600
#     t_long              1401..............2200
#                         |-shared-|
#
# Reads are placed strictly inside one of the three zones so their compatibility
# is a containment fact rather than an overlap-fraction judgement call.

T_SHORT = (1001, 1600)
T_LONG = (1401, 2200)

# 0-based starts; each read spans [start+1, start+READ_LEN] 1-based.
UNIQ_SHORT_STARTS = (1010, 1090, 1170, 1250)  # < T_LONG[0], t_short only
AMBIG_STARTS = (1400, 1430, 1460, 1490)  # inside both
UNIQ_LONG_STARTS = (1650, 1750, 1850, 1950)  # > T_SHORT[1], t_long only

# How many reads each priors bam holds. Both priors bams are the SAME SIZE on
# purpose (see the fourth test): equal cardinality means a difference in output
# can only have come from WHICH reads, not how many.
PRIORS_N = 3


def _header():
    return pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}],
        }
    )


def _write_bam(path, reads, xw=None):
    """``reads`` is (name, 0-based start, aligned length, is_reverse), sorted.

    Same shape as the helper in pylib/test_chunk_sg_bam.py. ``xw`` writes the
    coverage-normalization weight on every record, which is what makes a bam
    admissible as ``--bam_for_sg`` or ``--bam_for_priors``: LRAA's
    ``_require_thinning_weights`` (LRAA:166) exits on either one whose first
    aligned record has none, because an untagged read weighs 1 and the estimate
    would silently stop being weighted.
    """

    header = _header()
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for name, start, length, reverse in reads:
            aln = pysam.AlignedSegment(header)
            aln.query_name = name
            aln.reference_id = 0
            aln.reference_start = start
            aln.mapping_quality = 60
            aln.cigarstring = "{}M".format(length)
            aln.query_sequence = "A" * length
            aln.query_qualities = pysam.qualitystring_to_array("I" * length)
            aln.is_reverse = reverse
            if xw is not None:
                aln.set_tag("XW", xw)
            fh.write(aln)
    pysam.index(str(path))
    return Path(path)


def _reads(prefix, starts):
    return [
        ("{}{}".format(prefix, i), start, READ_LEN, False)
        for i, start in enumerate(starts)
    ]


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    """Genome, annotation, and the three bam ROLES.

    Module-scoped: four LRAA subprocesses share one immutable fixture, and each
    run writes under its own output prefix.
    """

    tmp_path = tmp_path_factory.mktemp("theta_priors")

    genome = tmp_path / "genome.fa"
    genome.write_text(">{}\n{}\n".format(CONTIG, "A" * CONTIG_LEN))
    pysam.faidx(str(genome))

    gtf = tmp_path / "annot.gtf"
    gtf.write_text(
        "\n".join(
            '{}\ttest\texon\t{}\t{}\t.\t+\t.\tgene_id "g1"; transcript_id "{}";'.format(
                CONTIG, lend, rend, name
            )
            for name, (lend, rend) in (("t_short", T_SHORT), ("t_long", T_LONG))
        )
        + "\n"
    )

    uniq_short = _reads("short", UNIQ_SHORT_STARTS)
    ambig = _reads("ambig", AMBIG_STARTS)
    uniq_long = _reads("long", UNIQ_LONG_STARTS)
    library = sorted(uniq_short + ambig + uniq_long, key=lambda r: r[1])

    # --bam: the library to quantify, UNTAGGED. LRAA refuses an XW-tagged --bam
    # outright (_require_no_thinning_weights, LRAA:2332) because a thinned
    # library invalidates the read denominator every reported count is scaled by.
    reads_bam = _write_bam(tmp_path / "reads.bam", library)

    # --bam_for_sg: the shared splice-graph evidence. Deliberately a STRICT
    # SUPERSET of either priors bam and 4x its size, so no reads_total can be
    # confused for the other bam's. In the cluster-guided shape this file is the
    # pool of every cluster's normalized reads and no single cluster's.
    sg_bam = _write_bam(tmp_path / "evidence.bam", library, xw=0.25)

    # --bam_for_priors: THIS caller's own normalized reads. Two of them, same
    # size, disjoint contents, each favouring a different isoform.
    priors_short = _write_bam(
        tmp_path / "priors_short.bam", uniq_short[:PRIORS_N], xw=0.5
    )
    priors_long = _write_bam(tmp_path / "priors_long.bam", uniq_long[:PRIORS_N], xw=0.5)

    return {
        "dir": tmp_path,
        "genome": genome,
        "gtf": gtf,
        "reads_bam": reads_bam,
        "sg_bam": sg_bam,
        "priors_short": priors_short,
        "priors_long": priors_long,
        "library_size": len(library),
    }


def _retained_count(bam_path):
    """The population pass 1 would iterate over this bam.

    Every record in these fixtures is a mapped, primary, forward, mapq-60,
    single-block alignment carrying no NM, so ``quant_discard_reason``
    (pylib/Util_funcs.py:331) keeps all of them and this simple count IS the
    retained count. Derived from the file rather than hard-coded so a fixture
    edit moves the expectation with it.
    """

    with pysam.AlignmentFile(str(bam_path), "rb") as fh:
        return sum(1 for aln in fh.fetch(until_eof=True) if not aln.is_unmapped)


def _quant_rows(expr_path):
    """{transcript_id: (all_reads, TPM)} from a .quant.expr.

    Parsed rather than compared as bytes on purpose. The file's second line is a
    ``# LRAA CMD:`` header that echoes --output_prefix and --bam_for_priors, so a
    raw text comparison between two runs differs no matter what LRAA computed --
    which would make the contents test below pass while the priors bam was
    ignored. This reads the numbers instead.
    """

    rows = {}
    header = None
    for line in expr_path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split("\t")
        if header is None:
            header = fields
            continue
        rows[fields[header.index("transcript_id")]] = (
            float(fields[header.index("all_reads")]),
            float(fields[header.index("TPM")]),
        )
    assert rows, "no quant rows in {}".format(expr_path)
    return rows


def _run_lraa(inputs, tag, priors_bam=None):
    """One quant-only LRAA run. Returns (reads_total, quant_rows)."""

    workdir = inputs["dir"] / tag
    workdir.mkdir(exist_ok=True)
    cmd = [
        sys.executable,
        str(LRAA),
        "--bam",
        str(inputs["reads_bam"]),
        "--genome",
        str(inputs["genome"]),
        "--gtf",
        str(inputs["gtf"]),
        "--quant_only",
        "--no_chunk",
        "--output_prefix",
        str(workdir / "out"),
        "--num_total_reads",
        str(inputs["library_size"]),
        "--cpu_budget",
        "2",
        "--bam_for_sg",
        str(inputs["sg_bam"]),
        # The sg bam arrives pre-normalized by contract; re-normalizing anything
        # here would put a second, run-local bam between the fixture and the
        # count under test.
        "--no_norm",
        # Rescue realigns unassigned reads against the transcriptome and folds
        # the survivors into reads_total, which is a SECOND population on top of
        # the one under test. Off, so reads_total is exactly the pass-1 bam's
        # retained reads.
        "--no_rescue_unassigned_reads_via_transcriptome_alignment",
        # An all-A genome has nothing to mask, and building the mask costs a
        # minimap2 pass per run.
        "--no_rdna_mask",
    ]
    if priors_bam is not None:
        cmd += ["--bam_for_priors", str(priors_bam)]

    result = subprocess.run(
        cmd, capture_output=True, text=True, cwd=str(workdir), timeout=1800
    )
    assert result.returncode == 0, (result.stdout + result.stderr)[-6000:]

    summary = workdir / "out.LRAA.quant-only.read_assignment.summary.tsv"
    assert summary.exists(), (result.stdout + result.stderr)[-4000:]
    rows = [line.split("\t") for line in summary.read_text().splitlines()]
    header, body = rows[0], rows[1:]
    total = [row for row in body if row[header.index("row_type")] == "TOTAL"]
    assert len(total) == 1, summary.read_text()
    reads_total = int(total[0][header.index("reads_total")])

    expr = workdir / "out.LRAA.quant-only.quant.expr"
    assert expr.exists(), (result.stdout + result.stderr)[-4000:]
    return reads_total, _quant_rows(expr)


@pytest.fixture(scope="module")
def runs(inputs):
    """The three runs, executed once and shared. Each is a fresh subprocess."""

    with_priors = _run_lraa(inputs, "with_priors", inputs["priors_short"])
    without_priors = _run_lraa(inputs, "without_priors", None)
    other_priors = _run_lraa(inputs, "other_priors", inputs["priors_long"])
    return {
        "with_priors": with_priors,
        "without_priors": without_priors,
        "other_priors": other_priors,
    }


# ------------------------------------------------------- the fixture's own shape


def test_fixture_makes_the_two_counts_unmistakable(inputs):
    """Guards the guard.

    If a later fixture edit collapsed the sg and priors retained counts, every
    assertion below would pass while pass 1 read the wrong file -- the exact
    failure mode that made the original defect invisible, since each cluster
    still reported plausible numbers. So the SEPARATION is asserted first, and
    the priors set is asserted to be a strict subset, which is what makes
    "several times larger" a fact about one population rather than two
    unrelated ones.
    """

    sg_n = _retained_count(inputs["sg_bam"])
    priors_n = _retained_count(inputs["priors_short"])
    other_n = _retained_count(inputs["priors_long"])

    assert sg_n != priors_n
    assert sg_n >= 4 * priors_n, "sg bam must be unmistakably larger, not merely bigger"
    # Same cardinality, disjoint contents: the fourth test's difference can then
    # only be attributable to WHICH reads, never to how many.
    assert other_n == priors_n
    assert _names(inputs["priors_short"]).isdisjoint(_names(inputs["priors_long"]))
    # Strict superset, so the priors population is a genuine sub-sample of the
    # shared evidence rather than a disjoint file that happens to be smaller.
    assert _names(inputs["priors_short"]) < _names(inputs["sg_bam"])
    assert _names(inputs["priors_long"]) < _names(inputs["sg_bam"])


def _names(bam_path):
    with pysam.AlignmentFile(str(bam_path), "rb") as fh:
        return {aln.query_name for aln in fh.fetch(until_eof=True) if not aln.is_unmapped}


# ----------------------------------------------------- which bam pass 1 iterated


def test_pass_one_iterates_the_priors_bam(inputs, runs):
    """WITH --bam_for_priors, reads_total is the PRIORS bam's retained count.

    Drop the priors branch of `_first_pass_assignment_bam` (LRAA:5694-5695) so the
    sg bam wins and this reports the sg bam's count instead -- 4x larger -- while
    `--bam_for_priors` is still parsed, still validated, and still names a real
    file. MEASURED against that patched tree: `assert 12 == 3`.
    """

    reads_total, _ = runs["with_priors"]
    assert reads_total == _retained_count(inputs["priors_short"])
    assert reads_total != _retained_count(inputs["sg_bam"])


def test_pass_one_falls_back_to_the_sg_bam_without_priors(inputs, runs):
    """WITHOUT the flag, today's behaviour verbatim: pass 1 reads the sg bam.

    Every existing caller passes no priors bam, so this is the compatibility
    half of the contract. `--stream_reads` is on by default, which is the branch
    that selects the sg bam over the quant bam (LRAA:5696-5698).
    """

    reads_total, _ = runs["without_priors"]
    assert reads_total == _retained_count(inputs["sg_bam"])
    assert reads_total != _retained_count(inputs["priors_short"])


def test_the_two_runs_disagree(runs):
    """So neither assertion above can pass by coincidence.

    The two runs differ in nothing but the presence of `--bam_for_priors`; same
    --bam, same --bam_for_sg, same annotation, same denominator.
    """

    with_priors, _ = runs["with_priors"]
    without_priors, _ = runs["without_priors"]
    assert with_priors != without_priors


# --------------------------------------------- the priors' CONTENTS are consumed


def test_a_different_priors_subset_moves_the_quantification(inputs, runs):
    """The flag must be OBEYED, not merely accepted.

    Two priors bams, same cardinality (asserted in
    test_fixture_makes_the_two_counts_unmistakable), disjoint reads, one
    favouring each isoform of the overlapping pair. The four AMBIG reads in
    --bam are compatible with both isoforms, so pass 2 apportions them by the
    theta pass 1 estimated -- and a theta estimated from t_short-only reads must
    apportion them differently from one estimated from t_long-only reads.

    This is the local analogue of the measurement in this module's docstring: a
    different cluster's priors moved 163 of 654 rows. A change that read the
    priors bam's PATH but not its CONTENTS -- e.g. counting it for reads_total
    while estimating theta elsewhere -- passes every test above and fails here.

    reads_total is deliberately NOT the readout here: both priors bams hold the
    same number of reads, so it cannot distinguish them. And the comparison is on
    PARSED numbers, not the file text -- the `# LRAA CMD:` header echoes the
    priors path, so a byte comparison would differ unconditionally.

    MEASURED on this fixture: t_short-favouring priors give t_short/t_long
    all_reads 8.0/0.0, t_long-favouring priors give 0.0/8.0, and no priors at all
    gives 7.8/4.2. Asserted as a relation, not as those literals, so the test
    survives a change to LRAA's EM without going vacuous.
    """

    _, quant_short = runs["with_priors"]
    _, quant_long = runs["other_priors"]

    # Same models quantified both times: the difference below is in the numbers,
    # not in which rows exist.
    assert set(quant_short) == set(quant_long) == {"t_short", "t_long"}
    assert quant_short != quant_long

    # And the difference has the SIGN the priors imply, which a coincidental
    # perturbation would not: each priors bam holds reads unique to one isoform,
    # so that isoform must come out ahead of where the other priors put it.
    assert quant_short["t_short"][0] > quant_long["t_short"][0]
    assert quant_long["t_long"][0] > quant_short["t_long"][0]
