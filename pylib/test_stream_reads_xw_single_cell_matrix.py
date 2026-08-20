#!/usr/bin/env python3

"""End-to-end proof: --stream_reads restores completeness, --use_XW_read_weights_for_quant
restores the split, for single-cell quant-only -- measured through the real, unmodified
util/sc/singlecell_tracking_to_sparse_matrix.py, not a hand-built stand-in for it.

Fixture: one gene (GENE_DEEP) with two isoforms sharing a first exon and diverging at a
second, plus a set of unspliced reads landing entirely inside the shared exon (ambiguous
between the two isoforms, which is what lets an EM theta bias actually move the matrix).
T1's second exon is deep (2000 reads); T2's is unique to it and stays below the
normalization target (40 reads), so normalize_bam_by_strand.py's per-window acceptance
probability provably thins them at different rates -- a scarce junction is kept whole
unconditionally (see normalize_bam_by_strand._acceptance_probability), while T1's deep
junction is thinned toward the target. A separate GENE_CTRL locus, well below the target,
is included as a regression guard that normalization must not touch. Reads carry CB/XM
tags across four barcodes.

Four runs, not three, so completeness (Run D) and split-correction (Run C) are
demonstrated separately rather than conflated -- see
XW_read_weight_stream_reads_plan.md, "Track 0" / "Track 1".

DEVIATION FROM THE PLAN'S LITERAL RUN B COMMAND -- confirmed empirically, not assumed:
Plain `--quant_only` (without --stream_reads) never reads the coverage-normalized bam for
quantification at all. `bam_file_for_quant = bam_filename` unconditionally (LRAA, near
"bam_file_for_quant = bam_filename"), and the auto-normalized `bam_file_for_sg` only feeds
the read-population pass when streaming is on: `bam_file_for_pass1 = bam_file_for_sg if
LRAA_Globals.config["stream_reads"] else bam_file_for_quant` (LRAA, inside
run_quant_only). Concretely: `LRAA --quant_only --library_type single_cell
--no_rescue_unassigned_reads_via_transcriptome_alignment --bam <bam>` on the SAME bam as
Run A was run and is BYTE-IDENTICAL to Run A -- no undercounting, because normalization's
only output (bam_file_for_sg) is never consulted. That contradicts this plan's own
description of Run B ("normalization on ... expect undercounting vs. Run A"), and is
exactly the class of finding the task instructions ask to report rather than silently
paper over.

The realistic topology where --use_XW_read_weights_for_quant has ever had an effect
without --stream_reads is the one implied by the flag's own existence: `--bam` is ITSELF
already a coverage-thinned bam (e.g. produced by an upstream normalize_bam_by_strand.py
pass, exactly what the chunked pipeline does per chunk before invoking a quant-only
stage) -- XW weighting reads the XW tag off whichever bam it is handed, and a bam LRAA
has not normalized carries no such tag. Run B below therefore supplies a bam pre-thinned
by util.normalize_bam_by_strand.sift_bam directly, using the SAME tool and the SAME
default parameters LRAA's own --normalize_max_cov_level uses internally (target 1000,
window 100, seed 42, origin 0, min_per_id 80 -- LRAA_Globals.config's
normalize_max_cov_level / chunk_depth_window / chunk_random_seed / chunk_grid_origin /
min_per_id defaults), rather than the plan's literal flag list on the unmodified bam,
which measurably makes no difference at all.

NOT COVERED, stated explicitly rather than silently treated as verified (per the plan's
own "Still open" item): this fixture's GENE_DEEP component contains exactly one gene (two
isoforms), so it covers the isoform-level split *within* a gene, not the gene-level split
within a multi-gene component sharing ambiguous reads (paralogs, overlapping loci). The
narrower claim this plan actually makes -- XW corrects the split within a shared
multi-gene component -- has no fixture coverage here and remains traced from code, not
tested.
"""

import csv
import os
import subprocess
import sys
from collections import defaultdict

import pysam
import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
LRAA = os.path.join(REPO_ROOT, "LRAA")
CONVERTER = os.path.join(REPO_ROOT, "util", "sc", "singlecell_tracking_to_sparse_matrix.py")

for _sub in ("util", "pylib"):
    _p = os.path.join(REPO_ROOT, _sub)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import normalize_bam_by_strand as norm  # noqa: E402
import LRAA_Globals  # noqa: E402

CONTIG = "chr1"
CONTIG_LEN = 20000
BARCODES = ["CELL1", "CELL2", "CELL3", "CELL4"]

# T1's second exon (3001-3200) is deep: normalization thins its junction toward the
# target. T2's second exon (8001-8200) is unique to T2 and stays below the target, so its
# junction is kept whole (normalize_bam_by_strand._acceptance_probability returns 1.0 for
# any junction supported below normalize_max_cov_level, unconditionally on window depth).
N_T1 = 2000
N_T2 = 40
# Unspliced, landing entirely inside the shared exon1 [1001,1200]: compatible with BOTH
# isoforms, so their EM split is theta-dependent -- this is what lets a theta bias move
# the matrix at all.
N_AMBIG = 400
# Control locus, well below the target: normalization must not touch it either way.
N_CTRL = 40


def _header():
    return pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"},
         "SQ": [{"SN": CONTIG, "LN": CONTIG_LEN}]}
    )


def _mk_read(hdr, name, ref_start, cigar, barcode, umi):
    a = pysam.AlignedSegment(hdr)
    a.query_name = name
    a.reference_id = 0
    a.reference_start = ref_start
    a.mapping_quality = 60
    a.cigarstring = cigar
    span = sum(l for op, l in a.cigartuples if op in (0, 1))
    a.query_sequence = "A" * span
    a.query_qualities = pysam.qualitystring_to_array("I" * span)
    a.set_tag("CB", barcode)
    a.set_tag("XM", umi)
    return a


def _build_bam(path):
    hdr = _header()
    reads = []
    # T1: exon1 (1001-1200) -> exon2 (3001-3200), intron length 1800.
    for i in range(N_T1):
        bc = BARCODES[i % len(BARCODES)]
        reads.append(_mk_read(hdr, f"t1_{i}", 1100, "100M1800N100M", bc, f"U{i}"))
    # T2: exon1 (1001-1200) -> exon3 (8001-8200), intron length 6800.
    for i in range(N_T2):
        bc = BARCODES[i % len(BARCODES)]
        reads.append(_mk_read(hdr, f"t2_{i}", 1100, "100M6800N100M", bc, f"U{i}"))
    # Ambiguous: unspliced, fully inside exon1 [1001,1200] (0-based [1000,1200)).
    for i in range(N_AMBIG):
        bc = BARCODES[i % len(BARCODES)]
        start = 1020 + (i % 30)
        reads.append(_mk_read(hdr, f"amb_{i}", start, "150M", bc, f"U{i}"))
    # Control gene: single exon 15001-15200, well below the normalization target.
    for i in range(N_CTRL):
        bc = BARCODES[i % len(BARCODES)]
        start = 15020 + (i % 30)
        reads.append(_mk_read(hdr, f"ctrl_{i}", start, "150M", bc, f"U{i}"))

    with pysam.AlignmentFile(str(path), "wb", header=hdr) as fh:
        for r in sorted(reads, key=lambda x: x.reference_start):
            fh.write(r)
    pysam.index(str(path))


def _build_gtf(path):
    lines = [
        f'{CONTIG}\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "GENE_DEEP"; transcript_id "T1";',
        f'{CONTIG}\ttest\texon\t3001\t3200\t.\t+\t.\tgene_id "GENE_DEEP"; transcript_id "T1";',
        f'{CONTIG}\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "GENE_DEEP"; transcript_id "T2";',
        f'{CONTIG}\ttest\texon\t8001\t8200\t.\t+\t.\tgene_id "GENE_DEEP"; transcript_id "T2";',
        f'{CONTIG}\ttest\texon\t15001\t15200\t.\t+\t.\tgene_id "GENE_CTRL"; transcript_id "TC";',
    ]
    path.write_text("\n".join(lines) + "\n")


def _build_genome(path):
    seq = "A" * CONTIG_LEN
    with open(path, "w") as f:
        f.write(f">{CONTIG}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i : i + 60] + "\n")


def _build_thinned_bam(true_bam, thinned_bam):
    """The bam LRAA's OWN --normalize_max_cov_level would produce internally, built with
    the same tool and the same defaults, for use directly as Run B's --bam (see module
    docstring for why plain --quant_only never reads such a bam on its own).
    """
    norm.sift_bam(
        str(true_bam),
        str(thinned_bam),
        normalize_max_cov_level=LRAA_Globals.config["normalize_max_cov_level"],
        depth_window=LRAA_Globals.config["chunk_depth_window"],
        random_seed=LRAA_Globals.config["chunk_random_seed"],
        window_origin=LRAA_Globals.config["chunk_grid_origin"],
        min_per_id=LRAA_Globals.config["min_per_id"],
        min_mapping_quality=LRAA_Globals.config["min_mapping_quality_for_final_quant"],
    )
    pysam.index(str(thinned_bam))


def _run_lraa(tmp_path, bam, gtf, genome, prefix, *extra):
    cmd = [sys.executable, LRAA, "--quant_only", "--bam", str(bam),
           "--gtf", str(gtf), "--genome", str(genome),
           "--output_prefix", str(tmp_path / prefix)] + list(extra)
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode == 0, r.stdout + r.stderr
    tracking = list(tmp_path.glob(f"{prefix}*quant.tracking.gz"))
    assert tracking, "expected a tracking file to be produced"
    return tracking[0]


def _convert(tmp_path, tracking, prefix):
    out_prefix = tmp_path / prefix
    cmd = [sys.executable, CONVERTER, "--tracking", str(tracking),
           "--output_prefix", str(out_prefix)]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    assert r.returncode == 0, r.stdout + r.stderr
    return out_prefix


def _feature_totals(out_prefix, level):
    """feature_id -> summed UMI_counts across every cell, from the real converter output."""
    totals = defaultdict(float)
    path = f"{out_prefix}.{level}_cell_counts.tsv"
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            totals[row["feature_id"]] += float(row["UMI_counts"])
    return dict(totals)


@pytest.fixture(scope="module")
def four_runs(tmp_path_factory):
    """Builds the fixture once and runs A/B/C/D + their conversions once each.

    Returns {run_label: {"gene": {...}, "isoform": {...}}}, every value a real number
    read back from singlecell_tracking_to_sparse_matrix.py's own output files.
    """
    tmp = tmp_path_factory.mktemp("xw_stream_matrix")
    true_bam = tmp / "true.bam"
    _build_bam(true_bam)
    gtf = tmp / "a.gtf"
    _build_gtf(gtf)
    genome = tmp / "g.fa"
    _build_genome(genome)
    thinned_bam = tmp / "thinned.bam"
    _build_thinned_bam(true_bam, thinned_bam)

    no_rescue = "--no_rescue_unassigned_reads_via_transcriptome_alignment"
    results = {}

    tracking_a = _run_lraa(tmp, true_bam, gtf, genome, "A", "--no_norm",
                            "--library_type", "single_cell", no_rescue,
                            "--no_chunk", "--no_stream_reads")
    prefix_a = _convert(tmp, tracking_a, "convA")
    results["A"] = {
        "gene": _feature_totals(prefix_a, "gene"),
        "isoform": _feature_totals(prefix_a, "isoform"),
    }

    # Run B: see module docstring -- pre-thinned bam stands in for "today's default",
    # since plain --quant_only on the unmodified bam is confirmed identical to Run A.
    tracking_b = _run_lraa(tmp, thinned_bam, gtf, genome, "B", "--no_norm",
                            "--library_type", "single_cell", no_rescue,
                            "--no_chunk", "--no_stream_reads")
    prefix_b = _convert(tmp, tracking_b, "convB")
    results["B"] = {
        "gene": _feature_totals(prefix_b, "gene"),
        "isoform": _feature_totals(prefix_b, "isoform"),
    }

    tracking_d = _run_lraa(tmp, true_bam, gtf, genome, "D", "--stream_reads",
                            "--library_type", "single_cell", no_rescue,
                            "--no_chunk")
    prefix_d = _convert(tmp, tracking_d, "convD")
    results["D"] = {
        "gene": _feature_totals(prefix_d, "gene"),
        "isoform": _feature_totals(prefix_d, "isoform"),
    }

    tracking_c = _run_lraa(tmp, true_bam, gtf, genome, "C", "--stream_reads",
                            "--use_XW_read_weights_for_quant",
                            "--library_type", "single_cell", no_rescue,
                            "--no_chunk")
    prefix_c = _convert(tmp, tracking_c, "convC")
    results["C"] = {
        "gene": _feature_totals(prefix_c, "gene"),
        "isoform": _feature_totals(prefix_c, "isoform"),
    }

    return results


def test_run_a_reference_accounts_for_every_read(four_runs):
    """Ground truth: --no_norm drops nothing, so the deep locus totals N_T1+N_T2+N_AMBIG
    and the control locus totals N_CTRL, exactly.
    """
    gene = four_runs["A"]["gene"]
    assert gene["GENE_DEEP"] == pytest.approx(N_T1 + N_T2 + N_AMBIG, abs=0.5)
    assert gene["GENE_CTRL"] == pytest.approx(N_CTRL, abs=0.01)


def test_gene_level_total_today_undercounts_stream_reads_restores_it(four_runs):
    """Deep-locus gene total: B (today's default) undercounts vs. A; D and C, both
    --stream_reads, track A closely -- completeness does not depend on theta or XW.
    """
    a = four_runs["A"]["gene"]["GENE_DEEP"]
    b = four_runs["B"]["gene"]["GENE_DEEP"]
    d = four_runs["D"]["gene"]["GENE_DEEP"]
    c = four_runs["C"]["gene"]["GENE_DEEP"]

    assert a == pytest.approx(2440.0, abs=1.0)
    # Measured: B totals 1257 of A's 2440, a 48.5% shortfall.
    assert b < 0.7 * a, f"expected a real undercount, got B={b} vs A={a}"
    assert d == pytest.approx(a, rel=0.01), f"D={d} should match A={a} (full accounting)"
    assert c == pytest.approx(a, rel=0.01), f"C={c} should match A={a} (full accounting)"


def test_isoform_level_split_today_and_stream_only_diverge_xw_restores_it(four_runs):
    """Deep-locus isoform split, measured as T2's share of the gene total (robust to the
    gene-level completeness differences already covered by the previous test): B and D
    both diverge from A's split; C, the only run with XW correction, tracks it closely.
    """

    def t2_share(run):
        gene = four_runs[run]["gene"]["GENE_DEEP"]
        t2 = four_runs[run]["isoform"]["T2"]
        return t2 / gene

    share_a = t2_share("A")
    share_b = t2_share("B")
    share_d = t2_share("D")
    share_c = t2_share("C")

    # Measured: A ~= 0.0199, B ~= 0.0428 (+115%), D ~= 0.0238 (+19%), C ~= 0.0199 (+0.09%).
    assert share_a == pytest.approx(0.0199, abs=0.002)
    assert share_b > 1.5 * share_a, (
        f"B's split should diverge sharply from A's (no completeness, no correction): "
        f"B={share_b} A={share_a}"
    )
    assert share_d > 1.1 * share_a, (
        f"D's split should diverge from A's (complete, but theta uncorrected): "
        f"D={share_d} A={share_a}"
    )
    assert share_c == pytest.approx(share_a, rel=0.05), (
        f"C's split should track A's closely (complete AND theta-corrected): "
        f"C={share_c} A={share_a}"
    )


def test_isoform_level_absolute_counts_xw_restores_them(four_runs):
    """The same claim in absolute terms: C's T1 and T2 totals track A's; D's do not."""
    a = four_runs["A"]["isoform"]
    c = four_runs["C"]["isoform"]
    d = four_runs["D"]["isoform"]

    assert c["T1"] == pytest.approx(a["T1"], rel=0.01)
    assert c["T2"] == pytest.approx(a["T2"], rel=0.05)

    # Measured: D's T2 (58.0) overshoots A's (48.6) by ~19%, well outside C's 5% band.
    assert d["T2"] > 1.1 * a["T2"]


def test_control_locus_is_a_regression_guard_across_all_four_runs(four_runs):
    """A locus normalization never touches must agree across every run, streamed,
    weighted, or neither -- otherwise a change here would be masking a real regression.
    """
    ctrl = {run: four_runs[run]["gene"]["GENE_CTRL"] for run in ("A", "B", "C", "D")}
    for run, total in ctrl.items():
        assert total == pytest.approx(N_CTRL, abs=0.01), f"{run}: {total}"
