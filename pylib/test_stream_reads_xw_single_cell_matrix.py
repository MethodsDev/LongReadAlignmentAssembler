#!/usr/bin/env python3

"""End-to-end proof that single-cell quant-only accounts for every read and lands on the
right within-gene isoform split -- measured through the real, unmodified
util/sc/singlecell_tracking_to_sparse_matrix.py, not a hand-built stand-in for it.

Fixture: one gene (GENE_DEEP) with two isoforms sharing a first exon and diverging at a
second, plus a set of unspliced reads landing entirely inside the shared exon (ambiguous
between the two isoforms, which is what lets an EM theta bias actually move the matrix).
T1's second exon is deep (2000 reads); T2's is unique to it and stays below the
normalization target (40 reads), so normalize_bam_by_strand.py's per-window acceptance
probability provably thins them at different rates -- a scarce junction is kept whole
unconditionally (see normalize_bam_by_strand._acceptance_probability), while T1's deep
junction is thinned toward the target. Measured: sift_bam at the stock target of 1000
retains 1297 of the 2480 records, i.e. roughly half the library disappears from thinned
evidence, and it disappears asymmetrically between the two isoforms. A separate
GENE_CTRL locus, well below the target, is included as a regression guard that
normalization must not touch. Reads carry CB/XM tags across four barcodes.

Three runs, each one a distinct supported way of feeding thinned evidence into quant:

  Run A -- --no_norm --no_stream_reads on the untagged full bam. Nothing is thinned
           anywhere, so this is ground truth: it must account for every read, and its
           isoform split is the number the other two runs are judged against.

  Run B -- the full untagged bam as --bam, the pre-thinned bam as --bam_for_sg, with
           --no_norm --no_stream_reads. This is "quantification against thinned
           evidence" expressed the way the inputs are now defined: --bam is the library
           everything is reported against and must NOT be thinned, while --bam_for_sg is
           the input that exists for coverage-normalized evidence and must BE thinned
           (LRAA:_require_no_thinning_weights / LRAA:_require_thinning_weights).
           _build_thinned_bam produces it through norm.sift_bam, which writes the XW tag,
           so the --bam_for_sg precondition is satisfied by construction. The splice
           graph is built from that thinned bam (LRAA, Splice_graph call inside
           run_quant_only, fed bam_file_for_sg) and divides each read's acceptance
           probability back out through XW; quantification reads bam_file_for_quant,
           which is the full bam. Nothing is re-normalized: --bam_for_sg is taken
           verbatim, whatever --normalize_max_cov_level says.

  Run C -- --stream_reads on the untagged full bam, with normalization left on, so LRAA
           thins the bam itself and hands the result to pass 1
           (bam_file_for_pass1 = bam_file_for_sg when streaming, LRAA inside
           run_quant_only) while the streaming second pass covers the full bam. This is a
           stock invocation: there is no weighting flag to pass. XW coverage-normalization
           weighting is unconditional -- _populate_read_multi_paths resolves
           use_XW_weights to True unless a caller passes weight_reads=False, which only
           discovery's pre-filter quant does.

Why the tight bands matter rather than a weighted-vs-unweighted contrast: weighting is no
longer something a run can decline, so there is no unweighted arm to diff against. The
claims are stated against Run A instead, and the bands are set tight enough to fail if
weighting regressed. Historically measured with weighting declined, T2 came out at 58.02
against a truth of 48.59 -- a +19% overshoot on the absolute count and T2's gene share at
0.0238 against 0.0199. Both are well outside the 5% bands asserted below, so a regression
that silently stopped dividing out XW would break these tests rather than slip through.

Run B needs the thinning to be real for its "matches Run A" claim to mean anything, so
the retained-record count is asserted directly rather than assumed: if sift_bam ever
stopped thinning this fixture, B would agree with A trivially and prove nothing.

NOT COVERED, stated explicitly rather than silently treated as verified: this fixture's
GENE_DEEP component contains exactly one gene (two isoforms), so it covers the isoform-
level split *within* a gene, not the gene-level split within a multi-gene component
sharing ambiguous reads (paralogs, overlapping loci). That narrower claim -- that XW
weighting corrects the split within a shared multi-gene component -- has no fixture
coverage here and remains traced from code, not tested.
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
    the same tool and the same defaults, for use as Run B's --bam_for_sg.

    sift_bam writes the XW acceptance weight on every record it keeps, which is exactly
    what --bam_for_sg requires (LRAA:_require_thinning_weights): the splice graph divides
    that weight back out, so an untagged thinned bam would be read as genuinely shallow
    evidence. Returns the retained-record count so the caller can assert the thinning was
    real rather than assume it.
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
    with pysam.AlignmentFile(str(thinned_bam)) as fh:
        return fh.count()


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
    # --emit_dense_counts because _feature_totals below reads the dense
    # *_cell_counts.tsv, which the converter no longer writes by default.
    cmd = [sys.executable, CONVERTER, "--tracking", str(tracking),
           "--output_prefix", str(out_prefix), "--emit_dense_counts"]
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
def three_runs(tmp_path_factory):
    """Builds the fixture once and runs A/B/C + their conversions once each.

    Returns {run_label: {"gene": {...}, "isoform": {...}}} plus a "records" entry holding
    the full and thinned record counts. Every quant value is a real number read back from
    singlecell_tracking_to_sparse_matrix.py's own output files.
    """
    tmp = tmp_path_factory.mktemp("xw_stream_matrix")
    true_bam = tmp / "true.bam"
    _build_bam(true_bam)
    gtf = tmp / "a.gtf"
    _build_gtf(gtf)
    genome = tmp / "g.fa"
    _build_genome(genome)
    thinned_bam = tmp / "thinned.bam"
    n_thinned = _build_thinned_bam(true_bam, thinned_bam)

    no_rescue = "--no_rescue_unassigned_reads_via_transcriptome_alignment"
    results = {"records": {"full": N_T1 + N_T2 + N_AMBIG + N_CTRL,
                           "thinned": n_thinned}}

    # Run A: nothing thinned anywhere. Ground truth for both claims below.
    tracking_a = _run_lraa(tmp, true_bam, gtf, genome, "A", "--no_norm",
                            "--library_type", "single_cell", no_rescue,
                            "--no_chunk", "--no_stream_reads")
    prefix_a = _convert(tmp, tracking_a, "convA")
    results["A"] = {
        "gene": _feature_totals(prefix_a, "gene"),
        "isoform": _feature_totals(prefix_a, "isoform"),
    }

    # Run B: thinned evidence supplied the supported way -- full library as --bam, the
    # XW-tagged thinned bam as --bam_for_sg. Handing the thinned bam to --bam is now
    # rejected outright (LRAA:_require_no_thinning_weights), because reported counts are
    # scaled by --bam and a thinned --bam silently reports against a fraction of the
    # library. --no_norm here proves --bam_for_sg is taken verbatim and never re-thinned.
    tracking_b = _run_lraa(tmp, true_bam, gtf, genome, "B", "--no_norm",
                            "--library_type", "single_cell", no_rescue,
                            "--no_chunk", "--no_stream_reads",
                            "--bam_for_sg", str(thinned_bam))
    prefix_b = _convert(tmp, tracking_b, "convB")
    results["B"] = {
        "gene": _feature_totals(prefix_b, "gene"),
        "isoform": _feature_totals(prefix_b, "isoform"),
    }

    # Run C: a stock streaming invocation. Normalization stays on, so LRAA thins the bam
    # itself for pass 1 and streams the full bam for pass 2. There is no weighting flag
    # to pass either way -- XW weighting is unconditional now -- so this run measures the
    # default rather than a configuration.
    tracking_c = _run_lraa(tmp, true_bam, gtf, genome, "C", "--stream_reads",
                            "--library_type", "single_cell", no_rescue,
                            "--no_chunk")
    prefix_c = _convert(tmp, tracking_c, "convC")
    results["C"] = {
        "gene": _feature_totals(prefix_c, "gene"),
        "isoform": _feature_totals(prefix_c, "isoform"),
    }

    return results


def test_run_a_reference_accounts_for_every_read(three_runs):
    """Ground truth: --no_norm drops nothing, so the deep locus totals N_T1+N_T2+N_AMBIG
    and the control locus totals N_CTRL, exactly.
    """
    gene = three_runs["A"]["gene"]
    assert gene["GENE_DEEP"] == pytest.approx(N_T1 + N_T2 + N_AMBIG, abs=0.5)
    assert gene["GENE_CTRL"] == pytest.approx(N_CTRL, abs=0.01)


def test_the_thinned_splice_graph_evidence_really_is_thinner(three_runs):
    """Run B's whole claim is that thinned evidence costs it nothing, which is only worth
    asserting if the evidence was actually thinned. If sift_bam's defaults ever stopped
    thinning this fixture -- a raised normalize_max_cov_level default, a changed window,
    a different acceptance rule -- Run B would agree with Run A for the uninteresting
    reason that the two bams are the same reads, and every B assertion below would go
    quietly vacuous. Measured: 1297 of 2480 records retained at the stock target of 1000.
    """
    records = three_runs["records"]
    assert records["full"] == 2480
    assert records["thinned"] == pytest.approx(1297, abs=60), (
        f"expected the stock target to drop roughly half the library, "
        f"retained {records['thinned']} of {records['full']}"
    )


def test_gene_level_totals_account_for_every_read_in_every_run(three_runs):
    """Deep-locus gene total: both supported ways of feeding thinned evidence still report
    against the whole library, so B and C match A.

    B matches because --bam_for_sg only replaces splice-graph evidence: quantification
    reads bam_file_for_quant, which is --bam. C matches because the streaming second pass
    covers the full bam even though pass 1 saw only the thinned one. The failure this
    guards is thinned evidence leaking into the reported denominator -- the pre-thinned
    --bam that is now rejected measured 1257 of A's 2440, a 48.5% shortfall.
    """
    a = three_runs["A"]["gene"]["GENE_DEEP"]
    b = three_runs["B"]["gene"]["GENE_DEEP"]
    c = three_runs["C"]["gene"]["GENE_DEEP"]

    assert a == pytest.approx(2440.0, abs=1.0)
    assert b == pytest.approx(a, rel=0.01), f"B={b} should match A={a} (full accounting)"
    assert c == pytest.approx(a, rel=0.01), f"C={c} should match A={a} (full accounting)"


def test_isoform_level_split_tracks_run_a_in_every_run(three_runs):
    """Deep-locus isoform split, measured as T2's share of the gene total (robust to any
    gene-level difference, which the previous test covers separately).

    T2's junction is scarce and survives thinning whole while T1's deep junction is thinned
    toward the target, so thinned evidence over-represents T2 unless the acceptance weight
    is divided back out. Both runs track A's split closely, which is the observable form of
    "XW weighting is applied unconditionally".

    Measured: A ~= 0.019916, B ~= 0.019916 (identical), C ~= 0.019934 (+0.09%). The 5%
    band is what makes this fail on a regression: with weighting declined the same
    fixture put T2's share at 0.0238, +19% over A.
    """

    def t2_share(run):
        gene = three_runs[run]["gene"]["GENE_DEEP"]
        t2 = three_runs[run]["isoform"]["T2"]
        return t2 / gene

    share_a = t2_share("A")
    share_b = t2_share("B")
    share_c = t2_share("C")

    assert share_a == pytest.approx(0.0199, abs=0.002)
    assert share_b == pytest.approx(share_a, rel=0.05), (
        f"B's split should track A's: thinned splice-graph evidence divides its own "
        f"acceptance weight back out. B={share_b} A={share_a}"
    )
    assert share_c == pytest.approx(share_a, rel=0.05), (
        f"C's split should track A's: complete AND theta-corrected. "
        f"C={share_c} A={share_a}"
    )


def test_isoform_level_absolute_counts_track_run_a(three_runs):
    """The same claim in absolute terms rather than as a share, so a proportional drift in
    both isoforms at once cannot hide inside the ratio.

    Measured against A's T1 2391.41 / T2 48.59: B is identical to the digit (its
    quantification bam is A's bam), C lands at T1 2391.36 / T2 48.64. The T2 band is the
    loose one because T2 is the scarce isoform the ambiguous reads get apportioned to;
    weighting declined put it at 58.02, outside this band.
    """
    a = three_runs["A"]["isoform"]
    b = three_runs["B"]["isoform"]
    c = three_runs["C"]["isoform"]

    assert a["T1"] == pytest.approx(2391.41, abs=1.0)
    assert a["T2"] == pytest.approx(48.59, abs=0.5)

    assert b["T1"] == pytest.approx(a["T1"], rel=0.01)
    assert b["T2"] == pytest.approx(a["T2"], rel=0.05)

    assert c["T1"] == pytest.approx(a["T1"], rel=0.01)
    assert c["T2"] == pytest.approx(a["T2"], rel=0.05)


def test_control_locus_is_a_regression_guard_across_every_run(three_runs):
    """A locus normalization never touches must agree across every run, streamed or not,
    with thinned splice-graph evidence or without -- otherwise a change here would be
    masking a real regression.
    """
    ctrl = {
        run: three_runs[run]["gene"]["GENE_CTRL"]
        for run in ("A", "B", "C")
    }
    for run, total in ctrl.items():
        assert total == pytest.approx(N_CTRL, abs=0.01), f"{run}: {total}"
